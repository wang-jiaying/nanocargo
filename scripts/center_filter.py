#!/usr/bin/env python3
import yaml
import os
from collections import defaultdict
import argparse

# -------------------------
# Argument parsing
# -------------------------
parser = argparse.ArgumentParser(
    description="Collapse reads and select representative sequences with IQR"
)
parser.add_argument("--config", "-c", required=True)
parser.add_argument("--paf", "-p", required=True)
parser.add_argument("--outdir", "-o", required=True)
parser.add_argument("--scov", "-s", type=float, help="Override SCOV_TH")
parser.add_argument("--ma", "-n", type=int, help="Override MA_TH")
parser.add_argument("--madelta", "-d", type=int, help="Override MA_DELTA")
parser.add_argument("--idth", "-m", type=float, help="Override ID_TH")
parser.add_argument("--minlen", type=int, help="Override MIN_READLEN")
args = parser.parse_args()

# -------------------------
# Prepare output directory
# -------------------------
outdir = os.path.abspath(args.outdir)
os.makedirs(outdir, exist_ok=True)

REP_STAT_FILE = os.path.join(outdir, "Initial_rep_stat.txt")

# -------------------------
# Load config
# -------------------------
with open(args.config) as f:
    cfg = yaml.safe_load(f)["collapse_reads"]

QCOV_TH = cfg.get("QCOV_TH", 0)
SCOV_TH = cfg.get("SCOV_TH", 0.99)
ID_TH = cfg.get("ID_TH", 0.8)
DELTA = cfg.get("DELTA", 60)
FA_COV_TH = cfg.get("FA_COV_TH", 0.99)
MA_TH = cfg.get("MA_TH", 0)
MA_DELTA = cfg.get("MA_DELTA", 1000)
SHORT_CENTER_LEN = cfg.get("SHORT_CENTER_LEN", 4000)
MIN_SHORT_CENTER_NUM = cfg.get("MIN_SHORT_CENTER_NUM", 5)
CENTER_MIN_SUPPORT = cfg.get("CENTER_MIN_SUPPORT", 3)
MIN_LA_RA = cfg.get("MIN_LA_RA", 5)

# IQR thresholds
LRA_MIN = cfg.get("LRA_MIN", 3)
IQR_MIN_SUPPORT = cfg.get("IQR_MIN_SUPPORT", 5) 
IQR_TH1 = cfg.get("IQR_TH1", 0.03)
IQR_TH2 = cfg.get("IQR_TH2", 0.08) 

# override by command line if provided
if args.scov is not None:
    SCOV_TH = args.scov
if args.ma is not None:
    MA_TH = args.ma
if args.madelta is not None:
    MA_DELTA = args.madelta
if args.idth is not None:
    ID_TH = args.idth
if args.minlen is not None and args.minlen > 0:
    MIN_READLEN = args.minlen
else:
    MIN_READLEN = 0
    if args.minlen is None:
        print("[WARNING] --minlen not provided. Setting MIN_READLEN=0; IQR normalization may be affected.")
    else:
        print(f"[WARNING] Invalid --minlen value ({args.minlen}). Setting MIN_READLEN=0; IQR normalization may be affected.")

# -------------------------
# Helper functions
# -------------------------
def classify_endpoint(qs, qe, qlen):
    if qs <= DELTA and qe < qlen - DELTA:
        return "LA"
    if qs > DELTA and qe >= qlen - DELTA:
        return "RA"
    if qs > MA_DELTA and qe < qlen - MA_DELTA:
        return "MA"
    return None

def iqr_from_pos(pos, qlen, min_len=None):
    if min_len is None:
        min_len = MIN_READLEN
    if len(pos) < IQR_MIN_SUPPORT:
        return 1.0
    pos = sorted(pos)
    n = len(pos)
    q25 = pos[int(0.25 * (n - 1))]
    q75 = pos[int(0.75 * (n - 1))]
    return (q75 - q25) / max(qlen - min_len, 1)

def compute_iqr(read, aln_by_query, qlen):
    LA_pos, RA_pos = [], []
    for qs, qe, _, _, _, identity, role, _ in aln_by_query.get(read, []):
        if role != "query" or identity <= ID_TH:
            continue
        at = classify_endpoint(qs, qe, qlen)
        if at == "LA":
            LA_pos.append(qe)
        elif at == "RA":
            RA_pos.append(qs)
    LA_iqr = iqr_from_pos(LA_pos, qlen)
    RA_iqr = iqr_from_pos(RA_pos, qlen)
    return min(LA_iqr, RA_iqr)

def endpoint_stat(read, aln_by_query):

    FAq = FAs = LA = RA = MA = 0
    LA_qcov, RA_qcov = [], []

    for qs, qe, qlen, qcov, scov, identity, role, mate in aln_by_query.get(read, []):

        if qcov > FA_COV_TH and scov > FA_COV_TH and identity > ID_TH:
            if role == "query":
                FAq += 1
            else:
                FAs += 1
            continue
        if role != "query":
            continue
        at = classify_endpoint(qs, qe, qlen)
        if at == "LA":
            LA += 1
            if identity > ID_TH:
                LA_qcov.append(qcov)
        elif at == "RA":
            RA += 1
            if identity > ID_TH:
                RA_qcov.append(qcov)
        elif at == "MA":
            MA += 1

    LA_qcov.sort(reverse=True)
    RA_qcov.sort(reverse=True)

    LRA = 0
    for i in range(min(len(LA_qcov), len(RA_qcov))):
        if LA_qcov[i] + RA_qcov[i] > 1:
            LRA += 1
        else:
            break

    return FAq, FAs, LRA, LA, RA, MA

def pass_center_cond(FAq, FAs, LRA, LA, RA, MA, IQR_value):
    if MA > MA_TH:
        return False
    if LRA <= 0:
        return False
    if FAq + FAs + LRA < CENTER_MIN_SUPPORT:
        return False
    if min(LA, RA) < MIN_LA_RA:
        return False
    if IQR_TH1 < IQR_value < IQR_TH2 and LRA < LRA_MIN:
        return False
    if IQR_value < IQR_TH1:
        return False
    return True


def get_all_members(center, contains):

    members = set()
    stack = [center]
    while stack:
        r = stack.pop()
        if r not in members:
            members.add(r)
            stack.extend(contains.get(r, []))
    return members

# -------------------------
# Main
# -------------------------
def main():
    contains = defaultdict(set)
    incoming = defaultdict(set)
    read_len = {}
    aln_by_query = defaultdict(list)

    # ---- Read PAF ----
    with open(args.paf) as f:
        for line in f:
            fs = line.rstrip().split("\t")
            qid, qlen, qs, qe = fs[0], int(fs[1]), int(fs[2]), int(fs[3])
            tid, tlen, ts, te = fs[5], int(fs[6]), int(fs[7]), int(fs[8])
            matches, aln_len = int(fs[9]), int(fs[10])

            identity = matches / aln_len
            qcov = (qe - qs) / qlen
            scov = (te - ts) / tlen

            read_len[qid] = qlen
            read_len[tid] = tlen

            aln_by_query[qid].append(
                (qs, qe, qlen, qcov, scov, identity, "query", tid)
            )
            aln_by_query[tid].append(
                (ts, te, tlen, scov, qcov, identity, "target", qid)
            )

            if qcov > QCOV_TH and scov > SCOV_TH and identity > ID_TH:
                contains[qid].add(tid)
                incoming[tid].add(qid)

    if not read_len:
        sys.exit("[ERROR] No alignments loaded from PAF file. The file may be empty or malformed.")

    # ---- Build cluster ----
    all_reads = set(read_len)
    star_centers = [r for r in all_reads if r not in incoming]

    temp_bins = {c: get_all_members(c, contains) for c in star_centers}

    def initial_center(members):
        return max(members, key=lambda r: (read_len[r], r))

    bins_sorted = sorted(
        temp_bins.items(),
        key=lambda x: (read_len[initial_center(x[1])], initial_center(x[1])),
        reverse=True
    )
    bins = [members for _, members in bins_sorted]

    # ---- Select representatives ----
    selected_bins = {}
    for bin_id, members in enumerate(bins):
        for r in sorted(members, key=lambda r: (-read_len[r], r)):
            FAq, FAs, LRA, LA, RA, MA = endpoint_stat(r, aln_by_query)
            IQR_value = compute_iqr(r, aln_by_query, read_len[r])

            if read_len[r] < SHORT_CENTER_LEN and (FAq + FAs) < MIN_SHORT_CENTER_NUM:
                continue

            if pass_center_cond(FAq, FAs, LRA, LA, RA, MA, IQR_value):
                selected_bins[bin_id] = (r, members)
                break

    # ---- Output ----
    try:
        with open(REP_STAT_FILE, "w") as osf:
            osf.write("ClusterID\trepID\tLength\tFA\tLRA\tLA\tRA\tMA\tIQR\n")
            for i in sorted(selected_bins):
                center, members = selected_bins[i]
                FAq, FAs, LRA, LA, RA, MA = endpoint_stat(center, aln_by_query)
                IQR_value = compute_iqr(center, aln_by_query, read_len[center])
                FA = FAq + FAs
                osf.write(
                    f"{i}\t{center}\t{read_len[center]}\t{FA}\t{LRA}\t{LA}\t{RA}\t{MA}\t{IQR_value:.4f}\n"
                )
    except IOError as e:
        sys.exit(f"[ERROR] Cannot write output file: {e}")

    if not selected_bins:
        print("[WARNING] No representative reads selected. Check alignment quality or thresholds.")


if __name__ == "__main__":
    main()
