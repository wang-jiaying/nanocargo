#!/usr/bin/env python3
import sys
import yaml
import os
import pandas as pd
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
DELTA = cfg.get("DELTA", 50)
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

# -------------------------
def endpoint_stat_and_iqr(read, aln_by_query, qlen):
    FAq = FAs = LA = RA = MA = 0
    LA_qcov, RA_qcov = [], []
    pos = []   

    for qs, qe, rqlen, qcov, scov, identity, role, mate in aln_by_query.get(read, []):

        if qcov > FA_COV_TH and scov > FA_COV_TH and identity > ID_TH:
            if role == "query":
                FAq += 1
            else:
                FAs += 1
            continue

        if role != "query":
            continue

        at = classify_endpoint(qs, qe, rqlen)

        if at == "LA":
            LA += 1
            if identity > ID_TH:
                LA_qcov.append(qcov)
                pos.append(qe)   
        elif at == "RA":
            RA += 1
            if identity > ID_TH:
                RA_qcov.append(qcov)
                pos.append(qs)   
        elif at == "MA":
            MA += 1

    # ---- LRA ----
    LA_qcov.sort(reverse=True)
    RA_qcov.sort(reverse=True)

    LRA = 0
    for i in range(min(len(LA_qcov), len(RA_qcov))):
        if LA_qcov[i] + RA_qcov[i] > 1:
            LRA += 1
        else:
            break

    # ---- IQR ----
    if len(pos) < IQR_MIN_SUPPORT:
        IQR_value = 1.0
    else:
        pos.sort()
        n = len(pos)
        q25 = pos[int(0.25 * (n - 1))]
        q75 = pos[int(0.75 * (n - 1))]
        IQR_value = (q75 - q25) / max(qlen - MIN_READLEN, 1)

    return FAq, FAs, LRA, LA, RA, MA, IQR_value

def pass_center_cond(FAq, FAs, LRA, LA, RA, MA, IQR_value):
    if MA > MA_TH:
        return False
    if LRA <= 0:
        return False
    if FAq + FAs + LRA < CENTER_MIN_SUPPORT:
        return False
    if min(LA, RA) < MIN_LA_RA:
        return False
    if IQR_TH1 <= IQR_value < IQR_TH2 and LRA < LRA_MIN:
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
    aln_by_query = defaultdict(list)

    # -------------------------
    try:
        df = pd.read_csv(
            args.paf,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 3, 5, 6, 7, 8, 9, 10],
            names=["qid", "qlen", "qs", "qe", "tid", "tlen", "ts", "te", "matches", "alen"],
            dtype={
                "qid": str, "qlen": int, "qs": int, "qe": int,
                "tid": str, "tlen": int, "ts": int, "te": int,
                "matches": int, "alen": int
            }
        )
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read PAF file: {e}")

    if df.empty:
        sys.exit("[ERROR] No alignments loaded from PAF file. The file may be empty or malformed.")

    df["identity"] = df["matches"] / df["alen"]
    df["qcov"]     = (df["qe"] - df["qs"]) / df["qlen"]
    df["scov"]     = (df["te"] - df["ts"]) / df["tlen"]

    read_len = {}
    for _, row in df[["qid", "qlen"]].drop_duplicates("qid").iterrows():
        read_len[row["qid"]] = row["qlen"]
    for _, row in df[["tid", "tlen"]].drop_duplicates("tid").iterrows():
        read_len[row["tid"]] = row["tlen"]

    for row in df.itertuples(index=False):
        aln_by_query[row.qid].append(
            (row.qs, row.qe, row.qlen, row.qcov, row.scov, row.identity, "query", row.tid)
        )
        aln_by_query[row.tid].append(
            (row.ts, row.te, row.tlen, row.scov, row.qcov, row.identity, "target", row.qid)
        )

        if row.qcov > QCOV_TH and row.scov > SCOV_TH and row.identity > ID_TH:
            contains[row.qid].add(row.tid)
            incoming[row.tid].add(row.qid)

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

    # ---- Select center reads ----
    selected_bins = {}
    for bin_id, members in enumerate(bins):
        for r in sorted(members, key=lambda r: (-read_len[r], r)):
            FAq, FAs, LRA, LA, RA, MA, IQR_value = endpoint_stat_and_iqr(
                r, aln_by_query, read_len[r]
            )

            if read_len[r] < SHORT_CENTER_LEN and (FAq + FAs) < MIN_SHORT_CENTER_NUM:
                continue
    
            if pass_center_cond(FAq, FAs, LRA, LA, RA, MA, IQR_value):
                selected_bins[bin_id] = (r, members)
                break


    # ---- Build FA candidate pool ----
    center_fa_pool = {}
    for bin_id, (center, members) in selected_bins.items():
        fa_reads = set()
        for _, _, _, qcov, scov, identity, role, mate in aln_by_query[center]:
            if qcov > FA_COV_TH and scov > FA_COV_TH and identity > ID_TH:
                fa_reads.add(mate)
        center_fa_pool[bin_id] = {center} | fa_reads


    # ---- Select representative ----
    rep_seen = {}
    rep_stat = {}

    for bin_id in sorted(selected_bins):
        pool = center_fa_pool[bin_id]
    
        def rep_key(r):
            FAq, FAs, LRA, LA, RA, MA, _ = endpoint_stat_and_iqr(
                r, aln_by_query, read_len[r]
            )
            return (-(FAq + FAs), -LRA, r)

        rep = min(pool, key=rep_key)

        if rep not in rep_seen:
            rep_seen[rep] = bin_id

            FAq, FAs, LRA, LA, RA, MA, IQR_value = endpoint_stat_and_iqr(
                rep, aln_by_query, read_len[rep]
            )

            rep_stat[bin_id] = (
                rep,
                FAq + FAs,
                LRA,
                LA,
                RA,
                MA,
                IQR_value
            )


    # ---- Output ----
    try:
        with open(REP_STAT_FILE, "w") as osf:
            osf.write("repID\tLength\tFA\tLRA\tLA\tRA\tMA\tIQR\n")
            for bin_id in sorted(rep_stat):
                rep, FA, LRA, LA, RA, MA, IQR_value = rep_stat[bin_id]
                osf.write(
                    f"{rep}\t{read_len[rep]}\t{FA}\t{LRA}\t{LA}\t{RA}\t{MA}\t{IQR_value:.4f}\n"
                )
    except IOError as e:
        sys.exit(f"[ERROR] Cannot write output file: {e}")

    if not rep_stat:
        print("[WARNING] No representative reads selected. Check alignment quality or thresholds.")

if __name__ == "__main__":
    main()