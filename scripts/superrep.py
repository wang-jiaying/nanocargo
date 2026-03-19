#!/usr/bin/env python3

from collections import defaultdict
import yaml
import argparse
import os
import sys

# =============================
# Argument parsing
# =============================
parser = argparse.ArgumentParser(description="SUPER representative selection")

parser.add_argument("--config", "-c", required=True)
parser.add_argument("--outdir", "-o", default=".")
parser.add_argument("--stdin", action="store_true")
parser.add_argument("--idth2", type=float, default=None,
                    help="Override ID_TH2 (graph identity threshold)")

args = parser.parse_args()

outdir = os.path.abspath(args.outdir)
os.makedirs(outdir, exist_ok=True)

# =============================
# Load config
# =============================
with open(args.config) as f:
    cfg = yaml.safe_load(f)["rep_superbin"]

REP_STAT_FILE = os.path.join(outdir, cfg["REP_STAT"])
SELFLIGATE_FILE = os.path.join(outdir, "selfligateid.txt")
OUT_SUPERREP = os.path.join(outdir, cfg["OUT_SUPERREP"])

# ---- Rep remove thresholds ----
ID_TH = cfg.get("ID_TH", 0.9)
SCOV_TH = cfg.get("SCOV_TH", 0.98)
ALEN_RATIO_LOW = cfg.get("ALEN_RATIO_LOW", 0.9)
ALEN_RATIO_HIGH = cfg.get("ALEN_RATIO_HIGH", 1.1)
QCOV_MAX = cfg.get("QCOV_MAX", 0.9)

SHORT_LEN = cfg.get("SHORT_LEN", 3000)
ID_TH_SHORT = cfg.get("ID_TH_SHORT", 0.9)
SCOV_TH_SHORT = cfg.get("SCOV_TH_SHORT", 0.9)
QCOV_MAX_SHORT = cfg.get("QCOV_MAX_SHORT", 0.7)

# ---- Graph thresholds ----
QCOV_TH = cfg.get("QCOV_TH", 0.98)
SCOV_TH2 = cfg.get("SCOV_TH2", 0.98)

ID_TH2 = args.idth2 if args.idth2 is not None else cfg.get("ID_TH2", 0.8)

LEN_RATIO_LOW = cfg.get("LEN_RATIO_LOW", 0.9)
LEN_RATIO_HIGH = cfg.get("LEN_RATIO_HIGH", 1.1)

# ---- Strong link thresholds ----
STRONG_ID_TH = cfg.get("STRONG_ID_TH", 0.9)
STRONG_SCOV_TH = cfg.get("STRONG_SCOV_TH", 0.99)
STRONG_QCOV_TH = cfg.get("STRONG_QCOV_TH", 0.95)

# =============================
# Helper functions
# =============================
def normalize_paf(fs):
    qid, qlen, qs, qe = fs[0], int(fs[1]), int(fs[2]), int(fs[3])
    strand = fs[4]
    tid, tlen, ts, te = fs[5], int(fs[6]), int(fs[7]), int(fs[8])
    matches, aln_len = int(fs[9]), int(fs[10])

    if qlen < tlen:
        qid, tid = tid, qid
        qlen, tlen = tlen, qlen
        qs, ts = ts, qs
        qe, te = te, qe

    return qid, qlen, qs, qe, tid, tlen, ts, te, matches, aln_len


def avg_identity(r):
    return identity_sum[r] / identity_cnt[r] if identity_cnt[r] else 0


# =============================
# Read PAF (stdin)
# =============================
if args.stdin:
    paf_lines = sys.stdin.readlines()
else:
    raise ValueError("PAF must come from stdin")

# =============================
# Step 1: detect reps to remove
# =============================
remove_reps = set()

for line in paf_lines:
    fs = line.rstrip().split("\t")
    qid, qlen, qs, qe, tid, tlen, ts, te, matches, aln_len = normalize_paf(fs)

    if qid == tid:
        continue

    identity = matches / aln_len
    qcov = (qe - qs) / qlen
    scov = (te - ts) / tlen
    alen_ratio = (qe - qs) / aln_len
    min_len = min(qlen, tlen)

    if (
        (
            identity > ID_TH
            and scov > SCOV_TH
            and ALEN_RATIO_LOW < alen_ratio < ALEN_RATIO_HIGH
            and qcov <= QCOV_MAX
        )
        or
        (
            min_len <= SHORT_LEN
            and identity > ID_TH_SHORT
            and scov > SCOV_TH_SHORT
            and ALEN_RATIO_LOW < alen_ratio < ALEN_RATIO_HIGH
            and qcov <= QCOV_MAX_SHORT
        )
    ):
        remove_reps.add(tid)

# =============================
# Step 2: read selfligate IDs
# =============================
selfligate_ids = set()

if os.path.exists(SELFLIGATE_FILE) and os.path.getsize(SELFLIGATE_FILE) > 0:
    with open(SELFLIGATE_FILE) as f:
        for line in f:
            selfligate_ids.add(line.strip())

# =============================
# Step 3: load REP_STAT
# =============================
rep_stat = {}
rep_len = {}
kept_reps = set()

if not os.path.exists(REP_STAT_FILE):
    sys.exit(f"[ERROR] Rep stat file not found: {REP_STAT_FILE}\nCheck if center_filter.py completed successfully.")

with open(REP_STAT_FILE) as fin:
    header = fin.readline()
    fields = header.rstrip().split("\t")
    repid_col = fields.index("repID")

    for line in fin:
        fs = line.rstrip().split("\t")
        r = fs[repid_col]

        if r in remove_reps or r in selfligate_ids:
            continue

        rec = dict(zip(fields, fs))

        kept_reps.add(r)
        rep_stat[r] = {
            "FA": int(rec["FA"]),
            "LRA": int(rec["LRA"]),
            "IQR": float(rec["IQR"]),
        }
        rep_len[r] = int(rec["Length"])

# =============================
# Step 4: build graph
# =============================
edges = defaultdict(set)
identity_sum = defaultdict(float)
identity_cnt = defaultdict(int)

for line in paf_lines:
    fs = line.rstrip().split("\t")
    q, qlen, qs, qe, t, tlen, ts, te, matches, aln_len = normalize_paf(fs)

    if q == t or q not in kept_reps or t not in kept_reps:
        continue

    qcov = (qe - qs) / qlen
    scov = (te - ts) / tlen
    identity = matches / aln_len
    lr = min(qlen, tlen) / max(qlen, tlen)

    if (
        (
            qcov >= QCOV_TH
            and scov >= SCOV_TH2
            and identity >= ID_TH2
        )
        or
        (
            identity > STRONG_ID_TH
            and scov > STRONG_SCOV_TH
            and qcov > STRONG_QCOV_TH
        )
    ) and (LEN_RATIO_LOW <= lr <= LEN_RATIO_HIGH):

        edges[q].add(t)
        edges[t].add(q)

        identity_sum[q] += identity
        identity_cnt[q] += 1
        identity_sum[t] += identity
        identity_cnt[t] += 1

# =============================
# Step 5: connected components
# =============================
visited = set()
components = []

for r in sorted(kept_reps, key=lambda x: (rep_len[x], x)):
    if r in visited:
        continue

    stack = [r]
    comp = set()

    while stack:
        x = stack.pop()
        if x in visited:
            continue
        visited.add(x)
        comp.add(x)
        stack.extend(edges.get(x, []))

    components.append(comp)

# =============================
# Step 6: select superrep
# =============================
with open(OUT_SUPERREP, "w") as out_super:
    out_super.write(
        "ComponentID\tComponentSize\tSuperRep\tLength\tFA\tLRA\tAvgIdentity\tIQR\tConnectNum\n"
    )

    for cid, comp in enumerate(components):
        best = None
        best_key = None

        for r in comp:
            FA = rep_stat[r]["FA"]
            LRA = rep_stat[r]["LRA"]
            IQR = rep_stat[r]["IQR"]
            connect = len(edges.get(r, []))
            avg_id = avg_identity(r)

            key = (FA, LRA, connect, avg_id, -rep_len[r], r)

            if best is None or key > best_key:
                best = r
                best_key = key

        out_super.write(
            f"{cid}\t{len(comp)}\t{best}\t{rep_len[best]}"
            f"\t{rep_stat[best]['FA']}"
            f"\t{rep_stat[best]['LRA']}\t{avg_identity(best):.4f}"
            f"\t{rep_stat[best]['IQR']}\t{len(edges.get(best, []))}\n"
        )

print("[INFO] SUPERrep finished successfully.")