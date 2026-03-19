#!/usr/bin/env python3

import sys
import os
import argparse
import yaml


def parse_paf(line):
    f = line.strip().split("\t")

    return {
        "queryid": f[0],
        "qlen": int(f[1]),
        "qs": int(f[2]),
        "qe": int(f[3]),
        "strand": f[4],
        "subjectid": f[5],
        "slen": int(f[6]),
        "ss": int(f[7]),
        "se": int(f[8]),
        "match": int(f[9]),
        "alen": int(f[10]),
        "mapq": int(f[11])
    }


def process_group(records, params, result_ids):

    if len(records) < 2:
        return

    min_align_len = params.get("min_align_len", 400)
    min_cov = params.get("min_cov", 0.05)
    max_cov = params.get("max_cov", 0.9)
    ratio_lower = params.get("ratio_lower", 0.93)
    ratio_upper = params.get("ratio_upper", 1.07)
    cov_sum = params.get("cov_sum", 0.95)
    query_chain_lower = params.get("query_chain_lower", 0.9)
    query_chain_upper = params.get("query_chain_upper", 1.11)
    subject_chain_lower = params.get("subject_chain_lower", 0.87)
    subject_chain_upper = params.get("subject_chain_upper", 1.15)
    for r in records:

        r["qlen_align"] = r["qe"] - r["qs"]
        r["slen_align"] = r["se"] - r["ss"]

        if r["qlen_align"] < min_align_len:
            r["valid"] = False
            continue

        if r["slen_align"] < min_align_len:
            r["valid"] = False
            continue

        if not (ratio_lower < r["qlen_align"] / r["alen"] < ratio_upper):
            r["valid"] = False
            continue

        if not (ratio_lower < r["slen_align"] / r["alen"] < ratio_upper):
            r["valid"] = False
            continue

        r["qcov"] = r["qlen_align"] / r["qlen"]
        r["scov"] = r["slen_align"] / r["slen"]

        r["valid"] = True

    records = [r for r in records if r["valid"]]

    if len(records) < 2:
        return

    qid = records[0]["queryid"]
    sid = records[0]["subjectid"]

    # query chain

    g1 = [r for r in records if min_cov < r["qcov"] <= max_cov]
    g1.sort(key=lambda x: x["qs"])

    for i in range(len(g1)-1):

        r1 = g1[i]
        r2 = g1[i+1]

        base = (
            query_chain_lower <= r1["qe"]/max(r2["qs"],1) <= query_chain_upper and
            r1["strand"] != r2["strand"] and
            (r1["qcov"] + r2["qcov"]) > cov_sum
        )

        if r1["strand"] == "-":

            condA = subject_chain_lower <= max(r1["ss"],1)/max(r2["ss"],1) <= subject_chain_upper
            condB = query_chain_lower <= r1["se"]/max(r2["se"],1) <= query_chain_upper

        else:

            condA = subject_chain_lower <= r1["se"]/max(r2["se"],1) <= subject_chain_upper
            condB = query_chain_lower <= max(r1["ss"],1)/max(r2["ss"],1) <= query_chain_upper

        cond = base and (condA or condB)

        if cond:
            result_ids.add(qid)
            return


    # subject chain

    g2 = [r for r in records if min_cov < r["scov"] <= max_cov]
    g2.sort(key=lambda x: x["ss"])

    for i in range(len(g2)-1):

        r1 = g2[i]
        r2 = g2[i+1]

        base = (
            query_chain_lower <= r1["se"]/max(r2["ss"],1) <= query_chain_upper and
            r1["strand"] != r2["strand"] and
            (r1["scov"] + r2["scov"]) > cov_sum
        )

        if r1["strand"] == "-":

            condA = subject_chain_lower <= max(r1["qs"],1)/max(r2["qs"],1) <= subject_chain_upper
            condB = query_chain_lower <= r1["qe"]/max(r2["qe"],1) <= query_chain_upper

        else:

            condA = subject_chain_lower <= r1["qe"]/max(r2["qe"],1) <= subject_chain_upper
            condB = query_chain_lower <= max(r1["qs"],1)/max(r2["qs"],1) <= query_chain_upper

        cond = base and (condA or condB)

        if cond:
            result_ids.add(sid)
            return


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--stdin", action="store_true")
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    params = cfg["self_ligation"]

    result_ids = set()

    current_pair = None
    buffer = []

    for line in sys.stdin:

        r = parse_paf(line)

        if r["queryid"] == r["subjectid"]:
            continue

        pair = (r["queryid"], r["subjectid"])

        if current_pair is None:
            current_pair = pair

        if pair != current_pair:

            process_group(buffer, params, result_ids)

            buffer = []
            current_pair = pair

        buffer.append(r)

    if buffer:
        process_group(buffer, params, result_ids)

    outfile = os.path.join(args.outdir, "selfligateid.txt")

    with open(outfile, "w") as f:
        for i in sorted(result_ids):
            f.write(i+"\n")

    print("Finished")
    print("Output:", outfile)


if __name__ == "__main__":
    main()