#!/usr/bin/env python3
"""
Build PR+RT and IN consensus sequences and read-depth stats directly from
fastq2codfreq .codfreq files.

    python codfreq_consensus.py -s SAMPLE -c SAMPLE.codfreq

Writes SAMPLE_regions.fasta and SAMPLE_depth.tsv (override with -o / -d).

Consensus = most frequent codon at each codon position (codon-aware, in frame).
Positions with depth < --min-depth, or whose top codon is below --min-freq,
or that are absent from the file, become NNN.
"""
import argparse
import csv
from collections import defaultdict

GENE_LEN = {"PR": 99, "RT": 560, "IN": 288}
REGIONS = {"PRRT": ["PR", "RT"], "IN": ["IN"]}


def find_col(header, *names):
    low = [h.strip().lower() for h in header]
    for n in names:
        if n in low:
            return low.index(n)
    raise SystemExit(f"Column {names} not found in header: {header}")


def read_codfreq(path):
    """Return {(gene, pos): {codon: count}} and {(gene, pos): total}."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        first = fh.readline()
        fh.seek(0)
        rdr = csv.reader(fh, delimiter="\t" if "\t" in first else ",")
        header = next(rdr)
        gi = find_col(header, "gene")
        pi = find_col(header, "position", "pos")
        ti = find_col(header, "total", "totalreads", "total_reads")
        ci = find_col(header, "codon")
        ni = find_col(header, "count", "reads", "codonreads")
        counts = defaultdict(dict)
        totals = {}
        for row in rdr:
            if not row or row[gi] not in GENE_LEN:
                continue
            key = (row[gi], int(row[pi]))
            counts[key][row[ci]] = counts[key].get(row[ci], 0) + int(float(row[ni]))
            totals[key] = int(float(row[ti]))
    return counts, totals


def call_codon(codon_counts, total, min_depth, min_freq):
    if total < min_depth or not codon_counts:
        return "NNN"
    codon, n = max(codon_counts.items(), key=lambda kv: kv[1])
    if n / total < min_freq:
        return "NNN"
    if "-" in codon or codon.lower() in ("del", ""):
        return ""                      # majority deletion: omit codon
    return codon.upper()               # >3 nt = codon plus in-frame insertion


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-s", "--sample", required=True, help="sample name used in FASTA headers and depth table")
    ap.add_argument("-c", "--codfreq", required=True, help="path to the sample's .codfreq file")
    ap.add_argument("-o", "--fasta", help="output FASTA (default: SAMPLE_regions.fasta)")
    ap.add_argument("-d", "--depth", help="output depth table (default: SAMPLE_depth.tsv)")
    ap.add_argument("--min-depth", type=int, default=50)
    ap.add_argument("--min-freq", type=float, default=0.5,
                    help="top codon must reach this fraction, else NNN (default 0.5)")
    a = ap.parse_args()

    s = a.sample
    fasta = a.fasta or f"{s}_regions.fasta"
    depth = a.depth or f"{s}_depth.tsv"
    counts, totals = read_codfreq(a.codfreq)

    with open(fasta, "w") as fa, open(depth, "w") as dt:
        dt.write("sample\tregion\tmean_depth\tmin_depth\tpct_ge_min\tcodons_covered\n")
        for region, genes in REGIONS.items():
            seq, depths = [], []
            for g in genes:
                for pos in range(1, GENE_LEN[g] + 1):
                    tot = totals.get((g, pos), 0)
                    depths.append(tot)
                    seq.append(call_codon(counts.get((g, pos), {}), tot, a.min_depth, a.min_freq))
            fa.write(f">{s}_{region}\n{''.join(seq)}\n")
            n = len(depths)
            dt.write(f"{s}\t{region}\t{sum(depths)/n:.1f}\t{min(depths)}\t"
                     f"{100*sum(d >= a.min_depth for d in depths)/n:.1f}\t"
                     f"{sum(d > 0 for d in depths)}/{n}\n")
    print(f"Wrote {fasta} and {depth}")


if __name__ == "__main__":
    main()