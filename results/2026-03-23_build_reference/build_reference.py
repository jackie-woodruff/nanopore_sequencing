#!/usr/bin/env python3
"""Build adaptive_sampling_reference.fasta.

New chr14 layout (coordinates documented in 1-based inclusive / pyfaidx 0-based half-open):

  Piece           | 1-based (samtools)          | 0-based (pyfaidx)          | Length
  ----------------|-----------------------------|-----------------------------|----------
  franken chr14   | chr14:1-105480000           | [0:105480000]               | 105,480,000
  franken ighc    | ighc:51364-401741           | [51363:401741]              |     350,378
  franken igh     | igh:3-1193129               | [2:1193129]                 |   1,193,127
  hg38 chr14 tail | chr14:106879845-107043718   | [106879844:107043718]       |     163,874
  ----------------|-----------------------------|-----------------------------|----------
  new chr14 total |                             |                             | 107,187,379

Note on ighc exclusion:
  ighc:1-51363 is a non-IGHC region that overlaps franken chr14:105428637-105480000.
  The chr14 version is retained; the ighc duplicate is removed.

All other franken contigs (except chr14, ighc, igh) are carried over unchanged.
chr14_hg19 = hg19 chr14 in its entirety, renamed.
"""

import argparse
import sys

from pyfaidx import Fasta

LINE_WIDTH = 60

EXPECTED = {
    "franken_chr14": 105_480_000,
    "ighc_portion":      350_378,   # ighc:51364-401741 (1-based) → [51363:401741]
    "igh_portion":     1_193_127,   # igh:3-1193129 (1-based)    → [2:1193129]
    "hg38_post_igh":     163_874,   # chr14 [106879844:107043718] (0-based half-open)
    "new_chr14":     107_187_379,
}


def check(label, got, exp):
    status = "OK" if got == exp else f"MISMATCH (expected {exp:,})"
    print(f"  [{status}] {label}: {got:,} bp", file=sys.stderr, flush=True)
    if got != exp:
        raise AssertionError(f"Length check failed for {label}: got {got:,}, expected {exp:,}")


def write_seq(fh, name, seq, lw=LINE_WIDTH):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), lw):
        fh.write(seq[i : i + lw] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--franken", required=True, help="Franken reference FASTA")
    parser.add_argument("--hg38",    required=True, help="hg38 reference FASTA")
    parser.add_argument("--hg19",    required=True, help="hg19 reference FASTA")
    parser.add_argument("--output",  required=True, help="Output FASTA path")
    args = parser.parse_args()

    print("Opening FASTA indexes ...", file=sys.stderr, flush=True)
    franken = Fasta(args.franken, rebuild=False)
    hg38    = Fasta(args.hg38,    rebuild=False)
    hg19    = Fasta(args.hg19,    rebuild=False)

    # ---- Build new chr14 --------------------------------------------------------
    print("\nExtracting new chr14 pieces ...", file=sys.stderr, flush=True)

    # Full franken chr14
    p1 = str(franken["chr14"])
    check("franken_chr14", len(p1), EXPECTED["franken_chr14"])

    # ighc:51364-401741 (1-based) → pyfaidx [51363:401741]
    # Excludes ighc:1-51363 which overlaps with chr14:105428637-105480000
    p2 = str(franken["ighc"][51363:401741])
    check("ighc_portion", len(p2), EXPECTED["ighc_portion"])

    # igh:3-1193129 (1-based) → pyfaidx [2:1193129]
    p3 = str(franken["igh"][2:1193129])
    check("igh_portion", len(p3), EXPECTED["igh_portion"])

    # hg38 chr14 remainder (0-based half-open): [106879844:107043718]
    p4 = str(hg38["chr14"][106_879_844:107_043_718])
    check("hg38_post_igh", len(p4), EXPECTED["hg38_post_igh"])

    new_chr14 = p1 + p2 + p3 + p4
    check("new_chr14", len(new_chr14), EXPECTED["new_chr14"])

    # ---- Write output -----------------------------------------------------------
    print(f"\nWriting {args.output} ...", file=sys.stderr, flush=True)

    skip = {"chr14", "ighc", "igh"}

    with open(args.output, "w") as fh:
        # 1. New assembled chr14
        print("  Writing new chr14 ...", file=sys.stderr, flush=True)
        write_seq(fh, "chr14", new_chr14)

        # 2. All remaining franken contigs (preserving order from .fai)
        for name in franken.keys():
            if name in skip:
                print(f"  Skipping franken {name} (folded into new chr14)", file=sys.stderr, flush=True)
                continue
            print(f"  Appending franken {name} ...", file=sys.stderr, flush=True)
            write_seq(fh, name, str(franken[name]))

        # 3. hg19 chr14 renamed to chr14_hg19
        print("  Appending chr14_hg19 from hg19 ...", file=sys.stderr, flush=True)
        write_seq(fh, "chr14_hg19", str(hg19["chr14"]))

    print("\nDone.", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
