#!/usr/bin/env python3
"""Pairwise identity over an aligned FASTA.

Used by 05_validate.sh to quantify the agreement between the de novo COI
(route A), the reference-guided consensus (route B) and the reference COX1.
Only columns at which both sequences carry an unambiguous base are compared,
so masked (N) and gap positions inflate neither identity nor divergence; they
are reported separately instead.

Emits one TSV row per pair:

    seq_a  seq_b  compared_sites  identical  differences  pct_identity  ungapped_overlap_frac

Usage:
    python lib/aln_identity.py results/05_validate/RHIN01_routes.aln.fasta
"""

from __future__ import annotations

import argparse
import itertools

from Bio import SeqIO

UNAMBIGUOUS = set("ACGT")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("alignment", help="aligned FASTA")
    parser.add_argument(
        "--no-header", action="store_true", help="omit the TSV header row"
    )
    args = parser.parse_args()

    records = list(SeqIO.parse(args.alignment, "fasta"))
    if len(records) < 2:
        raise SystemExit("ERROR: need at least two sequences to compare")

    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        raise SystemExit("ERROR: input is not aligned (unequal sequence lengths)")
    (aln_len,) = lengths

    if not args.no_header:
        print(
            "seq_a\tseq_b\tcompared_sites\tidentical\tdifferences\t"
            "pct_identity\tungapped_overlap_frac"
        )

    for a, b in itertools.combinations(records, 2):
        # MAFFT --adjustdirection prefixes reverse-complemented records with
        # "_R_"; strip it so the labels stay readable.
        name_a = a.id[3:] if a.id.startswith("_R_") else a.id
        name_b = b.id[3:] if b.id.startswith("_R_") else b.id

        compared = identical = 0
        for base_a, base_b in zip(str(a.seq).upper(), str(b.seq).upper()):
            if base_a in UNAMBIGUOUS and base_b in UNAMBIGUOUS:
                compared += 1
                identical += base_a == base_b

        pct = 100.0 * identical / compared if compared else 0.0
        overlap = compared / aln_len if aln_len else 0.0
        print(
            f"{name_a}\t{name_b}\t{compared}\t{identical}\t{compared - identical}\t"
            f"{pct:.2f}\t{overlap:.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
