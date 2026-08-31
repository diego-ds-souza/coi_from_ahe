#!/usr/bin/env python3
"""Reading-frame check for a recovered COI sequence.

A COI sequence pulled out of shotgun or capture data is only credible as a
mitochondrial barcode if it translates cleanly. Internal stop codons or a
length that is not a multiple of three under the best frame are the classic
signature of a nuclear mitochondrial pseudogene (NUMT), of an assembly
frameshift, or of the wrong genetic code.

Emits one TSV row per sequence:

    seq_id  length_bp  strand  frame  internal_stops  len_mod3  first_stop_aa  verdict

verdict is PASS when the best frame has no internal stops, FLAG otherwise.

Usage:
    python lib/check_orf.py results/coi_final/RHIN01_COI.fasta
"""

from __future__ import annotations

import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq

DEFAULT_TRANSL_TABLE = 5  # invertebrate mitochondrial


def evaluate(seq: Seq, table: int):
    """Return the best (fewest internal stops) of the six reading frames."""
    best = None
    for strand, nucleotides in (("+", seq), ("-", seq.reverse_complement())):
        for frame in (0, 1, 2):
            trimmed = nucleotides[frame:]
            trimmed = trimmed[: len(trimmed) - (len(trimmed) % 3)]
            if not trimmed:
                continue
            protein = str(trimmed.translate(table=table, to_stop=False))
            body = protein.rstrip("*")  # a single terminal stop is expected
            stops = body.count("*")
            first_stop = body.find("*") + 1 if stops else 0
            candidate = (stops, -len(protein), strand, frame + 1, first_stop)
            if best is None or candidate < best:
                best = candidate
    return best


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", nargs="+", help="one or more FASTA files")
    parser.add_argument("--transl-table", type=int, default=DEFAULT_TRANSL_TABLE)
    parser.add_argument(
        "--no-header", action="store_true", help="omit the TSV header row"
    )
    args = parser.parse_args()

    if not args.no_header:
        print(
            "seq_id\tlength_bp\tstrand\tframe\tinternal_stops\tlen_mod3\t"
            "first_stop_aa\tverdict"
        )

    flagged = 0
    for path in args.fasta:
        for record in SeqIO.parse(path, "fasta"):
            seq = Seq(str(record.seq).replace("-", "").upper())
            if not seq:
                continue
            best = evaluate(seq, args.transl_table)
            if best is None:
                continue
            stops, _, strand, frame, first_stop = best
            verdict = "PASS" if stops == 0 else "FLAG"
            flagged += stops > 0
            print(
                f"{record.id}\t{len(seq)}\t{strand}\t{frame}\t{stops}\t"
                f"{len(seq) % 3}\t{first_stop}\t{verdict}"
            )

    if flagged:
        print(
            f"WARNING: {flagged} sequence(s) carry internal stop codons under "
            "the best frame; treat them as candidate NUMTs or assembly errors.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
