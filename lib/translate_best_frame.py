#!/usr/bin/env python3
"""Translate a COI nucleotide sequence in its best reading frame.

Used wherever a COX1 protein is needed but no annotation supplies one: a COI cut
out of an assembly, or a COX1 interval located by blast in an unannotated
mitogenome. Neither starts reliably at the initiation codon, so the frame is
chosen as the one with the fewest internal stop codons rather than assumed.

Emits FASTA on stdout, or to --out, and a one-line report on stderr.

Usage:
    python lib/translate_best_frame.py ref/cox1.fna --name NC_059702_COX1 \\
        --out ref/cox1.faa
"""

from __future__ import annotations

import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq

DEFAULT_TRANSL_TABLE = 5  # invertebrate mitochondrial


def best_frame(seq: Seq, table: int):
    """Return (internal_stops, frame, protein) for the cleanest forward frame."""
    best = None
    for frame in (0, 1, 2):
        trimmed = seq[frame:]
        trimmed = trimmed[: len(trimmed) - (len(trimmed) % 3)]
        if not trimmed:
            continue
        protein = str(trimmed.translate(table=table, to_stop=False))
        stops = protein.rstrip("*").count("*")
        candidate = (stops, -len(protein), frame + 1, protein)
        if best is None or candidate < best:
            best = candidate
    if best is None:
        raise SystemExit("ERROR: sequence too short to translate")
    stops, _, frame, protein = best
    return stops, frame, protein


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="single-record nucleotide FASTA")
    parser.add_argument("--out", default="-", help="output FASTA (default: stdout)")
    parser.add_argument("--name", default=None, help="record name (default: input's)")
    parser.add_argument("--transl-table", type=int, default=DEFAULT_TRANSL_TABLE)
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if len(records) != 1:
        raise SystemExit(f"ERROR: expected one record in {args.fasta}, found {len(records)}")
    record = records[0]

    seq = Seq(str(record.seq).replace("-", "").upper())
    stops, frame, protein = best_frame(seq, args.transl_table)
    name = args.name or record.id

    # A trailing stop is expected and is dropped; internal ones are the signal.
    body = protein.rstrip("*")
    text = f">{name} frame={frame} table={args.transl_table}\n{body}\n"
    if args.out == "-":
        sys.stdout.write(text)
    else:
        with open(args.out, "w") as handle:
            handle.write(text)

    print(
        f"frame {frame}, {len(body)} aa, internal stops {stops}",
        file=sys.stderr,
    )
    if stops:
        print(
            "WARNING: internal stop codons in every frame. This sequence is a "
            "poor reference and may be a NUMT or a misassembly.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
