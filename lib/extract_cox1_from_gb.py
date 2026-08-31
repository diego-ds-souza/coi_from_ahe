#!/usr/bin/env python3
"""Extract the COX1 CDS from a mitochondrial GenBank flatfile.

Writes three files that the rest of the pipeline is anchored on:

  <prefix>.fna   spliced COX1 coding sequence, sense strand, nucleotides
  <prefix>.faa   its conceptual translation (genetic code from the record,
                 invertebrate mitochondrial / table 5 by default)
  <prefix>.bed   the COX1 interval on the reference contig, used by samtools
                 and bcftools to restrict depth and consensus calls

The CDS is identified from the /gene and /product qualifiers. Subunit I is
matched explicitly and subunits II and III are excluded, because a substring
search for "cytochrome c oxidase subunit I" also matches "... subunit III".

Usage:
    python lib/extract_cox1_from_gb.py --gb ref/NC_XXXXXX.gb --prefix ref/cox1
"""

from __future__ import annotations

import argparse
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq

# Accepted /gene qualifier spellings, normalized to lowercase.
GENE_ALIASES = {"cox1", "coi", "co1", "coxi", "mt-co1", "mtco1", "cox-1"}

# Default genetic code when the record carries no /transl_table qualifier.
DEFAULT_TRANSL_TABLE = 5  # invertebrate mitochondrial

# Distinct exit code for "this record has no COX1 annotation", so that a caller
# can fall back to a blast search for that case alone. Any other failure, a
# missing file or a broken interpreter, must abort loudly instead of being
# mistaken for an unannotated record.
NO_COX1_EXIT = 3


def _normalize(text: str) -> str:
    return re.sub(r"[\s_\-]+", " ", text.strip().lower())


def is_cox1(feature) -> bool:
    """True if a CDS feature is cytochrome c oxidase subunit I."""
    genes = [_normalize(g) for g in feature.qualifiers.get("gene", [])]
    if any(g in GENE_ALIASES for g in genes):
        return True

    for product in feature.qualifiers.get("product", []):
        norm = _normalize(product)
        if "oxidase" not in norm:
            continue
        # The subunit number is the final token; "i" and "1" are subunit I,
        # while "ii", "iii", "2" and "3" are not.
        if norm.split()[-1] in {"i", "1"}:
            return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gb", required=True, help="GenBank flatfile")
    parser.add_argument("--prefix", required=True, help="output path prefix")
    parser.add_argument(
        "--transl-table",
        type=int,
        default=None,
        help="override the genetic code (default: /transl_table, else 5)",
    )
    args = parser.parse_args()

    records = list(SeqIO.parse(args.gb, "genbank"))
    if len(records) != 1:
        sys.exit(f"ERROR: expected one record in {args.gb}, found {len(records)}")
    record = records[0]

    hits = [f for f in record.features if f.type == "CDS" and is_cox1(f)]
    if not hits:
        print(
            f"no COX1 CDS annotated in {args.gb} "
            f"({len(record.features)} features of any kind)",
            file=sys.stderr,
        )
        return NO_COX1_EXIT
    if len(hits) > 1:
        sys.exit(f"ERROR: {len(hits)} candidate COX1 CDS features in {args.gb}")
    feature = hits[0]

    table = args.transl_table or int(
        feature.qualifiers.get("transl_table", [DEFAULT_TRANSL_TABLE])[0]
    )
    nt: Seq = feature.extract(record.seq)
    # cds=False keeps the behaviour predictable for records whose COX1 uses an
    # incomplete stop codon, which is common in insect mitogenomes.
    aa: Seq = nt.translate(table=table, to_stop=False)

    start = int(feature.location.start)  # 0-based, BED convention
    end = int(feature.location.end)
    strand = "+" if feature.location.strand != -1 else "-"

    with open(f"{args.prefix}.fna", "w") as handle:
        handle.write(f">{record.id}_COX1 {len(nt)} bp table={table}\n{nt}\n")
    with open(f"{args.prefix}.faa", "w") as handle:
        handle.write(f">{record.id}_COX1 {len(aa)} aa table={table}\n{aa}\n")
    with open(f"{args.prefix}.bed", "w") as handle:
        handle.write(f"{record.id}\t{start}\t{end}\tCOX1\t0\t{strand}\n")

    internal_stops = str(aa).rstrip("*").count("*")
    print(
        f"COX1 {record.id}:{start + 1}-{end} strand {strand} "
        f"{len(nt)} bp / {len(aa)} aa, table {table}, "
        f"internal stops {internal_stops}",
        file=sys.stderr,
    )
    if internal_stops:
        print(
            "WARNING: the reference translation contains internal stop codons; "
            "the annotation or the genetic code may be wrong.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
