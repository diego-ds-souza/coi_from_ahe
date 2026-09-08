#!/usr/bin/env bash
#
# Reference preparation. Downloads a mitochondrial genome from NCBI (GenBank
# flatfile plus FASTA), establishes the COX1 CDS in nucleotide and amino-acid
# form together with its interval on the contig, and indexes the mitogenome for
# mapping. Everything downstream is anchored on these files.
#
# Two ways the COX1 interval is established:
#
#   1. From the record's own annotation, when it has one. Preferred.
#   2. By tblastn, when MITO_ACC is an unannotated record. Set ANNOT_ACC to any
#      annotated mitogenome, and its COX1 protein is used to locate COX1 in
#      MITO_ACC. This is what makes an UNVERIFIED GenBank record usable: those
#      often carry no features, but the sequence itself is fine, and a close
#      unannotated reference beats a distant annotated one for every step that
#      maps reads.
#
# Choosing MITO_ACC: the closest mitogenome available, annotated or not. To list
# candidates, filtering by length rather than by title, since many mitogenomes
# were deposited under titles that do not say "complete genome":
#
#   esearch -db nuccore -query "Cossoninae[Organism] AND mitochondrion[filter] \
#     AND 10000:20000[SLEN]" \
#   | efetch -format docsum \
#   | xtract -pattern DocumentSummary -element Caption,Organism,Slen,Title
#
# To check whether a candidate is annotated before committing to it:
#
#   efetch -db nuccore -id MH404140 -format gb | grep -c '^     CDS '
#
# Run from the project root:
#   conda activate coi_from_ahe
#   MITO_ACC=NC_038191 bash 00_reference.sh   # NC_038191 Cicindela anchoralis, https://www.ncbi.nlm.nih.gov/nuccore/NC_038191.1
#   MITO_ACC=MH404140 ANNOT_ACC=NC_059702 bash 00_reference.sh   # unannotated
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
MITO_ACC="${MITO_ACC:-}"                       # NCBI nuccore accession, required
ANNOT_ACC="${ANNOT_ACC:-}"                     # annotated donor, if MITO_ACC is not
REFDIR="${REFDIR:-ref}"                        # output directory
TRANSL_TABLE="${TRANSL_TABLE:-5}"              # invertebrate mitochondrial code
THREADS="${THREADS:-24}"                        # CPU threads, for tblastn
# -----------------------------------------------------------------------------

[[ -n "$MITO_ACC" ]] || die "set MITO_ACC to a reference mitogenome accession (see header)"
require_tools efetch bwa samtools seqkit
require_python

mkdir -p "$REFDIR"
GB="${REFDIR}/${MITO_ACC}.gb"
FA="${REFDIR}/${MITO_ACC}.fasta"

# 1) Download, with retries and a completeness check. Existing complete files
#    are reused, so the script is safe to re-run, does not hammer NCBI, and
#    works with no outbound network if the files are copied in by hand.
echo "[1/4] efetch ${MITO_ACC}"
fetch_ncbi "$MITO_ACC" gb    "$GB"
fetch_ncbi "$MITO_ACC" fasta "$FA"

# 2) Establish COX1, from annotation if possible.
echo "[2/4] locate COX1"
# Exit code 3 means specifically "no COX1 annotation in this record", which is
# the only condition that justifies the blast fallback. Any other nonzero exit
# is a real failure and must not be silently reinterpreted as an unannotated
# record: that turns a broken environment into a misleading message about the
# accession.
rc=0
"$PY" "${SCRIPT_DIR}/lib/extract_cox1_from_gb.py" \
  --gb "$GB" --prefix "${REFDIR}/cox1" --transl-table "$TRANSL_TABLE" || rc=$?

if [[ "$rc" -eq 0 ]]; then
  msg "COX1 taken from the record's own annotation"
elif [[ "$rc" -ne 3 ]]; then
  die "extract_cox1_from_gb.py failed with exit ${rc}. This is a tooling
       problem, not a problem with ${MITO_ACC}. The error above is the real one."
else
  msg "no COX1 annotation in ${MITO_ACC}; falling back to a tblastn search"
  [[ -n "$ANNOT_ACC" ]] || die "${MITO_ACC} carries no COX1 annotation.
       Set ANNOT_ACC to an annotated mitogenome (any weevil will do; its COX1
       protein is used only to locate the gene, not to define its sequence),
       for example:
         MITO_ACC=${MITO_ACC} ANNOT_ACC=NC_059702 bash 00_reference.sh"
  require_tools tblastn makeblastdb

  AGB="${REFDIR}/${ANNOT_ACC}.gb"
  fetch_ncbi "$ANNOT_ACC" gb "$AGB"

  drc=0
  "$PY" "${SCRIPT_DIR}/lib/extract_cox1_from_gb.py" \
    --gb "$AGB" --prefix "${REFDIR}/donor_cox1" --transl-table "$TRANSL_TABLE" || drc=$?
  if [[ "$drc" -eq 3 ]]; then
    die "${ANNOT_ACC} carries no COX1 annotation either; pick another donor.
       NC_059702 (Euwallacea fornicatus) is annotated and works."
  elif [[ "$drc" -ne 0 ]]; then
    die "extract_cox1_from_gb.py failed with exit ${drc} on the donor record.
       This is a tooling problem; the error above is the real one."
  fi

  mkdir -p "${REFDIR}/blastdb"
  makeblastdb -in "$FA" -dbtype nucl -out "${REFDIR}/blastdb/mito" \
    > "${REFDIR}/makeblastdb.log" 2>&1
  tblastn -query "${REFDIR}/donor_cox1.faa" -db "${REFDIR}/blastdb/mito" \
    -db_gencode "$TRANSL_TABLE" -evalue 1e-5 -num_threads "$THREADS" \
    -outfmt '6 sseqid sstart send pident length bitscore' \
    | sort -k6,6gr > "${REFDIR}/cox1_hits.tsv" || true
  [[ -s "${REFDIR}/cox1_hits.tsv" ]] \
    || die "no COX1 found in ${MITO_ACC} by tblastn against ${ANNOT_ACC}.
       Either MITO_ACC is not a mitogenome, or it does not span COX1."

  blast_span_bed "${REFDIR}/cox1_hits.tsv" COX1 > "${REFDIR}/cox1.bed"
  seqkit subseq --bed "${REFDIR}/cox1.bed" "$FA" \
    | seqkit replace -p '^.*$' -r "${MITO_ACC}_COX1" > "${REFDIR}/cox1.fna"
  "$PY" "${SCRIPT_DIR}/lib/translate_best_frame.py" "${REFDIR}/cox1.fna" \
    --name "${MITO_ACC}_COX1" --out "${REFDIR}/cox1.faa" \
    --transl-table "$TRANSL_TABLE"
fi

# 3) Index the mitogenome for bwa and samtools.
echo "[3/4] index"
bwa index "$FA" 2>/dev/null
samtools faidx "$FA"

# 4) Record the COX1 region string once, so later scripts need not re-derive it.
echo "[4/4] region"
region_from_bed "${REFDIR}/cox1.bed" > "${REFDIR}/cox1.region"

echo
echo "wrote ${REFDIR}/cox1.{fna,faa,bed,region} and indexed ${FA}"
echo "COX1 region: $(cat "${REFDIR}/cox1.region")  ($(seqkit fx2tab -nl "${REFDIR}/cox1.fna" | cut -f2) bp)"
echo
echo "A COX1 far from about 1530 bp means the interval is wrong. Check"
echo "${REFDIR}/cox1_hits.tsv if the tblastn fallback was used."
