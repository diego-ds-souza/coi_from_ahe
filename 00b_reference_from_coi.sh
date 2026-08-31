#!/usr/bin/env bash
#
# Build a reference directory from a COI sequence you already recovered, so that
# 02_bait.sh and 04_consensus.sh can be rerun against a conspecific target
# instead of a distant relative.
#
# Why this exists. With no Cossoninae mitogenome available, the reference in
# 00_reference.sh is subfamilies away from these samples. Route A tolerates that
# because it assembles the whole library and only uses the reference for a
# translated search. Route B does not: mapping recruits poorly at that distance,
# so the consensus comes back heavily masked and the route A versus route B
# identity in 05 becomes uninformative.
#
# Rerunning 02 and 04 against the route A COI fixes the recruitment problem and
# gives clean per-site depth. Read the result correctly, though:
#
#   Pass 1 (00_reference.sh, a related mitogenome): route B is INDEPENDENT of
#     route A. Their agreement is evidence the sequence is real.
#   Pass 2 (this script, route A's own COI): route B is DERIVED from route A.
#     Their agreement measures read support for route A, not corroboration of
#     it. It cannot detect a NUMT that assembled cleanly.
#
# Both are worth having, and they answer different questions. The independent
# second opinion on this platform is 06_mitogenome.sh, which reconstructs COX1
# through MitoFinder's own path.
#
# Run from the project root:
#   conda activate coi_from_ahe
#   COI_FASTA=results/coi_final/USNMENT01160338_COI.fasta \
#   REFDIR=ref_from_338 bash 00b_reference_from_coi.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
COI_FASTA="${COI_FASTA:-}"                     # a single-record COI FASTA, required
REFDIR="${REFDIR:-}"                           # output directory, required
REF_NAME="${REF_NAME:-COI_reference}"          # sequence ID written into the reference
TRANSL_TABLE="${TRANSL_TABLE:-5}"              # invertebrate mitochondrial code
# -----------------------------------------------------------------------------

# Argument checks first: a missing setting is the user's typo, a missing tool is
# the wrong environment, and the first message should name the right problem.
[[ -n "$COI_FASTA" ]] || die "set COI_FASTA to the COI sequence to use as reference"
[[ -n "$REFDIR" ]]    || die "set REFDIR to the output directory, e.g. ref_from_338"
require_files "$COI_FASTA"
require_tools bwa samtools seqkit
require_python

n_records=$(grep -c '^>' "$COI_FASTA")
[[ "$n_records" -eq 1 ]] || die "${COI_FASTA} holds ${n_records} records; expected exactly 1"

mkdir -p "$REFDIR"
FA="${REFDIR}/${REF_NAME}.fasta"

# 1) Normalize the header. Everything downstream keys on the sequence ID, and a
#    description containing spaces splits the ID differently in different tools.
echo "[1/4] normalize"
seqkit replace -p '^.*$' -r "$REF_NAME" "$COI_FASTA" | seqkit seq -u > "$FA"

# 2) The whole sequence is the COX1 interval, on the plus strand, because that
#    is how 03_assemble.sh wrote it out.
echo "[2/4] interval"
len=$(seqkit fx2tab -nl "$FA" | cut -f2)
printf '%s\t0\t%s\tCOX1\t0\t+\n' "$REF_NAME" "$len" > "${REFDIR}/cox1.bed"
region_from_bed "${REFDIR}/cox1.bed" > "${REFDIR}/cox1.region"
cp "$FA" "${REFDIR}/cox1.fna"

# 3) Protein form, for the tblastn search in 03. The frame is chosen as the one
#    with no internal stop codons rather than assumed to be frame 1, because a
#    contig-derived COI rarely starts exactly at the initiation codon.
echo "[3/4] translate"
"$PY" "${SCRIPT_DIR}/lib/translate_best_frame.py" "${REFDIR}/cox1.fna" \
  --name "${REF_NAME}_COX1" --out "${REFDIR}/cox1.faa" --transl-table "$TRANSL_TABLE"

# 4) Index.
echo "[4/4] index"
bwa index "$FA" 2>/dev/null
samtools faidx "$FA"

echo "wrote ${REFDIR}/ (${len} bp reference)"
echo
echo "Rerun the reference-dependent steps against it, into separate output dirs:"
echo "  SAMPLES=<sheet> REFDIR=${REFDIR} OUTDIR=results/02_bait_pass2 bash 02_bait.sh"
echo "  SAMPLES=<sheet> REFDIR=${REFDIR} BAITDIR=results/02_bait_pass2 \\"
echo "    OUTDIR=results/04_consensus_pass2 bash 04_consensus.sh"
