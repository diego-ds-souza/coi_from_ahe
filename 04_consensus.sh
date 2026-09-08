#!/usr/bin/env bash
#
# Route B: reference-guided consensus of COX1 from the alignment produced by
# 02_bait.sh. This is not a second opinion on the assembly graph so much as an
# independent reconstruction path: it uses the same reads but no assembler, and
# it fails in different ways from route A. Agreement between the two is the
# evidence that the sequence is real.
#
# Indels are excluded from calling (-V indels) so that consensus coordinates
# stay aligned to the reference and the COX1 interval can be cut out directly.
# Real indels in a protein-coding mitochondrial gene are rare between
# congeners, and any that exist are recovered by route A instead. Sites below
# MIN_DEPTH are masked to N rather than silently inherited from the reference,
# which is the most common way a reference-guided COI ends up subtly wrong.
#
# Run from the project root:
#   conda activate coi_from_ahe
#   bash 04_consensus.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
SAMPLES="${SAMPLES:-samples.tsv}"              # sample sheet
BAITDIR="${BAITDIR:-results/02_bait}"          # output of 02_bait.sh
REFDIR="${REFDIR:-ref}"                        # output of 00_reference.sh
OUTDIR="${OUTDIR:-results/04_consensus}"       # output directory
MIN_DEPTH="${MIN_DEPTH:-10}"                   # sites below this are masked to N
MIN_MQ="${MIN_MQ:-20}"                         # minimum mapping quality
MIN_BQ="${MIN_BQ:-20}"                         # minimum base quality
MAX_DEPTH="${MAX_DEPTH:-5000}"                 # mpileup depth cap
# -----------------------------------------------------------------------------

require_tools bcftools samtools seqkit
samples=$(read_samples "$SAMPLES")

REF=$(find "$REFDIR" -maxdepth 1 -name '*.fasta' | sort | head -1)
[[ -n "$REF" ]] || die "no reference FASTA in ${REFDIR}; run 00_reference.sh first"
require_files "$REF" "${REFDIR}/cox1.bed"

mkdir -p "$OUTDIR"
STATS="${OUTDIR}/consensus_stats.tsv"
printf 'sample\tcox1_length\tmasked_sites\tmasked_fraction\talt_calls\n' > "$STATS"

# process_sample <sample>
# Builds the reference-guided consensus and cuts out COX1 for one sample. An
# ordinary function, invoked by the loop below as "( process_sample "$sample" )"
# so that a failing command, or a die() inside it (e.g. a missing bam), exits
# only that subshell rather than the whole run.
process_sample() {
  local sample="$1" bam out len masked frac alts row
  bam="${BAITDIR}/${sample}.mito.bam"
  require_files "$bam"

  # 1) Pile up and call over the whole mitogenome. Duplicates are excluded by
  #    the mpileup default read filter.
  bcftools mpileup -f "$REF" -d "$MAX_DEPTH" -q "$MIN_MQ" -Q "$MIN_BQ" \
      -a FORMAT/AD -Ou "$bam" 2> "${OUTDIR}/${sample}.mpileup.log" \
    | bcftools call -m -V indels -Oz -o "${OUTDIR}/${sample}.calls.vcf.gz"
  bcftools index -f "${OUTDIR}/${sample}.calls.vcf.gz"

  # 2) Mask everything that is not supported by MIN_DEPTH reads.
  samtools depth -a "$bam" \
    | awk -F'\t' -v d="$MIN_DEPTH" 'BEGIN { OFS = "\t" } $3 < d { print $1, $2 - 1, $2 }' \
    > "${OUTDIR}/${sample}.lowcov.bed"

  bcftools consensus -f "$REF" -m "${OUTDIR}/${sample}.lowcov.bed" \
      "${OUTDIR}/${sample}.calls.vcf.gz" \
    > "${OUTDIR}/${sample}.mito_consensus.fasta" 2> "${OUTDIR}/${sample}.consensus.log"

  # 3) Cut out COX1. seqkit honors the strand column, so a minus-strand COX1
  #    comes back already reverse complemented.
  out="${OUTDIR}/${sample}_COI_routeB.fasta"
  seqkit subseq --bed "${REFDIR}/cox1.bed" "${OUTDIR}/${sample}.mito_consensus.fasta" \
    | seqkit replace -p '^.*$' -r "${sample}_COI_routeB" > "$out"

  len=$(seqkit fx2tab -nl "$out" | cut -f2)
  masked=$(seqkit fx2tab "$out" | cut -f2 | tr -cd 'Nn' | wc -c | tr -d ' ')
  frac=$(awk -v a="$masked" -v b="$len" 'BEGIN { printf "%.4f", (b ? a / b : 0) }')
  # GT="alt" is the documented predicate for "carries a non-reference allele";
  # a literal genotype string such as GT!="0/0" is not portable across
  # bcftools versions.
  alts=$(bcftools view -H -v snps -i 'GT="alt"' \
           -r "$(cat "${REFDIR}/cox1.region")" "${OUTDIR}/${sample}.calls.vcf.gz" | wc -l | tr -d ' ')

  row=$(printf '%s\t%s\t%s\t%s\t%s' "$sample" "$len" "$masked" "$frac" "$alts")
  printf '%s\n' "$row" >> "$STATS"
  printf '%s\n' "$row" > "${OUTDIR}/${sample}.consensus_stats.line"
}

FAILLOG="${OUTDIR}/failed_samples.tsv"
printf 'sample\tstage\terror_log\n' > "$FAILLOG"
nfail=0

n=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
i=0
while IFS=$'\t' read -r sample _ _; do
  i=$((i + 1))

  if [[ -s "${OUTDIR}/${sample}_COI_routeB.fasta" && -s "${OUTDIR}/${sample}.consensus_stats.line" ]]; then
    echo "[${i}/${n}] consensus ${sample}: already done, skipping"
    cat "${OUTDIR}/${sample}.consensus_stats.line" >> "$STATS"
    continue
  fi

  echo "[${i}/${n}] consensus ${sample}"
  errlog="${OUTDIR}/${sample}.consensus_error.log"

  if ( process_sample "$sample" ) 2> "$errlog"; then
    rm -f "$errlog"
  else
    nfail=$((nfail + 1))
    printf '%s\tconsensus\t%s\n' "$sample" "$errlog" >> "$FAILLOG"
    msg "  ${sample}: FAILED, continuing to the next sample. Detail:"
    sed 's/^/    /' "$errlog" >&2
  fi
done <<< "$samples"

echo
column -t "$STATS" 2>/dev/null || cat "$STATS"
echo
if [[ "$nfail" -gt 0 ]]; then
  echo "WARNING: ${nfail} of ${n} sample(s) failed; see ${FAILLOG} and each"
  echo "sample's <sample>.consensus_error.log for the reason. Rerunning this"
  echo "script will retry only those samples; the rest are already done."
  echo "05_validate.sh will skip these samples in turn rather than stopping."
fi
echo "wrote ${OUTDIR}/<sample>_COI_routeB.fasta"
