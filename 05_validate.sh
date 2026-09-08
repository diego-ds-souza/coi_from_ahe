#!/usr/bin/env bash
#
# Validation and final deliverable. Four independent checks are run on the
# route A sequence, and the two routes are compared:
#
#   1. Read support. The trimmed library is mapped back onto the de novo COI
#      and per-site depth is summarized. A contig with no read support over
#      part of its length is a chimera or a misassembly.
#   2. Competing alleles. Any site where a second allele reaches
#      MINOR_AF of a well covered position is reported. In a haploid,
#      maternally inherited marker these should be absent; a cluster of them is
#      the signature of a co-amplified NUMT, of two individuals in one library,
#      or of index hopping from a neighboring well.
#   3. Reading frame. Internal stop codons under genetic code 5 are the
#      standard NUMT and frameshift test.
#   4. Folmer window. The classic barcode primer binding sites are located in
#      the recovered sequence and the corresponding fragment is written out,
#      which is what makes the result comparable with BOLD and with published
#      COI barcodes.
#
# Route A and route B are aligned with the reference CDS and their pairwise
# identity is reported. Disagreement between routes at well covered sites means
# one of them is wrong; resolve it before using the sequence.
#
# Run from the project root:
#   conda activate coi_from_ahe
#   bash 05_validate.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
SAMPLES="${SAMPLES:-samples.tsv}"              # sample sheet
TRIMDIR="${TRIMDIR:-results/01_trim}"          # output of 01_trim.sh
DENOVODIR="${DENOVODIR:-results/03_denovo}"    # output of 03_assemble.sh
CONSDIR="${CONSDIR:-results/04_consensus}"     # output of 04_consensus.sh
REFDIR="${REFDIR:-ref}"                        # output of 00_reference.sh
OUTDIR="${OUTDIR:-results/05_validate}"        # diagnostics directory
FINALDIR="${FINALDIR:-results/coi_final}"      # deliverable directory
MIN_DEPTH="${MIN_DEPTH:-10}"                   # depth called sufficient per site
MINOR_AF="${MINOR_AF:-0.20}"                   # competing-allele threshold
MIN_SITE_DP="${MIN_SITE_DP:-20}"               # depth required to judge a site
GENETIC_CODE="${GENETIC_CODE:-5}"              # invertebrate mitochondrial
THREADS="${THREADS:-24}"                        # CPU threads
# Folmer et al. (1994) universal barcode primers. HCO2198 is the reverse
# primer, so it is found on the opposite strand; seqkit searches both.
LCO1490="${LCO1490:-GGTCAACAAATCATAAAGATATTGG}"
HCO2198="${HCO2198:-TAAACTTCAGGGTGACCAAAAAATCA}"
PRIMER_MM="${PRIMER_MM:-4}"                    # mismatches allowed per primer
# -----------------------------------------------------------------------------

require_tools bwa samtools bcftools seqkit mafft
require_python
samples=$(read_samples "$SAMPLES")
require_files "${REFDIR}/cox1.fna"

mkdir -p "$OUTDIR" "$FINALDIR"
SUMMARY="${OUTDIR}/validation_summary.tsv"
printf 'sample\tcoi_length\tmean_depth\tmedian_depth\tfrac_ge_%s\tcompeting_allele_sites\tinternal_stops\torf_verdict\tidentity_A_vs_B\tfolmer_bp\n' \
  "$MIN_DEPTH" > "$SUMMARY"

# process_sample <sample>
# Validates one sample. An ordinary function, invoked by the loop below as
# "( process_sample "$sample" )" so that a failing command, or a die() inside
# it (most commonly a missing route A or route B fasta from an earlier
# failed sample), exits only that subshell rather than the whole run.
process_sample() {
  local sample="$1" routeA routeB work r1 r2 contig coi_len dmean dmed dfrac
  local competing stops verdict ident folmer_bp range row

  routeA="${DENOVODIR}/${sample}_COI_routeA.fasta"
  routeB="${CONSDIR}/${sample}_COI_routeB.fasta"
  read -r r1 r2 < <(trimmed_reads "$TRIMDIR" "$sample")
  require_files "$routeA" "$routeB"

  work="${OUTDIR}/${sample}"
  mkdir -p "$work"

  # 1) Map the library back onto the de novo COI. Same species now, so default
  #    bwa stringency is appropriate and mismapped off-target reads are less
  #    likely to be retained.
  cp "$routeA" "${work}/target.fasta"
  bwa index "${work}/target.fasta" 2>/dev/null
  bwa mem -t "$THREADS" "${work}/target.fasta" "$r1" "$r2" 2> "${work}/bwa.log" \
    | samtools fixmate -u -m - - \
    | samtools sort -u -@ "$THREADS" -T "${work}/sorttmp" - \
    | samtools markdup -@ "$THREADS" -T "${work}/mdtmp" - "${work}/remap.bam"
  samtools index "${work}/remap.bam"

  contig=$(seqkit fx2tab -n "${work}/target.fasta" | head -1 | cut -f1)
  coi_len=$(seqkit fx2tab -nl "${work}/target.fasta" | cut -f2)
  read -r dmean dmed dfrac < <(depth_stats "${work}/remap.bam" "${contig}" "$MIN_DEPTH")

  # 2) Competing alleles at well covered sites.
  bcftools mpileup -f "${work}/target.fasta" -d 5000 -q 20 -Q 20 -a FORMAT/AD -Ou \
      "${work}/remap.bam" 2> "${work}/mpileup.log" \
    | bcftools call -m -Ov \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
    | awk -F'\t' -v thr="$MINOR_AF" -v mindp="$MIN_SITE_DP" 'BEGIN { OFS = "\t" }
        {
          split($5, counts, ",")
          total = 0; best = 0
          for (k in counts) { total += counts[k]; if (counts[k] + 0 > best) best = counts[k] + 0 }
          if (total >= mindp) {
            minor = (total - best) / total
            if (minor >= thr) print $1, $2, $3, $4, total, sprintf("%.3f", minor)
          }
        }' > "${work}/competing_alleles.tsv"
  competing=$(wc -l < "${work}/competing_alleles.tsv" | tr -d ' ')

  # 3) Reading frame under the invertebrate mitochondrial code.
  "$PY" "${SCRIPT_DIR}/lib/check_orf.py" --transl-table "$GENETIC_CODE" "$routeA" \
    > "${work}/orf.tsv"
  stops=$(awk -F'\t' 'NR == 2 { print $5 }' "${work}/orf.tsv")
  verdict=$(awk -F'\t' 'NR == 2 { print $8 }' "${work}/orf.tsv")

  # 4) Route A against route B against the reference CDS.
  cat "$routeA" "$routeB" "${REFDIR}/cox1.fna" > "${work}/routes.fasta"
  mafft --auto --adjustdirection --quiet "${work}/routes.fasta" > "${work}/routes.aln.fasta"
  "$PY" "${SCRIPT_DIR}/lib/aln_identity.py" "${work}/routes.aln.fasta" > "${work}/identity.tsv"
  ident=$(awk -F'\t' '$1 ~ /routeA/ && $2 ~ /routeB/ { print $6 }' "${work}/identity.tsv")
  ident="${ident:-NA}"

  # 5) Folmer window, and the final deliverables.
  seqkit locate -i -m "$PRIMER_MM" -p "$LCO1490" -p "$HCO2198" "$routeA" \
    > "${work}/folmer_hits.tsv"
  folmer_bp="NA"
  if [[ $(tail -n +2 "${work}/folmer_hits.tsv" | wc -l | tr -d ' ') -eq 2 ]]; then
    range=$(tail -n +2 "${work}/folmer_hits.tsv" \
      | awk -F'\t' 'NR == 1 { s = $5; e = $6 }
          { if ($5 < s) s = $5; if ($6 > e) e = $6 }
          END { print s ":" e }')
    seqkit subseq -r "$range" "$routeA" \
      | seqkit replace -p '^.*$' -r "${sample}_COI_folmer" \
      > "${FINALDIR}/${sample}_COI_folmer.fasta"
    folmer_bp=$(seqkit fx2tab -nl "${FINALDIR}/${sample}_COI_folmer.fasta" | cut -f2)
  else
    echo "  note: Folmer primer sites not both found in ${sample}; the recovered" >&2
    echo "        sequence may not span the standard barcode region." >&2
  fi

  seqkit replace -p '^.*$' -r "${sample}_COI" "$routeA" > "${FINALDIR}/${sample}_COI.fasta"

  row=$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$sample" "$coi_len" "$dmean" "$dmed" "$dfrac" "$competing" \
    "$stops" "$verdict" "$ident" "$folmer_bp")
  printf '%s\n' "$row" >> "$SUMMARY"
  printf '%s\n' "$row" > "${work}/validate_stats.line"
}

FAILLOG="${OUTDIR}/failed_samples.tsv"
printf 'sample\tstage\terror_log\n' > "$FAILLOG"
nfail=0

n=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
i=0
while IFS=$'\t' read -r sample _ _; do
  i=$((i + 1))

  if [[ -s "${FINALDIR}/${sample}_COI.fasta" && -s "${OUTDIR}/${sample}/validate_stats.line" ]]; then
    echo "[${i}/${n}] validate ${sample}: already done, skipping"
    cat "${OUTDIR}/${sample}/validate_stats.line" >> "$SUMMARY"
    continue
  fi

  echo "[${i}/${n}] validate ${sample}"
  errlog="${OUTDIR}/${sample}.validate_error.log"

  if ( process_sample "$sample" ) 2> "$errlog"; then
    rm -f "$errlog"
  else
    nfail=$((nfail + 1))
    printf '%s\tvalidate\t%s\n' "$sample" "$errlog" >> "$FAILLOG"
    msg "  ${sample}: FAILED, continuing to the next sample. Detail:"
    sed 's/^/    /' "$errlog" >&2
  fi
done <<< "$samples"

echo
column -t "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
echo
if [[ "$nfail" -gt 0 ]]; then
  echo "WARNING: ${nfail} of ${n} sample(s) failed validation; see ${FAILLOG} and"
  echo "each sample's <sample>.validate_error.log for the reason. They are absent"
  echo "from ${SUMMARY} and from ${FINALDIR}. Rerunning this script will retry"
  echo "only those samples; the rest are already done."
  echo
fi
echo "wrote ${FINALDIR}/<sample>_COI.fasta and ${SUMMARY}"
echo
echo "How to read the summary:"
echo "  frac_ge_${MIN_DEPTH}          should be 1.000; anything lower means part of the COI is thinly supported"
echo "  competing_allele_sites  should be 0; nonzero suggests a NUMT, a mixed library, or index hopping"
echo "  orf_verdict             must be PASS; FLAG means internal stop codons under code ${GENETIC_CODE}"
echo "  identity_A_vs_B         expect > 99.5 percent; lower means the two routes disagree"
