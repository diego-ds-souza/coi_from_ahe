#!/usr/bin/env bash
#
# Mapping of the trimmed reads onto the reference mitogenome. This step serves
# three purposes:
#
#   1. Diagnostics. The mitochondrial fraction of the library and the depth
#      across COX1 determine whether COI is recoverable at all. Read these
#      numbers before running anything else. Median COX1 depth above roughly
#      20x means any downstream approach will work; below about 10x, expect the
#      de novo route in 03 to fragment and rely more on 04.
#   2. The alignment consumed by the reference-guided consensus in 04.
#   3. A baited read subset, used by the SPAdes fallback in 03.
#
# The mapping parameters are deliberately permissive because the reference is
# a different species. Even so, sensitivity falls off above roughly 10 to 12
# percent COI divergence, and a low mitochondrial fraction reported here can
# mean a divergent reference rather than a poor library. The de novo route in
# 03 is the check on that, since it does not map reads to the reference.
#
# Duplicates are marked, not removed, so that the duplication rate stays
# visible in the statistics table; every downstream count excludes them.
#
# Run from the project root:
#   conda activate coi_from_ahe
#   bash 02_bait.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
SAMPLES="${SAMPLES:-samples.tsv}"              # sample sheet
TRIMDIR="${TRIMDIR:-results/01_trim}"          # output of 01_trim.sh
REFDIR="${REFDIR:-ref}"                        # output of 00_reference.sh
OUTDIR="${OUTDIR:-results/02_bait}"            # output directory
MIN_DEPTH="${MIN_DEPTH:-10}"                   # depth called sufficient per site
THREADS="${THREADS:-24}"                        # CPU threads
# Permissive alignment: shorter seeds, cheaper mismatches and gaps, lower score
# threshold. Tighten these if the reference is congeneric and the mitochondrial
# fraction looks implausibly high.
BWA_OPTS="${BWA_OPTS:--k 15 -B 3 -O 5,5 -E 2 -T 25}"
# -----------------------------------------------------------------------------

require_tools bwa samtools
samples=$(read_samples "$SAMPLES")

REF=$(find "$REFDIR" -maxdepth 1 -name '*.fasta' | sort | head -1)
[[ -n "$REF" ]] || die "no reference FASTA in ${REFDIR}; run 00_reference.sh first"
require_files "$REF" "${REFDIR}/cox1.region" "${REF}.bwt"
REGION=$(cat "${REFDIR}/cox1.region")

mkdir -p "$OUTDIR"
STATS="${OUTDIR}/bait_stats.tsv"
printf 'sample\ttotal_reads\tmito_mapped\tmito_mapped_nodup\tmito_fraction\tdup_rate\tcox1_mean_depth\tcox1_median_depth\tcox1_frac_ge_%s\n' \
  "$MIN_DEPTH" > "$STATS"

# process_sample <sample>
# Maps one sample, appends its bait_stats.tsv row, and extracts baited read
# pairs. An ordinary function, invoked by the loop below as
# "( process_sample "$sample" )" so that a failing command exits only that
# subshell rather than the whole run.
process_sample() {
  local sample="$1" bam r1 r2 total mapped mapped_nodup frac dup dmean dmed dfrac row
  bam="${OUTDIR}/${sample}.mito.bam"
  read -r r1 r2 < <(trimmed_reads "$TRIMDIR" "$sample")

  # 1) Map, fix mate information, sort, mark duplicates. Kept as one stream so
  #    that no intermediate BAM is written.
  # shellcheck disable=SC2086  # BWA_OPTS is intentionally word split
  bwa mem $BWA_OPTS -t "$THREADS" "$REF" "$r1" "$r2" 2> "${OUTDIR}/${sample}.bwa.log" \
    | samtools fixmate -u -m - - \
    | samtools sort -u -@ "$THREADS" -T "${OUTDIR}/${sample}.sorttmp" - \
    | samtools markdup -@ "$THREADS" -T "${OUTDIR}/${sample}.mdtmp" - "$bam"
  samtools index "$bam"

  # 2) Counts. 0x900 drops secondary and supplementary records, 0x4 unmapped,
  #    0x400 duplicates.
  total=$(samtools view -c -F 0x900 "$bam")
  mapped=$(samtools view -c -F 0x904 "$bam")
  mapped_nodup=$(samtools view -c -F 0xD04 "$bam")
  frac=$(awk -v a="$mapped_nodup" -v b="$total" 'BEGIN { printf "%.5f", (b ? a / b : 0) }')
  dup=$(awk -v a="$mapped_nodup" -v b="$mapped" 'BEGIN { printf "%.3f", (b ? 1 - a / b : 0) }')

  # 3) Depth across COX1 only, counting zero-coverage sites.
  read -r dmean dmed dfrac < <(depth_stats "$bam" "$REGION" "$MIN_DEPTH")
  row=$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$sample" "$total" "$mapped" "$mapped_nodup" "$frac" "$dup" "$dmean" "$dmed" "$dfrac")
  printf '%s\n' "$row" >> "$STATS"
  printf '%s\n' "$row" > "${OUTDIR}/${sample}.bait_stats.line"

  # 4) Baited read pairs for the SPAdes fallback: duplicates excluded, pairs
  #    kept when at least one mate maps.
  samtools view -u -F 0x400 -G 12 "$bam" \
    | samtools collate -u -O -@ "$THREADS" - \
    | samtools fastq -n \
        -1 "${OUTDIR}/${sample}.bait_R1.fq.gz" \
        -2 "${OUTDIR}/${sample}.bait_R2.fq.gz" \
        -0 /dev/null -s /dev/null - 2> "${OUTDIR}/${sample}.fastq.log"
}

FAILLOG="${OUTDIR}/failed_samples.tsv"
printf 'sample\tstage\terror_log\n' > "$FAILLOG"
nfail=0

n=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
i=0
while IFS=$'\t' read -r sample _ _; do
  i=$((i + 1))

  if [[ -s "${OUTDIR}/${sample}.mito.bam" && -s "${OUTDIR}/${sample}.mito.bam.bai" \
        && -s "${OUTDIR}/${sample}.bait_R1.fq.gz" && -s "${OUTDIR}/${sample}.bait_R2.fq.gz" \
        && -s "${OUTDIR}/${sample}.bait_stats.line" ]]; then
    echo "[${i}/${n}] bwa mem ${sample}: already baited, skipping"
    cat "${OUTDIR}/${sample}.bait_stats.line" >> "$STATS"
    continue
  fi

  echo "[${i}/${n}] bwa mem ${sample}"
  errlog="${OUTDIR}/${sample}.bait_error.log"

  if ( process_sample "$sample" ) 2> "$errlog"; then
    rm -f "$errlog"
  else
    nfail=$((nfail + 1))
    printf '%s\tbait\t%s\n' "$sample" "$errlog" >> "$FAILLOG"
    msg "  ${sample}: FAILED, continuing to the next sample. Detail:"
    sed 's/^/    /' "$errlog" >&2
  fi
done <<< "$samples"

echo
column -t "$STATS" 2>/dev/null || cat "$STATS"
echo
if [[ "$nfail" -gt 0 ]]; then
  echo "WARNING: ${nfail} of ${n} sample(s) failed; see ${FAILLOG} and each"
  echo "sample's <sample>.bait_error.log for the reason. Rerunning this script"
  echo "will retry only those samples; the rest are already baited."
  echo
fi
echo "wrote ${STATS}, ${OUTDIR}/<sample>.mito.bam and baited read pairs"
