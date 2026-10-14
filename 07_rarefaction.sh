#!/usr/bin/env bash
#
# Optional: rarefaction of mitochondrial coverage, to estimate how much
# sequencing each specimen would need for a complete mitogenome.
#
# For every sample the trimmed library is subsampled at a series of fractions
# and, at each one, duplicates are re-marked and the unique depth across COX1
# is measured. The curve of unique depth against sequencing effort shows both
# how fast mitochondrial coverage accrues and where the library runs out of
# unique molecules. lib/rarefaction_model.py fits that curve, extrapolates it,
# and converts a target depth into raw read pairs to order.
#
# Each library is mapped to the closest target available for that sample, in
# this order, recorded as ref_source:
#
#   own_mitogenome  the sample's own MitoFinder contigs from 06, when one of
#                   them carries COX1. No reference divergence, and long
#                   enough that both mates of nearly every pair map.
#   own_coi         the sample's own COI from 05. No divergence either, but
#                   at about 1.5 kb many pairs have one mate off the end, and
#                   duplicates among such reads are called on one coordinate
#                   only. At high depth unrelated reads then collide and are
#                   marked as duplicates, so the depth and the library ceiling
#                   are underestimated. That matters little at the low depths
#                   typical of samples that reach this fallback.
#   reference       the 02 alignment against the reference mitogenome (ref/),
#                   which for most samples is about 15 percent divergent in COI.
#                   Some reads fail to map, so depth is a lower bound.
#
# Why not the reference for everyone: at that divergence, beyond the 10 to 12
# percent at which 02_bait.sh says mapping sensitivity falls off, reference
# depth reflects divergence as well as library content.
#
# COX1 depth stands in for mitogenome depth: mitochondrial genes are present at
# the same copy number, so their read depth tracks one another closely, except
# in the AT-rich control region.
#
# Subsampling happens after mapping, on the reads that map, which is
# equivalent to subsampling the whole library: reads that do not map contribute
# nothing at any fraction. samtools view -s selects by read name, so mates stay
# together and duplicate marking stays valid in every subsample.
#
# Run from the project root:
#   conda activate coi_from_ahe
#   bash 07_rarefaction.sh
#   python lib/rarefaction_model.py
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
SAMPLES="${SAMPLES:-samples.tsv}"              # sample sheet
TRIMDIR="${TRIMDIR:-results/01_trim}"          # output of 01_trim.sh
BAITDIR="${BAITDIR:-results/02_bait}"          # output of 02_bait.sh (fallback)
REFDIR="${REFDIR:-ref}"                        # output of 00_reference.sh
FINALDIR="${FINALDIR:-results/coi_final}"      # output of 05_validate.sh
MITODIR="${MITODIR:-results/06_mitogenome}"    # output of 06_mitogenome.sh (read only)
OUTDIR="${OUTDIR:-results/07_rarefaction}"     # output directory
MIN_DEPTH="${MIN_DEPTH:-10}"                   # depth called sufficient per site
MIN_ALIGNED_FRAC="${MIN_ALIGNED_FRAC:-0.8}"    # minimum aligned fraction of a read
THREADS="${THREADS:-12}"                      # CPU threads per sample
JOBS="${JOBS:-4}"                              # samples processed concurrently
SEED="${SEED:-42}"                             # subsampling seed
FRACTIONS="${FRACTIONS:-0.005 0.01 0.02 0.05 0.1 0.2 0.35 0.5 0.75 1}"
# -----------------------------------------------------------------------------

require_tools bwa samtools blastn
samples=$(read_samples "$SAMPLES")
# Samples run concurrently, so a name listed twice would have two jobs
# writing the same files.
dups=$(printf '%s\n' "$samples" | cut -f1 | sort | uniq -d)
[[ -z "$dups" ]] || die "sample(s) listed more than once in ${SAMPLES}: $(echo $dups)"
mkdir -p "$OUTDIR"

REF_REGION=$(cat "${REFDIR}/cox1.region")

# rarefy <sample> <bam> <region> <ref_source>
# Writes one row per subsampling fraction to <sample>.rarefaction.tsv.
rarefy() {
  local sample="$1" bam="$2" region="$3" source="$4" f sub mapped uniq dmean dmed dfrac
  local out="${OUTDIR}/${sample}.rarefaction.tsv"
  printf 'sample\tref_source\tfraction\tmito_reads\tmito_unique\tcox1_mean_depth\tcox1_median_depth\tcox1_frac_ge_%s\n' \
    "$MIN_DEPTH" > "${out}.part"
  for f in $FRACTIONS; do
    sub="${OUTDIR}/${sample}.sub.bam"
    # -s SEED.FRACTION keeps a read name hash below FRACTION; at 1 the whole
    # file passes. Reads aligned over less than MIN_ALIGNED_FRAC of their
    # length are dropped: under the permissive settings of 02 most reads in
    # the COX1 window are short local hits, clipped to 25 to 60 bp, that are
    # not from this sample's mitogenome, and they would inflate depth several
    # fold. markdup -c clears the duplicate flags inherited from the full
    # library so that duplication is judged within the subsample only.
    local subsample=(-s "${SEED}${f#0}")
    [[ "$f" == "1" ]] && subsample=()
    samtools view -u "${subsample[@]}" -e "(qlen - sclen) >= ${MIN_ALIGNED_FRAC} * qlen" "$bam" \
      | samtools markdup -c -@ "$THREADS" -T "${OUTDIR}/${sample}.mdsub" - "$sub" 2> /dev/null
    samtools index "$sub"
    mapped=$(samtools view -c -F 0x904 "$sub" "$region")
    uniq=$(samtools view -c -F 0xD04 "$sub" "$region")
    read -r dmean dmed dfrac < <(depth_stats "$sub" "$region" "$MIN_DEPTH")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$sample" "$source" "$f" "$mapped" "$uniq" "$dmean" "$dmed" "$dfrac" >> "${out}.part"
  done
  rm -f "$sub" "${sub}.bai"
  mv "${out}.part" "$out"
}

# own_mitogenome_target <sample> <workdir>
# Writes <workdir>/target.fasta (the sample's MitoFinder contigs, renamed) and
# prints the COX1 region on it, located by blastn of MitoFinder's own COX1.
# Returns 1 when there are no contigs or none of them carries COX1.
own_mitogenome_target() {
  local sample="$1" work="$2" cox1 hit
  local -a contigs
  cox1="${MITODIR}/${sample}_COI_mitofinder.fasta"
  [[ -s "$cox1" ]] || return 1
  # Same selection as 06_mitogenome.sh: confined to *_Final_Results, where the
  # names cannot collide with the reference copies kept in _tmp.
  mapfile -t contigs < <(find "${MITODIR}/${sample}_mitofinder" -path '*_Final_Results/*' \
      \( -name "${sample}_mtDNA_contig.fasta" -o -name "${sample}_mtDNA_contig_[0-9]*.fasta" \) \
    | grep -v '_genes_' | sort)
  (( ${#contigs[@]} > 0 )) || return 1
  cat "${contigs[@]}" | awk -v s="$sample" '/^>/ { printf ">%s_c%d\n", s, ++i; next } { print }' \
    > "${work}/target.fasta"
  hit=$(blastn -task blastn -query "$cox1" -subject "${work}/target.fasta" \
          -outfmt '6 sseqid sstart send bitscore' 2> /dev/null \
        | sort -t$'\t' -k4,4gr | head -1)
  [[ -n "$hit" ]] || return 1
  awk -F'\t' '{ a = ($2 < $3) ? $2 : $3; b = ($2 < $3) ? $3 : $2; printf "%s:%d-%d\n", $1, a, b }' <<< "$hit"
}

# map_to <sample> <target.fasta> <bam>
# Maps the whole trimmed library to a per-sample target and keeps only pairs
# with at least one mate mapped.
map_to() {
  local sample="$1" target="$2" bam="$3" r1 r2
  read -r r1 r2 < <(trimmed_reads "$TRIMDIR" "$sample")
  bwa index "$target" 2> /dev/null
  # Default bwa settings: the target is the sample's own sequence, so the
  # permissive settings of 02, which exist for a divergent reference, would
  # only admit NUMT and contaminant reads. Pairs where neither mate maps are
  # dropped (-G 12) before duplicate tags are set up, which keeps the BAM to
  # the mitochondrial reads alone.
  bwa mem -t "$THREADS" "$target" "$r1" "$r2" 2> "${OUTDIR}/${sample}.bwa.log" \
    | samtools view -u -G 12 - \
    | samtools fixmate -u -m - - \
    | samtools sort -u -@ "$THREADS" -T "${OUTDIR}/${sample}.sorttmp" - \
    | samtools markdup -@ "$THREADS" -T "${OUTDIR}/${sample}.mdtmp" - "$bam" 2> /dev/null
  samtools index "$bam"
}

# process_sample <sample>
# Builds the per-sample alignment against the closest available target (see
# the header) and runs the rarefaction on it.
process_sample() {
  local sample="$1" coi work bam name region
  coi="${FINALDIR}/${sample}_COI.fasta"
  bam="${OUTDIR}/${sample}.cox1.bam"
  work="${OUTDIR}/${sample}_ref"
  mkdir -p "$work"

  if region=$(own_mitogenome_target "$sample" "$work"); then
    map_to "$sample" "${work}/target.fasta" "$bam"
    printf '%s\n' "$region" > "${work}/cox1.region"
    rarefy "$sample" "$bam" "$region" own_mitogenome
  elif [[ -s "$coi" ]]; then
    cp "$coi" "${work}/target.fasta"
    name=$(head -1 "$coi" | sed 's/^>//; s/[[:space:]].*//')
    map_to "$sample" "${work}/target.fasta" "$bam"
    rarefy "$sample" "$bam" "$name" own_coi
  else
    require_files "${BAITDIR}/${sample}.mito.bam"
    samtools view -b -G 12 -@ "$THREADS" "${BAITDIR}/${sample}.mito.bam" > "$bam"
    samtools index "$bam"
    rarefy "$sample" "$bam" "$REF_REGION" reference
  fi
}

FAILLOG="${OUTDIR}/failed_samples.tsv"
printf 'sample\tstage\terror_log\n' > "$FAILLOG"

# run_one <i> <n> <sample>
# Wrapper so that each sample can run in the background with its own error log.
run_one() {
  local i="$1" n="$2" sample="$3" errlog="${OUTDIR}/${3}.rarefaction_error.log"
  if [[ -s "${OUTDIR}/${sample}.rarefaction.tsv" ]]; then
    echo "[${i}/${n}] rarefaction ${sample}: already done, skipping"
    return 0
  fi
  echo "[${i}/${n}] rarefaction ${sample}"
  if ( process_sample "$sample" ) 2> "$errlog"; then
    rm -f "$errlog"
  else
    printf '%s\trarefaction\t%s\n' "$sample" "$errlog" >> "$FAILLOG"
    msg "  ${sample}: FAILED, continuing. Detail:"
    sed 's/^/    /' "$errlog" >&2
  fi
}

n=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
i=0
while IFS=$'\t' read -r sample _ _; do
  i=$((i + 1))
  while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || true; done
  run_one "$i" "$n" "$sample" &
done <<< "$samples"
wait

# Combined table, in sample-sheet order.
COMBINED="${OUTDIR}/rarefaction.tsv"
first=1
while IFS=$'\t' read -r sample _ _; do
  f="${OUTDIR}/${sample}.rarefaction.tsv"
  [[ -s "$f" ]] || continue
  if (( first )); then cat "$f"; first=0; else tail -n +2 "$f"; fi
done <<< "$samples" > "$COMBINED"

nfail=$(( $(wc -l < "$FAILLOG") - 1 ))
echo
echo "wrote ${COMBINED}"
if (( nfail > 0 )); then
  echo "WARNING: ${nfail} of ${n} sample(s) failed; see ${FAILLOG}. Rerunning"
  echo "this script retries only those samples."
fi
echo "next: python lib/rarefaction_model.py"
