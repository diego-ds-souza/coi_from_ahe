#!/usr/bin/env bash
#
# Route A: de novo reconstruction of COI.
#
# Two choices drive this script.
#
# METHOD selects the assembler:
#   megahit (default) memory efficient, the right choice for assembling whole
#           AHE libraries of 15 to 20 million read pairs.
#   spades  more accurate per base, considerably heavier. Reasonable on baited
#           reads, expensive on a full library.
#
# ASM_INPUT selects what is assembled:
#   full (default) the whole trimmed library. Independent of the reference
#           mitogenome, which is the point of having two routes: route B is
#           reference guided, so route A should not be.
#   baited  the read pairs recruited in 02. Fast, but it inherits 02's
#           sensitivity to how divergent the reference is. Use it for a quick
#           first look, not for the final answer.
#
# COX1 is then located in the contigs by tblastn against the reference COX1
# protein, falling back to blastn. The protein search is deliberate: across
# weevil tribes, nucleotide identity at COI drops far enough to weaken a blastn
# search while the translated search stays strong.
#
# MitoFinder is not used here. It wraps the same assemblers and adds whole
# mitogenome annotation, which is a separate deliverable; see 06_mitogenome.sh.
#
# Run from the project root:
#   conda activate coi_from_ahe
#   bash 03_assemble.sh
#   METHOD=spades ASM_INPUT=baited bash 03_assemble.sh    # quick first look
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
SAMPLES="${SAMPLES:-samples.tsv}"              # sample sheet
TRIMDIR="${TRIMDIR:-results/01_trim}"          # output of 01_trim.sh
BAITDIR="${BAITDIR:-results/02_bait}"          # output of 02_bait.sh
REFDIR="${REFDIR:-ref}"                        # output of 00_reference.sh
OUTDIR="${OUTDIR:-results/03_denovo}"          # output directory
METHOD="${METHOD:-megahit}"                    # megahit | spades
ASM_INPUT="${ASM_INPUT:-full}"                 # full | baited
SPADES_OPTS="${SPADES_OPTS:---careful}"        # set empty if your SPAdes rejects it
MIN_CONTIG="${MIN_CONTIG:-300}"                # megahit minimum contig length
GENETIC_CODE="${GENETIC_CODE:-5}"              # invertebrate mitochondrial
THREADS="${THREADS:-24}"                       # CPU threads
MEM_GB="${MEM_GB:-32}"                         # assembler memory ceiling
# -----------------------------------------------------------------------------

case "$METHOD" in
  megahit|spades) ;;
  mitofinder) die "MitoFinder is not part of route A here; see 06_mitogenome.sh" ;;
  *) die "METHOD must be megahit or spades, got: ${METHOD}" ;;
esac
case "$ASM_INPUT" in
  baited|full) ;;
  *) die "ASM_INPUT must be baited or full, got: ${ASM_INPUT}" ;;
esac
require_tools tblastn blastn makeblastdb seqkit
[[ "$METHOD" == "spades"  ]] && require_tools spades.py
[[ "$METHOD" == "megahit" ]] && require_tools megahit

samples=$(read_samples "$SAMPLES")
mkdir -p "$OUTDIR"

# assembly_input <sample> -> "r1<TAB>r2"
assembly_input() {
  local sample="$1"
  if [[ "$ASM_INPUT" == "baited" ]]; then
    local r1="${BAITDIR}/${sample}.bait_R1.fq.gz" r2="${BAITDIR}/${sample}.bait_R2.fq.gz"
    require_files "$r1" "$r2"
    printf '%s\t%s\n' "$r1" "$r2"
  else
    trimmed_reads "$TRIMDIR" "$sample"
  fi
}

# locate_cox1 <contigs> <sample> <out_fasta>
# Best-scoring contig by tblastn, then the span of all its HSPs. Falls back to
# blastn when the translated search finds nothing.
locate_cox1() {
  local contigs="$1" sample="$2" outfile="$3"
  local hits="${OUTDIR}/${sample}.cox1_hits.tsv"
  local db="${OUTDIR}/${sample}_blastdb"
  require_files "${REFDIR}/cox1.faa" "${REFDIR}/cox1.fna"

  # An indexed database rather than -subject. With a whole-library assembly the
  # contig set runs to hundreds of thousands of sequences, and -subject scans
  # them linearly with no index and no threading.
  mkdir -p "$db"
  makeblastdb -in "$contigs" -dbtype nucl -out "${db}/contigs" \
    > "${OUTDIR}/${sample}.makeblastdb.log" 2>&1

  tblastn -query "${REFDIR}/cox1.faa" -db "${db}/contigs" \
    -db_gencode "$GENETIC_CODE" -evalue 1e-5 -max_target_seqs 10 \
    -num_threads "$THREADS" \
    -outfmt '6 sseqid sstart send pident length bitscore' \
    | sort -k6,6gr > "$hits" || true

  if [[ ! -s "$hits" ]]; then
    msg "  tblastn found no COX1 in ${sample}, falling back to blastn"
    blastn -query "${REFDIR}/cox1.fna" -db "${db}/contigs" \
      -evalue 1e-10 -max_target_seqs 10 -num_threads "$THREADS" \
      -outfmt '6 sseqid sstart send pident length bitscore' \
      | sort -k6,6gr > "$hits" || true
  fi
  [[ -s "$hits" ]] || die "no COX1 hit in the ${sample} assembly.
       Either the mitochondrial fraction is too low, or, if ASM_INPUT=baited,
       the reference was too divergent to recruit reads in 02. Retry with
       METHOD=megahit ASM_INPUT=full."

  # Span every HSP on the best-scoring contig, so that a COX1 broken into
  # several HSPs by divergence is not truncated.
  blast_span_bed "$hits" COX1 > "${OUTDIR}/${sample}.cox1.bed"

  seqkit subseq --bed "${OUTDIR}/${sample}.cox1.bed" "$contigs" \
    | seqkit replace -p '^.*$' -r "${sample}_COI_routeA" > "$outfile"
}

# process_sample <sample>
# Runs the chosen assembler and locates COX1 for one sample. This is an
# ordinary function, not a subshell by itself: the loop below invokes it as
# "( process_sample "$sample" )" so that a failing command, or a die() inside
# it (e.g. no COX1 hit), exits only that subshell. set -e still applies inside
# it, so any step failing aborts the rest of this sample's work.
process_sample() {
  local sample="$1" out r1 r2 asmdir contigs
  out="${OUTDIR}/${sample}_COI_routeA.fasta"

  case "$METHOD" in

    spades)
      read -r r1 r2 < <(assembly_input "$sample")
      asmdir="${OUTDIR}/${sample}_spades"
      # shellcheck disable=SC2086  # SPADES_OPTS is intentionally word split
      spades.py -1 "$r1" -2 "$r2" -o "$asmdir" \
        -t "$THREADS" -m "$MEM_GB" -k 21,33,55,77 $SPADES_OPTS \
        > "${OUTDIR}/${sample}.spades.log" 2>&1 \
        || die "spades failed for ${sample}; see ${OUTDIR}/${sample}.spades.log"
      contigs="${asmdir}/contigs.fasta"
      require_files "$contigs"
      locate_cox1 "$contigs" "$sample" "$out"
      ;;

    megahit)
      read -r r1 r2 < <(assembly_input "$sample")
      asmdir="${OUTDIR}/${sample}_megahit"
      # megahit refuses to write into an existing directory.
      rm -rf "$asmdir"
      megahit -1 "$r1" -2 "$r2" -o "$asmdir" \
        -t "$THREADS" -m $(( MEM_GB * 1024 * 1024 * 1024 )) \
        --min-contig-len "$MIN_CONTIG" \
        > "${OUTDIR}/${sample}.megahit.log" 2>&1 \
        || die "megahit failed for ${sample}; see ${OUTDIR}/${sample}.megahit.log"
      contigs="${asmdir}/final.contigs.fa"
      require_files "$contigs"
      locate_cox1 "$contigs" "$sample" "$out"
      ;;

  esac

  printf '  %s: %s bp\n' "$sample" "$(seqkit fx2tab -nl "$out" | cut -f2)"
}

FAILLOG="${OUTDIR}/failed_samples.tsv"
printf 'sample\tstage\terror_log\n' > "$FAILLOG"
nfail=0

n=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
i=0
while IFS=$'\t' read -r sample _ _; do
  i=$((i + 1))
  out="${OUTDIR}/${sample}_COI_routeA.fasta"

  if [[ -s "$out" ]]; then
    echo "[${i}/${n}] ${METHOD} ${sample}: already assembled ($(seqkit fx2tab -nl "$out" | cut -f2) bp), skipping"
    continue
  fi

  echo "[${i}/${n}] ${METHOD} ${sample} (${ASM_INPUT} reads)"
  errlog="${OUTDIR}/${sample}.assemble_error.log"

  if ( process_sample "$sample" ) 2> "$errlog"; then
    rm -f "$errlog"
  else
    nfail=$((nfail + 1))
    printf '%s\tassemble\t%s\n' "$sample" "$errlog" >> "$FAILLOG"
    msg "  ${sample}: FAILED, continuing to the next sample. Detail:"
    sed 's/^/    /' "$errlog" >&2
  fi
done <<< "$samples"

echo
if [[ "$nfail" -gt 0 ]]; then
  echo "WARNING: ${nfail} of ${n} sample(s) failed; see ${FAILLOG} and each"
  echo "sample's <sample>.assemble_error.log for the reason and what to retry."
  echo "Rerunning this script will retry only those samples; the rest are"
  echo "already assembled. 04_consensus.sh is unaffected (it does not read 03's"
  echo "output); 05_validate.sh will skip any sample missing its route A fasta"
  echo "rather than stopping."
fi
echo "wrote ${OUTDIR}/<sample>_COI_routeA.fasta"
echo "A recovered length far from about 1530 bp means the HSP span is wrong;"
echo "check ${OUTDIR}/<sample>.cox1_hits.tsv before trusting it."
