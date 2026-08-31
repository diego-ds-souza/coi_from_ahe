#!/usr/bin/env bash
#
# Optional: full annotated mitogenome per sample, with MitoFinder.
#
# This is not part of the COI path. 03_assemble.sh already recovers COI, and
# MitoFinder wraps the same assemblers (MEGAHIT here), so it would not produce a
# better COI. What it adds is the rest of the mitogenome: the other twelve
# protein-coding genes, the rRNAs and the tRNAs, annotated, from reads you have
# already sequenced. For a phylogenomic program on wood-associated weevils that
# is a second dataset for the cost of one more run, and it is the reason to
# bother with MitoFinder on a machine that can install it.
#
# The script ends by aligning MitoFinder's COX1 against the COI that came out of
# 03 and 05. They should be effectively identical. A disagreement means one of
# the two assemblies picked up a NUMT or a chimera, which is worth knowing
# before either result is used.
#
# MitoFinder is Python 2.7 only, so it lives in its own environment:
#   conda activate coi_mitofinder
#   bash 06_mitogenome.sh
#
# The comparison at the end needs mafft and biopython, so rerun that part, or
# the whole script, under coi_from_ahe if coi_mitofinder lacks them. The script
# skips the comparison rather than failing when they are absent.
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
SAMPLES="${SAMPLES:-samples.tsv}"              # sample sheet
TRIMDIR="${TRIMDIR:-results/01_trim}"          # output of 01_trim.sh
REFDIR="${REFDIR:-ref}"                        # output of 00_reference.sh
FINALDIR="${FINALDIR:-results/coi_final}"      # output of 05_validate.sh
OUTDIR="${OUTDIR:-results/06_mitogenome}"      # output directory
GENETIC_CODE="${GENETIC_CODE:-5}"              # invertebrate mitochondrial
THREADS="${THREADS:-24}"                       # CPU threads
MEM_GB="${MEM_GB:-32}"                         # assembler memory ceiling
# -----------------------------------------------------------------------------

require_tools mitofinder seqkit readlink
samples=$(read_samples "$SAMPLES")

GB=""
for cand in "$REFDIR"/*.gb; do
  [[ -f "$cand" ]] || continue
  # MitoFinder needs an annotated reference: it builds its per-gene search
  # database from CDS/rRNA features. REFDIR can hold an unannotated MITO_ACC
  # record alongside an annotated ANNOT_ACC donor (00_reference.sh's tblastn
  # fallback), so pick by content, not by alphabetical order.
  grep -q '^     CDS ' "$cand" && { GB="$cand"; break; }
done
[[ -n "$GB" ]] || die "no annotated reference GenBank file (with CDS features) in ${REFDIR}.
     MitoFinder needs gene annotations to work from. Run 00_reference.sh with
     ANNOT_ACC set to an annotated donor if MITO_ACC itself is unannotated."
gb_abs=$(readlink -f "$GB")

mkdir -p "$OUTDIR"
STATS="${OUTDIR}/mitogenome_stats.tsv"
printf 'sample\tcontig_length\tgenes_annotated\tcox1_recovered\tcox1_vs_route_a\n' > "$STATS"

n=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
i=0
while IFS=$'\t' read -r sample _ _; do
  i=$((i + 1))
  read -r r1 r2 < <(trimmed_reads "$TRIMDIR" "$sample")
  r1_abs=$(readlink -f "$r1"); r2_abs=$(readlink -f "$r2")

  # MitoFinder writes into the current directory, so give it one of its own.
  rundir="${OUTDIR}/${sample}_mitofinder"
  mkdir -p "$rundir"

  echo "[${i}/${n}] mitofinder ${sample}"
  ( cd "$rundir" && mitofinder \
      -j "$sample" \
      -1 "$r1_abs" -2 "$r2_abs" \
      -r "$gb_abs" \
      -o "$GENETIC_CODE" \
      -p "$THREADS" \
      -m "$MEM_GB" \
      --tRNA-annotation arwen \
      --megahit ) > "${OUTDIR}/${sample}.mitofinder.log" 2>&1 \
    || die "mitofinder failed for ${sample}; see ${OUTDIR}/${sample}.mitofinder.log"

  # Output layout varies a little between MitoFinder releases. When more than
  # one contig matches the reference, MitoFinder numbers them
  # (*_mtDNA_contig_1.fasta, _2.fasta, ...) rather than writing a single
  # unnumbered file, so collect all candidates and keep the largest as "the"
  # mitogenome; smaller ones are usually NUMTs or unmerged fragments, and stay
  # on disk under rundir for inspection rather than being silently discarded.
  mapfile -t contig_candidates < <(find "$rundir" -name '*_mtDNA_contig_*.fasta' \
    | grep -v '_genes_' | grep -v '_ref')
  genes=$(find "$rundir" -name '*final_genes_NT.fasta' | sort | head -1)

  contig_len="NA"; n_genes="0"; cox1_len="NA"; ident="NA"

  contig=""; contig_len_bp=0
  for cand in "${contig_candidates[@]:-}"; do
    [[ -n "$cand" ]] || continue
    len=$(seqkit fx2tab -nl "$cand" | cut -f2 | head -1)
    if [[ -n "$len" && "$len" -gt "$contig_len_bp" ]]; then
      contig="$cand"; contig_len_bp="$len"
    fi
  done
  if (( ${#contig_candidates[@]} > 1 )); then
    msg "  ${sample}: ${#contig_candidates[@]} contigs matched the reference; using the largest (${contig_len_bp} bp)"
  fi

  if [[ -n "$contig" ]]; then
    cp "$contig" "${OUTDIR}/${sample}_mitogenome.fasta"
    contig_len=$(seqkit fx2tab -nl "${OUTDIR}/${sample}_mitogenome.fasta" | cut -f2 | head -1)
    gbk="${contig%.fasta}.gb"
    [[ -f "$gbk" ]] && cp "$gbk" "${OUTDIR}/${sample}_mitogenome.gb"
  fi

  if [[ -n "$genes" ]]; then
    cp "$genes" "${OUTDIR}/${sample}_genes_NT.fasta"
    n_genes=$(grep -c '^>' "${OUTDIR}/${sample}_genes_NT.fasta" || true)

    # \b prevents COI from also matching COII and COIII.
    seqkit grep -nri -p 'cox1\b' -p 'coi\b' -p 'co1\b' \
      "${OUTDIR}/${sample}_genes_NT.fasta" > "${OUTDIR}/${sample}_COI_mitofinder.fasta" || true
    if [[ -s "${OUTDIR}/${sample}_COI_mitofinder.fasta" ]]; then
      cox1_len=$(seqkit fx2tab -nl "${OUTDIR}/${sample}_COI_mitofinder.fasta" | cut -f2 | head -1)
    fi
  fi

  # Cross-check against the COI from route A, when both exist and the tools for
  # the comparison are present.
  routeA="${FINALDIR}/${sample}_COI.fasta"
  if [[ -s "${OUTDIR}/${sample}_COI_mitofinder.fasta" && -s "$routeA" ]] \
     && command -v mafft >/dev/null && PY=$(resolve_python_quiet); then
    cat "$routeA" "${OUTDIR}/${sample}_COI_mitofinder.fasta" > "${OUTDIR}/${sample}.compare.fasta"
    mafft --auto --adjustdirection --quiet "${OUTDIR}/${sample}.compare.fasta" \
      > "${OUTDIR}/${sample}.compare.aln.fasta"
    if "$PY" "${SCRIPT_DIR}/lib/aln_identity.py" "${OUTDIR}/${sample}.compare.aln.fasta" \
         > "${OUTDIR}/${sample}.compare.tsv" 2>/dev/null; then
      ident=$(awk -F'\t' 'NR == 2 { print $6 }' "${OUTDIR}/${sample}.compare.tsv")
    fi
  fi

  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$contig_len" "$n_genes" "$cox1_len" "${ident:-NA}" >> "$STATS"
done <<< "$samples"

echo
column -t "$STATS" 2>/dev/null || cat "$STATS"
echo
echo "wrote ${OUTDIR}/<sample>_mitogenome.fasta and .gb"
echo
echo "A complete beetle mitogenome is roughly 15 to 20 kb with 37 genes."
echo "cox1_vs_route_a below about 99.5 means MitoFinder and route A disagree;"
echo "resolve that before using either sequence."
