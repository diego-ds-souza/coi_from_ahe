#!/usr/bin/env bash
#
# Shared helpers for the coi_from_ahe pipeline. This file is sourced by every
# numbered script and is not meant to be executed on its own.

# msg <text>
# Timestamped progress line, written to stderr so that stdout stays clean for
# any script whose output is piped.
msg() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# die <text>
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# require_tools <tool>...
# Fails fast with an actionable message rather than part way through a run.
require_tools() {
  local t
  for t in "$@"; do
    command -v "$t" >/dev/null \
      || die "$t not found on PATH (did you 'conda activate coi_from_ahe'?)"
  done
}

# resolve_python_quiet
# Prints a Python 3.8+ interpreter that can import biopython, or returns 1
# without a message. For optional steps that should be skipped, not aborted.
resolve_python_quiet() {
  local cand
  for cand in ${PYTHON:-} python3 python; do
    command -v "$cand" >/dev/null 2>&1 || continue
    "$cand" -c 'import sys, Bio; sys.exit(0 if sys.version_info >= (3, 8) else 1)' \
      2>/dev/null || continue
    printf '%s\n' "$cand"
    return 0
  done
  return 1
}

# require_python
# Sets the global PY to a Python 3.8+ interpreter that can import biopython.
# Use this instead of require_tools python: on many servers "python" is still
# Python 2, which fails on type annotations with a SyntaxError several lines
# into a helper rather than at the point of invocation, and that error is easy
# to misread as a problem with the data.
require_python() {
  local cand
  PY=""
  for cand in ${PYTHON:-} python3 python; do
    command -v "$cand" >/dev/null 2>&1 || continue
    "$cand" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)' 2>/dev/null || continue
    PY="$cand"
    break
  done
  [[ -n "$PY" ]] || die "no Python 3.8 or newer on PATH.
       'python' currently resolves to: $(command -v python 2>/dev/null || echo 'nothing')
       Activate the pipeline environment (conda activate coi_from_ahe), or set
       PYTHON=/path/to/python3 for this run."
  "$PY" -c 'import Bio' 2>/dev/null || die "${PY} cannot import biopython.
       The environment is incomplete. Recreate it:
         conda env create -f environment.yml
         conda activate coi_from_ahe"
}

# validate_ncbi_file <path> <gb|fasta>
# A truncated download is not an empty file, and every later error it causes
# points somewhere unhelpful, so downloads are checked for their terminator.
validate_ncbi_file() {
  local f="$1" fmt="$2"
  [[ -s "$f" ]] || return 1
  case "$fmt" in
    gb)    grep -q '^LOCUS' "$f" && tail -5 "$f" | grep -q '^//' ;;
    fasta) head -1 "$f" | grep -q '^>' && [[ "$(wc -l < "$f")" -gt 2 ]] ;;
    *)     return 0 ;;
  esac
}

# fetch_ncbi <accession> <gb|fasta> <outfile>
# efetch with retries and a completeness check. Networks that terminate TLS
# connections early (a common corporate middlebox behavior, seen as
# "curl: (56) ... unexpected eof") truncate downloads silently.
fetch_ncbi() {
  local acc="$1" fmt="$2" out="$3" attempt
  if validate_ncbi_file "$out" "$fmt"; then
    msg "  using existing ${out}"
    return 0
  fi
  for attempt in 1 2 3; do
    efetch -db nuccore -id "$acc" -format "$fmt" > "${out}.part" 2> "${out}.efetch.log" || true
    if validate_ncbi_file "${out}.part" "$fmt"; then
      mv "${out}.part" "$out"
      rm -f "${out}.efetch.log"
      return 0
    fi
    msg "  ${acc} (${fmt}): attempt ${attempt} of 3 incomplete, retrying"
    sleep 3
  done
  rm -f "${out}.part"
  die "could not download ${acc} in ${fmt} format after 3 attempts.
       See ${out}.efetch.log. If this network keeps dropping TLS connections,
       fetch the file elsewhere and copy it to ${out}; the script picks it up."
}

# require_files <path>...
require_files() {
  local f
  for f in "$@"; do
    [[ -f "$f" ]] || die "no such file: $f"
  done
}

# read_samples <samples.tsv>
# Emits one "sample<TAB>r1<TAB>r2" line per data row. Blank lines, comment
# lines and the header row are skipped, and CR line endings are tolerated so
# that a sheet edited on Windows or in Excel does not silently produce paths
# with a trailing carriage return.
#
# The sheet is validated in a first pass and reported through die(), rather
# than by returning a nonzero status: bash does not reliably propagate the
# exit status of a command substitution under 'set -e', so a status-based
# check here would let a malformed sheet through.
read_samples() {
  local sheet="$1" bad rows
  [[ -f "$sheet" ]] || die "no such file: $sheet"

  bad=$(awk -F'\t' '
    { sub(/\r$/, "") }
    /^[[:space:]]*(#|$)/ { next }
    $1 == "sample"       { next }
    NF < 3               { print "  row " NR ": " $0 }' "$sheet")
  [[ -z "$bad" ]] || die "malformed row(s) in ${sheet}:"$'\n'"$bad"

  rows=$(awk -F'\t' 'BEGIN { OFS = "\t" }
    { sub(/\r$/, "") }
    /^[[:space:]]*(#|$)/ { next }
    $1 == "sample"       { next }
    { print $1, $2, $3 }' "$sheet")
  [[ -n "$rows" ]] || die "no data rows in ${sheet}"

  printf '%s\n' "$rows"
}

# trimmed_reads <trimdir> <sample>
# Emits "r1<TAB>r2" for the trimmed pair of one sample, accepting any of the
# usual extensions. This lets 01_trim.sh output be interchangeable with reads
# trimmed elsewhere and symlinked in, without forcing anyone to recompress a
# plain .fastq just to satisfy a filename pattern.
trimmed_reads() {
  local trimdir="$1" sample="$2" ext r1 r2
  for ext in fq.gz fastq.gz fq fastq; do
    r1="${trimdir}/${sample}_R1.${ext}"
    r2="${trimdir}/${sample}_R2.${ext}"
    if [[ -f "$r1" && -f "$r2" ]]; then
      printf '%s\t%s\n' "$r1" "$r2"
      return 0
    fi
  done
  die "no trimmed reads for ${sample} in ${trimdir}
       expected ${sample}_R1.<fq.gz|fastq.gz|fq|fastq> and the matching R2.
       Run 01_trim.sh, or symlink existing trimmed reads under those names."
}

# depth_stats <bam> <region> <min_depth>
# Emits "mean<TAB>median<TAB>fraction_of_sites_at_or_above_min_depth" over the
# region, counting zero-coverage sites (samtools depth -a).
depth_stats() {
  samtools depth -a -r "$2" "$1" \
    | cut -f3 \
    | sort -n \
    | awk -v m="$3" '
        { v[NR] = $1; s += $1; if ($1 >= m) c++ }
        END {
          if (NR == 0) { print "0\t0\t0"; exit }
          med = (NR % 2) ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2
          printf "%.1f\t%.1f\t%.3f\n", s / NR, med, c / NR
        }'
}

# blast_span_bed <hits_tsv> [feature_name]
# Reads a blast tabular file already sorted by descending bitscore, with columns
# sseqid sstart send pident length bitscore, and emits one BED line spanning
# every HSP on the best-scoring subject. Taking the span rather than the single
# best HSP matters at high divergence, where a translated search often breaks
# one gene into several HSPs and the top hit alone would truncate it.
blast_span_bed() {
  local hits="$1" name="${2:-COX1}" contig
  [[ -s "$hits" ]] || die "no blast hits in ${hits}"
  contig=$(head -1 "$hits" | cut -f1)
  awk -F'\t' -v c="$contig" -v n="$name" 'BEGIN { OFS = "\t" }
      $1 == c {
        s = ($2 < $3 ? $2 : $3); e = ($2 < $3 ? $3 : $2)
        if (min == "" || s < min) min = s
        if (e > max) max = e
        if (strand == "") strand = ($2 < $3 ? "+" : "-")
      }
      END { print c, min - 1, max, n, 0, strand }' "$hits"
}

# region_from_bed <bed>
# Converts the single-interval BED written by 00_reference.sh into the
# 1-based "chrom:start-end" string that samtools and bcftools expect.
region_from_bed() {
  [[ -f "$1" ]] || die "no such file: $1"
  awk -F'\t' 'NR == 1 { printf "%s:%d-%d\n", $1, $2 + 1, $3 }' "$1"
}
