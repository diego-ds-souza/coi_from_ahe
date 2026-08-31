#!/usr/bin/env bash
#
# Adapter and quality trimming of the raw AHE read pairs with fastp.
#
# Trimming is deliberately light. Aggressive quality trimming shortens reads
# and costs coverage on the off-target mitochondrial fraction, which is the
# scarce resource here. The fastp JSON reports are kept because their
# duplication estimate is the single most useful early indicator of whether the
# capture libraries have enough unique molecules to reconstruct COI.
#
# Run from the project root:
#   conda activate coi_from_ahe
#   bash 01_trim.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

# ------------------------------- settings ------------------------------------
SAMPLES="${SAMPLES:-samples.tsv}"              # sample sheet
OUTDIR="${OUTDIR:-results/01_trim}"            # output directory
MIN_LEN="${MIN_LEN:-50}"                       # discard shorter reads
MIN_QUAL="${MIN_QUAL:-20}"                     # per-base qualified threshold
MAX_UNQUAL="${MAX_UNQUAL:-30}"                 # max percent unqualified bases
THREADS="${THREADS:-24}"                        # fastp accepts at most 16
# -----------------------------------------------------------------------------

require_tools fastp
require_python
samples=$(read_samples "$SAMPLES")
mkdir -p "$OUTDIR"

# fastp rejects thread counts above 16, so a cluster-sized THREADS value set
# once for the whole pipeline does not have to be special cased by the caller.
FASTP_THREADS=$(( THREADS > 16 ? 16 : THREADS ))

n=$(printf '%s\n' "$samples" | wc -l | tr -d ' ')
i=0
while IFS=$'\t' read -r sample r1 r2; do
  i=$((i + 1))
  echo "[${i}/${n}] fastp ${sample}"
  require_files "$r1" "$r2"

  fastp \
    --in1 "$r1" --in2 "$r2" \
    --out1 "${OUTDIR}/${sample}_R1.fq.gz" \
    --out2 "${OUTDIR}/${sample}_R2.fq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred "$MIN_QUAL" \
    --unqualified_percent_limit "$MAX_UNQUAL" \
    --length_required "$MIN_LEN" \
    --thread "$FASTP_THREADS" \
    --json "${OUTDIR}/${sample}.fastp.json" \
    --html "${OUTDIR}/${sample}.fastp.html" \
    2> "${OUTDIR}/${sample}.fastp.log"
done <<< "$samples"

# Duplication rate is the number to read first: above roughly 60 percent the
# library has few unique molecules and de novo assembly of the off-target
# mitochondrial fraction will be depth limited whatever the raw read count.
echo
echo "duplication rate per sample (fastp estimate):"
for json in "${OUTDIR}"/*.fastp.json; do
  "$PY" - "$json" <<'PY'
import json, os, sys
with open(sys.argv[1]) as handle:
    report = json.load(handle)
name = os.path.basename(sys.argv[1]).replace(".fastp.json", "")
rate = report.get("duplication", {}).get("rate")
reads = report.get("summary", {}).get("after_filtering", {}).get("total_reads")
print(f"  {name}\t{'NA' if rate is None else format(rate, '.3f')}\t{reads} reads kept")
PY
done

echo "wrote ${OUTDIR}/<sample>_R{1,2}.fq.gz"
