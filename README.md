# coi_from_ahe

Recovery of the COI barcode from the off-target fraction of anchored hybrid
enrichment (AHE) libraries, by two independent routes that are reconstructed
separately and then compared.

AHE libraries are enriched shotgun libraries, not purified ones. Mitochondrial
DNA is present at high copy number in the input extraction, so the off-target
reads normally carry enough mitochondrial coverage to reconstruct COI. Whether
they do in a given library is an empirical question, and `02_bait.sh` answers it
before any assembly is attempted.

## Routes

| Route | Script | Method | Fails when |
| --- | --- | --- | --- |
| A | `03_assemble.sh` | de novo assembly of the whole trimmed library, COX1 located by tblastn | coverage is too low or too uneven to assemble; a NUMT assembles preferentially |
| B | `04_consensus.sh` | reference-guided consensus from reads mapped to a related mitogenome | the reference is too divergent to recruit reads; true indels are missed by design |

The two routes fail for different reasons, which is what makes their agreement
informative. `05_validate.sh` quantifies that agreement and applies the checks
that neither route makes on its own.

## Installation

```bash
conda env create -f environment.yml
conda activate coi_from_ahe
conda env export --no-builds > environment.lock.yml

conda env create -f environment_mitofinder.yml    # optional, for 06 only
```

`environment.yml` carries no version pins: bioconda builds for a given platform
only cover recent releases, so a version chosen elsewhere can make the solve
fail. Let conda resolve, then freeze. `environment.lock.yml` is what a methods
section cites and what reproduces a run exactly.

## Inputs

All scripts are run from the project root. Reads go in `data/raw_data/`, and
`samples.tsv` lists them:

```
sample	r1	r2
<sample>	data/raw_data/<sample>_R1.fastq.gz	data/raw_data/<sample>_R2.fastq.gz
```

Use whole-library reads, raw or lightly trimmed. Do not use reads that a
target-capture pipeline has already filtered, HybPiper output in particular:
those keep only reads matching the nuclear target file, which is precisely the
fraction that does not contain the mitogenome. Aggressive quality filtering is
also counterproductive here, because the off-target mitochondrial fraction is
the scarce resource; `01_trim.sh` trims lightly for that reason.

To reuse trimmed reads produced elsewhere, symlink them into `results/01_trim/`
as `<sample>_R1.<ext>` and `<sample>_R2.<ext>` and skip `01_trim.sh`. The
extensions `.fq.gz`, `.fastq.gz`, `.fq` and `.fastq` are all accepted.

## Reference mitogenome

`00_reference.sh` takes an NCBI accession in `MITO_ACC` and derives the COX1
CDS, its interval, and a protein form used by the tblastn search in `03`.

Choose the closest available mitogenome, annotated or not. To list candidates,
filtering by length rather than by title, since many mitogenomes were deposited
under titles that do not say "complete genome":

```bash
esearch -db nuccore -query "<Taxon>[Organism] AND mitochondrion[filter] \
  AND 10000:20000[SLEN]" \
| efetch -format docsum \
| xtract -pattern DocumentSummary -element Caption,Organism,Slen,Title
```

If the record carries a COX1 annotation, it is used directly. If it does not,
which is common for records flagged UNVERIFIED, set `ANNOT_ACC` to any annotated
mitogenome and its COX1 protein is used to locate COX1 by tblastn. The reference
sequence and interval still come from `MITO_ACC`; the donor contributes a query
and nothing else.

```bash
MITO_ACC=<accession> bash 00_reference.sh
MITO_ACC=<accession> ANNOT_ACC=<annotated accession> bash 00_reference.sh
```

The script reports the recovered COX1 length. A value far from about 1530 bp
means the interval is wrong.

Reference divergence affects the pipeline unevenly. Route A is barely affected,
since it assembles the whole library and uses the reference only for a
translated search. Route B and the `02` diagnostics are reference dependent, so
with a distant reference expect a low apparent mitochondrial fraction and a
heavily masked consensus, and do not read either as a verdict on the library.

## Run order

```bash
conda activate coi_from_ahe
export THREADS=24 MEM_GB=96                  # match the machine

MITO_ACC=<accession> bash 00_reference.sh    # reference, COX1 CDS, index
bash 01_trim.sh                              # fastp
bash 02_bait.sh                              # mapping, diagnostics, baited reads
```

Stop here and read `results/02_bait/bait_stats.tsv` before continuing. Then:

```bash
bash 03_assemble.sh                          # route A, de novo
bash 04_consensus.sh                         # route B, reference guided
bash 05_validate.sh                          # checks, comparison, deliverables
```

`03_assemble.sh` is the long step. Run it under `tmux` or `nohup`.

Every script reads its settings from environment variables with defaults in a
block at the top, so any parameter can be overridden without editing the file,
for example `THREADS=32 MEM_GB=64 bash 03_assemble.sh`. Redirecting `OUTDIR`
this way is what makes a second pass possible without disturbing the first.

## Reading the diagnostics

`results/02_bait/bait_stats.tsv` reports the mitochondrial read fraction, the
duplication rate, and the depth profile across COX1.

| Median COX1 depth | Expect |
| --- | --- |
| above about 20x | both routes work |
| 10x to 20x | route A may fragment; route B carries more weight |
| below about 10x | a partial COI |

High duplication with low unique depth means the library has few unique
molecules, and more sequencing of that library will not help.

## Optional steps

**Second pass against a closer reference.** Once route A has produced a COI,
`00b_reference_from_coi.sh` builds a reference directory from that sequence so
`02` and `04` can be rerun against a conspecific target, into separate output
directories. Read the result correctly: on the first pass route B is independent
of route A and their agreement is evidence; on a second pass route B is derived
from route A, so their agreement measures read support and cannot detect a NUMT
that assembled cleanly.

**Full annotated mitogenomes.** `06_mitogenome.sh` runs MitoFinder to recover
the rest of the mitogenome, the other protein-coding genes, rRNAs and tRNAs,
from the same reads. It is off the COI path: MitoFinder wraps the same
assemblers used in `03`, so it does not produce a better COI, and the script
exists for the additional genes. It also aligns its own COX1 against the COI
from routes A and B as a third opinion.

## Outputs

```
ref/                          reference mitogenome, COX1 CDS (nt and aa), interval
results/01_trim/              trimmed reads, fastp reports
results/02_bait/              alignments, baited read pairs, bait_stats.tsv
results/03_denovo/            <sample>_COI_routeA.fasta, assemblies, blast hits
results/04_consensus/         <sample>_COI_routeB.fasta, consensus_stats.tsv
results/05_validate/          per-sample diagnostics, validation_summary.tsv
results/coi_final/            <sample>_COI.fasta, <sample>_COI_folmer.fasta
results/06_mitogenome/        <sample>_mitogenome.fasta and .gb
```

`validation_summary.tsv` is the table to check before using any sequence:

| Column | Expected | Meaning if not |
| --- | --- | --- |
| `frac_ge_10` | 1.000 | part of the COI rests on fewer than 10 reads |
| `competing_allele_sites` | 0 | co-assembled NUMT, two individuals in one library, or index hopping |
| `orf_verdict` | PASS | internal stop codons under genetic code 5, so a pseudogene or a frameshift |
| `identity_A_vs_B` | > 99.5 | the two routes disagree; resolve before use |
| `folmer_bp` | about 658 | the sequence does not span the standard barcode region and is not directly comparable with BOLD |
