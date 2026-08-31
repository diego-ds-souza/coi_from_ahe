# coi_from_ahe (Linux server build)

Recovery of COI from the off-target fraction of anchored hybrid enrichment
(AHE) libraries, with two independent reconstruction routes and an explicit
validation step. This is the same pipeline as the macOS copy, with defaults set
for a machine that has cores and memory, plus one optional extra step.

AHE libraries are enriched shotgun libraries, not purified ones. Because
mitochondrial DNA is present at high copy number in the input extraction, the
off-target reads normally carry enough mitochondrial coverage to reconstruct
COI. Whether they do in a particular library is an empirical question, which
`02_bait.sh` answers before any assembly is attempted.

Two routes are run and compared:

| Route | Script | Method | Fails when |
| --- | --- | --- | --- |
| A | `03_assemble.sh` | MEGAHIT assembly of the whole trimmed library, COX1 located by tblastn | coverage is too low or too uneven to assemble; a NUMT assembles preferentially |
| B | `04_consensus.sh` | reference-guided consensus from reads mapped to a related mitogenome | the reference is too divergent to recruit reads; true indels are missed by design |

They fail for different reasons, which is why agreement between them is
informative. `05_validate.sh` quantifies that agreement and applies the checks
that catch the failure modes neither route detects on its own.

## What differs from the macOS copy

| | macOS | this build |
| --- | --- | --- |
| `METHOD` default | `spades` | `megahit` |
| `ASM_INPUT` default | `baited` | `full` |
| `THREADS` / `MEM_GB` defaults | 8 / 16 | 16 / 32 |
| MitoFinder | uninstallable (Python 2.7, no `osx-arm64` build) | available, in `06_mitogenome.sh` |

The important one is `ASM_INPUT=full`. Route B is reference guided by
construction, so route A should not be. Assembling the whole library rather
than the reads that a divergent reference managed to recruit makes the two
routes genuinely independent, which is the entire basis for treating their
agreement as evidence. On a laptop that was too expensive; on the server it is
the default.

## Installation

```bash
cd /path/to/coi_from_ahe_linux
conda env create -f environment.yml
conda activate coi_from_ahe
conda env export --no-builds > environment.lock.yml   # freeze what actually solved

conda env create -f environment_mitofinder.yml        # optional, for 06 only
```

`environment.yml` carries no version pins on purpose: pin after solving, not
before. `environment.lock.yml` is the file to cite in a methods section and to
reuse for an exact rerun, while `environment.yml` stays portable.

## Where MitoFinder fits, and where it does not

MitoFinder is a wrapper. It runs MEGAHIT, metaSPAdes or IDBA, identifies
mitochondrial contigs against a reference GenBank record, and annotates the
genes. The assembly is done by the same assemblers this pipeline calls directly,
so it cannot produce a better COI, and for a target this divergent its
nucleotide-level contig identification is arguably less sensitive than the
tblastn search in `03_assemble.sh`.

What it does add is the rest of the mitogenome: the other twelve protein-coding
genes, two rRNAs and 22 tRNAs, annotated, from reads you have already sequenced.
That is a second dataset for the cost of one more run, and a reasonable thing to
want from a phylogenomic program on wood-associated weevils. So it is here, as
`06_mitogenome.sh`, off the COI path and clearly optional.

It also cross-checks itself: the script aligns MitoFinder's COX1 against the COI
from routes A and B. Three independent reconstructions agreeing is a stronger
statement than two.

## If something fails

**`SyntaxError` in one of the `lib/*.py` helpers.** `python` on your PATH is
Python 2. The scripts now resolve a Python 3.8+ interpreter themselves and stop
with a clear message if none is found, so this should not recur, but you can
check with `python -V` and `python3 -V`, and override for one run with
`PYTHON=/path/to/python3 bash 00_reference.sh`.

**`curl: (56) ... unexpected eof while reading` during efetch.** Something on the
network is terminating TLS connections early. edirect usually retries past it,
but a download can end up truncated while still being a non-empty file.
`00_reference.sh` now retries up to three times and rejects any GenBank file
without a terminating `//` or any FASTA without a header, so a truncated
download fails loudly instead of silently. If it never completes, fetch the two
files elsewhere and copy them to `ref/<ACC>.gb` and `ref/<ACC>.fasta`; the
script reuses complete files and skips the download.

**"carries no COX1 annotation either; pick another donor".** Before this was
fixed, that message appeared whenever the extraction helper failed for any
reason, including a broken interpreter. It now appears only when the record
genuinely has no COX1 feature; a tooling failure aborts and says so instead.

## Which input to use

**Use the raw reads.** Not HybPiper output, and, for this dataset, not the
existing fastp output either.

HybPiper distributes reads to per-locus directories by mapping them against the
nuclear target file and keeps only what matches. The off-target fraction is what
it discards, and that fraction is where the mitogenome is.

The existing `data/fastp/` output is whole-library, but it was filtered hard:

| Sample | Reads in | Reads out | Retained | fastp duplication |
| --- | --- | --- | --- | --- |
| USNMENT01160338 | 33,030,796 | 13,819,187 | 42% | 0.55 |
| USNMENT01160339 | 41,168,424 | 14,046,013 | 34% | 0.62 |
| USNMENT01160347 | 28,524,486 | 4,175,986 | 15% | 0.80 |

Those settings are reasonable for HybPiper, where on-target coverage is
plentiful and stringency costs nothing. Here the off-target mitochondrial
fraction is the scarce resource, and discarding 58 to 85 percent of the library
discards most of it. `01_trim.sh` trims lightly by comparison, which is the
whole reason it exists as a separate step.

USNMENT01160347 is the one to watch: 15 percent retention with an 80 percent
duplication estimate points to a low-complexity library. If any sample fails to
yield COI, expect it to be that one, and note that more sequencing of the same
library would not fix it.

If you do want to reuse trimmed reads from elsewhere, symlink them into
`results/01_trim/` as `<sample>_R1.<ext>` and `<sample>_R2.<ext>`; `.fq.gz`,
`.fastq.gz`, `.fq` and `.fastq` are all accepted.

## Inputs

Everything runs from the project root, the directory holding these scripts.
The layout matches the macOS copy, so `samples.tsv` needs no editing once the
data are in place:

```
data/raw_data/<sample>_R{1,2}.fastq.gz      used by samples.tsv
data/fastp/<sample>/                        the earlier HybPiper-oriented run, unused
```

Only `data/raw_data/` is needed on the server, about 6.3 GB:

```bash
rsync -avP ~/Downloads/coi_from_ahe/data/raw_data/ user@server:/path/to/coi_from_ahe_linux/data/raw_data/
cd /path/to/coi_from_ahe_linux
md5sum -c data/raw_data/*.md5      # the .md5 files came with the reads; use them
```

## The reference mitogenome

Four cossonine mitogenomes exist in GenBank:

| Accession | Taxon | Length | Flag |
| --- | --- | --- | --- |
| JN163960 | *Brachytemnus porcatus* | 10,666 | partial genome |
| MH404131 | Cossoninae sp. 1 CG319 | 15,215 | UNVERIFIED |
| MH404140 | *Allopentarthrum elumbe* | 16,449 | UNVERIFIED |
| OQ716311 | *Amaurorhinus* sp. SCI190 | 15,205 | partial genome |

**The `NC_` prefix does not matter.** `NC_` marks RefSeq, which is a curated copy
of a GenBank record, not a quality tier the pipeline depends on. Three things
matter instead: does the record contain COX1, is COX1 annotated, and how close
is the taxon. All four of these are Cossoninae and therefore far closer to your
samples than any Scolytinae mitogenome, which makes them worth the extra care.

`UNVERIFIED` means GenBank staff could not verify the sequence and/or its
annotation. Such records are excluded from NCBI's BLAST databases and frequently
carry no feature table. That used to make them unusable here, because the COX1
interval was read from the annotation. It no longer does: `00_reference.sh`
falls back to locating COX1 by tblastn using the COX1 protein of any annotated
mitogenome you name in `ANNOT_ACC`.

```bash
# Is the candidate annotated, and does it contain COX1?
for acc in MH404140 MH404131 OQ716311 JN163960; do
  printf '%s\tCDS features: ' "$acc"
  efetch -db nuccore -id "$acc" -format gb | grep -c '^     CDS ' || true
done

# Annotated: use it directly.
MITO_ACC=MH404140 bash 00_reference.sh

# Unannotated: name an annotated donor for the protein query only.
MITO_ACC=MH404140 ANNOT_ACC=NC_059702 REFDIR=ref bash 00_reference.sh
```

The donor contributes a protein query, nothing else. The reference sequence,
the COX1 interval and the protein written to `ref/cox1.faa` all come from
`MITO_ACC`. Either way `00_reference.sh` prints the recovered COX1 length; a
value far from about 1530 bp means the interval is wrong, and
`ref/cox1_hits.tsv` shows why.

Two cautions on choosing among them. `JN163960` at 10,666 bp is partial and may
not span COX1 at all, so check it before relying on it. And confirm tribe
placements against a current catalogue rather than inferring them from genus
names: the Wikipedia Pentarthrini genus list includes *Microtrupis*, one of your
candidate genera for USNMENT01160338, but does not list *Allopentarthrum*.

Building a reference is cheap, so if two candidates look comparable, build both
into separate `REFDIR`s and compare what `02_bait.sh` recruits.

### The fallback for NC_059702

`MITO_ACC=NC_059702` (*Euwallacea fornicatus*) remains a valid choice and is the
right `ANNOT_ACC` donor regardless, since it is annotated. What it is and is
not: It is a genuine, complete, well-formed beetle
mitogenome: 15,745 bp, 73.2 percent AT, two ambiguous bases, with the Folmer
barcode window intact at positions 1348 to 2056 and a COX1 that translates
cleanly under genetic code 5. It is also Scolytinae (Xyleborini), a different
subfamily from your Cossoninae. Cossoninae sits near Scolytinae and Platypodinae
in the higher Curculionidae, so this is one of the better-placed options in the
absence of a cossonine mitogenome, but it is still a subfamily away.

The consequence is uneven across the pipeline:

| Step | Uses the reference for | Effect at this distance |
| --- | --- | --- |
| `03_assemble.sh` route A | recognizing COX1 in contigs, by tblastn | negligible; COX1 protein identity across weevil subfamilies stays high |
| `02_bait.sh` | recruiting reads by nucleotide mapping | substantial; expect a low apparent mitochondrial fraction |
| `04_consensus.sh` route B | the consensus scaffold | substantial; expect heavy masking |

So run it, and read `02` and `04` as pessimistic rather than as verdicts on the
libraries. Route A is the result; route B on this pass is a weak independent
check, not a second opinion of equal weight.

Once route A has produced a COI, `00b_reference_from_coi.sh` rebuilds a
reference directory from that sequence so `02` and `04` can be rerun against a
conspecific target. Note carefully what changes: on the first pass route B is
independent of route A and their agreement is evidence; on the second pass route
B is derived from route A, so their agreement measures read support and cannot
detect a NUMT that assembled cleanly. Both passes are worth having, and
`06_mitogenome.sh` is the genuinely independent reconstruction on this platform.

If the server has no outbound network, `00_reference.sh` skips the download when
the files already exist. Fetch them on a machine that does have network and copy
them to `ref/NC_059702.gb` and `ref/NC_059702.fasta`; the script picks up from
there. A FASTA alone is not enough, because the COX1 interval is read from the
GenBank annotation.

## Run order

```bash
cd /path/to/coi_from_ahe_linux
conda activate coi_from_ahe
export THREADS=24 MEM_GB=96               # raise to whatever the server allows

MITO_ACC=MH404140 ANNOT_ACC=NC_059702 \
  bash 00_reference.sh                    # reference mitogenome, COX1 CDS, index
bash 01_trim.sh                           # fastp
bash 02_bait.sh                           # mapping, diagnostics, baited reads
```

Stop here and read `results/02_bait/bait_stats.tsv` before going on. Then:

```bash
bash 03_assemble.sh                       # route A, MEGAHIT on the full library
bash 04_consensus.sh                      # route B, reference guided
bash 05_validate.sh                       # checks, comparison, deliverables
```

Optional, and independent of everything above:

```bash
conda activate coi_mitofinder
bash 06_mitogenome.sh                     # full annotated mitogenomes
```

Optional second pass, once route A has produced a COI, for clean read support
and a well-behaved route B. Separate output directories keep the first pass
intact so the two can be compared:

```bash
conda activate coi_from_ahe
COI_FASTA=results/coi_final/USNMENT01160338_COI.fasta \
  REFDIR=ref_from_338 bash 00b_reference_from_coi.sh

printf 'sample\tr1\tr2\nUSNMENT01160338\tdata/raw_data/USNMENT01160338_R1.fastq.gz\tdata/raw_data/USNMENT01160338_R2.fastq.gz\n' \
  > sample_338.tsv

SAMPLES=sample_338.tsv REFDIR=ref_from_338 OUTDIR=results/02_bait_pass2 bash 02_bait.sh
SAMPLES=sample_338.tsv REFDIR=ref_from_338 BAITDIR=results/02_bait_pass2 \
  OUTDIR=results/04_consensus_pass2 bash 04_consensus.sh
```

`03_assemble.sh` is the long step at 15 to 20 million read pairs per sample.
Run the whole thing detached rather than watching it:

```bash
tmux new -s coi
# or: nohup bash 03_assemble.sh > 03.log 2>&1 &
```

Rough expectations at 15 to 20 million read pairs per sample: `01_trim.sh` and
`02_bait.sh` are tens of minutes each, `03_assemble.sh` with MEGAHIT on the full
library is the long step, and `04` and `05` are minutes. Run the three samples
in one go rather than babysitting them; the scripts loop over `samples.tsv`.

## Stop and read the numbers after step 02

`results/02_bait/bait_stats.tsv` reports, per sample, the mitochondrial read
fraction, the duplication rate and the depth profile across COX1.

- Median COX1 depth above roughly 20x: both routes will work.
- Median depth 10x to 20x: route A may fragment; route B carries more weight.
- Median depth below 10x, or `cox1_frac_ge_10` well under 1: expect a partial
  COI. Consider sequencing more of the existing libraries rather than forcing an
  assembly, unless duplication is already high, in which case more sequencing of
  that library will not help.

## Outputs

```
ref/                          reference mitogenome, COX1 CDS (nt and aa), interval
results/01_trim/              trimmed reads, fastp reports
results/02_bait/              alignments, baited read pairs, bait_stats.tsv
results/03_denovo/            <sample>_COI_routeA.fasta, assemblies, blast hits
results/04_consensus/         <sample>_COI_routeB.fasta, consensus_stats.tsv
results/05_validate/          per-sample diagnostics, validation_summary.tsv
results/coi_final/            <sample>_COI.fasta, <sample>_COI_folmer.fasta
results/06_mitogenome/        <sample>_mitogenome.fasta and .gb, mitogenome_stats.tsv
```

`validation_summary.tsv` is the table to check before using anything:

| Column | Expected | Meaning if not |
| --- | --- | --- |
| `frac_ge_10` | 1.000 | part of the COI rests on fewer than 10 reads |
| `competing_allele_sites` | 0 | co-assembled NUMT, two individuals in one library, or index hopping from a neighboring well |
| `orf_verdict` | PASS | internal stop codons under genetic code 5, so a pseudogene or a frameshift |
| `identity_A_vs_B` | > 99.5 | the two routes disagree; resolve before use |
| `folmer_bp` | about 658 | the sequence does not span the standard barcode region and is not directly comparable with BOLD |

## Caveats worth stating in a manuscript

- **NUMTs.** The ORF check catches frameshifted and stop-containing pseudogenes.
  It does not catch an intact NUMT. Concordance between routes plus uniform read
  support is the practical evidence, not proof.
- **Index hopping.** Patterned flow cells hop indices between libraries on the
  same lane. With three samples this is easy to overlook: check that no sample's
  COI appears as another's minor allele.
- **Coverage reporting.** Report per-sample mean depth and the masking threshold
  alongside the sequences. Reviewers now expect this for bycatch-derived data.
- **Analysis.** COI is a single non-recombining linkage group. Keep it as its own
  partition rather than concatenating it with nuclear AHE loci, and examine
  mitonuclear discordance explicitly rather than averaging it away.
- **Reproducibility.** Record `MITO_ACC`, any overridden settings, and
  `environment.lock.yml`. Those are the only inputs not captured by the scripts.
