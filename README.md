# GenomeMed Kidney Transplant Genome-Wide Mismatch Study

This repository computes donor–recipient genomic mismatches and tests their
association with transplant outcomes.

## Requirements

- Slurm for the chromosome-array mismatch workflow
- `bcftools`, `plink2`, Python 3, `awk`, and standard Unix tools
- R with the `survival` package for Cox regression

## 1. Compute donor–recipient mismatches

`1.compute_MM.sh` accepts either one indexed multi-chromosome VCF or a directory
of per-chromosome VCFs. The pairs file has two whitespace-delimited columns,
recipient followed by donor, without a header. VCF sample IDs must match these
IDs.

With one VCF:

```bash
sbatch --array=1-22 1.compute_MM.sh \
  --vcf /path/to/genotypes.vcf.gz \
  --pairs /path/to/pairs.txt \
  --outdir /path/to/output
```

With per-chromosome VCFs:

```bash
sbatch --array=1-22 1.compute_MM.sh \
  --vcf-dir /path/to/vcfs \
  --pairs /path/to/pairs.txt \
  --outdir /path/to/output
```

The explicit `--array=1-22` is important because the script's current Slurm
header defaults to chromosome 22 only. Outputs include chromosome-level
mismatch VCF/PGEN files under `OUTDIR/MM` and the retained pair list under
`OUTDIR/PAIRS`.

For donor ALT dosage `x` and recipient ALT dosage `y`, this workflow computes:

```text
abs(y - 1) * abs(x - y)
```

This counts copies of a donor allele that the recipient does not carry.

## 2. Cox association analysis

`3.Cox_light.R` independently calculates a **binary** mismatch from paired
recipient and donor dosage matrices, then fits one adjusted Cox model per
variant. It does not directly read the PGEN output from step 1.

After rounding recipient and donor dosages to 0, 1, or 2, the R script uses:

```text
min(abs((recipient - donor) * (recipient != 1)), 1)
```

Thus, unlike step 1's copy count, the Cox predictor is always 0/1 and represents
absence/presence of donor-derived genetic material not carried by the recipient.

### Required inputs

The recipient and donor dosage tables must have identical variant rows and
these five named columns:

```text
CHROM  POS  ID  REF  ALT
```

Remaining columns contain dosages for donor–recipient pair IDs. Those column
names must match the phenotype pair-ID column, which defaults to `ID`.

The phenotype table must contain:

- a pair ID (`ID` by default);
- an event indicator (`FAILURE` by default, coded 0/1);
- follow-up time (`TIME_FAILURE` by default); and
- all requested covariates.

The legacy 62-column phenotype layout is detected automatically when the file
has 62 columns but the expected field names are absent. The file must still
contain a header row. Control this with `--legacy-schema auto|yes|no`.

### Basic usage

```bash
Rscript 3.Cox_light.R \
  --pheno /path/to/pheno-epi.txt \
  --recipient-dosage /path/to/recipients.dosage.failure.1e-3.txt \
  --donor-dosage /path/to/donors.dosage.failure.1e-3.txt \
  --output /path/to/results/FAILURE.hla-epi.txt \
  --event FAILURE \
  --hla epi \
  --cores 4
```

Run `Rscript 3.Cox_light.R --help` for the complete option list.

### Another outcome

Defaults for event and follow-up columns are derived from `--event`:

```bash
Rscript 3.Cox_light.R \
  --pheno /path/to/pheno.txt \
  --recipient-dosage /path/to/recipients.dosage.abmr.1e-3.txt \
  --donor-dosage /path/to/donors.dosage.abmr.1e-3.txt \
  --output /path/to/results/ABMR.hla-epi.txt \
  --event ABMR \
  --hla epi
```

This uses `ABMR` as the event column and `TIME_ABMR` as follow-up. Override
them with `--event-col` and `--time-col` when the phenotype uses other names.

### Custom covariates

Pass model terms as a comma-separated list:

```bash
Rscript 3.Cox_light.R \
  --pheno /path/to/pheno.txt \
  --recipient-dosage /path/to/recipient.tsv \
  --donor-dosage /path/to/donor.tsv \
  --output /path/to/results.tsv \
  --event-col graft_failure \
  --time-col follow_up_days \
  --pair-id-col pair_id \
  --hla none \
  --covariates 'PC1_R,PC2_R,PC1_D,PC2_D,recipient_age,donor_age'
```

In `--hla epi` mode, rows missing the column selected by `--epi-col` are
removed. The default covariates preserve the original analysis, including
`tt(AGE_D)` with the transformation `AGE_D * log(time)`.

### Output

The tab-delimited output contains:

```text
CHROM  POS  ID  REF  ALT  COEF  SE.COEF  P  FREQ
```

- `COEF` is the Cox log hazard ratio for carrying the binary mismatch.
- `exp(COEF)` is the mismatch hazard ratio.
- `FREQ` is the proportion of donor–recipient pairs with the mismatch.
- Failed or invariant variant models are retained with `NA` statistics.

`FREQ` is a pair-level mismatch-carrier frequency, not a conventional effect
allele frequency. The mismatch can reflect introduced ALT in some pairs and
introduced REF in others, so these results should not be labeled as an
ordinary allele-dosage GWAS without additional methodological guidance.
