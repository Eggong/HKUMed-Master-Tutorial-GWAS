---
layout: default
title: "Introduction to genome-wide association studies"
---

# Introduction to genome-wide association studies

### Genotype quality control, population structure, and association testing

**Suggested duration:** 1.5 - 2 hours 
**Software:** PLINK 1.9, R, and a Bash terminal

## Contents
{: .no_toc }

* Contents
{:toc}

## 1 Introduction

A genome-wide association study, or GWAS, tests whether genetic
variants are associated with a phenotype.

For a binary phenotype, we might ask whether a genetic variant is
associated with being a case rather than a control. For a quantitative
phenotype, we might investigate an association with a measurement such
as height.

Before testing associations, we need to examine the quality of the
genotype data and consider the design of the study.

In this practical, we will:

1. Inspect genotype and sample information.
2. Identify potential quality-control problems.
3. Apply documented filtering decisions.
4. Examine population structure.
5. Run an association analysis.
6. Interpret the resulting tables and plots.

The scientific background comes from Marees et al. (2018).
The accompanying repository provides the original practical materials.
See [References](#9-references).

<div class="warning" markdown="1">

### 1.1 Learning outcomes

By the end of this practical, you should be able to:

- Explain the purpose of the main GWAS quality-control steps.
- Recognize PLINK genotype, sample, and variant files.
- Distinguish diagnostic flags from justified exclusions.
- Explain why population structure matters.
- Run and interpret an additive logistic regression.
- Produce a reproducible record of your analysis.

## 2 Preparing the analysis

### 2.1 Software

The classroom environment should already contain:

- PLINK **1.9**, available as `plink`.
- R and `Rscript`.
- A Bash terminal.
- A text editor.

Check the installed programs:

```bash
./plink
./plink --version
Rscript
Rscript --version
```

Please record the output in your analysis notebook.

<div class="note" markdown="1">

**Which version of PLINK?**

The commands in this practical are written for PLINK 1.9 (https://www.cog-genomics.org/plink/1.9/).

Do not replace `plink` with `plink2` without adapting and checking the workflow.

</div>


### 2.2 Working directory

Open a terminal and move into `Course_Materials`.

Run all shell commands in this practical from that directory:

```bash
cd /path/to/Course_Materials

mkdir -p Results scripts

pwd
ls -la
```

Replace `/path/to/Course_Materials` with the real location.

Save shell scripts in `scripts/` and outputs in `Results/`.

### 2.3 Keep an analysis notebook

For each step, suggest to record:

- Input prefix.
- Command or script.
- Output prefix.
- Number of samples.
- Number of variants.
- Warnings.
- Reason for any exclusion.

Do not overwrite the original data!!!

## 3 Understanding PLINK files

### 3.1 The binary dataset

A PLINK binary dataset has three files sharing the same prefix:

| File | Contents |
|---|---|
| `study.bed` | Binary genotype data |
| `study.bim` | Variant information |
| `study.fam` | Sample information and phenotype |

The `.bed` file is binary and should not be opened as ordinary text.

Inspect the other files:

```bash
head Data/raw/study.bim
head Data/raw/study.fam
```

PLINK reads the three files together using their shared prefix:

./plink --bfile Data/raw/study


### 3.2 The sample information file

The six columns of a `.fam` file are:

| Column | Description |
|---|---|
| 1 | Family identifier: FID |
| 2 | Individual identifier: IID |
| 3 | Paternal identifier |
| 4 | Maternal identifier |
| 5 | Recorded sex code |
| 6 | Phenotype |

For the binary phenotype in this practical:

- `1` means control.
- `2` means case.
- `0` and `-9` indicate missing phenotype values.

These are PLINK input codes, not the mathematical coding used when
writing a logistic regression equation.

<div class="exercise" markdown="1">

**EXERCISE 1 — Inspect the input**

Find:

1. The number of samples.
2. The number of variants.
3. The counts of each phenotype code.

Complete the commands:

```bash
wc -l <<SAMPLE_FILE>>
wc -l <<VARIANT_FILE>>

awk '{print $6}' <<SAMPLE_FILE>> | sort | uniq -c
```

Write the results in your notebook.

</div>

<details markdown="1">
<summary>Hint</summary>

Each line of the `.fam` file represents a sample.

Each line of the `.bim` file represents a variant.

</details>

<details markdown="1">
<summary>Answer</summary>

```bash
wc -l Data/raw/study.fam
wc -l Data/raw/study.bim

awk '{print $6}' Data/raw/study.fam | sort | uniq -c
```

The numerical answers depend on the prepared teaching dataset.

</details>

## 4 Genotype quality control

Quality control involves both measurements and decisions.

For example, we can measure missingness automatically, but deciding
whether to exclude a sample may require considering several pieces
of evidence.

The thresholds below are classroom choices. They are not a universal
recipe for every GWAS.

### 4.1 Missing genotype calls

Missingness can be assessed at two levels:

- **Sample-level missingness:** the proportion of variants with missing genotype calls for each sample.
- **Variant-level missingness:** the proportion of samples with missing genotype calls for each variant.

Generate both reports:

```bash
./plink \
  --bfile Data/raw/study \
  --missing \
  --out Results/missingness
```

The main outputs are:

| Output | Contents |
|---|---|
| `missingness.imiss` | Missingness per sample |
| `missingness.lmiss` | Missingness per variant |
| `missingness.log` | Record of the PLINK run |

Inspect them:

```bash
head Results/missingness.imiss
head Results/missingness.lmiss
```

The `F_MISS` column contains the fraction of missing calls.

<div class="exercise" markdown="1">

**EXERCISE 2 — Apply missingness filters**

Create `scripts/01_missingness.sh`.

We will use a stringent threshold for illustration, however in real-world cases, you may try multiple threshold settings.

First exclude variants with more than 2% missing calls.
Then exclude samples with more than 2% missing calls.

Complete:

```bash
#!/usr/bin/env bash
set -euo pipefail

plink \
  --bfile Data/raw/study \
  --geno <<THRESHOLD>> \
  --make-bed \
  --out Results/01_variant_missingness

plink \
  --bfile Results/01_variant_missingness \
  --mind <<THRESHOLD>> \
  --make-bed \
  --out Results/02_missingness
```

Run the script and record the number of retained samples and variants.

</div>

<details markdown="1">
<summary>Hint</summary>

PLINK expects a proportion: 2% is written as `0.02`.

`--geno` filters variants; `--mind` filters samples.

</details>

<details markdown="1">
<summary>Answer</summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

plink \
  --bfile Data/raw/study \
  --geno 0.02 \
  --make-bed \
  --out Results/01_variant_missingness

plink \
  --bfile Results/01_variant_missingness \
  --mind 0.02 \
  --make-bed \
  --out Results/02_missingness
```

Run:

```bash
bash scripts/01_missingness.sh
```

Check:

```bash
wc -l Results/02_missingness.fam
wc -l Results/02_missingness.bim
```

</details>

**Discuss:** Could changing the order of sample and variant filtering
change the retained dataset? Explain your reasoning.

### 4.2 Minor allele frequency

Minor allele frequency, or MAF, describes the frequency of the less
common allele at a biallelic variant.

Calculate allele frequencies:

```bash
plink \
  --bfile Results/02_missingness \
  --freq \
  --out Results/allele_frequency
```

Inspect the report:

```bash
head Results/allele_frequency.frq
```

For this small teaching dataset, retain variants with MAF of at least 5%:

```bash
plink \
  --bfile Results/02_missingness \
  --maf 0.05 \
  --make-bed \
  --out Results/03_maf
```

<div class="exercise" markdown="1">

**EXERCISE 3 — Think about allele counts**

Suppose a variant has MAF = 0.01 in 200 completely genotyped
diploid individuals.

Approximately how many copies of the minor allele are present?

Why could an association estimate based on very few minor alleles
be difficult to interpret?

</div>

<details markdown="1">
<summary>Answer</summary>

There are 400 allele copies:

```text
2 × 200 = 400
```

The approximate minor allele count is:

```text
400 × 0.01 = 4
```

Ask whether four allele copies provide enough information for a stable
estimate, particularly after dividing individuals into cases and controls.

</details>

### 4.3 Hardy–Weinberg equilibrium

For a biallelic variant with allele frequencies \(p\) and \(q\),
Hardy–Weinberg equilibrium gives expected genotype frequencies:

```text
p², 2pq, q²
```

Generate a diagnostic report:

```bash
plink \
  --bfile Results/03_maf \
  --hardy \
  --out Results/hwe
```

For the classroom filtering exercise:

```bash
plink \
  --bfile Results/03_maf \
  --hwe 1e-6 \
  --make-bed \
  --out Results/04_variant_qc
```

For a recognized binary phenotype, PLINK 1.9 normally uses controls for
this filter and excludes nonfounders unless instructed otherwise.

Check the phenotype coding and read the log before interpreting the
result.

**Discuss:**

- Why examine controls rather than combining cases and controls?
- Does a small HWE p-value prove that a variant is technically unreliable?
- Why should filtering decisions be specified before examining GWAS hits?

### 4.4 Sample-level review

We will now examine three additional areas:

1. Sex-chromosome consistency.
2. Heterozygosity.
3. Relatedness.

For these steps, use the evidence and review notes supplied by the
instructor alongside your own diagnostic output.


#### 4.4.1 Sex-chromosome consistency

Check for sex discrepancy.

Subjects who were a priori determined as females are expected to have an X-chromosome inbreeding coefficient, \(F\), of less than 0.2. Subjects who were a priori determined as males are expected to have an \(F\) value greater than 0.8. This \(F\) value is based on the X-chromosome inbreeding, or homozygosity, estimate.

Subjects who do not meet these criteria are flagged as `PROBLEM` by PLINK.

```bash
plink \
  --bfile Results/04_variant_qc \
  --check-sex \
  --out Results/sex_check
```

The main output is:

| Output | Contents |
|---|---|
| `Results/sex_check.sexcheck` | Sex-check results for each sample |
| `Results/sex_check.log` | Record of the PLINK run |

Inspect the output:

```bash
head Results/sex_check.sexcheck
```

<div class="note" markdown="1">

**Important**

A sex-discrepancy flag requires investigation. It does not automatically justify removing an individual.

In this classroom workflow, individuals flagged as `PROBLEM` are removed to generate the next binary PLINK fileset, `Results/05_sexcheck`.

</div>

<div class="exercise" markdown="1">

**EXERCISE 4 — Remove individuals with sex discrepancy**

Generate a list of individuals with the `PROBLEM` status:

```bash
grep "PROBLEM" Results/sex_check.sexcheck | awk '{print $1, $2}' > Results/sex_discrepancy.txt
```

Remove those individuals and create a new binary PLINK fileset:

```bash
plink \
  --bfile Results/04_variant_qc \
  --remove Results/sex_discrepancy.txt \
  --make-bed \
  --out Results/05_sexcheck
```

This command creates:

| Output | Contents |
|---|---|
| `Results/05_sexcheck.bed` | Filtered binary genotype data |
| `Results/05_sexcheck.bim` | Filtered variant information |
| `Results/05_sexcheck.fam` | Filtered sample information |
| `Results/05_sexcheck.log` | Record of the PLINK run |

After running the command, record the number of retained samples and variants.

</div>

Check the resulting fileset:

```bash
wc -l Results/05_sexcheck.fam
wc -l Results/05_sexcheck.bim
```


#### 4.4.2 LD pruning

Nearby variants can carry correlated genotype information.

We will create a reduced marker set for heterozygosity, relatedness,
and population-structure calculations:

```bash
plink \
  --bfile Results/04_variant_qc \
  --autosome \
  --indep-pairwise 50 5 0.2 \
  --out Results/qc_pruning
```

The arguments specify:

- A window of 50 variants.
- A step of 5 variants.
- An \(r^2\) threshold of 0.2.

The file `qc_pruning.prune.in` contains retained marker identifiers.

<div class="note" markdown="1">

**Pruning is not the final association filter**

The reduced marker set is used for diagnostic calculations.

We will still test the retained full set of autosomal variants in the
association analysis.

</div>

#### 4.4.3 Heterozygosity

Generate a heterozygosity report:

```bash
plink \
  --bfile Results/04_variant_qc \
  --autosome \
  --extract Results/qc_pruning.prune.in \
  --het \
  --out Results/heterozygosity
```

The report includes observed and expected homozygous genotype counts
and an inbreeding-coefficient estimate.

In R, calculate the observed heterozygosity fraction:

```r
het <- read.table(
  "Results/heterozygosity.het",
  header = TRUE,
  check.names = FALSE
)

het$heterozygosity <- (
  het[["N(NM)"]] - het[["O(HOM)"]]
) / het[["N(NM)"]]

pdf("Results/heterozygosity.pdf")
hist(
  het$heterozygosity,
  breaks = 20,
  main = "Observed heterozygosity",
  xlab = "Fraction of heterozygous genotypes"
)
dev.off()
```

Compare unusual observations with the instructor’s review notes.

#### 4.4.4 Relatedness

Estimate identity by descent using the pruned markers:

```bash
plink \
  --bfile Results/04_variant_qc \
  --autosome \
  --extract Results/qc_pruning.prune.in \
  --genome \
  --out Results/relatedness
```

Inspect:

```bash
head Results/relatedness.genome
```

The `PI_HAT` column summarizes estimated sharing identical by descent.

Interpretation depends on the marker set and allele-frequency estimates,
so do not treat the output as an infallible family classification.

<div class="exercise" markdown="1">

**EXERCISE 4 — Review a sample-exclusion decision**

Choose one entry in `Data/review/review_notes.tsv`.

Explain:

1. What evidence triggered the review?
2. Whether other QC measurements support the decision.
3. Why the proposed action is appropriate.
4. What additional information would help.

For a related pair, explain why one sample was retained rather than
the other.

</div>

<details markdown="1">
<summary>Answer</summary>

There is no single numerical answer.

A good response distinguishes a diagnostic flag from an exclusion
decision, considers multiple measurements, and documents uncertainty.

For related samples, check whether the decision also resolves other
relationships in the group, rather than considering each pair in
isolation.

</details>

### 4.5 Apply reviewed sample exclusions

Use the reviewed exclusion list:

```bash
if [ -s Data/review/qc_exclusions.txt ]; then
  plink \
    --bfile Results/04_variant_qc \
    --remove Data/review/qc_exclusions.txt \
    --make-bed \
    --out Results/05_sample_qc
else
  plink \
    --bfile Results/04_variant_qc \
    --make-bed \
    --out Results/05_sample_qc
fi
```

The conditional handles an empty exclusion list.

Create a QC table in your notebook:

| Stage | Samples retained | Variants retained | Reason |
|---|---:|---:|---|
| Input | | | |
| Missingness | | | |
| MAF | | | |
| HWE | | | |
| Sample review | | | |

## 5 Population structure

### 5.1 Why does population structure matter?

Consider a study in which cases and controls were recruited from
populations with different allele frequencies.

A variant could appear associated with the phenotype because of those
population differences rather than an effect on the phenotype.

We will investigate structure in two ways:

1. Review the reference-based analysis.
2. Calculate covariates within the retained study sample.

### 5.2 Reference-based analysis

The second module of the original repository uses 1000 Genomes
reference data.

Your instructor will provide the reference-based plot and the
corresponding analysis notes, or guide you through
`2_Main_script_MDS.txt`.

Before combining reference and study data, check:

- Genome build.
- Variant positions and identifiers.
- Allele compatibility.
- Sample identifiers.
- The marker set used for the analysis.

Do not resolve merge failures by indiscriminately flipping alleles.

<div class="exercise" markdown="1">

**EXERCISE 5 — Interpret the reference plot**

Using the reference-based plot, answer:

1. Where do study samples lie relative to the reference samples?
2. Are any samples unusual?
3. Does the plot alone justify excluding them?
4. How would exclusion change the population addressed by the study?

Compare your interpretation with `review_notes.tsv`.

</div>

<details markdown="1">
<summary>Answer</summary>

A strong answer describes the observed pattern before proposing a
decision.

Exclusions should follow the study’s analysis plan and supporting
evidence. A sample should not be labelled poor quality merely because
it has a different ancestry profile.

Reference individuals are used here to interpret structure; they must
not be added as extra controls.

</details>

### 5.3 Define the analysis sample

Apply only the reviewed population exclusions and restrict this
introductory association analysis to autosomes:

```bash
if [ -s Data/review/population_exclusions.txt ]; then
  plink \
    --bfile Results/05_sample_qc \
    --remove Data/review/population_exclusions.txt \
    --autosome \
    --make-bed \
    --out Results/06_analysis
else
  plink \
    --bfile Results/05_sample_qc \
    --autosome \
    --make-bed \
    --out Results/06_analysis
fi
```

This is a study-only dataset.

### 5.4 Calculate MDS dimensions

We will use multidimensional scaling, or MDS, to summarize structure
within the retained study samples.

Recalculate the pruned marker set:

```bash
plink \
  --bfile Results/06_analysis \
  --indep-pairwise 50 5 0.2 \
  --out Results/structure_pruning
```

Calculate ten dimensions:

```bash
plink \
  --bfile Results/06_analysis \
  --extract Results/structure_pruning.prune.in \
  --cluster \
  --mds-plot 10 \
  --out Results/structure
```

MDS and principal component analysis are related approaches to
dimension reduction, but they are not identical procedures.

### 5.5 Plot the dimensions and prepare covariates

Create `scripts/02_structure.R`:

```r
mds <- read.table(
  "Results/structure.mds",
  header = TRUE,
  colClasses = c(FID = "character", IID = "character")
)

pdf("Results/structure.pdf")
plot(
  mds$C1,
  mds$C2,
  pch = 19,
  xlab = "MDS dimension 1",
  ylab = "MDS dimension 2",
  main = "Structure within the study sample"
)
dev.off()

covariates <- mds[
  , c("FID", "IID", paste0("C", 1:10))
]

stopifnot(!anyNA(covariates))

write.table(
  covariates,
  "Results/covariates.tsv",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
```

Run:

```bash
Rscript scripts/02_structure.R
```

Inspect `Results/structure.pdf`.

<div class="note" markdown="1">

**How many dimensions?**

Ten dimensions are a classroom choice, not a universal recommendation.

As a sensitivity exercise, compare models using two, five, and ten
dimensions. Consider sample size and model stability as well as changes
in the results.

</div>

## 6 Association testing

### 6.1 The model

For the binary teaching phenotype, we will use logistic regression.

Each variant is tested using an additive allele-count model:

```text
log odds of being a case =
    intercept
  + genotype effect × allele count
  + covariate effects
```

The allele count is 0, 1, or 2 copies of the tested allele.

The exponentiated genotype coefficient is the odds ratio.

An odds ratio above 1 indicates higher estimated odds per additional
copy of that allele, conditional on the model covariates.

### 6.2 Unadjusted analysis

Run:

```bash
plink \
  --bfile Results/06_analysis \
  --logistic hide-covar \
  --ci 0.95 \
  --out Results/gwas_unadjusted
```

This model does not include the MDS covariates.

### 6.3 Structure-adjusted analysis

<div class="exercise" markdown="1">

**EXERCISE 6 — Add population-structure covariates**

Create `scripts/03_association.sh`.

Use the study dataset and the covariates generated in Section 5.

Complete:

```bash
#!/usr/bin/env bash
set -euo pipefail

plink \
  --bfile <<DATASET_PREFIX>> \
  --covar <<COVARIATE_FILE>> \
  --covar-name C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 \
  --logistic hide-covar \
  --ci 0.95 \
  --out Results/gwas_adjusted
```

</div>

<details markdown="1">
<summary>Hint</summary>

Use the full retained autosomal dataset, not the pruned marker list.

PLINK matches covariates to genotypes using `FID` and `IID`.

</details>

<details markdown="1">
<summary>Answer</summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

plink \
  --bfile Results/06_analysis \
  --covar Results/covariates.tsv \
  --covar-name C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 \
  --logistic hide-covar \
  --ci 0.95 \
  --out Results/gwas_adjusted
```

Run:

```bash
bash scripts/03_association.sh
```

`hide-covar` hides covariate-specific result rows. It does not remove
covariates from the model.

</details>

### 6.4 Inspect the results

```bash
head Results/gwas_adjusted.assoc.logistic
less Results/gwas_adjusted.log
```

Press `q` to leave `less`.

Important result fields include:

| Field | Meaning |
|---|---|
| `CHR` | Chromosome |
| `SNP` | Variant identifier |
| `BP` | Genomic position |
| `A1` | Tested allele |
| `TEST` | Model term; `ADD` is the additive genotype term |
| `NMISS` | Number of observations contributing to the test |
| `OR` | Odds ratio |
| `P` | Association p-value |

Before interpreting a result:

- Check warnings.
- Check missing or failed tests.
- Compare analyzed sample counts.
- Confirm the tested allele.
- Check the confidence interval.

<div class="exercise" markdown="1">

**EXERCISE 7 — Interpret an odds ratio**

A variant has an estimated odds ratio of 1.20.

Does this mean:

1. A 20-percentage-point increase in disease probability?
2. A 20% increase in odds per additional tested allele?
3. Proof that the variant causes the phenotype?

Explain your answer.

</div>

<details markdown="1">
<summary>Answer</summary>

The correct interpretation is the second, conditional on the model
covariates.

Odds are not the same as probabilities.

An association estimate alone does not establish causality.

</details>

## 7 Visualizing and interpreting results

### 7.1 Multiple testing

A GWAS performs many statistical tests.

We will compare:

- A dataset-specific Bonferroni threshold: `0.05 / number of tests`.
- The conventional common-variant GWAS reference line: `5e-8`.

The second is a reference convention, not a threshold derived from the
number of variants in this small teaching dataset.

### 7.2 Manhattan and Q–Q plots

The following R code uses base R and does not require additional
plotting packages.

Create `scripts/04_plots.R`:

```r
plot_gwas <- function(prefix) {

  result <- read.table(
    paste0(prefix, ".assoc.logistic"),
    header = TRUE,
    stringsAsFactors = FALSE
  )

  additive <- subset(result, TEST == "ADD")

  valid <- with(
    additive,
    is.finite(P) & P > 0 & P <= 1
  )

  message("Additive tests: ", nrow(additive))
  message("Tests omitted from plots: ", sum(!valid))

  gwas <- additive[valid, ]

  if (nrow(gwas) == 0) {
    stop("No valid p-values available.")
  }

  gwas <- gwas[order(gwas$CHR, gwas$BP), ]

  chromosomes <- sort(unique(gwas$CHR))

  offset <- 0
  ticks <- numeric(length(chromosomes))
  gwas$cumulative_position <- NA_real_

  for (j in seq_along(chromosomes)) {
    chr <- chromosomes[j]
    selected <- gwas$CHR == chr

    gwas$cumulative_position[selected] <-
      gwas$BP[selected] + offset

    ticks[j] <- mean(
      range(gwas$cumulative_position[selected])
    )

    offset <- offset + max(gwas$BP[selected])
  }

  colours <- c("#176b87", "#75a9ba")[
    1 + match(gwas$CHR, chromosomes) %% 2
  ]

  pdf(paste0(prefix, "_manhattan.pdf"), width = 11, height = 5)

  plot(
    gwas$cumulative_position,
    -log10(gwas$P),
    pch = 20,
    cex = 0.5,
    col = colours,
    xaxt = "n",
    xlab = "Chromosome",
    ylab = expression(-log[10](p)),
    main = "Genome-wide association results"
  )

  axis(1, at = ticks, labels = chromosomes)
  abline(h = -log10(5e-8), col = "red", lty = 2)

  dev.off()

  observed <- -log10(sort(gwas$P))
  expected <- -log10(ppoints(nrow(gwas)))

  pdf(paste0(prefix, "_qq.pdf"), width = 6, height = 6)

  plot(
    expected,
    observed,
    pch = 20,
    cex = 0.6,
    xlab = expression(Expected~~-log[10](p)),
    ylab = expression(Observed~~-log[10](p)),
    main = "Q–Q plot"
  )

  abline(0, 1, col = "red", lty = 2)

  dev.off()

  write.table(
    head(gwas[order(gwas$P), ], 10),
    paste0(prefix, "_top10.tsv"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  message("Valid tests plotted: ", nrow(gwas))
  message(
    "Bonferroni threshold for these valid tests: ",
    0.05 / nrow(gwas)
  )
}

plot_gwas("Results/gwas_unadjusted")
plot_gwas("Results/gwas_adjusted")
```

Run:

```bash
Rscript scripts/04_plots.R
```

<div class="warning" markdown="1">

**Do not silently discard failed tests**

The script reports how many additive tests were omitted from plotting.

Investigate missing p-values, numerical zeros, and model warnings.
A plot containing only valid results is not evidence that every model
fitted successfully.

</div>

### 7.3 Interpreting the Manhattan plot

The Manhattan plot shows:

- Genomic location on the horizontal axis.
- Association strength, expressed as `−log10(p)`, on the vertical axis.

Inspect the most prominent regions.

Do not assume that several nearby significant variants represent
several independent discoveries.

### 7.4 Interpreting the Q–Q plot

The Q–Q plot compares observed p-values with the distribution expected
under the null.

Use it as a diagnostic rather than a pass/fail certificate.

Consider the pattern alongside:

- Population structure.
- Relatedness.
- Technical quality.
- Model fit.
- The simulated nature of the phenotype.

<div class="exercise" markdown="1">

**EXERCISE 8 — Compare adjusted and unadjusted results**

Compare both Manhattan plots and both Q–Q plots.

Then select five variants and complete:

| Variant | Tested allele | Unadjusted OR | Adjusted OR | Unadjusted p | Adjusted p |
|---|---|---:|---:|---:|---:|
| | | | | | |

Answer:

1. Did the odds ratios change?
2. Did the p-values change?
3. Were the analyzed sample counts comparable?
4. Did every adjusted p-value become larger?
5. What can you conclude—and what remains uncertain?

</div>

<details markdown="1">
<summary>Answer</summary>

There is no fixed numerical answer.

Adjustment does not have to increase every p-value, and it may have
little visible effect in this teaching dataset.

A good interpretation compares effect estimates, uncertainty, sample
counts, and diagnostics rather than declaring one analysis better
because it produces more significant results.

</details>

## 8 Outputs and submission

Your output directory should now contain:

| Output | Purpose |
|---|---|
| `missingness.imiss` | Sample missingness |
| `missingness.lmiss` | Variant missingness |
| `heterozygosity.het` | Heterozygosity diagnostics |
| `relatedness.genome` | Pairwise relatedness estimates |
| `06_analysis.bed/.bim/.fam` | Retained analysis dataset |
| `structure.mds` | MDS coordinates |
| `covariates.tsv` | Association covariates |
| `gwas_adjusted.assoc.logistic` | Adjusted association results |
| `gwas_adjusted_manhattan.pdf` | Manhattan plot |
| `gwas_adjusted_qq.pdf` | Q–Q plot |

Submit:

1. Your analysis scripts.
2. A table of sample and variant counts at each QC stage.
3. A short justification of reviewed exclusions.
4. The population-structure plot.
5. Adjusted and unadjusted association plots.
6. A results table with tested alleles and effect estimates.
7. A short discussion of limitations.
8. Software versions and the teaching-data provenance record.

**Final question:** Which conclusions are supported by this practical,
and which would require additional data or a different study design?

## References

### Scientific background

Marees AT, de Kluiver H, Stringer S, et al. (2018).
*A tutorial on conducting genome-wide association studies:
Quality control and statistical analysis.*
International Journal of Methods in Psychiatric Research, 27, e1608.

- [Paper](https://doi.org/10.1002/mpr.1608)
- [Paper](https://www.nature.com/articles/s43586-021-00056-9)
- [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC1950838/)
- [Source](https://github.com/MareesAT/GWA_tutorial/)
- [Source](https://cloufield.github.io/GWASTutorial/)
- [Source](https://pbreheny.github.io/adv-gwas-tutorial/index.html)
- [Source](https://gwaslab.org)

### Software documentation

- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
- [PLINK input filtering](https://www.cog-genomics.org/plink/1.9/filter)
- [PLINK basic statistics](https://www.cog-genomics.org/plink/1.9/basic_stats)
- [PLINK file formats](https://www.cog-genomics.org/plink/1.9/formats)
- [PLINK LD pruning](https://www.cog-genomics.org/plink/1.9/ld)
- [PLINK population structure](https://www.cog-genomics.org/plink/1.9/strat)
- [PLINK association analysis](https://www.cog-genomics.org/plink/1.9/assoc)
