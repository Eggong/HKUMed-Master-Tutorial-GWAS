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

These are PLINK input codes, not the mathematical coding used when writing a logistic regression equation.

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

For example, we can measure missingness automatically, but deciding whether to exclude a sample may require considering several pieces of evidence.

The thresholds below are choices for this tutorial. They are not a universal recipe for every GWAS.

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

./plink \
  --bfile Data/raw/study \
  --geno <<THRESHOLD>> \
  --make-bed \
  --out Results/01_variant_missingness

./plink \
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

./plink \
  --bfile Data/raw/study \
  --geno 0.02 \
  --make-bed \
  --out Results/01_variant_missingness

./plink \
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
./plink \
  --bfile Results/02_missingness \
  --freq \
  --out Results/allele_frequency
```

Inspect the report:

```bash
head Results/allele_frequency.frq
```

For this teaching dataset, retain variants with MAF of at least 5%:

```bash
./plink \
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

### 4.3 Hardy–Weinberg equilibrium (HWE)

For a biallelic variant with allele frequencies \(p\) and \(q\),
Hardy–Weinberg equilibrium gives expected genotype frequencies:

```text
p², 2pq, q²
```

Generate a diagnostic report:

```bash
./plink \
  --bfile Results/03_maf \
  --hardy \
  --out Results/hwe
```

For the classroom filtering exercise:

```bash
./plink \
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

Check for gender discrepancy.

Subjects who were a priori determined as females are expected to have an X-chromosome inbreeding coefficient, \(F\), of less than 0.2. Subjects who were a priori determined as males are expected to have an \(F\) value greater than 0.8. This \(F\) value is based on the X-chromosome inbreeding, or homozygosity, estimate.

Subjects who do not meet these criteria are flagged as `PROBLEM` by PLINK.

```bash
./plink \
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

In this tutorial, individuals flagged as `PROBLEM` are removed to generate the next binary PLINK fileset, `Results/05_sexcheck`.

</div>

<div class="exercise" markdown="1">

**EXERCISE 4 — Remove individuals with gender discrepancy**

Generate a list of individuals with the `PROBLEM` status:

```bash
grep "PROBLEM" Results/sex_check.sexcheck | awk '{print $1, $2}' > Results/sex_discrepancy.txt
```

Remove those individuals and create a new binary PLINK fileset:

```bash
./plink \
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

In genetics, LD occurs when alleles at nearby loci are inherited together more often than expected by chance. If two SNPs are strongly correlated, they provide largely redundant information. LD pruning keeps a subset of relatively independent variants.

LD pruning means reducing a set of genetic variants by removing SNPs that are highly correlated due to linkage disequilibrium, or LD.

We will next create a reduced marker set for heterozygosity, relatedness, and population-structure calculations:

```bash
./plink \
  --bfile Results/05_sexcheck \
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
Heterozygosity is an individual-level quality-control metric used to identify samples with unusually high or low levels of heterozygous genotypes. In diploid organisms, a genotype is heterozygous when the two alleles at a locus differ. Across many SNPs, each individual is expected to have a broadly similar heterozygosity rate, assuming comparable ancestry and good genotype quality.

Unusual heterozygosity can indicate sample or data problems.

Generate a heterozygosity report:

```bash
./plink \
  --bfile Results/05_sexcheck \
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

This plot shows the distribution of observed heterozygosity across individuals in a GWAS QC step.

#### 4.4.4 Relatedness

Relatedness is used to identify pairs or groups of individuals who are genetically more similar than expected for unrelated samples. In many genome-wide association studies and population genetic analyses, the assumption is that individuals are unrelated. If close relatives are included without accounting for relatedness, test statistics can become inflated and association results may be biased.

Relatedness checks are therefore commonly used to detect duplicate samples, monozygotic twins, parent–offspring pairs, siblings, and other close relatives.

Estimate identity by descent using the pruned markers:

```bash
./plink \
  --bfile Results/05_sexcheck \
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
  ./plink \
    --bfile Results/05_sexcheck \
    --remove Data/review/qc_exclusions.txt \
    --make-bed \
    --out Results/06_sample_qc
else
  ./plink \
    --bfile Results/05_sexcheck \
    --make-bed \
    --out Results/06_sample_qc
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

1. Examine population structure using a PCA-like approach.
2. Calculate covariates within the retained study sample.


### 5.2 Define the analysis sample

Apply only the reviewed population exclusions and restrict this introductory association analysis to autosomes:

```bash
if [ -s Data/review/population_exclusions.txt ]; then
  ./plink \
    --bfile Results/06_sample_qc \
    --remove Data/review/population_exclusions.txt \
    --autosome \
    --make-bed \
    --out Results/07_analysis
else
  ./plink \
    --bfile Results/06_sample_qc \
    --autosome \
    --make-bed \
    --out Results/07_analysis
fi
```

### 5.4 Calculate MDS dimensions

We will use multidimensional scaling (MDS), which is widely used to summarize genetic similarity among individuals and to visualize population structure, to examine structure within the retained study samples. In GWAS quality control, MDS dimensions are often calculated from genome-wide identity-by-state or genetic distance estimates and used to identify ancestry outliers or as covariates in association analyses.

Recalculate the pruned marker set:

```bash
./plink \
  --bfile Results/07_analysis \
  --indep-pairwise 50 5 0.2 \
  --out Results/structure_pruning
```

Calculate 50 dimensions:

```bash
./plink \
  --bfile Results/07_analysis \
  --extract Results/structure_pruning.prune.in \
  --cluster \
  --mds-plot 50 \
  --out Results/structure
```

MDS and principal component analysis (PCA) are related approaches to dimension reduction, but they are not identical procedures. Briefly speaking, MDS starts from a distance or dissimilarity matrix between individuals and finds low-dimensional coordinates that preserve those distances, whereas PCA starts from the original genotype or feature matrix and finds axes that explain the greatest variance in the data.

### 5.5 Plot the dimensions and prepare covariates

Create `scripts/02_structure.R` and run:

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

Inspect `Results/structure.pdf`.

<div class="note" markdown="1">

**How many dimensions?**

Ten dimensions are a classroom choice, not a universal recommendation.

As a sensitivity exercise, compare models using two, five, and ten dimensions. Consider sample size and model stability as well as changes in the results.

**Using an elbow plot**

An elbow plot can help choose the number of dimensions to retain. Plot the variance explained, eigenvalues, or another measure of information captured by each MDS dimension. The point where the curve begins to flatten, often called the “elbow,” suggests that additional dimensions contribute relatively little extra information.

In practice, the elbow plot should be used as a guide rather than a strict rule. If the elbow occurs around dimension \(k\), consider using approximately \(k\) dimensions and compare results with nearby choices, such as \(k - 1\) and \(k + 1\), to assess sensitivity.

Create `scripts/03_mds_elbow.R` and run:

```r
mds <- read.table(
  "Results/structure.mds",
  header = TRUE,
  colClasses = c(FID = "character", IID = "character")
)

dimensions <- grep("^C[0-9]+$", names(mds), value = TRUE)

# Calculate the variance of each MDS dimension
dimension_variance <- apply(
  mds[, dimensions],
  2,
  var,
  na.rm = TRUE
)

elbow <- data.frame(
  dimension = seq_along(dimension_variance),
  variance = dimension_variance,
  proportion = dimension_variance / sum(dimension_variance)
)

pdf("Results/mds_elbow.pdf")
plot(
  elbow$dimension,
  elbow$proportion,
  type = "b",
  pch = 19,
  xlab = "MDS dimension",
  ylab = "Proportion of variance",
  main = "Elbow plot for MDS dimensions"
)
dev.off()
```

Inspect `Results/mds_elbow.pdf`.

</div>

## 6 Association testing

### 6.1 The model

After preparing phenotype, genotype, and covariate files, run association testing using logistic regression. Logistic regression is appropriate when the phenotype is binary, such as case/control disease status.

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

First, run an unadjusted association analysis without population-structure covariates. This provides a baseline result that can be compared with the adjusted analysis.

Run:

```bash
./plink \
  --bfile Results/07_analysis \
  --logistic hide-covar \
  --ci 0.95 \
  --out Results/gwas_unadjusted
```

The unadjusted analysis tests the association between each SNP and case/control status without accounting for ancestry or other covariates.

This analysis is useful as a comparison, but it should not usually be treated as the final result because population structure, sex, age, batch effects, or other confounders may inflate or bias associations.

### 6.3 Structure-adjusted analysis

Next, run a structure-adjusted association analysis using the MDS dimensions as covariates. This helps account for population structure within the study sample.

<div class="exercise" markdown="1">

**EXERCISE 6 — Add population-structure covariates**

Use the study dataset and the covariates generated in Section 5. In this example, we use the top ten dimensions as covariates. However, you are encouraged to vary the number of dimensions and compare the results to examine how adjustment for population structure affects the GWAS analysis.

Complete:

```bash
./plink \
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

Run:

```bash
./plink \
  --bfile Results/07_analysis \
  --covar Results/covariates.tsv \
  --covar-name C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 \
  --logistic hide-covar \
  --ci 0.95 \
  --out Results/gwas_adjusted
```

`hide-covar` hides covariate-specific result rows. It does not remove covariates from the model.

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

The second is a reference convention, not a threshold derived from the number of variants in this small teaching dataset.

### 7.2 Manhattan and Q–Q plots

The following R code uses base R and does not require additional plotting packages.

Create `scripts/04_plots.R` and run:

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
    main = "Q-Q plot"
  )
  
  abline(0, 1, col = "red", lty = 2)
  
  dev.off()
  
  ## Calculate lambda GC
  chisq <- qchisq(
    p = 1 - gwas$P,
    df = 1
  )
  
  lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
  
  message("Lambda GC: ", lambda_gc)
  
  write.table(
    head(gwas[order(gwas$P), ], 10),
    paste0(prefix, "_top10.tsv"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  write.table(
    data.frame(
      prefix = prefix,
      n_valid_tests = nrow(gwas),
      lambda_gc = lambda_gc,
      bonferroni_threshold = 0.05 / nrow(gwas)
    ),
    paste0(prefix, "_summary.tsv"),
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

<div class="warning" markdown="1">

**Do not silently discard failed tests**

The script reports how many additive tests were omitted from plotting.

Investigate missing p-values, numerical zeros, and model warnings. A plot containing only valid results is not evidence that every model fitted successfully.

</div>

### 7.3 Interpreting the Manhattan plot

The Manhattan plot shows:

- Genomic location on the horizontal axis.
- Association strength, expressed as `−log10(p)`, on the vertical axis.

Inspect the most prominent regions.

Do not assume that several nearby significant variants represent several independent discoveries.

### 7.4 Interpreting the Q–Q plot

The Q–Q plot compares observed p-values with the distribution expected under the null.

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

Adjustment does not have to increase every p-value, and it may have little visible effect in this teaching dataset.

A good interpretation compares effect estimates, uncertainty, sample counts, and diagnostics rather than declaring one analysis better because it produces more significant results.

</details>

## 8 Expected outputs

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
