# Bayesian multiple regression methods for genomic prediction and genome-wide association studies


## Introduction

Advance in technologies of next-generation sequencing and SNP arrays has stimulated the development of whole-genome regression methods for genome-wide prediction and association studies in the last several decades. These whole-genome regression methods are able to impose a wide range of assumptions about the distributions on the effect sizes of causal variants and allow for inference of the causal variants themselves. JWAS is an open-source software that can perform different Bayesian regression methods with genomic data for both genome-wide prediction and association studies. The Bayesian regression methods implemented in JWAS for genome-enable analysis include BayesA, BayesB, BayesC, BayesCpi, Bayesian Lasso, RR-BLUP and GBLUP. Furthermore, JWAS has been developed to perform pedigree-based analysis, including, but not limited to, animal model, animal model with maternal effects, reduced animal model. By combining the genomic and pedigree information, JWAS    


## Installation
* Required software:
    * Julia  
* Install JWAS package inside Julia
    ```julia
    using Pkg; Pkg.add("JWAS")
    ```
## Input Data
The publicaly available pig genotypes, phenotypes, and pedigree data in [^fn1] were used. For demonstration purpose, 5000 SNPs were randomly sampled from raw data.


## Procedure
The workflow of JWAS is shown in below figure. 


Below is an example of multi-trait genomic prediction and GWAS analysis.
```julia
# Step 1: Load packages
using JWAS,DataFrames,CSV,Statistics,Random
Random.seed!(1)

# Step 2: Read data
phenotypes = CSV.read("phenotypes.txt",DataFrame,delim = ',',header=true,missingstrings=["."])
pedigree   = get_pedigree("pedigree.txt",separator=",",header=true)
genotypes  = get_genotypes("genotypes_5k.txt",separator=',',method="BayesC")
first(phenotypes,5)

# Step 3: Build Model Equations
model_equation  ="t1 = intercept + ID + genotypes
                  t2 = intercept + ID + genotypes
                  t3 = intercept + ID + genotypes";
model = build_model(model_equation);

# Step 4: Set Factors or Covariates
# No Factors or Covariates for this data

# Step 5: Set Random or Fixed Effects
set_random(model,"ID",pedigree);

# Step 6: Run Analysis
out=runMCMC(model,phenotypes);

# Check Accuruacy
results    = innerjoin(out["EBV_t3"], phenotypes, on = :ID)
ind_id     = findall(x -> !ismissing(x), results[!,:t3])  #find individuals with phenotypes
accuruacy  = cor(results[ind_id,:EBV],results[ind_id,:t3])

# Compute the model frequency for each marker
marker_effects_file="results/MCMC_samples_marker_effects_genotypes_t1.txt"
GWAS(marker_effects_file,header=true)

```


[^fn1]: Matthew A Cleveland, John M Hickey, Selma Forni, A Common Dataset for Genomic Analysis of Livestock Populations, G3 Genes|Genomes|Genetics, Volume 2, Issue 4, 1 April 2012, Pages 429–435, https://doi.org/10.1534/g3.111.001453
