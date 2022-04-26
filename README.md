# Bayesian multiple regression methods for genomic prediction and genome-wide association studies


## Introduction

Advancement in technologies of next-generation sequencing and SNP arrays has stimulated the development of whole-genome regression methods for genome-wide prediction and association studies in the last several decades. These whole-genome regression methods are able to impose a wide range of assumptions about the distributions on the effect sizes of causal variants and allow for inference of the causal variants themselves. JWAS [^fn2] is an open-source software that can perform different Bayesian regression methods with genomic data for both genome-wide prediction and association studies. The Bayesian regression methods implemented in JWAS for genome-enable analysis include BayesA, BayesB, BayesC, BayesCpi, Bayesian Lasso, RR-BLUP and GBLUP. Furthermore, JWAS has been developed to perform pedigree-based analysis, including, but not limited to, animal model, animal model with maternal effects, reduced animal model. By providing both genomic and pedigree data, JWAS can also perform single-step Bayesian regression models to integrate both sources of information in genome-enabled analysis. Lastly, all the aforementioned methods that are implemented in JWAS allow for both single-trait and multi-trait analysis.   


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
The workflow of JWAS is shown in below figure. The major steps are (1) reading data, (2) building model, and (3) running analysis.

![Screen Shot 2022-04-26 at 9 59 35 AM](https://user-images.githubusercontent.com/18593116/165353767-65da93ba-2b24-4b79-82d4-007f34637b8d.png)



Below is an example of multi-trait genomic prediction and GWAS analysis. Estimates, standard deviations, and MCMC samples for variables of interest can be found in the output folder.

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
[^fn2]: Cheng, Hao, Rohan Fernando, and Dorian Garrick. "JWAS: Julia implementation of whole-genome analysis software." Proceedings of the world congress on genetics applied to livestock production. Vol. 11. 2018.
