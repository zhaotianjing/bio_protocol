# A. install packages
# 1. download Julia programming language (https://julialang.org/)
# 2. open Julia and install JWAS package by running: using Pkg; Pkg.add("JWAS");
# 3. repeat A.2. to install other required packages in Julia:
#    using Pkg; Pkg.add("DataFrames"); Pkg.add("CSV"); Pkg.add("Statistics");


# B. Run Bayesian multiple regression methods in JWAS
# 1. Load packages
using JWAS,DataFrames,CSV,Statistics,Random
Random.seed!(1)

# 2. Read data
phenotypes = CSV.read("phenotypes.txt",DataFrame,delim = ',',header=true,missingstrings=["."])
pedigree   = get_pedigree("pedigree.txt",separator=",",header=true)
genotypes  = get_genotypes("genotypes_5k.txt",separator=',',method="BayesC")
first(phenotypes,5)

# 3. Build Model Equations
model_equation  ="t1 = intercept + genotypes
                  t2 = intercept + ID + genotypes
                  t3 = intercept + ID + genotypes";
model = build_model(model_equation);

# 4. Set Factors or Covariates
# No Factors or Covariates for this data

# 5. Set Random or Fixed Effects
set_random(model,"ID",pedigree);

# 6. Run Analysis
out=runMCMC(model,phenotypes);

# C. Genomic prediction or GWAS analysis
# 1. Calculate accuracy
results    = innerjoin(out["EBV_t3"], phenotypes, on = :ID)
ind_id     = findall(x -> !ismissing(x), results[!,:t3])
accuruacy  = cor(results[ind_id,:EBV],results[ind_id,:t3])

# 2. Compute the model frequency for each marker
marker_effects_file="results/MCMC_samples_marker_effects_genotypes_t1.txt"
GWAS(marker_effects_file,header=true)
