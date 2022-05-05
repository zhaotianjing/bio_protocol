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


# 4. Set Factors or Covariates (no covariates for this protocol)
# example:
#  set_covariate(model,"x1");


# 5. Set Random or Fixed Effects
set_random(model,"ID",pedigree);
# example:
#  set_random(model,"x2");
#  set_random(model,"ID dam",pedigree);


# 6. Run Analysis
out=runMCMC(model,phenotypes,chain_length=1000);
# example:
#  out=runMCMC(model, phenotypes, chain_length=50000, burnin=10000, output_samples_frequency= 100);
#  out=runMCMC(model, phenotypes, chain_length=50000, single_step_analysis=true, pedigree=pedigree);

# 7. Check Results
# The estimates (i.e., posterior mean), standard error, and MCMC samples for parameters of interest 
# will be saved in the output folder after the analysis. 
# example:
# -------------------------------------------------------------------------------------------------
# EBV_t1.txt            | estimated breeding values for trait named "t1" 
# EBV_t2.txt            | estimated breeding values for trait named "t2" 
# EBV_t3.txt            | estimated breeding values for trait named "t3" 
# genetic_variance.txt  | estimated genetic variance-covariance for all traits
# residual_variance.txt | estimated residual variance-covariance for all traits
# heritability.txt      | estimated heritability
# pi_genotypes.txt      | estimated pi
# marker_effects_genotypes.txt            | estimated marker effects for all traits
# polygenic_effects_covariance_matrix.txt | estimated variance-covariance between polygenic effects 
# --------------------------------------------------------------------------------------------------



# C. Genomic prediction or GWAS analysis
# 1. Calculate accuracy
results    = innerjoin(out["EBV_t3"], phenotypes, on = :ID)
ind_id     = findall(x -> !ismissing(x), results[!,:t3])
accuruacy  = cor(results[ind_id,:EBV],results[ind_id,:t3]) #cor(phenotype,ebv)

# 2. GWAS
# (a) Compute the model frequency for each marker
marker_effects_file = "results/MCMC_samples_marker_effects_genotypes_t1.txt"
GWAS(marker_effects_file,header=true)

# (b) Compute the posterior probability of association of the genomic window that explains a large proportion of the total genetic variance
# example1:
#  map_file = "map.csv"
#  marker_effects_file = "results/MCMC_samples_marker_effects_genotypes_t1.txt"
#  out=GWAS(model,map_file, marker_effects_file, header=true, window_size="1 Mb",threshold=0.01)

# example2:
#  map_file="map.csv"
#  marker_effects_file1="results/MCMC_samples_marker_effects_genotypes_t1.txt"
#  marker_effects_file2="results/MCMC_samples_marker_effects_genotypes_t2.txt"
#  marker_effects_file3="results/MCMC_samples_marker_effects_genotypes_t3.txt"
#  out=GWAS(model,map_file,marker_effects_file1,marker_effects_file2,marker_effects_file3,header=true,window_size="1 Mb")

# (c) Compute the genetic correlation between two correlated traits for each genomic window
# example:
#  map_file="map.csv"
#  marker_effects_file1 ="results/MCMC_samples_marker_effects_genotypes_y1.txt"
#  marker_effects_file2 ="results/MCMC_samples_marker_effects_genotypes_y2.txt"
#  out=GWAS(model,map_file,marker_effects_file1,marker_effects_file2,genetic_correlation=true,header=true,window_size="1 Mb")



