# This is a script to run corncob on a large dataset
# I will run this on hipergator
# Make sure all files are in the directory of interest (/blue/juliemeyer/cynthiabecker/R_data_analysis/corncob)

# load packages
library(phyloseq)
library(corncob)

# save version numbers of all packages and software
writeLines(capture.output(sessionInfo()), "sessionInfo_9.10.24.txt")

# load in the phyloseq object - only annotated genes and their counts across all files
ps.count.annotated <- readRDS("FLK_coassembly_phyloseq_count.annotated.rds")

# Set the levels so that ALDEx2 will choose unaffected as baseline
sample_data(ps.count.annotated)$Disease_on_colony_when_coring_binary <- factor(sample_data(ps.count.annotated)$Disease_on_colony_when_coring_binary,
                                                                     levels = c("Unaffected", "Diseased"))

set.seed(10) # set seed for reproducibility
corncobDA.dis <- differentialTest(formula = ~ Disease_on_colony_when_coring_binary, 
                                   phi.formula = ~ Disease_on_colony_when_coring_binary,
                                   formula_null = ~ 1,
                                   phi.formula_null = ~ Disease_on_colony_when_coring_binary,
                                   test = "Wald", boot = FALSE,
                                   data = ps.count.annotated,
                                   fdr_cutoff = 0.05)

# Save the corncob phyloseq object so you can manipulate it in R
saveRDS(corncobDA.dis, "corncob_output_disease_at_sampling.rds")

## Try controlling for the differences in the symbiont community
set.seed(10) # set seed for reproducibility
corncobDA.dis.symb <- differentialTest(formula = ~ Disease_on_colony_when_coring_binary + Zoox_background, 
                                  phi.formula = ~ Disease_on_colony_when_coring_binary + Zoox_background,
                                  formula_null = ~ Zoox_background,
                                  phi.formula_null = ~ Disease_on_colony_when_coring_binary + Zoox_background,
                                  test = "Wald", boot = FALSE,
                                  data = ps.count.annotated,
                                  fdr_cutoff = 0.05)

# Save the corncob phyloseq object so you can manipulate it in R
saveRDS(corncobDA.dis.symb, "corncob_output_disease_at_sampling_zoox.rds")
