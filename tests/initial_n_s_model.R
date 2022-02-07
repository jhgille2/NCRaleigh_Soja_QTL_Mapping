load(here("data", "all_refined.RData"))
load(here("data", "all_stepwise.RData"))
load(here("data", "filled_genos.RData"))

# Start with the full model as fit by stepwise qtl 
fit_init <- fitqtl(filled_genos[[5]], 
                   pheno.col = 2, 
                   qtl = all_refined[[5]], 
                   formula = attr(all_refined[[5]], "formula"))

# Lookup table to map formula qtl names to long qtl names
qtl_name_lookup <- tibble(qtl_shortname = all_refined[[5]]$name, 
                          qtl_longname = all_refined[[5]]$altname)

# Many of the interaction terms do not seems to be significant, remove these from the model and refine/refit with the new formula
model_reduced_1  <- "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10 + Q11 + Q12 + Q13 + Q14 + Q15 + Q16 + Q17 + Q18 + Q19 + Q20 + Q1:Q17 + Q9:Q14 + Q18:Q19"
refine_reduced_1 <- refineqtl(filled_genos[[5]], pheno.col = 2, qtl = all_refined[[5]], formula = model_reduced_1)


qtl_name_lookup <- tibble(qtl_shortname = refine_reduced_1$name, 
                          qtl_longname = refine_reduced_1$altname)

model_reduced_2  <- "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q7 + Q8 + Q9 + Q11 + Q13 + Q14 + Q15 + Q17 + Q18 + Q19 + Q20 + Q1:Q17 + Q9:Q14 + Q18:Q19"
refine_reduced_2 <- refineqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_1, formula = model_reduced_2)

qtl_name_lookup <- tibble(qtl_shortname = refine_reduced_2$name, 
                          qtl_longname = refine_reduced_2$altname)

fitqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_2, formula = model_reduced_2)

# Marker regression for comparison
n_s_regression_imp <- scanone(filled_genos[[5]], pheno.col = 2, method = "imp")
n_s_regression_mr  <- scanone(LinkageMap, pheno.col = 2, method = "mr")
n_s_regression_em  <- scanone(calc.genoprob(LinkageMap), pheno.col = 2, method = "em")

# Permutations
job::job({
  
  n_s_regression_imp_perms <- scanone(filled_genos[[5]], pheno.col = 2, method = "imp", n.perm = 1000, n.cluster = 10)
  n_s_regression_mr_perms  <- scanone(LinkageMap, pheno.col = 2, method = "mr", n.perm = 1000, n.cluster = 10)
  n_s_regression_em_perms  <- scanone(calc.genoprob(LinkageMap), pheno.col = 2, method = "em", n.perm = 1000, n.cluster = 10)
  
})

plot(n_s_regression_imp, n_s_regression_mr, n_s_regression_em)

cross_1000_imp <- filled_genos[[5]]



model_reduced_3  <- "y ~ Q1 + Q2 + Q4 + Q7 + Q8 + Q9 + Q11 + Q13 + Q14 + Q15 + Q17 + Q18 + Q19 + Q20 + Q1:Q17 + Q9:Q14 + Q18:Q19"
refine_reduced_3 <- refineqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_2, formula = model_reduced_3)
fitqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_3, formula = model_reduced_3)




model_reduced_4  <- "y ~ Q1 + Q4 + Q7 + Q8 + Q9 + Q11 + Q13 + Q14 + Q15 + Q17 + Q18 + Q19 + Q20 + Q9:Q14 + Q18:Q19"
refine_reduced_4 <- refineqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_3, formula = model_reduced_4)
fitqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_4, formula = model_reduced_4)

qtl_name_lookup <- tibble(qtl_shortname = refine_reduced_4$name, 
                          qtl_longname = refine_reduced_4$altname)

model_reduced_5  <- "y ~ Q1 + Q4 + Q7 + Q8 + Q9 + Q11 + Q13 + Q14 + Q15 + Q18 + Q19 + Q20 + Q9:Q14 + Q18:Q19"
refine_reduced_5 <- refineqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_4, formula = model_reduced_5)
fitqtl(filled_genos[[5]], pheno.col = 2, qtl = refine_reduced_5, formula = model_reduced_5)


save(n_s_regression_imp_perms, file = here("data", "mapping_example", "n_s_regression_imp_perms.RData"))
save(n_s_regression_imp, file = here("data", "mapping_example", "n_s_regression_imp.RData"))
save(cross_1000_imp, file = here("data", "mapping_example", "filled_genos_1000.RData"))
save(refine_reduced_5, file = here("data", "mapping_example", "final_qtl.RData"))
save(model_reduced_5, file = here("data", "mapping_example", "final_model.RData"))
