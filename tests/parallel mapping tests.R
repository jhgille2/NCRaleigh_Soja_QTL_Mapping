tar_load(LinkageMap)

map_simgeno <- sim.geno(LinkageMap, n.draws = 500, map.function = "kosambi")

scantwo_perms <- scantwo(LinkageMap, pheno.col = 2, method = "hk", n.perm = 1000, n.cluster = 10)
scantwo_pens <- calc.penalties(scantwo_perms)

LinkageMap <- calc.genoprob(LinkageMap)
stepwise_n_s <- stepwiseqtl(LinkageMap, 
            pheno.col = 2, 
            max.qtl = 20, 
            method = "hk", 
            penalties = calc.penalties(scantwo_perms))


LinkageMap_simgeno <- sim.geno(LinkageMap, n.draws = 500)
stepwise_n_s_imp <- stepwiseqtl(LinkageMap_simgeno, 
                                pheno.col = 2, 
                                max.qtl = 20, 
                                method = "imp", 
                                penalties = calc.penalties(scantwo_perms))


drawnums <- c(100, 250, 500, 750, 1000)

future::plan(multisession, workers = 5)

filled_genos <- future_map(drawnums, function(x) sim.geno(LinkageMap, n.draws = x, map.function = "kosambi"), .options = furrr_options(seed = TRUE))
all_stepwise <- future_map(filled_genos, function(x) stepwiseqtl(x, pheno.col = 2, max.qtl = 20, method = "imp", penalties = scantwo_pens), .options = furrr_options(seed = TRUE))


test <- fitqtl(filled_genos[[1]], pheno.col = 2, qtl = all_stepwise[[1]])
test2 <- dropfromqtl(all_stepwise[[1]], which(test$result.drop[, 7] > 0.05))

# A function to drop non-significant (pvalue(F) > 0.05) from the fit model and 
# return the new model minus these terms
drop_non_significant_qtl <- function(Cross, qtl_obj)
{
  # First, fit the full model
  full_model    <- fitqtl(cross = Cross, pheno.col = 2, qtl_obj)
  
  # Drop qtl with a f value less than 0.05
  reduced_model <- dropfromqtl(qtl_obj, which(full_model$result.drop[, 7] > 0.05))
  
  final_model <- fitqtl(Cross, pheno.col = 2, reduced_model)
  
  return(list(qtl = reduced_model, fit = final_model))
}

# All the reduced models
all_reduced <- map2(filled_genos, all_stepwise, drop_non_significant_qtl)

hk_stepwise <- stepwiseqtl(calc.genoprob(LinkageMap), 
                           pheno.col = 2, 
                           max.qtl = 20,
                           method = "hk", 
                           keeplodprofile = TRUE, 
                           penalties = calc.penalties(scantwo_perms))

drop_non_significant_qtl(calc.genoprob(LinkageMap), hk_stepwise)$result.drop

test_nitrogen_hk_stepwise <- stepwiseqtl(calc.genoprob(LinkageMap), 
                                         pheno.col = 5, 
                                         max.qtl = 20, 
                                         method = "hk", 
                                         keeplodprofile = TRUE, 
                                         penalties = calc.penalties(scantwo_perms))
