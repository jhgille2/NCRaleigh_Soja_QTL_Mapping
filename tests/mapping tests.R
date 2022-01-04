tar_load(LinkageMap)

genetic_clones <- genClones(LinkageMap, id = "id")
LinkageMap_noClones <- fixClones(LinkageMap, gc = genetic_clones$cgd, id = "id", consensus = FALSE)

gen_stats <- statGen(LinkageMap_noClones, id = "id", bychr = FALSE)
which(gen_stats$miss > sum(nmar(LinkageMap_noClones)) * 0.1)

Genotypes_to_keep <-  which(gen_stats$miss < sum(nmar(LinkageMap_noClones)) * 0.05)

LinkageMap_noClones_NoMissing <- subset(LinkageMap_noClones, ind = Genotypes_to_keep)


test_scan <- scanone(LinkageMap_noClones_NoMissing, 
                     method = "em", 
                     pheno.col = 2:5)

test_cross <- calc.genoprob(LinkageMap_noClones_NoMissing)

test_cross_imputed <- sim.geno(LinkageMap_noClones_NoMissing, n.draws = 1000, map.function = "kosambi")
test_cross_imputed_small <- sim.geno(LinkageMap_noClones_NoMissing, n.draws = 100, map.function = "kosambi")

scantwo_perms_1000 <- scantwo(test_cross, 
                         pheno.col = 2, 
                         method = "hk", 
                         n.perm = 1000, 
                         n.cluster = 5)

test_stepwise <- stepwiseqtl(test_cross, 
                             pheno.col = 2, 
                             keeptrace = TRUE, 
                             penalties = calc.penalties(scantwo_perms), 
                             max.qtl = 20)

test_stepwise_10 <- stepwiseqtl(test_cross, 
                             pheno.col = 2, 
                             keeptrace = TRUE, 
                             penalties = calc.penalties(scantwo_perms_1000), 
                             max.qtl = 15)

test_stepwise_imputed <- stepwiseqtl(test_cross_imputed, 
                                     pheno.col = 2, 
                                     keeptrace = TRUE, 
                                     penalties = calc.penalties(scantwo_perms_1000), 
                                     max.qtl = 15)

test_stepwise_imputed_small <- stepwiseqtl(test_cross_imputed_small, 
                                           pheno.col = 2, 
                                           keeptrace = TRUE, 
                                           penalties = calc.penalties(scantwo_perms_1000), 
                                           max.qtl = 15)

stepwise_10_fitqtl            <- fitqtl(test_cross, pheno.col = 2, test_stepwise_10)
stepwise_imputed_fitqtl       <- fitqtl(test_cross_imputed, pheno.col = 2, test_stepwise_imputed)
stepwise_imputed_fitqtl_small <- fitqtl(test_cross_imputed_small, pheno.col = 2, test_stepwise_imputed_small)

test_stepwise_10_reduced      <- dropfromqtl(test_stepwise_10, which(stepwise_10_fitqtl$result.drop[, 7] >= 0.05))
test_stepwise_imputed_reduced <- dropfromqtl(test_stepwise_imputed, which(stepwise_imputed_fitqtl$result.drop[, 7] >= 0.05))
test_stepwise_imputed_reduced_small <- dropfromqtl(test_stepwise_imputed_small, which(stepwise_imputed_fitqtl_small$result.drop[, 7] >= 0.05))

fitqtl(test_cross, pheno.col = 2, test_stepwise_10_reduced)$result.drop
fitqtl(test_cross_imputed, pheno.col = 2, test_stepwise_imputed_reduced)$result.drop
fitqtl(test_cross_imputed_small, pheno.col = 2, test_stepwise_imputed_reduced_small)$result.drop

plan(multisession, workers = 2)

cross_df <- tibble(n_draws = c(100, 250, 500, 750)) %>% 
  mutate(cross_imputed = future_map(n_draws, function(x) sim.geno(LinkageMap_noClones_NoMissing, n.draws = x, map.function = "kosambi"), .options = furrr_options(seed = TRUE)), 
         pheno_index   = 2)

# A function to drop qtl with a f score below 0.05 from a fitqtl object
drop_non_sig_qtl <- function(fitted_qtl, sig.level = 0.05)
{
  dropfromqtl(fitted_qtl, which(fitted_qtl$result.drop[, 7] >= 0.05))
}

stepwise_df <- cross_df %>% 
  mutate(stepwise_fit = future_map(cross_imputed, function(x) stepwiseqtl(x, 
                                                                          pheno.col = 2, 
                                                                          penalties = calc.penalties(scantwo_perms_1000), 
                                                                          max.qtl = 15)))


# fit_summary = future_map(stepwise_fit, function(x) fitqtl(cross_imputed, pheno.col = 2, stepwise_fit)), 
# sig_qtl = map(fit_summary, drop_non_sig_qtl)

test_tbl %>% 
mutate(sum = pmap_dbl(., test_fn, z = 2))

job::job({
  
})