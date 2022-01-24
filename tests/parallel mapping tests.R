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

future::plan(remote, workers = paste0("node", str_pad(0:5, 3, pad = "0")), verbose = TRUE)

filled_genos <- future_map(drawnums, function(x) sim.geno(LinkageMap, n.draws = x, map.function = "kosambi"))
all_stepwise <- future_map(drawnums, function(x) stepwiseqtl(x, pheno.col = 2, max.qtl = 20, method = "imp", penalties = scantwo_pens))
