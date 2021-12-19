
genofile1 <- here("data", "old_samples_AB.csv")
genofile2 <- here("data", "new_samples_AB.csv")

genodata1 <- vroom::vroom(genofile1, skip = 9)
genodata2 <- vroom::vroom(genofile2, skip = 9)


snpData <- left_join(genodata1, genodata2, by = "...1") %>% 
  column_to_rownames("...1") %>% 
  relocate(c("2104", "2105"))


# Functions to get the chromosome and bp position from a snp name
get_lg_from_snp_name <- function(snp_name)
{
  snp_name %>% 
    str_split(., pattern = "_") %>% 
    pluck(1, 1) %>%
    parse_number() %>% 
    as.numeric()
}

get_bp_from_snp_name <- function(snp_name)
{
  snp_name %>% 
    str_split(., pattern = "_") %>% 
    pluck(1, 2) %>%
    parse_number() %>% 
    as.numeric()
}

mp_samples <- sort(as.numeric(colnames(snpData))) %>% as.character()

snpData <- snpData[, mp_samples] %>% 
  relocate(c("2104", "2105"))

apply(snpData, 1, function(x) table(unlist(x)))

snpData[snpData == snpData[, "2104"]] <- "A"
snpData[snpData == "AB"] <- "H"
snpData[snpData == "--"] <- "-"
snpData[snpData == "BB"] <- "B"
snpData[snpData == "AA"] <- "B"

snpData %<>% 
  mutate(chr = map_dbl(rownames(.), get_lg_from_snp_name), 
         pos = map_dbl(rownames(.), get_bp_from_snp_name)) %>% 
  relocate(chr, pos) %>% 
  rownames_to_column(var = "id")

colnames(snpData)[2:3] <- ""

write.table(snpData, file = here("data", "fullgenodata.csv"), row.names = FALSE, col.names = TRUE, sep = ",")

phenofile <- read.csv(here("data", ))








snpData_cross <- read.cross(format = "csvsr", 
                            genfile = here("data", "fullgenodata.csv"), 
                            phefile = here("data", "rqtl_phenotypes.csv"), 
                            F.gen = 4)

heatMap(snpData_cross)


snpData_cross <- pullCross(snpData_cross, type = "co.located")
snpData_cross_segdist <- pullCross(snpData_cross, type = "seg.distortion", pars = list(seg.thresh = 0.001))
snpData_cross_noMiss_5 <- pullCross(snpData_cross_segdist, type = "missing", pars = list("miss.thresh" = 0.05))

snpData_mst_pass1 <- mstmap.cross(snpData_cross_noMiss_5, id = "id", anchor = TRUE)
