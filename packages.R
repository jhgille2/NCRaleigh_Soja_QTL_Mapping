# Install pacman if it does not already exist
if(!require(pacman)){
  install.packages(pacman)
}

# Use pacman to load/install packages
pacman::p_load(conflicted, 
               dotenv, 
               targets, 
               tarchetypes, 
               tidyverse, 
               rmarkdown, 
               qtl, 
               here, 
               readxl, 
               janitor, 
               magrittr, 
               ggthemes, 
               lattice, 
               reactable,
               emmeans, 
               vroom, 
               vcfR, 
               ASMap)

# Conflict preferences
conflict_prefer("filter", "dplyr")

