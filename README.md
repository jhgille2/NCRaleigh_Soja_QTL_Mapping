# Raleigh Soja Mapping
Code for analysis of a Raleigh x Soja soybean mapping population.  

[writeup here](https://github.com/jhgille2/Raleigh_Soja_Mapping/blob/main/docs/Notebook.pdf)  
[journal here](https://github.com/jhgille2/NCRaleigh_Soja_QTL_Mapping/blob/main/docs/Journal.pdf)  

## Workflow overview

I'm using the [targets](https://github.com/ropensci/targets) package to orchestrate this project, and the [tflow](https://github.com/MilesMcBain/tflow) package to organize the file structure. Because the organization might seem a bit wonky, I'll give a quick "roadmap" of how the repo is organized. A good starting point is the **\_targets.R** file in the main directory. This is the file that tells the targets package how to run the workflow, essentially it is a list of named steps where most of the steps call a custom function. These custom functions in turn are found in the R folder. The **packages.R** file holds the names of the packages used in the workflow. Input/output data is found in the **data** folder, and a couple writeups are in the **docs** folder. The **\_targets** folder is used only by the targets package to manage the workflow. 
