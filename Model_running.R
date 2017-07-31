
### Clear the workspace
rm(list=ls())

### load the libraries
library(ape)
library(BioGeoBEARS)
library(optimx)
library(FD)     # for FD::maxent() (make sure this is up-to-date)
library(snow)
library(parallel)
library(plotrix)

### load the updates recommended by the developer
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)  # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)


###GET FILES
#my_newick<-"MrBayes/Tree_files/MrBayes_treesets/run1_MCC.newick"

#### Rename tree tips
#import the MCC tree
#code.tree<-read.tree(file = "MrBayes/Tree_files/MrBayes_treesets/run1_MCC.newick")
#import names
#names<-read.csv(file="Metadata/Multilocus_names.csv", stringsAsFactors = FALSE)

#name.tree<-code.tree
#name.tree$tip.label<-names$Alt_label[match(name.tree$tip.label,names$Name)]
#write.tree(name.tree, file="MrBayes/Tree_files/MrBayes_treesets/run1_MCC_named.newick")


##########PREPARE NEWICK FORMATTED TREE FROM NEXUS
### read in the MCC Cratopus tree from tree annotator
my.nexus.trees<-read.nexus(file = "MrBayes/Tree_files/MrBayes_treesets/Run1_dated.t")

for(i in sample(1:length(my.nexus.trees), 10, replace = FALSE)){
  sampled.tree<-my.nexus.trees[i]
  write.tree(sampled.tree, file=paste("MrBayes/Tree_files/MrBayes_treesets/Run1_",i,"_newick.tree", sep=""))
}

###Run DEC+J+X

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001  ###Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE  ###set to FALSE for e.g. DEC* model, DEC*+J, etc.
###(For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
### Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
### the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
### DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
###Also: search script on "include_null_range" for other places to change

###Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = dist.mat
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = FALSE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE  ###if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 6
BioGeoBEARS_run_object$force_sparse = FALSE  ###force_sparse=TRUE causes pathology & isn't much faster at this scale

###This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
###(It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

###Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
###The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

###Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE  ###get ancestral states from optim run

###Set up DEC+J+X model
###Get the ML parameter values from the 2-parameter nested model
###(this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

###Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DEC+J+X.Rdata"
runslow = TRUE
if (runslow)
{
  
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECjx = res
} else {
  ###Loads to "res"
  load(resfn)
  resDECjx = res
}