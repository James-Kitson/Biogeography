

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
my.nexus.tree<-read.nexus(file = "MrBayes/Tree_files/MrBayes_treesets/run1_out.t.MCC")

### BioGeoBEARS doesn't like nexus files so we export this as .newick format and reimport it.
write.tree(my.nexus.trees, file="MrBayes/Tree_files/MrBayes_treesets/Run1_newick.tree")


###GET FILES
my_newick<-"MrBayes/Tree_files/MrBayes_treesets/Run1_newick.tree"

###"trfn" = "tree file name"
trfn <- my_newick
  
###Look at the raw Newick file:
moref(trfn)
  
###Look at your phylogeny:
my.tree <- read.tree(trfn)

#####provide a path to the geography file
geogfn <- "Metadata/BioGeoBEARS_islands_old.txt"
  
#####rescale your distance matrix and export it for BioGeoBEARS
distances<-as.matrix(read.table("Metadata/BioGeoBEARS_dist.txt", sep="\t", header = TRUE))
distances2<-rescale(distances,c(0,1))
write.table(distances2, file="Metadata/BioGeoBEARS_dist_scale.txt", quote=FALSE, row.names = FALSE, sep="\t")

### add an end of file line to our rescalled distance matrix to satify BioGeoBEARS formatting
line1<-""
line2<-"END"
write(c(line1,line2),file="Metadata/BioGeoBEARS_dist_scale.txt",append=TRUE)

### provide a path to the distance matrix
dist.mat<-"Metadata/BioGeoBEARS_dist_scale_old.txt"
  
###Look at the raw geography text file:
moref(geogfn)
  
###Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
  
###Set the maximum number of areas any species may occupy; this cannot be larger
###than the number of areas you set up, but it can be smaller.
max_range_size = 2
  



###DEC, DEC+J AND DEC+J+X ANALYSIS


###NOTE: The BioGeoBEARS "DEC" model is identical with 
###the Lagrange DEC model, and should return identical
###ML estimates of parameters, and the same 
###log-likelihoods, for the same datasets.
#
###Ancestral state probabilities at nodes will be slightly 
###different, since BioGeoBEARS is reporting the 
###ancestral state probabilities under the global ML
###model, and Lagrange is reporting ancestral state
###probabilities after re-optimizing the likelihood
###after fixing the state at each node. These will 
###be similar, but not identical. See Matzke (2014),
###Systematic Biology, for discussion.
#
###Also see Matzke (2014) for presentation of the 
###DEC+J model.







###Run DEC


###Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

###Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

###Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

###Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001  ###Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE  ###set to FALSE for e.g. DEC* model, DEC*+J, etc.
###(For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
### Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
### the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
### DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
###Also: search script on "include_null_range" for other places to change

###Set up a time-stratified analysis:
###1. Here, un-comment ONLY the files you want to use.
###2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
###3. For example files see (a) extdata_dir, 
### or (b) http://phylo.wikidot.com/biogeobears#files
### and BioGeoBEARS Google Group posts for further hints)
#
###Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE  ###if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 6
###(use more cores to speed it up; this requires
###library(parallel) and/or library(snow). The package "parallel" 
###is now default on Macs in R 3.0+, but apparently still 
###has to be typed on some Windows machines. Note: apparently 
###parallel works on Mac command-line R, but not R.app.
###BioGeoBEARS checks for this and resets to 1
###core with R.app)

###Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
###I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
###but the results are imprecise and so I haven't explored it further.
###In a Bayesian analysis, it might work OK, but the ML point estimates are
###not identical.
###Also, I have not implemented all functions to work with force_sparse=TRUE.
###Volunteers are welcome to work on it!!
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

###Set up DEC model
###(nothing to do; defaults)

###Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

###This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

###This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

###Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#####These analyses spit out some R.data files so make a subdirectory od /Data to save them in
dir.create("Metadata/BioGeoBEARS_save_files", showWarnings = FALSE)

###For a slow analysis, run once, then set runslow=FALSE to just
###load the saved result.
runslow = TRUE
resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DEC.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
###Loads to "res"
  load(resfn)
  resDEC = res
}


###Run DEC+J

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
#BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
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

###Set up DEC+J model
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

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DEC+J.Rdata"
runslow = TRUE
if (runslow)
{

  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
###Loads to "res"
  load(resfn)
  resDECj = res
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
xstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DEC+X.Rdata"
runslow = TRUE
if (runslow)
{
  
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECx = res
} else {
  ###Loads to "res"
  load(resfn)
  resDECx = res
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




###DIVALIKE, DIVALIKE+J, DIVALIKE+X AND DIVALIKE+J+X ANALYSIS


###NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
###Ronquist (1997)'s parsimony DIVA. It is a likelihood
###interpretation of DIVA, constructed by modelling DIVA's
###processes the way DEC does, but only allowing the 
###processes DIVA allows (widespread vicariance: yes; subset
###sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
###DIVALIKE is a likelihood interpretation of parsimony
###DIVA, and it is "like DIVA" -- similar to, but not
###identical to, parsimony DIVA.
#
###I thus now call the model "DIVALIKE", and you should also. ;-)




###Run DIVALIKE

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
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

###Set up DIVALIKE model
###Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

###Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

###No jump dispersal/founder-event speciation
###BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
###BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
###BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DIVALIKE.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
###Loads to "res"
  load(resfn)
  resDIVALIKE = res
}


###Run DIVALIKE+J

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
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

###Set up DIVALIKE+J model
###Get the ML parameter values from the 2-parameter nested model
###(this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

###Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

###Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

###Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DIVALIKE+J.Rdata"
runslow = TRUE
if (runslow)
{

  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
###Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}


###Run DIVALIKE+X

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
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

###Set up DIVALIKE+X model
###Get the ML parameter values from the 2-parameter nested model
###(this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

###Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

###Add distance based penalties
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DIVALIKE+X.Rdata"
runslow = TRUE
if (runslow)
{
  
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEx = res
} else {
  ###Loads to "res"
  load(resfn)
  resDIVALIKEx = res
}


###Run DIVALIKE+J+X

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
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

###Set up DIVALIKE+J model
###Get the ML parameter values from the 2-parameter nested model
###(this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

###Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

###Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

###Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

###Add distance based penalties
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_DIVALIKE+J+X.Rdata"
runslow = TRUE
if (runslow)
{
  
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEjx = res
} else {
  ###Loads to "res"
  load(resfn)
  resDIVALIKEjx = res
}



###BAYAREALIKE, BAYAREALIKE+J, BAYAREALIKE+X AND BAYAREALIKE+J+X ANALYSIS


###NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
###not identical with the full Bayesian model implemented 
###in the "BayArea" program of Landis et al. (2013). 
#
###Instead, this is a simplified likelihood interpretation
###of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
###"d" and "e" work like they do in the DEC model of Lagrange 
###(and BioGeoBEARS), and then BayArea's cladogenesis assumption
###(which is that nothing in particular happens at cladogenesis) is 
###replicated by BioGeoBEARS.
#
###This leaves out 3 important things that are in BayArea:
###1. Distance dependence (you can add this with a distances 
###   matrix + the "x" parameter in BioGeoBEARS, however)
###2. A correction for disallowing "e" events that drive
###   a species extinct (a null geographic range)
###3. The neat Bayesian sampling of histories, which allows
###   analyses on large numbers of areas.
#
###The main purpose of having a "BAYAREALIKE" model is 
###to test the importance of the cladogenesis model on 
###particular datasets. Does it help or hurt the data 
###likelihood if there is no special cladogenesis process?
###
###BAYAREALIKE is a likelihood interpretation of BayArea,
###and it is "like BayArea" -- similar to, but not
###identical to, Bayesian BayArea.
###I thus now call the model "BAYAREALIKE", and you should also. ;-)




###Run BAYAREALIKE

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
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

###Set up BAYAREALIKE model
###No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

###No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

###No jump dispersal/founder-event speciation
###BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
###BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
###BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

###Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

###Only sympatric/range-copying (y) events allowed, and with 
###exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

###Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_BAYAREALIKE.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
###Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}


###Run BAYAREALIKE+J

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
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

###Set up BAYAREALIKE+J model
###Get the ML parameter values from the 2-parameter nested model
###(this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

###No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

###*DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

###Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

###Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

###Only sympatric/range-copying (y) events allowed, and with 
###exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

###NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
###machines. I can't replicate this on my Mac machines, but it is almost certainly
###just some precision under-run issue, when optim/optimx tries some parameter value 
###just below zero.  The "min" and "max" options on each parameter are supposed to
###prevent this, but apparently optim/optimx sometimes go slightly beyond 
###these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
###slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_BAYAREALIKE+J.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
###Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}


###Run BAYAREALIKE+X

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
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

###Set up BAYAREALIKE+X model
###Get the ML parameter values from the 2-parameter nested model
###(this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
xstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

###No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

###Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart

###Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

###Only sympatric/range-copying (y) events allowed, and with 
###exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

###NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
###machines. I can't replicate this on my Mac machines, but it is almost certainly
###just some precision under-run issue, when optim/optimx tries some parameter value 
###just below zero.  The "min" and "max" options on each parameter are supposed to
###prevent this, but apparently optim/optimx sometimes go slightly beyond 
###these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
###slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_BAYAREALIKE+X.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEx = res
} else {
  ###Loads to "res"
  load(resfn)
  resBAYAREALIKEx = res
}


###Run BAYAREALIKE+J+X

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
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
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

###Set up BAYAREALIKE+J+X model
###Get the ML parameter values from the 2-parameter nested model
###(this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

###Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

###No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

###No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

###*DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

###Under BAYAREALIKE+J+X, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

###Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = jstart

###Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

###Only sympatric/range-copying (y) events allowed, and with 
###exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

###NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
###machines. I can't replicate this on my Mac machines, but it is almost certainly
###just some precision under-run issue, when optim/optimx tries some parameter value 
###just below zero.  The "min" and "max" options on each parameter are supposed to
###prevent this, but apparently optim/optimx sometimes go slightly beyond 
###these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
###slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Metadata/BioGeoBEARS_save_files/Cratopus_BAYAREALIKE+J+X.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEjx = res
} else {
  ###Loads to "res"
  load(resfn)
  resBAYAREALIKEjx = res
}

###Plot the graphical output of all models as a set of pdfs

###PDF plots for DEC based models
pdffn = "Metadata/BioGeoBEARS_save_files/Cratopus_DEC_vs_DEC+J_vs_DEC+X_vs_DEC+J+X.pdf"
pdf(pdffn, width=10, height=16)

###DEC
###Setup
analysis_titletxt ="BioGeoBEARS DEC on Cratopus Beast MCC tree"
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
###States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)
###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###DECJ
analysis_titletxt ="BioGeoBEARS DEC+J on Cratopus Beast MCC tree"
###Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)
###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###DECX
analysis_titletxt ="BioGeoBEARS DEC+X on Cratopus Beast MCC tree"
###Setup
results_object = resDECx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)
###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###DECJX
analysis_titletxt ="BioGeoBEARS DEC+J+X on Cratopus Beast MCC tree"

###Setup
results_object = resDECjx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), 
                                plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7,
                                splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                                include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)


dev.off()###Turn off PDF




###PDF plots for DIVALIKE based models


pdffn = "Metadata/BioGeoBEARS_save_files/Cratopus_DIVALIKE_vs_DIVALIKE+J_vs_DIVALIKE+X_vs_DIVALIKE+J+X.pdf"
pdf(pdffn, width=6, height=12)

analysis_titletxt ="BioGeoBEARS DIVALIKE on Cratopus Beast MCC tree"

###Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)


###DIVALIKE+J

analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Cratopus Beast MCC tree"

###Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)


###DIVALIKE+X

analysis_titletxt ="BioGeoBEARS DIVALIKE+X on Cratopus Beast MCC tree"

###Setup
results_object = resDIVALIKEx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)


###DIVALIKE+J+X

analysis_titletxt ="BioGeoBEARS DIVALIKE+J+X on Cratopus Beast MCC tree"

###Setup
results_object = resDIVALIKEjx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

dev.off()


###PDF plots for BAYAREALIKE based models

pdffn = "Metadata/BioGeoBEARS_save_files/Cratopus_BAYAREALIKE_vs_BAYAREALIKE+J_vs_BAYAREALIKE+X_vs_BAYAREALIKE+J+X.pdf"
pdf(pdffn, width=6, height=12)

analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Cratopus Beast MCC tree"

###Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)


###BAYAREALIKE+J

analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Cratopus Beast MCC tree"

###Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)


###BAYAREALIKE+X

analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Cratopus Beast MCC tree"

###Setup
results_object = resBAYAREALIKEx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=my.tree, tipranges=tipranges)


###BAYAREALIKE+J+X

analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J+X on Cratopus Beast MCC tree"

###Setup
results_object = resBAYAREALIKEjx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

###States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), 
                                plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, 
                                splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                                include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

###Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), 
                         plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, 
                         splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                         include_null_range=TRUE, tr=my.tree, tipranges=tipranges)

dev.off()

###CALCULATE SUMMARY STATISTICS TO COMPARE
###DEC, DEC+J, DEC+X, DEC+J+X, DIVALIKE, DIVALIKE+J, DIVALIKE+X, DIVALIKE+J+X, BAYAREALIKE, BAYAREALIKE+J, BAYAREALIKE+X, BAYAREALIKE+J+X

###Set up empty results tables and extract data from results objects

###Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

###We have to extract the log-likelihood differently, depending on the 
###version of optim/optimx
LnL_DEC = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_DECJ = get_LnL_from_BioGeoBEARS_results_object(resDECj)
LnL_DECX = get_LnL_from_BioGeoBEARS_results_object(resDECx)
LnL_DECJX = get_LnL_from_BioGeoBEARS_results_object(resDECjx)

LnL_DIVALIKE = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_DIVALIKEJ = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
LnL_DIVALIKEX = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEx)
LnL_DIVALIKEJX = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEjx)

LnL_BAYAREALIKE = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_BAYAREALIKEJ = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
LnL_BAYAREALIKEX = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEx)
LnL_BAYAREALIKEJX = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEjx)


###Statistics -- DEC vs DEC+J vs DEC+X vs DEC+J+X


numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_DECJ, LnL_DEC, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_DECX, LnL_DEC, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_DECJX, LnL_DECX, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams0 = 4
numparams1 = 3
stats = AICstats_2models(LnL_DECJX, LnL_DECJ, numparams0, numparams1)
stats

###DEC, null model for Likelihood Ratio Test (LRT)
DECparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###DEC+J, alternative model for Likelihood Ratio Test (LRT)
DECJparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###DEC+X, alternative model for Likelihood Ratio Test (LRT)
DECXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDECx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###DEC+J+X, bonus model for additional Likelihood Ratio Tests (LRT)
DECJXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDECjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

###The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
###confer the same likelihood on the data. See: Brian O'Meara's webpage:
###http://www.brianomeara.info/tutorials/aic
###...for an intro to LRT, AIC, and AICc

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, DECparams, DECJparams, DECXparams, DECJXparams)
teststable = rbind(teststable, tmp_tests)


###Statistics -- DIVALIKE vs DIVALIKE+J vs DIVALIKE+X vs DIVALIKE+J+X


numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_DIVALIKEJ, LnL_DIVALIKE, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_DIVALIKEX, LnL_DIVALIKE, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_DIVALIKEJX, LnL_DIVALIKEX, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams0 = 4
numparams1 = 3
stats = AICstats_2models(LnL_DIVALIKEJX, LnL_DIVALIKEJ, numparams0, numparams1)
stats

###DIVALIKE, null model for Likelihood Ratio Test (LRT)
DIVALIKEparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
DIVALIKEJparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###DIVALIKE+X, alternative model for Likelihood Ratio Test (LRT)
DIVALIKEXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###DIVALIKE+J+X, bonus model for additional Likelihood Ratio Tests (LRT)
DIVALIKEJXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

###The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
###confer the same likelihood on the data. See: Brian O'Meara's webpage:
###http://www.brianomeara.info/tutorials/aic
###...for an intro to LRT, AIC, and AICc

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, DIVALIKEparams, DIVALIKEJparams, DIVALIKEXparams, DIVALIKEJXparams)
teststable = rbind(teststable, tmp_tests)


###Statistics -- BAYAREALIKE vs BAYAREALIKE+J vs BAYAREALIKE+X vs BAYAREALIKE+J+X

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_BAYAREALIKEJ, LnL_BAYAREALIKE, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_BAYAREALIKEX, LnL_BAYAREALIKE, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_BAYAREALIKEJX, LnL_BAYAREALIKEX, numparams1, numparams2)
stats

tmp_tests = conditional_format_table(stats)

teststable = rbind(teststable, tmp_tests)

numparams0 = 4
numparams1 = 3
stats = AICstats_2models(LnL_BAYAREALIKEJX, LnL_BAYAREALIKEJ, numparams0, numparams1)
stats

###BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
BAYAREALIKEparams = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
BAYAREALIKEJparams = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###BAYAREALIKE+X, alternative model for Likelihood Ratio Test (LRT)
BAYAREALIKEXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
###BAYAREALIKE+J+X, bonus model for additional Likelihood Ratio Tests (LRT)
BAYAREALIKEJXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

###The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
###confer the same likelihood on the data. See: Brian O'Meara's webpage:
###http://www.brianomeara.info/tutorials/aic
###...for an intro to LRT, AIC, and AICc

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, BAYAREALIKEparams, BAYAREALIKEJparams, BAYAREALIKEXparams, BAYAREALIKEJXparams)
teststable = rbind(teststable, tmp_tests)


###ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J

teststable$alt = c("DEC+J", "DEC+X", "DEC+J+X", "DEC+J+X", "DIVALIKE+J", "DIVALIKE+X", "DIVALIKE+J+X", "DIVALIKE+J+X","BAYAREALIKE+J", "BAYAREALIKE+X", "BAYAREALIKE+J+X", "BAYAREALIKE+J+X")
teststable$null = c("DEC", "DEC", "DEC+X", "DEC+J", "DIVALIKE", "DIVALIKE", "DIVALIKE+X", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE", "BAYAREALIKE+X", "BAYAREALIKE+J")
row.names(restable) = c("DEC", "DEC+J", "DEC+X", "DEC+J+X", "DIVALIKE", "DIVALIKE+J", "DIVALIKE+X", "DIVALIKE+J+X", "BAYAREALIKE", "BAYAREALIKE+J", "BAYAREALIKE+X", "BAYAREALIKE+J+X")
restable = put_jcol_after_ecol(restable)


###Save the results tables for later -- check for e.g.
###convergence issues


###Loads to "restable"
save(restable, file="Metadata/BioGeoBEARS_save_files/restable_v1.Rdata")
load(file="Metadata/BioGeoBEARS_save_files/restable_v1.Rdata")

###Loads to "teststable"
save(teststable, file="Metadata/BioGeoBEARS_save_files/teststable_v1.Rdata")
load(file="Metadata/BioGeoBEARS_save_files/teststable_v1.Rdata")

###Also save to text files
write.csv(restable, file="Metadata/BioGeoBEARS_save_files/restable.csv", quote=FALSE)
write.csv(unlist_df(teststable), file="Metadata/BioGeoBEARS_save_files/teststable.csv", quote=FALSE)


###Model weights of all six models

restable2 = restable

###With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

###With AICcs -- factors in sample size
samplesize = length(my.tree$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike


###Calculating relative model probabilities for AIC and AICc model selection

restable_AIC_rellike$deltaAIC<-min(restable_AIC_rellike$AIC)-restable_AIC_rellike$AIC
restable_AIC_rellike$rel_mod_prob<-exp(restable_AIC_rellike$deltaAIC/2)
#Order the table by decreasing model probabilility
restable_AIC_rellike<-restable_AIC_rellike[with(restable_AIC_rellike, order(-rel_mod_prob)), ]

restable_AICc_rellike$deltaAICc<-min(restable_AICc_rellike$AICc)-restable_AICc_rellike$AICc
restable_AICc_rellike$rel_mod_prob<-exp(restable_AICc_rellike$deltaAICc/2)
#Order the table by decreasing model probabilility
restable_AICc_rellike<-restable_AICc_rellike[with(restable_AICc_rellike, order(-rel_mod_prob)), ]

###Export tables of results

###Also save to text files
write.csv(restable_AIC_rellike, file="Metadata/BioGeoBEARS_save_files/restable_AIC_rellike.csv", quote=FALSE)
write.csv(restable_AICc_rellike, file="Metadata/BioGeoBEARS_save_files/restable_AICc_rellike.csv", quote=FALSE)

###Save with nice conditional formatting
write.csv(conditional_format_table(restable_AIC_rellike), file="Metadata/BioGeoBEARS_save_files/restable_AIC_rellike_formatted.csv", quote=FALSE)
write.csv(conditional_format_table(restable_AICc_rellike), file="Metadata/BioGeoBEARS_save_files/restable_AICc_rellike_formatted.csv", quote=FALSE)

###print all the results to screen for convenience
teststable
restable_AIC_rellike
restable_AICc_rellike
