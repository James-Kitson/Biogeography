
### Clear the workspace
rm(list=ls())


# Load the package (after installation, see above).
library(optimx)         # You need to have some version of optimx available
# as it is a BioGeoBEARS dependency; however, if you
# don't want to use optimx, and use optim() (from R core) 
# you can set:
# BioGeoBEARS_run_object$use_optimx = FALSE
# ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
library(stringr)

########################################################
# TO GET THE OPTIMX/OPTIM FIX, AND THE UPPASS FIX, 
# SOURCE THE REVISED FUNCTIONS WITH THESE COMMANDS
#
# CRUCIAL CRUCIAL CRUCIAL: 
# YOU HAVE TO RUN THE SOURCE COMMANDS AFTER 
# *EVERY TIME* YOU DO library(BioGeoBEARS). THE CHANGES ARE NOT "PERMANENT", 
# THEY HAVE TO BE MADE EACH TIME.  IF YOU ARE GOING TO BE OFFLINE, 
# YOU CAN DOWNLOAD EACH .R FILE TO YOUR HARD DRIVE AND REFER THE source()
# COMMANDS TO THE FULL PATH AND FILENAME OF EACH FILE ON YOUR
# LOCAL SYSTEM INSTEAD.
########################################################
library(BioGeoBEARS)
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
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully

#######################################################
# Local source()-ing method -- uses BioGeoBEARS sourceall() function 
# on a directory of .R files, so you don't have to type them out.
# The directories here are on my machine, you would have to make a 
# directory, save the .R files there, and refer to them.
#
# NOTE: it's best to source the "cladoRcpp.R" update first, to avoid warnings like this:
##
## Note: possible error in 'rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs = tmpca_1, ': 
##         unused arguments (m = m, m_null_range = include_null_range, jts_matrix = jts_matrix) 
##
#
# TO USE: Delete or comment out the 'source("http://...")' commands above, and un-comment
#              the below...
########################################################################
# Un-comment (and fix directory paths) to use:
#library(BioGeoBEARS)
#source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
#sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
#calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
#calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
########################################################################

#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files"

#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory 
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# 
# Each computer operating system might install BioGeoBEARS in a different place, 
# depending on your OS and settings. 
# 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# Time-stratified files are here
#extdata_dir2 = np(slashslash(paste(extdata_dir, "examples/Psychotria_M3strat/BGB/", sep="/")))

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

##########PREPARE NEWICK FORMATTED TREE FROM NEXUS
### read in all the trees from BEAST that have had the burnin removed and the branch lengths rescaled to mya.
#my.nexus.trees<-list.files("../Trees/", pattern = "*dated*")

### set how many trees per BEAST replicated we want trees_per_file=10 = 10*10 trees to analyse
#trees_per_file<-10

### Create some output folders
dir.create("Output", showWarnings = FALSE)
dir.create("BioGeoBEARs_save_files_final_model", showWarnings = FALSE)
### Create a place to store individual trees incase we want to plot them later
dir.create("BioGeoBEARs_save_files_final_model/Working_trees")

###Run DEC+J+X

###Set the maximum number of areas any species may occupy; this cannot be larger
###than the number of areas you set up, but it can be smaller.
max_range_size = 2

### set up a dataframe to fill with colonisation ages
colonisation_means<-NULL

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################

### empty the colonisation means so that you don't just keep adding analyses together if you rerun the script
colonisation_means<-NULL

trfn<-"BioGeoBEARs_save_files_final_model/Working_trees/Relaxed_clock_MCC_newick.tree"
my.tree<-read.tree(trfn)

# Example files are given below. To run your own data,
# make the below lines point to your own files, e.g.
geogfn = "../BioGeoBEARS_islands.txt"
### provide a path to the distance matrix
dist.mat<-"../BioGeoBEARS_dist_scale.txt"
### provide a path to the time periods file
#eras<-"../timeperiods_BEAST.txt"
### provide a path to the allowed areas file
#allowed<-"../areas_allowed.txt"

tr=read.tree(trfn)
#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################
# This is the example Newick file for Hawaiian JuliaDupin
# (from Ree & Smith 2008)

#######################################################
# Geography file
# Notes:
# 1. This is a PHLYIP-formatted file. This means that in the 
#    first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
# 2. This is the same format used for C++ LAGRANGE geography files.
# 3. All names in the geography file must match names in the phylogeny file.
# 4. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
# 5. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 2

#######################################################
# Run DEC+J+X
#######################################################


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
#BioGeoBEARS_run_object$timesfn = eras
#BioGeoBEARS_run_object$dispersal_multipliers_fn = dist.mat
#BioGeoBEARS_run_object$areas_allowed_fn = allowed
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = dist.mat
###See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

###Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = FALSE        ###shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE  ###if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 28
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
dstart = 0.0001
estart = 0.0001
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
# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = paste("BioGeoBEARs_save_files_final_model/Cratopus_DEC+J+X_BSM.Rdata", sep="")
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
} else {
  # Loads to "res"
  load(resfn)
}

################################################################################
### This first section performs the DEC+J+X analysis on a subsample of trees ###
### We will use the save objects generated here to perform the BSM analyses ####
################################################################################

restable = NULL

pdffn = "Output/Cratopus_DEC+J+X.pdf"
pdf(pdffn, width=10, height=10)

    #my.results.objects<-list.files("BioGeoBEARs_save_files_final_model/")
    #for(l in my.results.objects){
    #  resfn<-paste("BioGeoBEARs_save_files_final_model/",l, sep="")
    #  load(resfn)
    resDECjx = res
    DECJXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDECjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
    restable = rbind(restable, DECJXparams)
    
    ###Setup
    results_object = resDECjx
    scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
    
    analysis_titletxt ="BioGeoBEARS DEC+J+X on Cratopus BEAST tree"
    
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
    

  dev.off()###Turn off PDF
  
  ####################################################################################
  ############ Extract the colonisation nodes from the BioGeoBEARS analysis ##########
  ################ and write their ages to a file for later analysis #################
  ####################################################################################
  
  ### first of all you need to generate a list of area states to annotate the table generated by prt()
  # Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
  max_range_size = 2
  ### this must be a range of the number of areas you have in your data set - i.e. not STATES
  areas = LETTERS[1:12]
  
  # This is the list of states/ranges, where each state/range
  # is a list of areas, counting from 0
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
  
  # Make the list of ranges
  ranges_list = NULL
  for (i in 1:length(states_list_0based))
  {    
    if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
    {
      tmprange = "_"
    } else {
      tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
    }
    ranges_list = c(ranges_list, tmprange)
  }
  
  # Look at the ranges list
  ranges_list
  
  ####################################################################################
  ################ This bit refers to the standard analysis not BSM #####################
  ####################################################################################
  AT_node<-results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
  max<- apply(AT_node, 1, which.max)
  
  ### annotate the above data onto the tree _summary from prt() using a match function
  tree_summary<-prt(my.tree)
  tree_summary$most_prob_state<-max
  islands<-as.data.frame(cbind(seq(1,length(ranges_list),1),ranges_list),stringsAsFactors = FALSE)
  colnames(islands)<-c("State","Island")
  tree_summary$most_prob_island<-islands$Island[match(tree_summary$most_prob_state,islands$State)]
  
  #########################
  ###extract most probable state changes
  ### for each node, find it's ancestor and it's ancestors most probable state, write that to a data frame
  
  tree_iter="MCC_tree"
  
  state_changes<-NULL
  for(i in 1:nrow(tree_summary)){
    daughter_node<-tree_summary$node[i]
    ancestor<-tree_summary$ancestor[i]
    daughter_node_island<-tree_summary$most_prob_island[i]
    daughter_node_age<-tree_summary$time_bp[i]
    daughter_node_label<-tree_summary$label[i]
    new_line<-cbind(ancestor,daughter_node,daughter_node_island,daughter_node_age,daughter_node_label)
    state_changes<-rbind(state_changes,new_line)
  }
  
  state_changes<-as.data.frame(state_changes, stringsAsFactors=FALSE)
  state_changes$daughter_node_age<-as.numeric(state_changes$daughter_node_age)
  state_changes$ancestor_node_island<-tree_summary$most_prob_island[match(state_changes$ancestor,tree_summary$node)]
  state_changes$ancestor_node_age<-tree_summary$time_bp[match(state_changes$ancestor,tree_summary$node)]
  
  # mean age of colonisations from B->C i.e. Reunion to Mauritius
  Reu_Mau<-mean(as.numeric(state_changes$ancestor_node_age[which(state_changes$ancestor_node_island=="B" & state_changes$daughter_node_island=="C")]))
  # mean age of colonisations from C->B i.e. Mauritius to Reunion
  Mau_Reu<-mean(as.numeric(state_changes$ancestor_node_age[which(state_changes$ancestor_node_island=="C" & state_changes$daughter_node_island=="B")]))
  
  # count the number of colonisations from B->C and C->B
  Reu_Mau_count<-length(which(state_changes$ancestor_node_island=="B" & state_changes$daughter_node_island=="C"))
  Mau_Reu_count<-length(which(state_changes$ancestor_node_island=="C" & state_changes$daughter_node_island=="B"))
  
  # List them for convenience
  Reu_Mau_list<-which(state_changes$ancestor_node_island=="B" & state_changes$daughter_node_island=="C")
  Mau_Reu_list<-which(state_changes$ancestor_node_island=="C" & state_changes$daughter_node_island=="B")
  
  Reu_Mau_list
  Mau_Reu_list
  
  # mean age of in situ speciation events on B (Reunion)
  Reu_insitu<-mean(as.numeric(state_changes$ancestor_node_age[which(state_changes$ancestor_node_island=="B" & state_changes$daughter_node_island=="B")]))
  # mean age of in situ speciation events on C (Mauritius)
  Mau_insitu<-mean(as.numeric(state_changes$ancestor_node_age[which(state_changes$ancestor_node_island=="C" & state_changes$daughter_node_island=="C")]))
  
  ### Make a vector of the flightless taxa to pull out
  fl_tips<-c("E07_Mau3707","B10_Reu2320","Reu_mur","E09_Reu268","A06_Mau3268","A11_Rod56")
  ### extract the ancestor nodes for each flightles tip, we know that one pair are sister species so reduce this to unique values
  fl_tip_ancestors<-unique(tree_summary$ancestor[tree_summary$label %in% fl_tips])
  ### calculate the mean age of all these flight loss ancestors, all ancestor nodes are represented twice in tree_summary as they are the ancestor of two nodes, so again we take unique values
  fl_age<-mean(unique(tree_summary$time_bp[tree_summary$node %in% fl_tip_ancestors]))
  ### extract all other node ages
  all_other_nodes<-unique(tree_summary$ancestor[!(tree_summary$label %in% fl_tips)])
  ### calculate the mean age of all the other nodes, as above, reduce this to unique values as each node is represented twice
  not_fl_age<-mean(unique(tree_summary$time_bp[tree_summary$node %in% all_other_nodes]))
  
  # put the values above into a data frame with the tree number
  colonisation_means<-as.data.frame(cbind(tree_iter,Reu_Mau,Mau_Reu,Reu_insitu,Mau_insitu,fl_age,not_fl_age,Reu_Mau_count,Mau_Reu_count))
  
  ### write the data out to a csv for convenience
  write.csv(colonisation_means, paste(file ="Output/colonisation_means_unstratified_notBSM.csv",sep=""))  
  
  
  ################################################################################
  ########### This mark starts the BSM analysis the MCC tree #####
  ################################################################################
  
  ###########################################
  # Pick your model name:
  ###########################################
  model_name = "DEC+J+X"
  
  analysis_titletxt = paste0(model_name, " on Cratopus")
  
  # Setup
  results_object = res
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  dev.off()  # Turn off PDF
  
  #######################################################
  # Stochastic mapping on DEC M3b stratified with islands coming up
  #######################################################
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  #######################################################
  # Get the inputs for Biogeographical Stochastic Mapping
  # Note: this can be slow for large state spaces and trees, since 
  # the independent likelihoods for each branch are being pre-calculated
  # E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
  # for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
  # for storage of "BSM_inputs_file.Rdata".
  # Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
  # the same settings will be used for get_inputs_for_stochastic_mapping().
  #######################################################
  BSM_inputs_fn = "BioGeoBEARs_save_files_final_model/BSM_inputs_file.Rdata"
  runInputsSlow = TRUE
  if (runInputsSlow)
  {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr)
  names(stochastic_mapping_inputs_list)
  stochastic_mapping_inputs_list$phy2
  stochastic_mapping_inputs_list$COO_weights_columnar
  stochastic_mapping_inputs_list$unconstr
  set.seed(seed=as.numeric(Sys.time()))
  
  runBSMslow = TRUE
  if (runBSMslow == TRUE)
  {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=1000, nummaps_goal=500, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
  } else {
    # Load previously saved...
    
    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
  } # END if (runBSMslow == TRUE)
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  head(clado_events_tables[[1]])
  head(ana_events_tables[[1]])
  length(clado_events_tables)
  length(ana_events_tables)
  
  #######################################################
  # Plot one stochastic map, manual method
  #######################################################
  # (we have to convert the stochastic maps into event
  #  maps for plotting)
  
  ######################
  # Get the color scheme
  ######################
  include_null_range = TRUE
  areanames = names(tipranges@df)
  areas = areanames
  max_range_size = 2
  
  # Note: If you did something to change the states_list from the default given the number of areas, you would
  # have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  
  ############################################
  # Setup for painting a single stochastic map
  ############################################
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified=FALSE
  clado_events_table = clado_events_tables[[1]]
  ana_events_table = ana_events_tables[[1]]
  
  # cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
  # colnums = match(cols_to_get, names(ana_events_table))
  # ana_events_table_cols_to_add = ana_events_table[,colnums]
  # anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
  # ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
  # rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
  # master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)
  
  ############################################
  # Open a PDF
  ############################################
  pdffn = paste0("Output/",model_name, "_single_stochastic_map_n1.pdf")
  pdf(file=pdffn, width=6, height=6)
  
  # Convert the BSM into a modified res object
  master_table_cladogenetic_events = clado_events_tables[[1]]
  resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)
  
  plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
  
  # Paint on the branch states
  paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)
  
  plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
  
  ############################################
  # Close PDF
  ############################################
  dev.off()
  #cmdstr = paste("open ", pdffn, sep="")
  #system(cmdstr)
  
  #######################################################
  # Plot all 50 stochastic maps to PDF
  #######################################################
  # Setup
  include_null_range = include_null_range
  areanames = areanames
  areas = areanames
  max_range_size = max_range_size
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified = stratified
  
  # Loop through the maps and plot to PDF
  pdffn = paste0("Output/",model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
  pdf(file=pdffn, width=6, height=6)
  
  nummaps_goal = 50
  for (i in 1:nummaps_goal)
  {
    clado_events_table = clado_events_tables[[i]]
    analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
    plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
  } # END for (i in 1:nummaps_goal)
  
  dev.off()
  #cmdstr = paste("open ", pdffn, sep="")
  #system(cmdstr)
  
  #######################################################
  # Summarize stochastic map tables
  #######################################################
  length(clado_events_tables)
  length(ana_events_tables)
  
  head(clado_events_tables[[1]][,-20])
  tail(clado_events_tables[[1]][,-20])
  
  head(ana_events_tables[[1]])
  tail(ana_events_tables[[1]])
  
  areanames = names(tipranges@df)
  actual_names = areanames
  actual_names
  
  # Get the dmat and times (if any)
  dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
  dmat_times
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  
  # Simulate the source areas
  BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
  clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
  ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
  
  # Count all anagenetic and cladogenetic events
  counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)
  
  summary_counts_BSMs = counts_list$summary_counts_BSMs
  print(conditional_format_table(summary_counts_BSMs))
  
  # Histogram of event counts
  hist_event_counts(counts_list, pdffn=paste0("Output/",model_name, "_histograms_of_event_counts.pdf"))
  
  #######################################################
  # Print counts to files
  #######################################################
  tmpnames = names(counts_list)
  cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
  for (i in 1:length(tmpnames))
  {
    cmdtxt = paste0("item = counts_list$", tmpnames[i])
    eval(parse(text=cmdtxt))
    
    # Skip cubes
    if (length(dim(item)) != 2)
    {
      next()
    }
    
    outfn = paste0(tmpnames[i], ".txt")
    if (length(item) == 0)
    {
      cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
      cat("\n")
    } else {
      cat(outfn)
      cat("\n")
      write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    } # END if (length(item) == 0)
  } # END for (i in 1:length(tmpnames))
  cat("...done.\n")
  
  #######################################################
  # Check that ML ancestral state/range probabilities and
  # the mean of the BSMs approximately line up
  #######################################################
  #library(MultinomialCI)    # For 95% CIs on BSM counts
  #check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)
  
  ########################################################################################################################################
  ### BSM is complete, the next stage extracts the colonisation events and writes them to a file for analysis
  ########################################################################################################################################
  ### For each BSM build a data table of the columns we want
  for(i in 1:length(clado_events_tables)){
    BSM_iter=i
    ### Pull ages of colonisation nodes from each BSM
    BSM_nodes<-as.data.frame(cbind(clado_events_tables[[i]]$node,
                                   clado_events_tables[[i]]$clado_event_txt,
                                   rep(NA, nrow(clado_events_tables[[i]])),
                                   clado_events_tables[[i]]$time_bp,
                                   clado_events_tables[[i]]$node.type))
    BSM_nodes[,4]<-as.numeric(as.character(BSM_nodes[,4]))
    
    ### Drop the tip data as we need internal nodes
    processed_nodes<-subset(BSM_nodes, V5 != "tip")
    
    # split the vector describing the events into before and after vectors
    processed_nodes[,c(2,3)]<-str_split_fixed(processed_nodes$V2, "->", 2)
    
    # mean age of colonisations from B->C
    Reu_Mau<-mean(as.numeric(processed_nodes$V4[which(processed_nodes$V2=="B" & grepl("C",processed_nodes$V3))]))
    # mean age of colonisations from C->B
    Mau_Reu<-mean(as.numeric(processed_nodes$V4[which(processed_nodes$V2=="C" & grepl("B",processed_nodes$V3))]))
    # mean age of in situ speciation on B
    Reu_insitu<-mean(as.numeric(processed_nodes$V4[which(processed_nodes$V2=="B" & grepl("B,B",processed_nodes$V3))]))
    # mean age of in situ speciation on C
    Mau_insitu<-mean(as.numeric(processed_nodes$V4[which(processed_nodes$V2=="C" & grepl("C,C",processed_nodes$V3))]))
    
    ### Make a vector of the flightless taxa to pull out
    fl_tips<-c("E07_Mau3707","B10_Reu2320","Reu_mur","E09_Reu268","A06_Mau3268","A11_Rod56")
    ### extract the ancestor nodes for each flightles tip, we know that one pair are sister species so reduce this to unique values
    fl_tip_ancestors<-unique(clado_events_table$ancestor[clado_events_table$label %in% fl_tips])
    ### calculate the mean age of all these flight loss ancestors, all ancestor nodes are represented twice in clado_events_table as they are the ancestor of two nodes, so again we take unique values
    fl_age<-mean(unique(clado_events_table$time_bp[clado_events_table$node %in% fl_tip_ancestors]))
    ### extract all other node ages
    all_other_nodes<-unique(clado_events_table$ancestor[!(clado_events_table$label %in% fl_tips)])
    ### calculate the mean age of all the other nodes, as above, reduce this to unique values as each node is represented twice
    not_fl_age<-mean(unique(clado_events_table$time_bp[clado_events_table$node %in% all_other_nodes]))
    
    # put the values above into a data frame with the BSM number
    next_line<-cbind(BSM_iter,Reu_Mau,Mau_Reu,Reu_insitu,Mau_insitu,fl_age,not_fl_age)
    
    # add the next line to the end of the colonisation means
    colonisation_means<-rbind(colonisation_means, next_line)
    
    ### make the colonisation_means list into a dataframe
    colonisation_means<-as.data.frame(colonisation_means)
    
    ### write the data out to a csv for convenience
    write.csv(colonisation_means, file = paste("Output/colonisation_means_Unstratified_BSM.csv", sep=""))
  }
  
  ### test whether colonisations of Mauritius are significantly younger than colonisations of Reunion
  t.test(colonisation_means$Reu_Mau,colonisation_means$Mau_Reu, paired=TRUE)
  
  ### test whether colonisation of Mauritius are younger than in situ speciation on Mauritius
  t.test(colonisation_means$Reu_Mau,colonisation_means$Mau_insitu, paired=TRUE)
  
  ### test whether in situ speciation events on Reunion are significantly younger than colonisations of Reunion
  t.test(colonisation_means$Reu_insitu,colonisation_means$Mau_Reu, paired=TRUE)
  
  ### test whether flightloss nodes are yonger than other nodes
  t.test(colonisation_means$fl_age,colonisation_means$not_fl_age, paired=TRUE)
  
