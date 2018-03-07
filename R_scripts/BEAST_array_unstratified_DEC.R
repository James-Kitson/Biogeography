
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

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

####################################################################################
##################### Start setting up the analysis space ##########################
####################################################################################

### Create some output folders, usually easier to create these manually first
dir.create("Output", showWarnings = FALSE)
dir.create("BioGeoBEARS_save_files_final_model", showWarnings = FALSE)
### Create a place to store individual trees incase we want to plot them later
dir.create("BioGeoBEARS_save_files_final_model/Working_trees", showWarnings = FALSE)

#This line ensures that the intermediate files written below don't simply grow in size each time you run the analysis but are instead created from scratch each time.
file.remove("Output/colonisation_means.csv")


###Run DEC+J+X

### read the  list of trees back in as we'll iterate across this in a bit
my.newick.trees<-read.tree(file="BioGeoBEARS_save_files_final_model/Working_trees/BEAST_newick.trees")

#####provide a path to the geography file
geogfn <- "../BioGeoBEARS_islands.txt"
###Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
### provide a path to the distance matrix
dist.mat<-"../BioGeoBEARS_dist_scale.txt"
### provide a path to the time periods file
#eras<-"../timeperiods_BEAST.txt"
### provide a path to the allowed areas file
#allowed<-"../areas_allowed.txt"

###Set the maximum number of areas any species may occupy; this cannot be larger
###than the number of areas you set up, but it can be smaller.
max_range_size = 2

### set up a dataframe to fill with colonisation ages
colonisation_means<-NULL

################################################################################
###### This  section performs the DEC+J+X analysis on a subsample of trees #####
### We will use the save objects generated here to perform the BSM analyses ####
################################################################################

restable = NULL

### normally you would iterate across a series of values of k (as in the hashed
### out code below) but for this analysis we will set all analyses running simultaneously 
### on the cluster using an array job. This requires a recognisable tag ("array_tag") in the text for a 
### find and replace to change for each job. In this case, sbatch --array=1-100 will set 100 jobs running.
### in each case k will be set to a different number. Hash out the line below and unhash the loop to run normally.
k= array_tag

### create a tree tag to label the output .csv later
tree_iter=k

#pdffn = paste("Output/Cratopus_DEC+J+X",k,".pdf",sep="")
#pdf(pdffn, width=10, height=10)

### for each tree in my.newick.trees, perform the DEC+J+X analysis - unhash this loop (both ends) to run the analysis outside a cluster
### THIS WILL TAKE A LONG TIME!!!!

  #for(k in 1:100){
    
    write.tree(my.newick.trees[k], file=paste("BioGeoBEARS_save_files_final_model/Working_trees/Tree_",k,".tree",sep=""), append = FALSE)
    
    trfn<-paste("BioGeoBEARS_save_files_final_model/Working_trees/Tree_",k,".tree",sep="")
    my.tree<-read.tree(trfn)
    
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
    #BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    #BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    #BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    ###Add x as a free parameter
    #BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
    #BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = jstart
    #BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = jstart
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste("BioGeoBEARS_save_files_final_model/Cratopus_DEC+J+X_",k,".Rdata", sep="")
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
    
    DECJXparams = extract_params_from_BioGeoBEARS_results_object(results_object=resDECjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
    restable = rbind(restable, DECJXparams)
    
    ###Setup
    results_object = resDECjx
    scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
    
    analysis_titletxt ="BioGeoBEARS DEC+J+X on Cratopus BEAST tree"
    
    ###States
  #  res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), 
     #                               plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7,
    #                                splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
   #                                 include_null_range=TRUE, tr=my.tree, tipranges=tipranges)
    
    ###Pie chart
  #  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")),
      #                       plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7,
       #                      splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
        #                     include_null_range=TRUE, tr=my.tree, tipranges=tipranges)
  
    
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
    ################ Identify the most probable state at each node #####################
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
    write.csv(colonisation_means, paste(file ="Output/colonisation_means_unstratified_",k,".csv",sep=""))  
# }
# dev.off() ###Turn off PDF

