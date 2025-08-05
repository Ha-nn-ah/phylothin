## subsampling-PhyloThin: R-script for subsampling large phylogenetic trees and applying PhyloThin

# by Hannah Götsch

# Compile this code using:
# Rscript phylothin_subsampling.r path_to_folder input_tree (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)

#### TUNING PARAMETER #############################################################################################

# !! TO DO: some theoretical computations (add them to "Supplementary Note 1: PhyloThin - Implementation" of manuscript)

# proportion on how often a sample has to be in a oversampling-clade such that it gets classified as oversampled
alpha_5 <- 0.9 # !! TO DO: what is the right threshold here?

#### INPUT ########################################################################################################

library(ape) # package for the analysis of phylogenetics and evolution

# extract command-line arguments:
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 7) {
  stop("Wrong command line input\n
       Usage is \"Rscript phylothin_subsampling.r path_to_folder input_tree 
       (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)\"", call.=FALSE)
}
default <- 1
runs <- NULL
subsamplesize <- NULL
basepath <- NULL
input_tree_file <- NULL
pathd8 <- 1
i <- 1
while (i <= length(args)) {
  if (args[i] == "-r") {
    runs <- as.numeric(args[i + 1]) # number of subsamples
    default <- 0
    i <- i + 2
  } else if (args[i] == "-s") {
    subsamplesize <- as.numeric(args[i + 1]) # size of subsamples
    default <- 0
    i <- i + 2
  } else if (args[i] == "no_PATHd8") {
    pathd8 <- 0 # skip PATHd8
    i <- i + 1
  } else {
    if (is.null(basepath)) { # first positional argument assumed to be basepath
      basepath <- args[i]
      if (substr(basepath, nchar(basepath), nchar(basepath)) == "/") { # path should NOT end with /
        basepath <- substr(basepath, 1, nchar(basepath)-1)
      }
    } else if (is.null(input_tree_file)){ # second positional argument assumed to be input tree
      input_tree_file <- args[i]
      if (substr(input_tree_file, nchar(input_tree_file)-3, nchar(input_tree_file)) == ".nwk") {
        tree_name <- substr(input_tree_file, 1, nchar(input_tree_file)-4)
      } else {
        tree_name <- input_tree_file
      }
    }
    i <- i + 1
  }
}

# sanity check:
if (is.null(basepath)) {
  stop("You need to specify the path to the data folder\n
       Usage is \"Rscript phylothin_subsampling.r path_to_folder input_tree 
       (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)\"", call.=FALSE)
} 
if (is.null(input_tree_file)) {
  stop("Missing argument input_tree\n
       Usage is \"Rscript phylothin_subsampling.r path_to_folder input_tree 
       (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)\"", call.=FALSE)
} 

# load (ultrametric) tree:
input_tree <- read.tree(file.path(basepath, input_tree_file))

# sanity check:
if (input_tree$Nnode < 2){ # "cintervals" can not be computed
  stop("Your tree has only one node. PhyloThin can not be used.")
}

num_sample <- length(um_tree$tip.label) # sample size
# default setting !! TO DO: runs and subsamplesize should depend on sample size
if (is.null(subsamplesize)) {subsamplesize <- 10}
if (is.null(runs)) {runs <- 10}

if (default == 1) {
  print(paste("The default parameter setting is used for subsampling:", 
              runs, "drawings with a subsample size of", subsamplesize, "."))
} else {
  print(paste("The following parameter setting has been specified for subsampling:", 
              runs, "drawings with a subsample size of", subsamplesize, "."))
  # sanity check:
  if (subsamplesize > num_sample){
    stop("Error: subsample size has to be smaller than sample size.")
  }
}

if (!file.exists(file.path(basepath, "phylothinoutput"))){ # output-folder for PhyloThin
  dir.create(file.path(basepath, "phylothinoutput"))
}
if (!file.exists(file.path(basepath, "phylothinoutput/subsampling"))){ # folder for sub-sampling output
  dir.create(file.path(basepath, "phylothinoutput/subsampling"))
}

if (pathd8==0){ # skip PATHd8 since requested by input
  print("Making the given tree ultrametric with PATHd8 is skipped.")
  um_tree <- input_tree
} else{ # make ultrametric tree with PATHd8 (Britton et al 2007)
  print("start calculating ultrametric tree with PATHd8.")
  # define paths needed in bash-script
  Sys.setenv(PATHd8 = paste(basepath, "/PATHd8", sep = ""))
  Sys.setenv(tree_file = file.path(basepath, input_tree_file))
  Sys.setenv(PATHd8_tree_file = paste(basepath, "/pathd8_", input_tree_file, sep = ""))
  Sys.setenv(um_tree_file = paste(basepath, "/um_", input_tree_file, sep = ""))
  # bash-script for PATHd8
  system('chmod u+x $tree_file') # to get permission
  system('$PATHd8 $tree_file $PATHd8_tree_file')
  system('sed -n "/d8.*;/p" $PATHd8_tree_file > $um_tree_file')
  # ultrametric tree
  um_tree <- read.tree(paste(basepath, "/um_", input_tree_file, sep = ""))
  if(class(um_tree)=="multiPhylo"){um_tree <- um_tree$`d8tree:`}
}

###################################################################################################################
#### SUBSAMPLING ##################################################################################################
###################################################################################################################

subsample_table <- data.frame(sample = um_tree$tip.label, 
                              num_sampled = rep(0, num_sample), 
                              num_in_cluster = rep(0, num_sample))
write.csv(subsample_table, file = paste(basepath, "/phylothinoutput/subsample_table.csv", sep = "" ),
          quote = F, row.names = F)

for (i in 1:runs) {
  print(paste("%%%%%%% START Subsampling", i, "%%%%%%%"))
  set.seed(i)
  nosample_tips <- um_tree$tip.label[sample(1:num_sample, num_sample-subsamplesize)]
  sample_tree <- drop.tip(um_tree, nosample_tips)
  write.tree(sample_tree, file = paste(basepath, "/phylothinoutput/subsampling/um_", tree_name, "_", i, ".nwk", sep = ""))

###################################################################################################################
#### PHYLOTHIN ####################################################################################################
###################################################################################################################
  
  Sys.setenv(phylothin = paste(basepath, "/phylothin.r", sep = ""))
  Sys.setenv(path_to_folder = paste(basepath, "/phylothinoutput/subsampling", sep = ""))
  Sys.setenv(input_tree = paste("um_", tree_name, "_", i, ".nwk", sep = ""))
  system('Rscript $phylothin $path_to_folder $input_tree no_PATHd8 no_clade')

###################################################################################################################
#### OUTPUT #######################################################################################################
###################################################################################################################

  subsample_table <- read.csv(paste(basepath, "/phylothinoutput/subsample_table.csv", sep = "" ), header=T)
  for (subsample_sample in setdiff(um_tree$tip.label, nosample_tips)) {
    subsample_table[subsample_table$sample == subsample_sample,]$num_sampled <- 1 + 
      subsample_table[subsample_table$sample == subsample_sample,]$num_sampled
  }
  subtree_name <- paste("um_", tree_name, "_", i, sep = "")
  if (file.exists(paste(basepath, "/phylothinoutput/subsampling/clades_", subtree_name, ".csv", sep = "" ))) {
    clades_file <- read.csv(paste(basepath, "/phylothinoutput/subsampling/clades_", subtree_name, ".csv", sep = "" ),
                            header=T)
    clades_tips <- clades_file[!is.na(clades_file$clade),]$samples
    for (clades_sample in clades_tips) {
     subsample_table[subsample_table$sample == clades_sample,]$num_in_cluster <- 1 + 
        subsample_table[subsample_table$sample == clades_sample,]$num_in_cluster
    }
  }
  write.csv(subsample_table, file = paste(basepath, "/phylothinoutput/subsample_table.csv", sep = "" ),
            quote = F, row.names = F)
}

subsample_table <- read.csv(paste(basepath, "/phylothinoutput/subsample_table.csv", sep = "" ), header=T)
subsample_table$oversampled <- subsample_table$num_in_cluster/subsample_table$num_sampled
write.csv(subsample_table, file = paste(basepath, "/phylothinoutput/subsample_table.csv", sep = "" ),
          quote = F, row.names = F)

# samples classified as oversampled
removed_labels <- subsample_table[subsample_table$oversampled >= alpha_5,]$sample

# removed sample-ids
if(length(removed_labels)>0){
  write.table(removed_labels, 
              file = paste(basepath, "/phylothinoutput/removed_ids_", tree_name, "_subsampling.txt" ,sep = "" ),
              quote = F, row.names = F, col.names = F)
}

# kept sample-ids
write.table(setdiff(um_tree$tip.label,removed_labels), 
            file = paste(basepath, "/phylothinoutput/kept_ids_", tree_name, "_subsampling.txt" ,sep = "" ),
            quote = F, row.names = F, col.names = F)

# tree with only kept samples
dropped_tree <- drop.tip(input_tree, tip = removed_labels) # tree where tips removed
write.tree(dropped_tree, 
           file = paste(basepath, "/phylothinoutput/reduced_tree_", tree_name, "_subsampling.nwk", sep = ""))

###################################################################################################################