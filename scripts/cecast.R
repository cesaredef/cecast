#!/usr/bin/env Rscript

## Cesare de Filippo
## Maximum likelihood Co-Estimate of Contamination And Split Time (cecast) in low-coverage samples.

library("argparse", quietly=TRUE)
library("parallel", quietly=TRUE)
library("RColorBrewer", quietly=TRUE)
library("jsonlite", quietly=TRUE)
suppressMessages(library(fields, warn.conflict = FALSE, quietly = TRUE))

# Create a parser
parser <- ArgumentParser(description="Maximum likelihood Co-Estimates of human Contamination And Split Times (cecast).")

# Add command line arguments
parser$add_argument("input_file",
help="input file", nargs='?', type="character")

parser$add_argument("-c", dest="contamination",
help="Lower and upper bound of human contamination to be explored like '0,1' (the default value).", type="character", default='0,1')

parser$add_argument("-t", dest="split_time_range",
help="Range of split times to be searched in the form of 'min,max'.", type="character")

parser$add_argument("-n", dest="number_boot",
help="Number of bootstraps to compute the confidence intervals (CI). Default is 100.", type="integer", default=100)

parser$add_argument("-b", dest="boot_type",
help="Type of bootstraps to be performed: 'd' data (the default), 'h' history, and 'b' both. The uncertainity of demographic 'history' can use at most 100 bootstraps (-n flag).", type="character", default='d')

parser$add_argument("-q", dest="quantile",
help="Quantile for the confidence intervals (CI). Default is 0.95, which is 0.025 and 0.975 for the two tails.", type="double", default=0.95)

parser$add_argument("-k", dest="chromosomes",
help="Chromosome(s) to be considered. 'A|a' is for all autosomes, 'X|x' for X and 1..22 for any other chromosome.", default='1:22', type="character")

parser$add_argument("-p", dest="population",
help="Human population as source of contaminantion. Default is Mbuti, and the ohers are: Dai, French, Han, Mandenka, Papuan, San, Sardinian, Yoruba, Karitiana", default="Mbuti", type="character")

parser$add_argument("-r", dest="human_reference",
help="The human sample to be used: y for Yoruba (default) and m for Mbuti. Notice that Yoruba is the same individual as the one for source of contamination, while the Mbuti is a differenct individual. Therefore, another individual should be specified with -p. ", default="m", type="character")

parser$add_argument("-e", dest="exclude_class",
help="Remove a given class of sites (lineages) for the analyses. Use the commas to separate the different classes.", type="character")

parser$add_argument("-a", dest="ascertainment",
help="The type of ascertainment as one of the following snps arrays: (b) big_steffi; (p) archaic_plus; (v) archaic_var. It's possible to specify just one of the tree letters (b, p, v), the full or partial name.",  type="character")

parser$add_argument("-o", dest="outgroup",
help="The type of outgroup or better ancestry to be used: CBGO, CBG, CGO, CBO, BGO. Default is to use all.",  type="character", default=NULL)

parser$add_argument("-s", dest="source", help="Input file with the count of the source of  contamination.", type="character")

parser$add_argument("-x", dest="expectations", help="The file with the expecations in .RData format", type="character")

parser$add_argument("-W", dest="window",
help="Run the estimate for each window in the input file.", action="store_true")

parser$add_argument("-O", dest="output",
help="Output the results of all estimates and bootstraps. The estimates that differ less than 3.4 from the max log(likelihood) are printed.", action="store_true")

parser$add_argument("-B", dest="Binomial", help="Use the binomial to calculate the likelihood.", action="store_true")

parser$add_argument("-H", dest="Header", help="Print only the header.", action="store_true")


## Parse the command line arguments and exit if something is not fine.
argv <- parser$parse_args()
if(!is.null(argv$help)) {
	parser$print_help()
	q(status=1)
}

if(argv$Header) {
	cat("t\tt_low\tt_high\tpop\tc\tc_low\tc_high\tlogLike\tc_source\tnseq\tchr\tboot\tcomments\n")
	q(status=1)
}

if(is.null(argv$input_file)) {
	cat("Input file required.\n")
	parser$print_usage()
	q(status=1)
}


#################################################################################
##-------------------------------------------------------------------------------

## Load configuration 
config <- fromJSON("~/.config_cecast.json")

## First load the functions 

source(paste(config$scripts_dir, config$functions_file, sep="/"))

# Functions loaded
# LikeCon: calculate the likelihood integrating human contamination.
# LikeFun: calulate the likelihood over the entire set of parameters. It uses the previous function 'LikeCon'.
# estimate: estimate the parameters and report the result in various formats to be processed later.
# proc_est: process the likehoods after running the function 'estimate'.


## Define some default values

if(is.null(argv$Binomial)) {
	argv$Binomial = FALSE
}

if(!is.null(argv$outgroup)) {
	argv$outgroup = toupper(unlist(strsplit(argv$outgroup, split=",")))
	outgroups = c("CBGO", "CBG", "CGO", "CBO", "BGO")
	if(sum(argv$outgroup %in% outgroups) == 0) {
		cat("Outgroup must be at least one of the following:", outgroups, "\n")
		q(status=1)
	}
}

Warnings = c()

## The default counts from simulations and some humans to be used to calculate the likelihood.
data_dir <- config$data_dir
counts_files <- file.path(data_dir, config$data_files)
for(f in counts_files) {
	load(f)
}

if(!is.null(argv$expectations)) {
	if(file.exists(argv$expectations)) {
		# still need to check that the file contains the proper objects.
		load(argv$expectations)
	}
}


## SNPs arrays
array_list <- list(b=c("afr_combined", "big_steffi"), v="var")
array_list_s <- sapply(array_list, function(x) strsplit(x, split=""))
array_names <- c("big_steffi", "var")

array_sims <- sapply(c("archaics_afr_combined", "archaics_var"), function(snp) paste(config$data_dir, "/expected_counts_fine_grid_", snp,".RData" ))
names(array_sims) <- names(array_list)

## Load the specified expectation for the spescified ascertainment.
if(!is.null(argv$ascertainment)) {
	ids <- which(tolower(argv$ascertainment) == names(array_sims))
	if(length(ids) < 1) {
		ids <- which(sapply(array_list, function(x) sum(tolower(argv$ascertainment) == x) ) == 1)
		if(length(ids) < 1) {
			a <- tolower(strsplit(argv$ascertainment, split="")[[1]])
			ids <- which.max(sapply(array_list_s, function(x) max(sapply(x, function(z) sum(z %in% a)/length(z)))))
		}
	}
	load(array_sims[ids])
	rm(ids)
}

names(Lineages) <- c("vin", "cha", "alt", "den", "v-c", "nea", "arc")

## -c: Human contamination to be explored
a <-  strsplit(argv$contamination, split=",")[[1]]
z <- sum(is.na(sapply(a, as.numeric))) # This is in case some non-numeric values have been given.
if (length(a) == 2 & z == 0) {
	Contaminations = sort(as.numeric(a))
} else {
	Contaminations = c(0, 1)
}

rm(a, z, array_list, array_list_s, array_sims)

## -s: The split time ranges to be explored
if (!is.null(argv[["split_time_range"]])) {
	a <-  as.numeric(strsplit(argv[["split_time_range"]], split=",")[[1]])
	if (length(a) == 2) {
		argv[["split_time_range"]] = sort(a)
	} else {
		argv[["split_time_range"]] = c(0,3000)
		Warnings = c(Warnings, "(-t) split_time_range: only two time values in Kya are accepted.")
	}
	rm(a)
}


# The name of the human samples as a source of contamination
human_sources <- names(human_frequencies) 
names(human_sources) <- c("Dai", "Fre", "Han", "Man",  "Mbu", "Pap", "San", "Sar", "Yor", "Kar")


## Get the most similar sample or population specified with -p option.
## This is to avoid to write the entire name of the population.
if(sum(!is.null(argv$population) & nchar(argv$population) > 2) == 1) {
	if (sum(toupper(argv$population) %in% c("SIMS", "SIM","SIMULATIONS", "SIMULATION")) > 0) {
		argv$population = NULL
		c_source = "simulation"
	} else {
		perfect_match <- c_source <- human_sources[toupper(human_sources) == toupper(argv$population)]
		if(length(perfect_match) == 0) {
			best_match <- seq_match(argv$population, human_sources)
			if(is.na(best_match)) {
				Warnings <- c(Warnings, paste("(-p) population: there is no match with", argv$population, "specified with -p option, Yoruba is taken."))
				argv$population <- NULL
			} else {
				Warnings <- c(Warnings, paste("(-p) population: the best match for", argv$population, "is", best_match))
				argv$population <- c_source <- best_match
			}
		} else {
			argv$population <- c_source <- perfect_match
		}
		rm(perfect_match)
	}
}

# the human as a reference, this prevent the error in case is not specified as 'y' or 'm'. 
human_reference <- seq_match(argv$human_reference, c("y","m"))
if(is.na(human_reference)) {
	human_reference = "y"
	Warnings <- c(Warnings, paste("(-r) reference: no match for human reference", argv$population, ", y (Yoruba) is chosen as by default."))
}

# Make sure that type and number of bootstraps are compatible.
argv$boot_type <- tolower(strsplit(argv$boot_type, split="")[[1]][1])
if(argv$number_boot > 100) {
	if(argv$boot_type == "b" | argv$boot_type == "h" ) {
		argv$number_boot <- 100
	}
}
names(argv$boot_type) = c(d="data", h="history", b="both")[argv$boot_type]



## Other objects:
## I should change these with one letter abbreviation like 'a', 'a-c', 'a-c-v-d', etc...
ClassOfSites <- c("altai", "altai-chagyrskaya", "altai-chagyrskaya-denisova", "altai-denisova", "altai-vindija", "altai-vindija-denisova", "chagyrskaya", "chagyrskaya-denisova", "denisova", "human", "human-altai", "human-altai-chagyrskaya", "human-altai-chagyrskaya-denisova", "human-altai-denisova", "human-altai-vindija", "human-altai-vindija-denisova", "human-chagyrskaya", "human-chagyrskaya-denisova", "human-denisova", "human-neandertal", "human-vindija", "human-vindija-chagyrskaya", "human-vindija-chagyrskaya-denisova", "human-vindija-denisova", "neandertal", "neandertal-denisova", "vindija", "vindija-chagyrskaya", "vindija-chagyrskaya-denisova", "vindija-denisova")

names(ClassOfSites) <- c("a", "a-c", "a-c-d", "a-d", "a-v", "a-v-d", "c", "c-d", "d", "x", "x-a", "x-a-c", "x-a-c-d", "x-a-d", "x-a-v", "x-a-v-d", "x-c", "x-c-d", "x-d", "x-a-v-c", "x-v", "x-v-c", "x-v-c-d", "x-v-d", "a-v-c", "a-v-c-d", "v", "v-c", "v-c-d", "v-d")

if(!is.null(argv$exclude_class)) {
	l = unlist(strsplit(argv$exclude_class, split=","))
	idx = which(ClassOfSites %in% l | names(ClassOfSites) %in% l )
	ClassOfSites = ClassOfSites[-idx]
	n1=names(Simulations)
	n3=names(human_frequencies)
	if(length(idx) > 0) {
		Simulations = lapply(Simulations, function(x) x[-idx,])
		Simulations_bootstraps = lapply(Simulations_bootstraps, function(z) {n2=names(z); w=lapply(z, function(x) x[-idx,]); names(w) = n2; w })
		human_frequencies = lapply(human_frequencies, function(x) x[-idx,])
		names(Simulations) = n1
		names(human_frequencies) = n3
	}
}


######################################################
### Read the data
if (!is.null(argv$input_file) & file.exists(argv$input_file) & file.info(argv$input_file)$size) {
	xdata = read_lineages(argv$input_file, anc=argv$outgroup, hum=human_reference)
	# Check that the xdata is alligned with the Simulations. This will crash if some matching will not work.
	# the column within xdata of the human as source of contamination
	h_id = grep(names(c_source), colnames(xdata[[1]]))
	# I will drop this later on.
	if(sum(rownames(xdata[[1]]) %in% ClassOfSites) == 0 ) {
		lin0 = strsplit(gsub("-","",names(ClassOfSites)), split="")
		lin1 = strsplit(gsub("-","",rownames(xdata[[1]])), split="")
		l0 = sapply(lin0, function(x) paste(sort(x),collapse="-"))
		l1 = sapply(lin1, function(x) paste(sort(x),collapse="-"))
		ids = match(l0,l1)
		# In case there is some lineage missing in the data
		if(sum(is.na(ids))) {
			idx = l0[is.na(ids)]
			m = matrix(0, ncol=ncol(xdata[[1]]),nrow=length(idx), dimnames=list(idx,colnames(xdata[[1]])))
			xdata = lapply(xdata, function(x) rbind(x,m))
			l1 = c(l1,idx)
		}
		ids = match(l0,l1)
		xdata = lapply(xdata, function(x) {x=x[ids, c(1:2,h_id[1:2])]; rownames(x) = ClassOfSites; x})
	} else {
		xdata = lapply(xdata, function(x) x[rownames(x) != "all", c(1:2, h_id[1:2])])
	}

	if(!is.null(argv$exclude_class)) {
		CoS = unlist(strsplit((argv$exclude_class), split="," ))
		if ( sum(CoS %in% ClassOfSites) > 0) {
			xdata = lapply(xdata, function(x) x[!rownames(x) %in% CoS, ])
			Simulations = lapply(Simulations, function(x) x[!rownames(x) %in% CoS,])
		}
	}
	# Unify the ClassOfSites in xdata, Simulations, Simulations_bootstraps, human_frequencies
	s1 = sub("africa", "human", sub("neandertal", "altai-vindija-chagyrskaya", rownames(Simulations[[1]])))
	for (i in 1:length(xdata)) {
		rownames(xdata[[i]]) = sub("africa", "human", sub("neandertal", "altai-vindija-chagyrskaya", rownames(xdata[[i]])))
	}
	for (i in 1:length(Simulations)) {
		rownames(Simulations[[i]]) = s1
	}
	for(j in 1:length(Simulations_bootstraps)) {
		for(i in 1:length(Simulations_bootstraps[[j]])) {
			rownames(Simulations_bootstraps[[j]][[i]]) = s1
		}
	}
	for (i in 1:length(human_frequencies)) {
		rownames(human_frequencies[[i]]) = s1
	}
	n = sapply(xdata, function(x) sum(x[,1:2]))
	## Remove blocks without observations in the observed data and also in the source of contaminations
	xdata = xdata[n > 0]
	if (argv$chromosomes %in% c("a", "all")) {
		argv$chromosomes = "1:22"
	}
	chrs = strsplit(gsub(":","-", toupper(argv$chromosomes)), split="-")[[1]]
	if (length(chrs) == 2) {
		chrs=seq(as.numeric(chrs[1]), as.numeric(chrs[2]), 1)
	} else {
		if (length(chrs) > 2) {
			Warnings = c(Warnings,"(-k) chromosomes: not conform to format. All autosomes were used.")
			chrs = 1:22
		}
	}
	ids = unique(unlist(sapply(chrs, function(k) grep(paste("^",k,":",sep=""), names(xdata), perl=TRUE))))
	if(length(ids) > 0) {
		xdata = xdata[ids]
	}
	n = length(xdata)
	n_sequences = do.call("sum", lapply(xdata, function(x) x[,1:2]))
} else {
	stop("Check input file or ask for help.")
	quit()
}

## Do the analyses for each window of the input file.
if(argv$window & n == 1) {
	stop("Input file is not conform for estimates per window (-w). Check input file or ask for help.\n")
	quit()
}

######################################################
### Run the estimates and return the output(s).

main = function() {
	if(argv$window) {
		e1 = mclapply(xdata, function(z) proc_est(estimate(list(z), return_for_plot=FALSE, expectations=Simulations, return_LLratio=TRUE, Binomial=argv$Binomial)), mc.cores=detectCores()-1 )
		output = do.call("rbind",lapply(e1, function(x) x[["Estimates"]]))
		comments = sapply(e1, function(x) if(is.null(x[["Comments"]])) { NA } else {x[["Comments"]]})
		nseq = sapply(xdata, sum)
		output=cbind('#window'=rownames(output), output, nseq, comments)
		output[,"pop"]=names(Lineages)[output[,"pop"]]
		cat(paste(colnames(output), collapse="\t"), "\n")
		apply(output,1, function(x) cat(x,"\n",sep="\t"))
	} else {
		e1 = estimate(xdata, return_for_plot=TRUE, expectations=Simulations, return_LLratio=TRUE, Binomial=argv$Binomial);
		res_to_plot = e1[["Results"]]
		res_to_out = proc_est(e1[["Estimates"]])
		comments = ""
		if(!is.null(res_to_out[["Comments"]])) {
			comments = paste(nrow(e1[["Estimates"]]), "_estimates_with_logLike<3.4:", res_to_out[["Comments"]], sep="")
		}
		res_to_out = res_to_out[["Estimates"]][,-1]
		## The confidence intervals...
		if (length(xdata) == 1 | argv[["boot_type"]] == 'h') {
			# If the input is only one window with anc/der counts, it is only possible to bootstrap the demongraphy/history (h).
			boot = do.call("rbind", lapply(sample(1:length(Simulations_bootstraps), argv[["number_boot"]]), function(x) estimate(xdata, expectations=Simulations_bootstraps[[x]], Binomial=argv$Binomial)))
		} else {
			if (argv[["boot_type"]] == 'd') {
				boot = do.call("rbind", lapply(1:argv[["number_boot"]], function(x) estimate(xdata[sample(1:n, n, replace=T)], Binomial=argv$Binomial)))
			}
			if (argv[["boot_type"]] == 'b') {
			boot = do.call("rbind", lapply(sample(1:length(Simulations_bootstraps), argv[["number_boot"]]), function(x) estimate(xdata[sample(1:n, n, replace=T)], expectations=Simulations_bootstraps[[x]], Binomial=argv$Binomial)))
			}
		}
		boot1 = proc_est(boot, best_pop=res_to_out[1,"pop"], quantile.prob=argv$quantile)
		bbtype = names(argv[["boot_type"]]) # To show the type of bootstrap used.
		b = boot1[["Estimates"]]
		res_to_out[,c("t_low", "t_high", "c_low","c_high")] = b[,c("t_low", "t_high", "c_low","c_high")]
		if(!is.null(boot1[["Comments"]])) {
			bbtype = paste(bbtype, boot1[["Comments"]], sep=":")
		}
		if (argv$output) { ## Print all values: estimate(s) and bootstraps
			output = data.frame(rbind(cbind(e1[["Estimates"]], Type=rep("estimate", nrow(e1[["Estimates"]]))), cbind(boot, Type=rep(bbtype, nrow(boot)))))
			output[,"pop"]=names(Lineages)[output[,"pop"]]
			cat("#t_x\tt_ND\tpop\tc\tlogLike\tType\n")
			apply(output,1, function(x) cat(x,"\n",sep="\t"))
		} else { ## Print the normal output.
			res_to_out[,"pop"] = names(Lineages)[res_to_out[,"pop"]]
			cat("#", paste(c(names(res_to_out), "c_source\tnseq\tchr\tboot\tcomments"),collpase="\t"), "\n", sep="")
			cat(paste(c(res_to_out, c_source, n_sequences, argv$chromosomes, bbtype, comments),collpase="\t",sep=""), "\n", sep="")
		}
		if (length(Warnings) > 0) {
			warning(paste(Warnings, collapse="\n"), call.=FALSE)
		}
	}
}

if(n_sequences > 10) {
	tryCatch(main(), finally=quit(status=1))
} else {
	stop(paste("There are only", n_sequences, "sequences in the input."))
	quit()
}
