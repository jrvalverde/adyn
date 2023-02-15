#!/usr/bin/env Rscript

#
#	This is an R script file
#
suppressPackageStartupMessages(require(methods))
suppressPackageStartupMessages(require(tools))

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(sybil))
suppressPackageStartupMessages(library(sybilSBML))
suppressPackageStartupMessages(library(sybilDynFBA))

library(adfba)
library('banneR')

SYBIL_SETTINGS('SOLVER', 'glpkAPI') # should work but fails on AML in Ilumina  
#SYBIL_SETTINGS('METHOD', 'simplex')
#SYBIL_SETTINGS('METHOD', 'interior')
#SYBIL_SETTINGS('METHOD', 'exact')
#SYBIL_SETTINGS('METHOD', 'mip')

# lpSolve should work but is slow
#SYBIL_SETTINGS('SOLVER', 'lpSolveAPI') # should work but gives different results
#SYBIL_SETTINGS('METHOD', 'lp_solve')

# free version refuses to work with >1000 equations
#SYBIL_SETTINGS('SOLVER', 'cplexAPI') 
#SYBIL_SETTINGS('METHOD', 'lpopt')
#SYBIL_SETTINGS('METHOD', 'primopt')
#SYBIL_SETTINGS('METHOD', 'dualopt')
#SYBIL_SETTINGS('METHOD', 'baropt')
#SYBIL_SETTINGS('METHOD', 'hybbaropt')
#SYBIL_SETTINGS('METHOD', 'hybnetopt')
#SYBIL_SETTINGS('METHOD', 'siftopt')
#SYBIL_SETTINGS('METHOD', 'mipopt')
#SYBIL_SETTINGS('METHOD', 'qpopt')


# gives results similar to glpk
#SYBIL_SETTINGS('SOLVER', 'clpAPI')
#SYBIL_SETTINGS('METHOD', 'general_solve')
#SYBIL_SETTINGS('METHOD', 'inidual')
#SYBIL_SETTINGS('METHOD', 'iniprimal')
#SYBIL_SETTINGS('METHOD', 'inibarrier')
#SYBIL_SETTINGS('METHOD', 'inibarriernoc')
#SYBIL_SETTINGS('METHOD', 'idiot')
#SYBIL_SETTINGS('METHOD', 'dual')
#SYBIL_SETTINGS('METHOD', 'primal')

#SYBIL_SETTINGS('TOLERANCE', '1E-6')
#SYBIL_SETTINGS('TOLERANCE', '1E-5')

# Usage: R --slave -f dynFBA.R --args <model name> <nSteps> 
#
# Or, for many files
#    for i in xxx*_react.tsv ; do 
#	echo doing $i ; 
#	R --slave -e 'source("dynFBA.R")' --args `basename $i _react.tsv` $constrains $medium 60 ; 
#	echo $i done ; 
#    done
#
# for i in xxx*_react,tsv ; do R --slave -e 'source("dynFBA.R")' --args `basename $i _react.tsv` basename $contrains .tab` $medium 60 ; echo $i ; done
#
# for i in *.tab ; do R --slave -e 'source("dynFBA.R")' --args Slp_aml `basename $i .tab` NMMP 60 ; echo Slp_aml $i ; done


#
# DEFINE COMPUTATION PARAMETERS
# =============================
#
# --------------------------------------------
# This model uses as units mmol per gDW per hr
# gDW (grams Dry Weight)
# --------------------------------------------
#
# Objective function(s)
# ---------------------
# the objective coefficients specify the weights given to the
# respective reactions in the rxnNameList
#
# according to the TNF XML file, the default objective is R701,
#	EX_BIOMASS (exchange of Biomass)
#
	biomass <- 'Biomass_SLI';
	obj <- biomass
#obj <- 'EX_BIOMASS';
#obj <- 'BIOMASS_SCO';	# For S. coelicolor
#obj <- 'TNF'
#obj <- 'AML'
#ibj <- 'DAG'
#obj <- 'SEC_TNF'
#obj <- 'SEC_AML'
#obj <- 'SEC_DAG'
#obj <- 'EX_tnf(e)'
#obj <- 'EX_aml(e)';
#obj <- 'EX_dag(e)';
	ci <- 1;
#obj <- c('EX_aml(e)', 'Biomass_SLI');
#obj <- c('Biomass_SLI', 'EX_aml(e)');
#ci <- c(1, 1);


# Default control for DFBA
# ---------------------
	# time step and length of simulation
	ts <- 1;	# hr
	nSteps <- 60	# 60 * 1 hours
	# experimental value of inoculum size
	inoculum <- 0.01;	# grams / L
 	substrateName <- ""

# select defaults for calculations
# --------------------------------
        doFBA  <- FALSE
        doMTF  <- FALSE
        doFVA  <- FALSE
        doRA   <- FALSE
        doPPPA <- FALSE
        doDFBA <- FALSE
        doADFBA<- TRUE

# whether to save output by default
# ---------------------------------
	log=1


# to be used instead of library(): this function ensures
# that the package is installed if not present in the system
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE)
    require(p, character.only = TRUE)
}


# source all Q/R/S files in a directory
#	this may be convenient to source more than one file at once
#	i.e. when getting dynamic constraints as a function,
#	although it would likely make more sense if the file defining
#	the function did include (source) itself all other required
#	files.
sourceDir <- function(path, trace = TRUE, ...) {
   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
   }
}


openLogFile <- function(logName) {

    if (file.exists(logName)) file.remove(logName)

    logFile <- file(logName, open="wt");
    sink(logFile, type=c("output", "message"), split=TRUE);
    return (logFile)
}


closeLogFile <- function(logFile)  {
    close(logFile)
    while (sink.number() > 0) { sink() }
}


file.extension <- function (path) 
{
    parts <- strsplit(path, ".", fixed=T)[[1]]
    last <- parts[length(parts)]
    last
}


# Load model
# ==========
#
# first try SBML and if that fails, try TSV
#
# NOTE: we should also check if the full filename has been given by testing
# its extension(s)
#
loadModel <- function(modelName) {
    xmlName <- paste(modelName, ".xml", sep="");
    tsvName <- modelName;	# files must be called tsvName_{react|met|desc}.tsv


    if (file.exists(xmlName)) {
	    model <- readSBMLmod(xmlName);
    } else if (file.exists(paste(tsvName, "_react.tsv", sep=""))) {
	    # we only need the reaction table as a minimum
	    model <- readTSVmod(tsvName);
    } else if (file.exists(modelName)) {
    	    # as a last resort try to load it as an SBML file
            xmlName <- modelName
            model <- readSBMLmod(xmlName);
    } else {
 	    stop('No suitable model (.xml or .tsv) found');
    }
    return( model )
}


saveModel <- function(model, tsv=TRUE, sbml=TRUE, overwrite=TRUE, outDir=".") {
    modelName <- model@mod_name
    if (outDir == "") { outDir <- "." }
    tsvdir <- paste(outDir, "/", "s_tsv/", sep="")
    sbmldir <- paste(outDir, "/", "s_sbml/", sep="")
    dir.create(tsvdir, showWarnings = FALSE)
    dir.create(sbmldir, showWarnings = FALSE)
    tsvName   <- paste( tsvdir, modelName, sep="" );	# files will be named modelName_{desc|met|react}.tsv
    sbmlName  <- paste( sbmldir, modelName, "_L2v1.xml", sep="" );

    if (overwrite == FALSE) {
	if ((tsv == TRUE)
           && (! file.exists(paste(tsvName, '_react.tsv', sep='')))) {
	    cat('Writing file ', tsvName, '_{desc|met|react}.tsv\n')
	    modelorg2tsv(model, 
   			prefix=tsvName, suffix="tsv",
			extMetFlag="b", makeClosedNetwork=TRUE)
	}
	# and in SBML format
	if ((sbml == TRUE) && (! file.exists(sbmlName))) {
	    cat('Writing file ', sbmlName, '\n')
	    writeSBML(model, filename=sbmlName, level=2, version=1)
	}
    }
    else {
	if (tsv == TRUE) {
	    cat('Writing file ', tsvName, '_{desc|met|react}.tsv\n')
            modelorg2tsv(model, 
   		    prefix=tsvName, suffix="tsv",
		    extMetFlag="b", makeClosedNetwork=TRUE)
	}
        if (sbml == TRUE) {
	    cat('Writing file ', sbmlName, '\n')
	    writeSBML(model, filename=sbmlName, level=2, version=1)
	}
    }
}

# This one could be simplified by ignoring the extension and smply
# trying all known formats in order
#	That'd be safer but "dirtier"
loadTableFile <- function(file="", 
                          stringsAsFactors=FALSE, 
                          header=FALSE,
                          check.names=FALSE)
{
    if (file == "") {
        return(NULL)
    }
    
    options(stringsAsFactors=stringsAsFactors)
    
    # make a vector with all acceptable suffixes 
    suff <- c('tsv', 'csv', 'tab', 'dat', 'txt')
    if ((is.element('xlsx', installed.packages()[,1]))) {
        suff <- c(suff, 'xlsx')
    }
    if ((is.element('foreign', installed.packages()[,1]))) {
        suff <- c(suff, 'arff', 'dbf', 'dta', 'rec', 'mtp', 'mat', 
                  'sav', 'ssd', 'sas', 'sys', 'syd', 'xpt', 'dump')
    }
    fnames <- c(file)
    fnames <- c(fnames, paste(file, suff, sep='.'))
    avail <- file.exists(fnames)

    if (all(avail == FALSE)) {
        file <- file.choose()
        fnames <- c(file)
        avail <- c(TRUE)
	suff <- c(file_ext(file))
    }
        
    table <- data.frame()
    for (i in 1:length(fnames)) {
        if (! avail[i]) {
            next
        }
        # get suffix
        ext <- file_ext(fnames[i])
        if (ext == 'tsv') {
            table <- read.table(fnames[i], sep="\t", 
                                header=header, check.names=check.names)
            if (! empty(table)) { return(table) }
        } else if (ext == 'tab') {
            table <- read.table(fnames[i], sep="\t", 
                                header=header, check.names=check.names)
            if (! empty(table)) { return(table) }
        } else if (ext == 'csv') {
            table <- read.table(fnames[i], sep=",", 
                                header=header, check.names=check.names)
            if (! empty(table)) { return(table) }
        } else if (ext == 'txt') {
            table <- read.table(fnames[i], 
                                header=header, check.names=check.names)
            if (! empty(table)) { return(table) }
        } else if (ext == 'dat') {
            table <- read.table(fnames[i], header=header, check.names=check.names)
            if (! empty(table)) { return(table) }
        } else if ((is.element('xlsx', installed.packages()[,1]))) {
	    if (ext == 'xlsx') {
                library(xlsx)
                # use xlsx2 to convert dates to POSIXct
                table <- read.xlsx2(fnames[i], 
                                    header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            }            
        } else if ((is.element('foreign', installed.packages()[,1]))) {
            library(foreign)
            if (ext == 'arff') {
                table <- read.arff(fnames[i], 
                                   header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'dbf') {
                table <- read.dbf(fnames[i], 
                                  as.is=stringsAsFactors, 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'dta') {
                table <- read.dta(fnames[i], 
                         convert.factors=stringsAsFactors,
                         convert.dates=TRUE,
                         convert.underscore=TRUE, 
                         header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'rec') {
                table <- read.epiinfo(fnames[i], 
                                      get.broken.dates=TRUE,
                                      header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'mtp') {
                table <- read.mtp(fnames[i], 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'mat') {
                table <- read.mtp(fnames[i], 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if (ext == 'sav') {
                table <- read.spss(fnames[i], 
                                   use.value.abels=stringsAsFactors,
                                   to.data.frame=TRUE,
                                   header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if ((ext == 'ssd') || (ext == 'sas')) {
                table <- read.ssd(fnames[i], 
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if ((ext == 'sys') || (ext == 'syd')) {
                table <- read.systat(fnames[i], 
                                  to.data.frame=TRUE,
                                  header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if ((ext == 'xpt')) {
                table <- read.xport(fnames[i], 
                                    header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            } else if ((ext == 'dump')) {
                #table <- data.restore(fnames[i]) 
                table <- read.S(fnames[i], 
                                header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            }
        } else {
            # At this point we have tried all known suffixes. As a last resort, 
            # we will try to read it as a text file
            table <- read.table(fnames[i], header=header, check.names=check.names)
	    # and leave it as default if nothing else works
        }
    }
    return(table)    
}


loadConstraints <- function(constraints="")
{
    if (constraints == "") {
        return (NULL)
    }
    # Read mods. file (if any) and apply it
    options(stringsAsFactors=FALSE)		# we want strings as strings
    # by default, # will be used as comment.char
    if (file.exists(paste(constraints, ".tsv", sep="")))
    {
    	limits <- read.table(paste(constraints, ".tsv", sep=""), sep="\t", header=FALSE)
    } else if (file.exists(paste(constraints, ".tab", sep=""))) {
     	limits <- read.table(paste(constraints, ".tab", sep=""), sep="\t", header=FALSE)
    } else if (file.exists(paste(constraints, ".csv", sep=""))) {
    	limits <- read.table(paste(constraints, ".csv", sep=""), sep=",", header=FALSE)
    } else if (file.exists(paste(constraints, ".txt", sep=""))) {
    	limits <- read.table(paste(constraints, ".txt", sep=""), header=FALSE)
    } else if (file.exists(constraints)) {
        # last resort, try to read it as a TAB file
   	limits <- read.table(constraints, sep="\t", header=FALSE)
    } else if ((is.element('xlsx', installed.packages()[,1])) &&
             (file.exists(paste(constraints, ".xlsx", sep="")))) {
        limits <- read.xlsx(paste(constraints, ".xlsx", sep=""), sheet.index=1)
    } else {
        cat("File ", constraints, "{.txt, .tsv, .tab, .csv} not found\n")
    	limits <- read.table(file.choose(), header=FALSE)
    }
    return(limits)
}

applyConstraints <- function(model, obj=NULL, ci=NULL, limits=NULL, verbose=0)
{
    # set objective function
    # ----------------------
    if (verbose > 0) {
        cat("\n\nCustomizing model\n")
        cat(    "-----------------\n\n")
    }
    
    if (obj != NULL) {
        if (ci == NULL) {
            ci = rep(1, length(obj))
        }
        cat("Setting Objective Function(s) to ", obj, ci, "\n\n")
        model <- changeObjFunc(model, obj, ci)
    }
    if (limits == NULL) {
        return (model)
    }
    
    # loop over all the limits and apply them
    if (verbose > 0) {
        cat('Setting model constraints:\n')
        cat('    REACTION\tLB\tUB\n')
        cat('    --------------------------------------------------\n')
    }
    for ( i in 1:dim(limits)[1] ) {
        reactName <- limits[i,1]
        lb <- limits[i,2]
        ub <- limits[i,3]
        cat('   ', reactName, '\t', lb, '\t', ub, '\n')
        lowbnd(model)[react_id(model)== reactName] <- lb;
	uppbnd(model)[react_id(model)== reactName] <- ub;
    }
    if (verbose > 0) {
        cat('    --------------------------------------------------\n')
        cat('\n')
    }        
    return (model)

}

# Unused, redefined below. Work in progress.
customizeModel <- function(model, obj, ci, constraints="", verbose=0)
{
    if (constraints != '') {
        # Read constraints from file
        limits <- loadConstraints(constraints)
        if (limits == NULL) {
            cat('customizeModel: no constraints obtained from', constraints, '\n\n')
        } else {
            if (verbose > 0) {
                cat('Setting constraints from file ', constraints, ':\n')
            }
        }
    }

    # Apply required changes
    applyConstraints(model, obj,ci, limits)
    
    # reflect changes in model name
    if (verbose > 0) {
        cat("Changing model name from", model@mod_name, "to", 
    	    paste(model@mod_name, "_", constraints, sep=""), "\n\n")
    }
    model@mod_name <- paste(model@mod_name, "_", constraints, sep="")
        
    return (model)

}


# MODIFY MODEL TO REFLECT OUR CUSTOMIZATION
# =========================================
# Limit uptake rates
# ------------------
# Reactions with upper and lower bounds set to zero are not functional.
# Exchange reactions with lower bound = 0 and upper bound > 0 only
#	allow secretion.
# To enable metabolite uptake, set lower bound < 0
#
#	THIS FUNCTION MODIFIES THE MODEL
#
customizeModel <- function(model, obj, ci, constraints="") 
{
    # set objective function
    # ----------------------
    cat("\n\nCustomizing model\n")
    cat(    "-----------------\n\n")
    cat("Setting Objective Function(s) to ", obj, ci, "\n\n")
    model <- changeObjFunc(model, obj, ci)

    if (constraints == "") {
        return (model)
    }
    # Read mods. file (if any) and apply it
    options(stringsAsFactors=FALSE)		# we want strings as strings
    # by default, # will be used as comment.char
    if (file.exists(paste(constraints, ".tsv", sep="")))
    {
    	limits <- read.table(paste(constraints, ".tsv", sep=""), sep="\t", header=FALSE)
    } else if (file.exists(paste(constraints, ".tab", sep=""))) {
     	limits <- read.table(paste(constraints, ".tab", sep=""), sep="\t", header=FALSE)
    } else if (file.exists(paste(constraints, ".csv", sep=""))) {
    	limits <- read.table(paste(constraints, ".csv", sep=""), sep=",", header=FALSE)
    } else if (file.exists(paste(constraints, ".txt", sep=""))) {
    	limits <- read.table(paste(constraints, ".txt", sep=""), header=FALSE)
    } else if (file.exists(constraints)) {
        # last resort, try to read it as a TAB file
   	limits <- read.table(constraints, sep="\t", header=FALSE)
    } else {
        cat("File ", constraints, "{.txt, .tsv, .tab, .csv} not found\n")
    	limits <- read.table(file.choose(), header=FALSE)
    }
    
    # loop over all the limits and apply them
    cat('Setting constraints from file ', constraints, ':\n')
    cat('    REACTION\tLB\tUB\n')
    cat('    --------------------------------------------------\n')
    for ( i in 1:dim(limits)[1] ) {
        reactName <- limits[i,1]
        lb <- limits[i,2]
        ub <- limits[i,3]
        cat('   ', reactName, '\t', lb, '\t', ub, '\n')
        lowbnd(model)[react_id(model)== reactName] <- lb;
	uppbnd(model)[react_id(model)== reactName] <- ub;
    }
    cat('    --------------------------------------------------\n')
    cat('\n')
    cat("Changing model name from", model@mod_name, "to", 
    	paste(model@mod_name, "_", constraints, sep=""), "\n\n")
    model@mod_name <- paste(model@mod_name, "_", constraints, sep="")
        
    return (model)
}


# Check the model
#----------------
# inspect model: which reactions have l/u bounds within the limits
#
checkModel <- function(model)  {
    xch <- findExchReact(model);
    upt <- uptReact(xch);

    cat('\n\nExchange reactions\n\n')
    print(xch)
    cat('\n\nUptake reactions\n\n')
    print(upt)
    
    # add additional checks here
}


# Setup experimental data
# -----------------------
#
# This function reads a file with experimental values measured at specific
# times, and interpolates the values for the times that we are going to
# simulate, returning a data frame that we can use to steer the dynamically
# adjusted FBA calculations.
#
# Data is read from a file in TAB-separated format
#	Data MUST be organized as follows
#	First line contains 'Time' and reaction names
#	Second value containes observed values for each reaction at each time
#	Missing values are indicated by NA
#
# Values are interpolated from the data in the file to fill in
# a table with nSteps values separated tStep hours each.
#
# The Biomass reaction may be specified so it can be approximated
# to a Sigmoid function.
#
# Any non-Biomass reaction will be approximated by linear interpolation
# though other interpolation methods may be specified.
#
# If you also want Biomass interpolated using the same method as all
# other reactions, simply do not specify it.
#

# The sigmoid function that we will model. We'll use this function
# to generate values using the predicted model.
sigm <- function(p, x) {
    p[1] + ((p[2] - p[1]) / (1 + exp(-p[3] * (x - p[4]))))
}

# A function to convert "observed" values (note that they might as well be
# the result of interpolation between a few key points). There must be
# one value per time point, and that a time point may be repeated.
#
# In effect this produces the rate that will change a concentration at a 
# given time into the concentration at the next time step.
obs2rates <- function(time, conc) {
    n <- length(time)
    t <- c()
    delta <- c()

    for (i in 1:(n-1)) {
        # AVOID USING THIS "FEATURE" UNTIL WE CAN SEE IF IT MAKES SENSE
    	# This should allow us to specify jumps in concentrations
        # e.g. during fed-batch: we would obtain a calculated concentration
        # for a time point, then after the calculation, we would add more
        # nutrient: the next time point will start from the new concentration,
        # not the one calculated after the last step
        # 
        # If we didn't do this in two steps, then we would compute a rate
        # that led to the modified concentration, without considering that
        # there has been an external contribution. This way we define a
        # "virtual" time step of length zero that does not require a rate
        # calculation.
        if (time[i] == time[i+1]) next        
        
        delta <- c(delta, (conc[i+1] - conc[i]) / (time[i+1] - time[i]) );
        t <- c(t, time[i])
    }
    delta <- c(delta, delta[n-1]) # maintain last value; this should not be 
    t <- c(t, t[n])		  # used, as it cannot be computed, we only
    			  	  # provide it to keep the number of rows!
    
    return( data.frame(time=t, rate=delta) )
}


# A function to convert "observed" growth values (note that they might as
# well be the result of interpolation) to rates that can be used to drive 
# DFBA. This assumes we start at time zero
obs2mu <- function(x, y) {
    n <- length(x)
    mu <- 1:n
    
    # we cannot compute mu[n] because there is no subsequent value
    for (i in 1:n-1) {
    	# this should allow us to specify jumps
        if (x[i] == x[i+1]) next
	mu[i] = (ln(y[i+1]) - ln(y[i])) / (x[i+1] - x[i])
    }
    mu[n] = mu[n-1]	# end with a non-zero value

    return( data.frame(x, y=mu) )
}


loadDynamicRates <- function(data='', verbose=0) {
    if (data == '') {
    	return (NULL)
    }
    if (file.exists(data)) {
    	dat <- read.table(data, sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.dat', sep=""))) {
    	dat <- read.table(paste(data, '.dat', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.tsv', sep=""))) {
    	dat <- read.table(paste(data, '.tsv', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.tab', sep=""))) {
    	dat <- read.table(paste(data, '.tab', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.csv', sep=""))) {
        dat <- read.table(paste(data, '.csv', sep=""), sep=',', header=TRUE, check.names='FALSE')
    } else if ((is.element('xlsx', installed.packages()[,1])) &&
             (file.exists(paste(data, ".xlsx", sep="")))) {
        dat <- read.xlsx(paste(data, ".xlsx", sep=""), sheet.index=1)
    } else {
        cat("\n\ngetRateChanges: File", data, "does not exist!\n\n\n")
    	dat <- read.table(file.choose(), header=FALSE)
    }
    
    # provide feedback
    if (verbose > 0) {
        cat("\nloadRateChanges: reading", data, "\n\n")
        print(dat)
    }
    return(dat)
}

interpolateDynamicRates <- function(dynamic.changes, tStep=1, nSteps=60, 
		interpol='linear',
                biomassInterpol='linear',
		biomassRxn='',
                verbose=0) {
    dat <- dynamic.changes
    
    # get list of reaction names
    reactions <- names(dat)
    reactions <- reactions[ -which(reactions %in% 'Time') ]

    tStep = as.numeric(tStep)	# in case it is passed as a string
    nSteps = as.numeric(nSteps)	# in case it is passed as a string

    # create output data frame with target times
    out = data.frame( Time = seq(1, nSteps * tStep, by=tStep))

    x <- dat['Time'][,1]	# Prepare a vector with time values
    # We can
    #	get all as rates (use linear/spline interpolation)
    #	get biomass as observed (sigmoid) and others as rates
    #	get all as observed (should convert all to rates)
    for (r in reactions) {
	# These give pairs x,y (time, r)
	
        y <- dat[r][,1]		# Prepare a vector with reaction values
        # Check if this is the Biomass reaction (which we can expect to
        # follow a well-defined sigmoid or logistic curve).
        #    We need to check for biomass, biomass.low. and biomass.upp.
        #    since the three are valid.
        #
        # Note that if biomassRxn is not found, then the general interpolation
        # scheme will be used.
        if ((r == biomassRxn || 
             r == paste(biomassRxn, '.low.', sep="") ||
             r == paste(biomassRxn, '.upp.', sep=""))) {
	    if (biomasssInterpol == 'sigmoid') {
                # obsdata <- interpol_sigmoid(x, y, time)
                # interdata <- obs2mu(obsdata['x'][,1], obsdata['y'][,1])
                #
                # do sigmoid interpolation followed by conversion to rate
	        mdl <- nls(y ~ a + ((b - a) / (1 + exp(-c * (x - d)))), 
	            start=list(a = min(y), b=max(y), c= 1, d=median(x)), 
		    trace=FALSE, algorithm="port")
		p <- coef(mdl)
		# Sigmoid interpolation only makes sense if Biomass
                # has been given as observed values
		obsdata <- data.frame(
			x = seq(1, nSteps * tStep, by=tStep),
			y = sigm(p, seq(1, nSteps * tStep, by=tStep))
		)
                # convert to growth rates
		interdata <- obs2mu(obsdata['x'][,1], obsdata['y'][,1])
                
	    } else if (biomassInterpol == 'lognormal') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE
                library(MASS)
                # obtain fit 
                fit <- fitdistr(y, "lognormal")
                if (verbose > 0) {
                    cat('\ngetRateChanges: lognormal fit for Biomass\n\n')
                    print(summary(fit))
                }
                
            } else if (biomassInterpol == 'lognormal2') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE 
                library(fitdistrplus)

                # obtain fit
                fit <- fitdist(y, "lnorm")
                if (verbose > 0) {
                    cat('\ngetRateChanges: goodness of lognormal2 fit for Biomass\n\n')
                    print(gofstat(fit))
                    print(summary(fit))
                }
                
            } else if (biomassInterpol == 'lognormal3') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE 
                library(vcd)
                library(VGAM)
                fit <- vglm(y ~ 1, lognormal3, trace=TRUE, crit="c")
                if (verbose > 0) {
                    cat('\ngetRateChanges: lognormal3 fit for Biomass\n\n')
		    coef(fit)
                    Coef(fit)
                    summary(fit)
                    
                }
                
            } else if (biomassInterpol == 'poisson-lognormal') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE 
                # fit to a Poisson log-normal using ML
                library(MASS)
                library(permute)
                library(vegan)

    	        # NEED TO FIND OUT BEST WAY TO INPUT DATA XXX JR XXX
	        #data <- y
    	        #data <- as.fisher(y)

    	        # fisherfit fits Fisher's logseries to abundance data. Function
	        # prestonfit groups species into double octave classes and fits
	        # Preston's lognormal model, and function prestondistr fits the
	        # truncated lognormal model without pooling data into octaves
		mod <- fisherfit(data)
                mod
	        # prestonfit needs large samples
	        mod.oct <- prestonfit(ColSums(data))
	        mod.ll  <- prestondistr(ColSums(data))
	        mod.oct
	        mod.ll
	        plot(mod.oct)
	        lines(mod.ll, line.col="blue3") # Different
    	        ## Smoothed density
    	        den <- density(log2(ColSums(data)))
    	        lines(den$x, ncol(data)*den$y, lwd=2) # Fairly similar to mod.oct
    	        ## Extrapolated richness
    	        veiledspec(mod.oct)
    	        veiledspec(mod.ll)

	        summary(mod)
            } else if (biomassInterpol == 'logistic') {
                # obsdata <- interpol_logistic(x, y, time)
                # interdata <- obs2mu(obsdata$x, obsdata$y)
            
                # WORK IN PROGRESS!!!
                # Method of Hall et al. 2014 (it is possible to modify the
                # default parameters h and quota, but we won't
                # UNTESTED INCOMPLETE DO NOT USE 
                library(growthrates)
		
                # NOT GOOD ENOUGH FOR NOW. A linear fit of a logarithmic scale
                # fit <- fit_easylinear(x, y)
                # par(mfrow = c(1, 2))
                # plot(fit, log = "y")
                # plot(fit)

                
                # A parametric growth model consists of a mathematical formula
                # that describes the growth of a population (e.g.
                # grow_logistic) and its parameters (e.g. y0, mumax, and K,).
                # Fitting a parametric model is the process of estimating an
                # optimal parameter set that minimizes a given
                # quality criterion. Here we use the method of least squares,
                # also known as ordinary least squares regression (OLS). As most
                # of the growth models are non-linear, we need always a good set
                # of start parameters p. It is wise to choose values for start
                # parameters carefully by considering the main properties of
                # the selected growth model (e.g. that the carrying capacity K
                # should be around the observed maximum of the data), or by
                # experimentation, i.e. plotting the model together with the
                # data. In order to prevent unrealistic (e.g.
                # negative) parameter values, it is optionally possible to
                # specify box-constraints (upper and lower). For difficult
                # problems one may consider to change the involved model fitting
                # algorithm from Marquardt ("Marq") to something else, e.g.
                # to"L-BFGS-B". Details can be found on the ?modFit help page.
                
                #p     <- c(y0 = 0.01, mumax = 0.2, K = 0.1)
                #lower <- c(y0 = 1e-6, mumax = 0,   K = 0)
                #upper <- c(y0 = 0.05, mumax = 5,   K = 0.5)

                p     <- c(y0 = min(y), K=max(y), mumax=((max(y)-min(y))/2))
                lower <- c(y0 = min(y)/5, K=max(y)/5, mumax=((max(y)-min(y))/20))
                upper <- c(y0 = min(y)*5, K=max(y)*5, mumax=((max(y)-min(y))*10))
                
                fit <- fit_growthmodel(FUN = grow_logistic, p = p, x, y,
                                        lower = lower, upper = upper)
		if (verbose > 0) {
                    print(coef(fit))
                }
		# use formula to get Y values
                #y0 <- coef(fit)["y0"]
                #mumax <- coef(fit)["mumax"]
                #K <- coef(fit)["K"]
                time <- seq(1, nSteps * tStep, by=tStep)
                #y = (K * y0) / (y0 + (K - y0) * exp(-mumax * time))
		interdata <- obs2mu(grow_logistic(time, coef(fit)))

            } else if (biomassInterpol == 'twostep') {
                # obsdata <- interpol_twostep(x, y, time)
                # interdata <- obs2mu(obsdata$x, obsdata$y)

                # WORK IN PROGRESS!!!
                # Method of Hall et al. 2014 (it is possible to modify the
                # default parameters h and quota, but we won't
                # UNTESTED INCOMPLETE DO NOT USE 
                library(growthrates)
		

                p     <- c(yi = 0.02, ya = 0.001, kw = 0.1, mumax = 0.2, K = 0.1)
                lower <- c(yi = 1e-6, ya = 1e-6, kw = 0,    mumax = 0,   K = 0)
                upper <- c(yi = 0.05, ya = 0.05, kw = 10,   mumax = 5,   K = 0.5)

                fit <- fit_growthmodel(FUN = grow_twostep, p = p, time = x, y = y,
                                        lower = lower, upper = upper)

                interdata <- obs2mu(grow_twostep(time, coef(fit))[ ,'y'])

		if (verbose > 0) {
                    print(coef(fit))
		    par(mfrow = c(1, 2))

                    plot(fit)
                    lines(fit, col = "red")

                    plot(fit, log = "y")
                    lines(fit, col = "red")
		}
                # Despite the fact that the above model is solved as a 
                # differential equation, the relatively high number of
                # parameters may need special care, too. In such cases, package
                # growthrates allows to fit subsets of parameters while setting
                # the others to fixed values. In the following, this is done by
                # specifying a subset without initial abundances y_a and y_i
		# in which:

		fit <- fit_growthmodel(FUN = grow_twostep, p = p, 
                	time = x, y = y,
                        lower = lower, upper = upper, 
                        which = c("kw", "mumax", "K"))

		interdata <- obs2mu(grow_twostep(time, coef(fit))[ ,'y'])

		if (verbose > 0) {		
                    summary(fit)
		    coef(fit)
                    plot(fit)
                }
                
                # Smoothing splines are a quick method to estimate maximum 
                # growth. The method is called nonparametric, because the
                # growth rate is directly estimated from the smoothed data
                # without being restricted to a specific model formula.

                ## automatic smoothing with cv
                #res <- fit_spline(x, y)

                #par(mfrow = c(1, 2))
                #plot(res, log = "y")
                #plot(res)
                
                # also check grow_baranyi(), grow_gompertz(), 
                # grow_exponential(), grow_richards(), ode_genlogistic(),
                # ode_twostep(); and growthmodel() to allow user-specified
                # growth model functions.
                
            } else if (biomassInterpol == 'spline') {
	        interdata <- data.frame(
		    spline(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
            } else if (biomassInterpol == 'linear') {
 	        interdata <- data.frame(
		    approx(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
            }
        } else {
	    # For general rates, we cannot assume a sigmoid distribution
	    # use either linear or spline interpolation
	    # In principle, linear is recommended as it makes no assumptions.
	    if (interpol == 'linear') {
	        interdata <- data.frame(
		    approx(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
	    } else if (interpol == 'spline') {
	        interdata <- data.frame(
		    spline(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
	    }
	}
	# add the corresponding column to the output data frame
	out[r] <- interdata['y']
    }
    cat("\n\nWE WILL USE THE FOLLOWING EXTENDED DYNAMIC CONSTRAINTS:\n\n\n")
    print(out) 
    #plot(out)

    # return the data frame
    return( out )
}

# Get the changes to constraints that have to be applied to the model
# during the simulation. This is used when we know in advance how reaction
# rates will change at each time point.
# For ease of use, we do not require ALL time points to be specified, but
# only key values, and will interpolate the rest.
# NOTE that we do not have a default interpolation (last else) in the if
# cascades: if an "invalid" interpolation (e.g. "none") is specified, then
# no interpolation will occur.
# Valid options are "linear" and "spline" for interpol and "sigmoid",
# "linear" and "spline" for biomassInterpol. Note that "sigmoid" only makes
# sense for Biomass if we are given observed biomass values, not growth 
# rates.
#
# Currently we only get rates and use linear interpolation for everything.
#
getRateChanges <- function(data, tStep=1, nSteps=60, 
		interpol='linear',
                biomassInterpol='linear',
		biomassRxn='',
                verbose=0) {
    if (missing(biomassRxn)) {
        biomassRxn=''
    }
    if (data == '') {
    	return (NULL)
    }
    if (file.exists(data)) {
    	dat <- read.table(data, sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.dat', sep=""))) {
    	dat <- read.table(paste(data, '.dat', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.tsv', sep=""))) {
    	dat <- read.table(paste(data, '.tsv', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.tab', sep=""))) {
    	dat <- read.table(paste(data, '.tab', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.csv', sep=""))) {
        dat <- read.table(paste(data, '.csv', sep=""), sep=',', header=TRUE, check.names='FALSE')
    } else {
        cat("\n\ngetRateChanges: File", data, "does not exist!\n\n\n")
	return(NULL)
    }
    
    # provide feedback
    cat("\ngetRateChanges: processing", data, "\n\n")
    print(dat)
    
    # get list of reaction names
    reactions <- names(dat)
    reactions <- reactions[ -which(reactions %in% 'Time') ]

    tStep = as.numeric(tStep)	# in case it is passed as a string
    nSteps = as.numeric(nSteps)	# in case it is passed as a string

    # create output data frame with target times
    out = data.frame( Time = seq(1, nSteps * tStep, by=tStep))

    x <- dat['Time'][,1]	# Prepare a vector with time values
    # We can
    #	get all as rates (use linear/spline interpolation)
    #	get biomass as observed (sigmoid) and others as rates
    #	get all as observed (should convert all to rates)
    for (r in reactions) {
	# These give pairs x,y (time, r)
	
        y <- dat[r][,1]		# Prepare a vector with reaction values
        # Check if this is the Biomass reaction (which we can expect to
        # follow a well-defined sigmoid or logistic curve).
        #    We need to check for biomass, biomass.low. and biomass.upp.
        #    since the three are valid.
        #
        # Note that if biomassRxn is not found, then the general interpolation
        # scheme will be used.
        if ((r == biomassRxn || 
             r == paste(biomassRxn, '.low.', sep="") ||
             r == paste(biomassRxn, '.upp.', sep=""))) {
	    if (biomasssInterpol == 'sigmoid') {
                # do sigmoid interpolation followed by conversion to rate
	        mdl <- nls(y ~ a + ((b - a) / (1 + exp(-c * (x - d)))), 
	            start=list(a = min(y), b=max(y), c= 1, d=median(x)), 
		    trace=FALSE, algorithm="port")
		p <- coef(mdl)
		# Sigmoid interpolation only makes sense if Biomass
                # has been given as observed values
		obsdata <- data.frame(
			x = seq(1, nSteps * tStep, by=tStep),
			y = sigm(p, seq(1, nSteps * tStep, by=tStep))
		)
                # convert to growth rates
		interdata <- obs2mu(obsdata['x'][,1], obsdata['y'][,1])
                
	    } else if (biomassInterpol == 'lognormal') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE
                library(MASS)
                # obtain fit 
                fit <- fitdistr(y, "lognormal")
                if (verbose > 0) {
                    cat('\ngetRateChanges: lognormal fit for Biomass\n\n')
                    print(summary(fit))
                }
                
            } else if (biomassInterpol == 'lognormal2') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE 
                library(fitdistrplus)

                # obtain fit
                fit <- fitdist(y, "lnorm")
                if (verbose > 0) {
                    cat('\ngetRateChanges: goodness of lognormal2 fit for Biomass\n\n')
                    print(gofstat(fit))
                    print(summary(fit))
                }
                
            } else if (biomassInterpol == 'lognormal3') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE 
                library(vcd)
                library(VGAM)
                fit <- vglm(y ~ 1, lognormal3, trace=TRUE, crit="c")
                if (verbose > 0) {
                    cat('\ngetRateChanges: lognormal3 fit for Biomass\n\n')
		    coef(fit)
                    Coef(fit)
                    summary(fit)
                    
                }
                
            } else if (biomassInterpol == 'poisson-lognormal') {
                # WORK IN PROGRESS!!!
                # UNTESTED INCOMPLETE DO NOT USE 
                # fit to a Poisson log-normal using ML
                library(MASS)
                library(permute)
                library(vegan)

    	        # NEED TO FIND OUT BEST WAY TO INPUT DATA XXX JR XXX
	        #data <- y
    	        #data <- as.fisher(y)

    	        # fisherfit fits Fisher's logseries to abundance data. Function
	        # prestonfit groups species into double octave classes and fits
	        # Preston's lognormal model, and function prestondistr fits the
	        # truncated lognormal model without pooling data into octaves
		mod <- fisherfit(data)
                mod
	        # prestonfit needs large samples
	        mod.oct <- prestonfit(ColSums(data))
	        mod.ll  <- prestondistr(ColSums(data))
	        mod.oct
	        mod.ll
	        plot(mod.oct)
	        lines(mod.ll, line.col="blue3") # Different
    	        ## Smoothed density
    	        den <- density(log2(ColSums(data)))
    	        lines(den$x, ncol(data)*den$y, lwd=2) # Fairly similar to mod.oct
    	        ## Extrapolated richness
    	        veiledspec(mod.oct)
    	        veiledspec(mod.ll)

	        summary(mod)
            } else if (biomassInterpol == 'logistic') {
                # WORK IN PROGRESS!!!
                # Method of Hall et al. 2014 (it is possible to modify the
                # default parameters h and quota, but we won't
                # UNTESTED INCOMPLETE DO NOT USE 
                library(growthrates)
		
                # NOT GOOD ENOUGH FOR NOW. A linear fit of a logarithmic scale
                # fit <- fit_easylinear(x, y)
                # par(mfrow = c(1, 2))
                # plot(fit, log = "y")
                # plot(fit)

                
                # A parametric growth model consists of a mathematical formula
                # that describes the growth of a population (e.g.
                # grow_logistic) and its parameters (e.g. y0, mumax, and K,).
                # Fitting a parametric model is the process of estimating an
                # optimal parameter set that minimizes a given
                # quality criterion. Here we use the method of least squares,
                # also known as ordinary least squares regression (OLS). As most
                # of the growth models are non-linear, we need always a good set
                # of start parameters p. It is wise to choose values for start
                # parameters carefully by considering the main properties of
                # the selected growth model (e.g. that the carrying capacity K
                # should be around the observed maximum of the data), or by
                # experimentation, i.e. plotting the model together with the
                # data. In order to prevent unrealistic (e.g.
                # negative) parameter values, it is optionally possible to
                # specify box-constraints (upper and lower). For difficult
                # problems one may consider to change the involved model fitting
                # algorithm from Marquardt ("Marq") to something else, e.g.
                # to"L-BFGS-B". Details can be found on the ?modFit help page.
                
                #p     <- c(y0 = 0.01, mumax = 0.2, K = 0.1)
                #lower <- c(y0 = 1e-6, mumax = 0,   K = 0)
                #upper <- c(y0 = 0.05, mumax = 5,   K = 0.5)

                p     <- c(y0 = min(y), K=max(y), mumax=((max(y)-min(y))/2))
                lower <- c(y0 = min(y)/5, K=max(y)/5, mumax=((max(y)-min(y))/20))
                upper <- c(y0 = min(y)*5, K=max(y)*5, mumax=((max(y)-min(y))*10))
                
                fit <- fit_growthmodel(FUN = grow_logistic, p = p, x, y,
                                        lower = lower, upper = upper)
		if (verbose > 0) {
                    print(coef(fit))
                }
		# use formula to get Y values
                #y0 <- coef(fit)["y0"]
                #mumax <- coef(fit)["mumax"]
                #K <- coef(fit)["K"]
                #time <- seq(1, nSteps * tStep, by=tStep)
                #y = (K * y0) / (y0 + (K - y0) * exp(-mumax * time))
		interdata <- obs2mu(grow_logistic(time, coef(fit)))

            } else if (biomassInterpol == 'twostep') {
                # WORK IN PROGRESS!!!
                # Method of Hall et al. 2014 (it is possible to modify the
                # default parameters h and quota, but we won't
                # UNTESTED INCOMPLETE DO NOT USE 
                library(growthrates)
		

                p     <- c(yi = 0.02, ya = 0.001, kw = 0.1, mumax = 0.2, K = 0.1)
                lower <- c(yi = 1e-6, ya = 1e-6, kw = 0,    mumax = 0,   K = 0)
                upper <- c(yi = 0.05, ya = 0.05, kw = 10,   mumax = 5,   K = 0.5)

                fit <- fit_growthmodel(FUN = grow_twostep, p = p, time = x, y = y,
                                        lower = lower, upper = upper)

                interdata <- obs2mu(grow_twostep(time, coef(fit))[ ,'y'])

		if (verbose > 0) {
                    print(coef(fit))
		    par(mfrow = c(1, 2))

                    plot(fit)
                    lines(fit, col = "red")

                    plot(fit, log = "y")
                    lines(fit, col = "red")
		}                
                # Despite the fact that the above model is solved as a 
                # differential equation, the relatively high number of
                # parameters may need special care, too. In such cases, package
                # growthrates allows to fit subsets of parameters while setting
                # the others to fixed values. In the following, this is done by
                # specifying a subset without initial abundances y_a and y_i
		# in which:

		fit <- fit_growthmodel(FUN = grow_twostep, p = p, 
                	time = x, y = y,
                        lower = lower, upper = upper, 
                        which = c("kw", "mumax", "K"))

		interdata <- obs2mu(grow_twostep(time, coef(fit))[ ,'y'])

		if (verbose > 0) {		
                    summary(fit)
		    coef(fit)
                    plot(fit)
                }
                # Smoothing splines are a quick method to estimate maximum 
                # growth. The method is called nonparametric, because the
                # growth rate is directly estimated from the smoothed data
                # without being restricted to a specific model formula.

                ## automatic smoothing with cv
                #res <- fit_spline(x, y)

                #par(mfrow = c(1, 2))
                #plot(res, log = "y")
                #plot(res)
                
                # also check grow_baranyi(), grow_gompertz(), 
                # grow_exponential(), grow_richards(), ode_genlogistic(),
                # ode_twostep(); and growthmodel() to allow user-specified
                # growth model functions.
                
            } else if (biomassInterpol == 'spline') {
	        interdata <- data.frame(
		    spline(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
            } else if (biomassInterpol == 'linear') {
 	        interdata <- data.frame(
		    approx(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
            }
        } else {
	    # For general rates, we cannot assume a sigmoid distribution
	    # use either linear or spline interpolation
	    # In principle, linear is recommended as it makes no assumptions.
	    if (interpol == 'linear') {
	        interdata <- data.frame(
		    approx(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
	    } else if (interpol == 'spline') {
	        interdata <- data.frame(
		    spline(x, y, xout=seq(1, nSteps * tStep, by=tStep))
		)
	    }
	}
	# add the corresponding column to the output data frame
	out[r] <- interdata['y']
    }
    cat("\n\nWE WILL USE THE FOLLOWING EXTENDED DYNAMIC CONSTRAINTS:\n\n\n")
    print(out) 
    #plot(out)

    # return the data frame
    return( out )
}

loadNutrientDeltas <- function() {
    return()
}

# Get modifications to nutrient concentrations as a statically defined
# function of time (i.e. we know in advance at which time the concentration
# will be modified and by how much
#
getNutrientChanges <- function(data, tStep=1, nSteps=60) {
    if (data == '') {
    	return (NULL)
    }
    if (file.exists(data)) {
    	dat <- read.table(data, sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.dat', sep=""))) {
    	dat <- read.table(paste(data, '.dat', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.tsv', sep=""))) {
    	dat <- read.table(paste(data, '.tsv', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.tab', sep=""))) {
    	dat <- read.table(paste(data, '.tab', sep=""), sep='\t', header=TRUE, check.names='FALSE')
    } else if (file.exists(paste(data, '.csv', sep=""))) {
        dat <- read.table(paste(data, '.csv', sep=""), sep=',', header=TRUE, check.names='FALSE')
    } else {
        cat("\n\ngetNutrientChanges: File", data, "does not exist!\n\n\n")
	return(NULL)
    }
    
    # provide feedback
    cat("\ngetNutrientChanges: processing", data, "\n\n")
    print(dat)
    
    # get list of nutrients
    nutrients <- names(dat)
    nutrients <- nutrients[ -which(nutrients %in% 'Time') ]

    tStep = as.numeric(tStep)	# in case it is passed as a string
    nSteps = as.numeric(nSteps)	# in case it is passed as a string

    # create output data frame with target times
    # starting at first time step
    out = data.frame( Time = seq(tStep, nSteps * tStep, by=tStep))
    #out = data.frame( Time = seq(1, nSteps * tStep, by=tStep))

    x <- dat['Time'][,1]	# Prepare a vector with time values

    for (n in nutrients) {
        # create an empty vector of nSteps zeros
        changes=c(rep(0, nSteps))
        # set the values we got in 'dat'
        changes[x]=dat[,n]
        # join the vector with the output data frame
        out[ paste('EX_', n, '(e)', sep='') ] <- changes
    }
    
    cat('\n\nWE WILL USE THE FOLLOWING NUTRIENT CHANGES:\n\n\n')
    print(out)
    
    return(out)
}

# 
# Standard FBA 
# ------------ 
# 
# Perform  ux-balance analysis (FBA) by using method optimizeProb of 
# class modelorg . Method optimizeProb performs flux-balance analysis 
# [Edwards et al., 2002a, Orth et al., 2010b]. 
#
#
fba <- function(model, outDir='.', verbose=0)
{

    cat('\nStandard Flux Balance Analysis(FBA)')
    cat('\n===================================\n')

    #optL <- optimizeProb(model, algorithm = "fba", retOptSol = FALSE);
    #print(optL)
    #print(optL$fluxes)
    #print(optL$fluxes[react_id(model) == "Biomass_SLI"])
    #
    opt <- optimizeProb(model, algorithm = "fba", retOptSol = TRUE);
    #
    cat('\nResult\n')
    print( checkOptSol(opt) )

    cat('\nValue of the objective function:\n')
    print( lp_obj(opt) )

    # The function summaryOptsol() returns an object of class optsolSummary  
    # and needs the object of class optsol (results of simulations) and the 
    # corresponding object of class modelorg (the entire metabolic network).  
    # The generated object of class summaryOptsol contains some information 
    # about the flux distribution, substrates, products and limiting reactions. 
    #
    # The method printExchange() prints a subset of the flux distribution for 
    # the exchange reactions in the model.  Each column represents the 
    # environment of one optimization. The symbol "-" indicates
    # that the corresponding metabolite is imported (is a substrate); 
    # the symbol "+" indicates, that the corresponding metabolite is 
    # excreted (is a product).
    #
    cat('\nSummary:\n')
    sum <- summaryOptsol(opt, model)
    print(sum)
    # or only exchange reactions:
    #printExchange(sum, dense = TRUE)

    # Save FBA result to file

    # return value of the objective function
    return (lp_obj(opt))
}

# Minimize total flux
# -------------------
#
# Usually, an FBA solution is not unique.  There can be many equivalent flux
# distributions supporting the same objective value.  A method to decide for 
# one out of these solutions is to compute the flux distribution minimizing 
# the total absolute flux (MTF) but still supporting the objective value of 
# the FBA solution.  At first, an objective value, for example calculated 
# via FBA, is required
#
mtf <- function(model, outDir='.', verbose=0) {
    modelName <- model@mod_name
    dir.create(outDir, showWarnings = FALSE)
    mtfDir=paste(outDir, "/mtf", sep="")
    dir.create(mtfDir, showWarnings = FALSE)
    tabMTFvalues <- paste(mtfDir, "/", modelName, "_s_MTF.tab", sep="");

    cat('\nMinimize total flux')
    cat('\n===================\n')

    cat('\nObtain an initial optimized value using FBA\n');
    fba <- optimizeProb(model, algorithm="fba");
    #
    cat('\nOptimized value of the objective function\n')
    print( mod_obj(fba) )
    #
    cat('\nNumber of variables (fba)\n')
    print( nvar(fluxdist(fba)) )
    #
    cat('\nUsing FBA value to obtain MTF\n')
    mtf_sol <- optimizeProb(model, algorithm="mtf", wtobj = mod_obj(fba));
    #
    cat('\nValue of the objective function for the MTF:\n')
    print( lp_obj(mtf_sol) )
    #
    cat('\nNumber of variables (mtf)\n')
    print( nvar(fluxdist(mtf_sol)) )
    #
    cat('\nFlux distribution of the MTF solution\n')
    fl <- getFluxDist(mtf_sol);

    cat('\nLength of the flux distribution:\n')
    print( length(fl) )

    if (verbose > 2) {
        print(fl)
    }

    # SAVE MTF VALUES TO FILE

    cat(paste('\nSaving optimized output from MTF to "', tabMTFvalues,'"\n', sep=""))
    #print( lp_obj(fl) )
    nreact=length(fl)
    for (i in 1:nreact) {
        # output as CSV to log file
        #print(paste('MTF,', i, ',', react_id(model)[i], ',', 
    	#	    fl[i]))
        out <- paste(i, '	', react_id(model)[i], '	', 
    		    fl[i])
        cat(out, file=tabMTFvalues, sep="\n", append=TRUE)
      }
    # END SAVE MTF VALUES

    if (verbose > 2) {
	cat('\nFlux distribution of ALL the reactions\n')
	print(getFluxDist(mtf_sol, checkReactId(model, react_id(model))))
    }
    #
    exchngRxns <- findExchReact(model)
    uptakeRxns <- uptReact(exchngRxns)

    if (verbose > 1) {
	cat('\nFlux distribution of the exchange reactions\n')
	cat(  '-------------------------------------------\n')
	fd <- getFluxDist(mtf_sol, exchngRxns);
	print(fd)
    }
    # 
    if (verbose > 0) {
	cat('\nNet flux of the exchange reactions\n');
	cat(  '----------------------------------\n');
	print( getNetFlux(fd) )
    }
    #
    cat('\nValue of the objective function in the MTF model:\n')
    print( mod_obj(mtf_sol) )
    #
    cat('\nMTF Summary:\n')
    cat(  '------------\n');
    sum <- summaryOptsol(mtf_sol, model)
    if (verbose > 2) {
	print(sum)
    }
    # or only exchange reactions:
    cat('\nSummary of the exchange reactions\n')
    printExchange(sum, dense = TRUE)

    # return the value of the objective function
    #return ( mod_obj(mtf_sol) )
    # return the solution
    return (mtf_sol)
}

# Flux Variability Analysis (FVA)
# -------------------------------
#
# FBA only returns a single flux distribution that corresponds to maximal growth
# under given growth conditions. However, alternate optimal solutions may exist
# which correspond to maximal growth. FVA calculates the full range of numerical
# values for each reaction flux within the network
# 
# The function fluxVar performs a flux variability analysis with a given model 
# [Mahadevan and Schilling, 2003].  The minimum and maximum flux values for 
# each reaction in the model are calculated, which still support a certain 
# percentage of a given optimal functional state Z_opt
#
fva <- function(model, outDir='.', verbose=0) {
    modelName <- model@mod_name
    dir.create(outDir, showWarnings = FALSE)
    fvaDir=paste(outDir, "/fva", sep="")
    dir.create(fvaDir, showWarnings = FALSE)
    pngFVAplot <- paste(fvaDir, '/', modelName, "_FVA.png", sep="")
    tabFVAvalues <- paste(fvaDir, '/', modelName, "_s_FVA.tab", sep="");


    cat('\nFlux Variability Analysis (FVA)')
    cat('\n===============================\n')

    ranges <- fluxVar(model, percentage=80, verboseMode=verbose)
    #
    cat('\nSee plot to visualize minimum and maximum flux values for each reaction\n')
    png(pngFVAplot)
    plot(ranges)
    dev.off()

    # SAVE FVA VALUES
    cat('\nSummary:\n')
    sum <- summaryOptsol(ranges, model)
    print(sum)
    # or only exchange reactions:
    #printExchange(sum, dense = TRUE)

    if (verbose > 3) {
        #capture.output(print(sum), file="FVA.log", split=TRUE)
        print(ranges)
    }
    cat(paste('\nSaving optimized output from FVA to "', tabFVAvalues,'"\n', sep=""))
    # print( lp_obj(ranges) )
    cat(paste('\nTotal length: ', length((ranges)), 
    	      '\nNo. of reactions: ', length(ranges)/2, '\n\n'))

    nreact=length(ranges)/2
    for (i in 1:nreact) {
        if (verbose > 2) {
            #print(paste(i, react_id(model)[i], ': min: ', lp_obj(ranges)[i]))
            #print(paste(i, react_id(model)[i],': max: ', lp_obj(ranges)[i+nreact]))
            print(paste('FVA,', i, ',', react_id(model)[i], ',', 
    			lp_obj(ranges)[i], ',', lp_obj(ranges)[i+nreact]))
        }
        out <- paste(i, '	', react_id(model)[i], '	', 
    		    lp_obj(ranges)[i], '	', lp_obj(ranges)[i+nreact])
        cat(out, file=tabFVAvalues, sep="\n", append=TRUE)
      }
    # END SAVE FVA VALUES

}

# Robustness analysis
# -------------------
#
# The function robAna performs a robustness analysis with a given model.  
# The flux of a control reaction will be varied stepwise between the maximum 
# and minimum value the flux of the control reaction can reach [Palsson, 2006].
#
#
ra <- function(model, reaction, outDir='.', verbose=0) {
    modelName <- model@mod_name
    dir.create(outDir, showWarnings = FALSE)
    raDir=paste(outDir, "/ra", sep="")
    dir.create(raDir, showWarnings = FALSE)
    pngRAplot <- paste(raDir, "/", modelName, "_RA.png", sep="")

    cat('\nRobustness Analysis')
    cat('\n===================\n')

    cat('\nSee plot to visualize RA of ', reaction, '\n')
    # ideally this would be a 'for' loop over all substrateRxns
    opt <- robAna(model, ctrlreact=reaction, verboseMode=verbose)

    png(pngRAplot)
    plot(opt)
    dev.off()

}


# Phenotypic phase plane analysis
# -------------------------------
#
# The function phpp performs a phenotypic phase plane analysis [Edwards et al.,
# 2001, 2002b] with a given model.  The flux of two control reactions will be 
# varied stepwise between a given maximum and minimum value.  
#
# We should select a meaningful pair of reactions here (e.g. BIOMASS/SECRETION)
pppa <- function(model, rxns, outDir='.', verbose=0) {
    modelName <- model@mod_name
    
    dir.create(outDir, showWarnings = FALSE)
    pppaDir=paste(outDir, "/pppa", sep="")
    dir.create(pppaDir, showWarnings = FALSE)
    pngPPPAplot <- paste(raDir, "/", modelName, "_PPPA.png", sep="")
    
    cat('\nPhenotypic Phase Plane Analysis')
    cat('\n===============================\n')
    cat(' ', rxns[1], ' vs. ', rxns[2])
    cat('\n-------------------------------\n')

    opt <- phpp(model, 
		  ctrlreact = rxns,
		  redCosts = TRUE,
        	  numP = 25,
        	  verboseMode = verbose)
    #
    png(pngPPPAplot)

    plot(opt)
    plot(opt, rxn[1])
    plot(opt, rxn[2])
    #plot(opt, 'EX_aml(e)')

    dev.off()

}


readSubstrate <- function(substrate="") 
{
    if (substrate == "") {
        return (list(
    	    substrateRxns=c(""), 
	    initConcentrations=c(""),
    	    exclUptakeRxns=c(""))
        )
    }
    # Read subst. file (if any) and apply it
    #
    # by default, # will be used as comment.char, i.e. lines
    # starting with # will be ignored.
    #
    options(stringsAsFactors=FALSE)		# we want strings as strings
    if (file.exists(paste(substrate, ".tsv", sep=""))) {
    	subs <- read.table(paste(substrate, ".tsv", sep=""), sep="\t", header=FALSE)
    } else if (file.exists(paste(substrate, ".tab", sep=""))) {
     	subs <- read.table(paste(substrate, ".tab", sep=""), sep="\t", header=FALSE)
    } else if (file.exists(paste(substrate, ".csv", sep=""))) {
    	subs <- read.table(paste(substrate, ".csv", sep=""), sep=",", header=FALSE)
    } else if (file.exists(paste(substrate, ".txt", sep=""))) {
    	subs <- read.table(paste(substrate, ".txt", sep=""), header=FALSE)
    } else if (file.exists(substrate)) {
        # last resort, try to read it as a TAB file
        subs <- read.table(substrate, sep="\t", header=FALSE)
    } else {
        cat("File ", substrate, "{.txt, .tsv, .tab, .csv} not found\n")
    	subs <- read.table(file.choose(), header=FALSE)
    }
    
    # loop over all the substrates and apply them
    # We'll read in
    #    Metabolite	Mol.Weight	mg/L
    # and compute the concentration in mmol/L
    # When g/L <= 0 we'll take it as an excluded uptake reaction
    #
    substrateRxns <- vector()
    initConcentrations <- vector()
    exclUptakeRxns <- vector()
    
    cat('\nSetting up substrate from file ', substrate, ':\n\n')
    cat('    Met.\tMWT\tmg/L\texchange\tmmol/L\n')
    cat('    --------------------------------------------------\n')
    for ( i in 1:dim(subs)[1] ) {
        met <- subs[i,1]
        mwt <- subs[i,2]
        mg_per_l <- subs[i,3]
        exch <- met	# use this to use exchange reactions instead
        		# of metabolite names
        exch <- paste('EX_', met, '(e)', sep="")   # use this to use metabolite
        		# names instead of exchange reactions
                        # NOTE that this assumes a specific naming convention!!!
        mmol_per_l <- mg_per_l/ mwt
        cat('   ', met, '\t', mwt, '\t', mg_per_l, '\t', exch, '\t', mmol_per_l, '\n')
        if (mg_per_l >= 0) {
            substrateRxns <- c(substrateRxns, exch)
            initConcentrations <- c(initConcentrations, mmol_per_l)
        } else {
            exclUptakeRxns <- c(exclUptakeRxns, exch)
        }
    }
    cat('    --------------------------------------------------\n')
    cat('\n')
    return (list(
    	substrateRxns=substrateRxns, 
	initConcentrations=initConcentrations,
    	exclUptakeRxns=exclUptakeRxns)
    )
}


# Dynamic Flux Balance Analysis (DFBA)
# ------------------------------------
#
# Calculate concentrations of metabolites of exchange reactions at defined 
# time points given the initial concentrations. To accomplish this task 
# this function calls optimizeProb function to get the fluxes then update 
# the concentrations and the reaction boundaries ..etc
#
#

dfba <- function(model, dfbaProblem, outDir=".", verbose=0) {
    modelName <- model@mod_name
    # prepare output file names

    dir.create(outDir, showWarnings = FALSE)
    dfbaDir=paste(outDir, "/dfba", sep="")
    dir.create(dfbaDir, showWarnings = FALSE)

    pngDFBAplot <- paste(dfbaDir, '/', modelName, "_DFBA.png", sep="")
    pngDFBAlplot <- paste(dfbaDir, '/', modelName, "_DFBA-l.png", sep="")
    pngMetsAplot <- paste(dfbaDir, '/', modelName, "_MetsA.png", sep="")
    pngMetsBplot <- paste(dfbaDir, '/', modelName, "_MetsB.png", sep="")
    pngAaAplot <- paste(dfbaDir, '/', modelName, "_AaA.png", sep="")
    pngAaBplot <- paste(dfbaDir, '/', modelName, "_AaB.png", sep="")
    pngAaCplot <- paste(dfbaDir, '/', modelName, "_AaC.png", sep="")
    pngAaDplot <- paste(dfbaDir, '/', modelName, "_AaD.png", sep="")
    pngExchMetsplot <- paste(dfbaDir, '/', modelName, "_ExchMets.png", sep="")
    pngExchMetsNormalizedplot <- paste(dfbaDir, '/', modelName, "_ExchMetsNormalized.png", sep="")


    # extract control parameters
    #	XXX j XXX We should check they exist before assignment!
    substrateRxns <- dfbaProblem$substrateRxns
    initConcentrations <- dfbaProblem$initConcentrations
    initBiomass <- dfbaProblem$initBiomass
    timeStep <- dfbaProblem$timeStep
    nSteps <- dfbaProblem$nSteps
    plotRxns <- dfbaProblem$plotRxns
    exclUptakeRxns <- dfbaProblem$exclUptakeRxns
    dynamicConstraints <- dfbaProblem$dynamicConstraints
    biomassRxn <- dfbaProblem$biomassRxn


    cat('\n------------------------------------------------------------------------\n')
    cat('\n====================================\n')
    cat('\nDynamic Flux Balance Analysis (DFBA)')
    cat('\n====================================\n')

    if (verbose > 0) {
        cat('\nSubstrate Concentration\n')
        for (i in 1:length(substrateRxns))
    	    cat(substrateRxns[i], "\t", initConcentrations[i], "\n")
        cat('\nBiomass reaction', biomass, '\n')
        cat('\nInitial Biomass', initBiomass, '\n')
        cat('\nTime Step', timeStep, '\n')
        cat('\nNumber of Steps', nSteps, '\n')
        cat('\nExcluded uptake reactions\n')
        print(exclUptakeRxns)
        cat('\nPlot reactions\n')
        print(plotRxns)
    }
    
    cat('\n------------------------------------------------------------------------\n')
    cat('\nDoing DFBA\n')
    #
    df_sol <- dynamicFBA(model, substrateRxns=substrateRxns, 
		    initConcentrations=initConcentrations,
		    initBiomass = initBiomass,
		    timeStep = timeStep,
		    nSteps = nSteps,
		    exclUptakeRxns=exclUptakeRxns,
		    retOptSol=TRUE,
		    fld=TRUE,
### j		    biomassRxn=biomassRxn,
### j		    dynamicConstraints=dynamicConstraints,
		    verboseMode=verbose);

    if (verbose > 2) {
	# Plot concentrations measured by D'Huys in 2011
	cat("Plot concentrations measured by D'Huys in 2011: ")
	plotMetsA=c('EX_mnl(e)', 'EX_glc(e)', 'EX_nh4(e)')
	plotMetsB=c('EX_pyr(e)', 'EX_lac_D(e)', 'EX_akg(e)', 'EX_succ(e)')
	plotAaA=c('EX_glu_L(e)', 'EX_asp_L(e)', 'EX_ala_L(e)', 'EX_pro_L(e)')
	plotAaB=c('EX_leu_L(e)', 'EX_ile_L(e)', 'EX_val_L(e)', 'EX_met_L(e)')
	plotAaC=c('EX_ser_L(e)', 'EX_thr_L(e)', 'EX_gly(e)')
	plotAaD=c('EX_lys_L(e)', 'EX_his_L(e)', 'EX_tyr_L(e)', 'EX_phe_L(e)')
	cat("metsA ")
	png(pngMetsAplot, width=1000, height=1000)
	plot(df_sol, plotRxns=plotMetsA)
	dev.off()
	cat("metsB ")
	png(pngMetsBplot, width=1000, height=1000)
	plot(df_sol, plotRxns=plotMetsB)
	dev.off()
	cat("aaA ")
	png(pngAaAplot, width=1000, height=1000)
	plot(df_sol, plotRxns=plotAaA)
	dev.off()
	cat("aaB ")
	png(pngAaBplot, width=1000, height=1000)
	plot(df_sol, plotRxns=plotAaB)
	dev.off()
	cat("aaC ")
	png(pngAaCplot, width=1000, height=1000)
	plot(df_sol, plotRxns=plotAaC)
	dev.off()
	cat("aaD\n")
	png(pngAaDplot, width=1000, height=1000)
	plot(df_sol, plotRxns=plotAaD)
	dev.off()
    }
        
    ## Plot all exchange metabolites
    #png(pngExtraMetsplot, width=1000, height=1000)
    #plot.new()
    #plotAllExchangeMetabolites(df_sol)
    #dev.off()

    if (verbose > 1) {
	# plot active exchange metabolites
	cat("Plot active exchange metabolites\n")
	x11()
	plotExchMets(df_sol)
	png(pngExchMetsplot, width=1000, height=1000)
	plotExchMets(df_sol)
	dev.off()
    }
    
    # make a normalized plot of all exchange metabolites superposed
    # to facilitate visualization of consumption rate change 
    # relationships
    cat("Plot all exchange metabolites normalized and superposed\n")
    x11()
    plotExchMetsNormalized(df_sol)
    png(pngExchMetsNormalizedplot, width=1000, height=1000)
    plotExchMetsNormalized(df_sol)
    dev.off()

    # Create typical COBRA-Style DFBA plot
    #png(pngDFBAplot)
    x11()
    plot(df_sol, plotRxns=plotRxns);

    #dev.copy(png, file=pngDFBAplot, width=1000, height=1000)
    png(pngDFBAplot, width=1000, height=1000)
    plot(df_sol, plotRxns=plotRxns);
    dev.off()

    # Create typical COBRA-Style DFBA plot
    #png(pngDFBAplot)
    x11()
    lplotOptSolDFBA(df_sol, plotRxns=plotRxns);

    #dev.copy(png, file=pngDFBAplot, width=1000, height=1000)
    png(pngDFBAlplot, width=1000, height=1000)
    lplotOptSolDFBA(df_sol, plotRxns=plotRxns);
    dev.off()

    if (verbose > 2) {
	# plot active exchange metabolites, each one in its own file.
	cat("Plot active exchange metabolites, each one in its own file\n")
	cat("Please, wait...\n")
        if (verbose > 4)
	    # This function saves images in files named by itself!!!
	    #	That is why it wants to know te model's name to use as a prefix
	    plotExchRxns1by1(df_sol, outDir=dfbaDir, prefix=modelName, all=TRUE)
	else
	    plotExchRxns1by1(df_sol, outDir=dfbaDir, prefix=modelName)
    }

    return(df_sol)
}



# Adaptive Dynamic Flux Balance Analysis (ADFBA)
# ----------------------------------------------
#
# Calculate concentrations of metabolites of exchange reactions at defined 
# time points given the initial concentrations. To accomplish this task 
# this function calls optimizeProb function to get the fluxes then update 
# the concentrations and the reaction boundaries adapting to changing
# conditions ..etc
#
#

adfba <- function(model, adfbaProblem, outDir=".", verbose=0) {
    modelName <- model@mod_name
    dir.create(outDir, showWarnings = FALSE)

    # prepare output file names
    if ((exists("adfbaProblem$dynamicConstraints")) ||
        (! is.null(adfbaProblem$dynamicConstraints))) {
        # default output dir is 'adfba'
        adfbaDir=paste(outDir, "/adfba", sep="")   
    } else {
        # we'll behave as normal DFBA
        adfbaDir=paste(outDir, "/dfba", sep="")
        # if it didn't exist
        adfbaProblem$dynamicConstraints <- NULL
    }
    dir.create(adfbaDir, showWarnings = FALSE)

    # prepare plot filenames
    pngDFBAplot <- paste(adfbaDir, '/', modelName, "_ADFBA.png", sep="")
    pngDFBAlplot <- paste(adfbaDir, '/', modelName, "_ADFBA-l.png", sep="")
    pngMetsAplot <- paste(adfbaDir, '/', modelName, "_MetsA.png", sep="")
    pngMetsBplot <- paste(adfbaDir, '/', modelName, "_MetsB.png", sep="")
    pngAaAplot <- paste(adfbaDir, '/', modelName, "_AaA.png", sep="")
    pngAaBplot <- paste(adfbaDir, '/', modelName, "_AaB.png", sep="")
    pngAaCplot <- paste(adfbaDir, '/', modelName, "_AaC.png", sep="")
    pngAaDplot <- paste(adfbaDir, '/', modelName, "_AaD.png", sep="")
    pngExchMetsplot <- paste(adfbaDir, '/', modelName, "_ExchMets.png", sep="")
    pngExchMetsNormalizedplot <- paste(adfbaDir, '/', modelName, "_ExchMetsNormalized.png", sep="")
    concentrationsFile <- paste(adfbaDir, '/', modelName, "_concs.tab", sep="")
    fluxesFile <- paste(adfbaDir, '/', modelName, "_fluxes.tab", sep="")


    # extract control parameters
    #	XXX j XXX We should check they exist before assignment!
    substrateRxns <- adfbaProblem$substrateRxns
    initConcentrations <- adfbaProblem$initConcentrations
    initBiomass <- adfbaProblem$initBiomass
    timeStep <- adfbaProblem$timeStep
    nSteps <- adfbaProblem$nSteps
    plotRxns <- adfbaProblem$plotRxns
    exclUptakeRxns <- adfbaProblem$exclUptakeRxns
    dynamicConstraints <- adfbaProblem$dynamicConstraints
    nutrientChanges <- adfbaProblem$nutrientChanges
    biomassRxn <- adfbaProblem$biomassRxn
    method <- adfbaProblem$method


    cat('\n------------------------------------------------------------------------\n')
    cat('==============================================\n')
    cat('Adaptive Dynamic Flux Balance Analysis (ADFBA)\n')
    cat('==============================================\n')

    if (verbose > 0) {
        cat('\nSubstrate Concentration\n')
        for (i in 1:length(substrateRxns))
    	    cat(substrateRxns[i], "\t", initConcentrations[i], "\n")
        cat('\nBiomass reaction:', biomass, '\n')
        cat('\nInitial Biomass:', initBiomass, '\n')
        cat('\nTime Step:', timeStep, '\n')
        cat('\nNumber of Steps:', nSteps, '\n')
        cat('\nExcluded uptake reactions\n')
        print(exclUptakeRxns)
        cat('\nPlot reactions\n')
        print(plotRxns)
        #print(timeStep);
        #print(nSteps);
        #print(initBiomass)
    }
    
    cat('\n------------------------------------------------------------------------\n')
    cat('\nDoing ADFBA\n')
    #
# NOTE: this is only for publication, so we can provide sample datasets
#	to reproduce our results
#    Slividans <- model
#    SlividansWT <- model
#    Slividans_pIJ486 <- model
#    SlividansTNF <- model
#    SlividansAML <- model
#    SlividansDAG <- model
#    SlividansCEL <- model
#    save(model, substrateRxns, initConcentrations, initBiomass,
#    	timeStep, nSteps, exclUptakeRxns, biomassRxn, dynamicConstraints,
#        plotRxns, verbose, Slividans, file="Slividans.RData")
    af_sol <- adaptiveDFBA(model, substrateRxns=substrateRxns, 
		    initConcentrations=initConcentrations,
		    initBiomass = initBiomass,
		    timeStep = timeStep,
		    nSteps = nSteps,
		    exclUptakeRxns=exclUptakeRxns,
		    retOptSol=TRUE,
		    fld=TRUE,
		    biomassRxn=biomassRxn,
		    dynamicConstraints=dynamicConstraints,
                    nutrientChanges=nutrientChanges,
		    verboseMode=verbose, 
                    method=method);

    ### JR ### we should check if adaptiveDFBA ended correctly.

    if (verbose > 3) {
        saveConcentrations(af_sol, concentrationsFile)
        saveFluxes(af_sol, fluxesFile)
    }

    if (verbose > 2) {
	# Plot concentrations measured by D'Huys in 2011
	cat("Plot concentrations measured by D'Huys in 2011: ")
	plotMetsA=c('EX_mnl(e)', 'EX_glc(e)', 'EX_nh4(e)')
	plotMetsB=c('EX_pyr(e)', 'EX_lac_D(e)', 'EX_akg(e)', 'EX_succ(e)')
	plotAaA=c('EX_glu_L(e)', 'EX_asp_L(e)', 'EX_ala_L(e)', 'EX_pro_L(e)')
	plotAaB=c('EX_leu_L(e)', 'EX_ile_L(e)', 'EX_val_L(e)', 'EX_met_L(e)')
	plotAaC=c('EX_ser_L(e)', 'EX_thr_L(e)', 'EX_gly(e)')
	plotAaD=c('EX_lys_L(e)', 'EX_his_L(e)', 'EX_tyr_L(e)', 'EX_phe_L(e)')
	cat("metsA ")
	png(pngMetsAplot, width=1000, height=1000)
	plot(af_sol, plotRxns=plotMetsA)
	dev.off()
	cat("metsB ")
	png(pngMetsBplot, width=1000, height=1000)
	plot(af_sol, plotRxns=plotMetsB)
	dev.off()
	cat("aaA ")
	png(pngAaAplot, width=1000, height=1000)
	plot(af_sol, plotRxns=plotAaA)
	dev.off()
	cat("aaB ")
	png(pngAaBplot, width=1000, height=1000)
	plot(af_sol, plotRxns=plotAaB)
	dev.off()
	cat("aaC ")
	png(pngAaCplot, width=1000, height=1000)
	plot(af_sol, plotRxns=plotAaC)
	dev.off()
	cat("aaD\n")
	png(pngAaDplot, width=1000, height=1000)
	plot(af_sol, plotRxns=plotAaD)
	dev.off()
    }
    
    ## Plot all exchange metabolites
    #png(pngExtraMetsplot, width=1000, height=1000)
    #plot.new()
    #plotAllExchangeMetabolites(af_sol)
    #dev.off()

    if (verbose > 1) {
	# plot active exchange metabolites
	cat("Plot active exchange metabolites\n")
	x11()
	plotExchMets(af_sol)
	png(pngExchMetsplot, width=1000, height=1000)
	plotExchMets(af_sol)
	dev.off()
    }

    # make a normalized plot of all exchange metabolites superposed
    # to facilitate visualization of consumption rate change 
    # relationships
    cat("Plot all exchange metabolites normalized and superposed\n")
    x11()
    plotExchMetsNormalized(af_sol)
    png(pngExchMetsNormalizedplot, width=1000, height=1000)
    plotExchMetsNormalized(af_sol)
    dev.off()

    # Create typical COBRA-Style DFBA plot
    #png(pngDFBAplot)
    x11()
    plot(af_sol, plotRxns=plotRxns);

    #dev.copy(png, file=pngDFBAplot, width=1000, height=1000)
    png(pngDFBAplot, width=1000, height=1000)
    plot(af_sol, plotRxns=plotRxns);
    dev.off()

    # Create typical COBRA-Style DFBA plot but this time using lines 
    # instead of splines
    #png(pngDFBAplot)
    x11()
    lplotOptSolDFBA(af_sol, plotRxns=plotRxns);

    #dev.copy(png, file=pngDFBAplot, width=1000, height=1000)
    png(pngDFBAlplot, width=1000, height=1000)
    lplotOptSolDFBA(af_sol, plotRxns=plotRxns);
    dev.off()

    if (verbose > 2) {
	# plot active exchange metabolites, each one in its own file.
	cat("Plot active exchange metabolites, each one in its own file\n")
	cat("Please, wait...\n")
        if (verbose > 4)
	    # This function saves images in files named by itself!!!
	    #	That is why it wants to know te model's name to use as a prefix
	    plotExchRxns1by1(af_sol, outDir=adfbaDir, prefix=modelName, all=TRUE)
	else
	    plotExchRxns1by1(af_sol, outDir=adfbaDir, prefix=modelName)
    }

    return(af_sol)
}


saveConcentrations <- function(dfba, file="dynConcs.tab", all=FALSE)
{
    feats <- attributes(dfba)
    concs <- feats$concentrationMatrix

    if (all != TRUE) {
#        concs <- feats$concentrationMatrix[
#        	     apply(feats$concentrationMatrix[, -1], 1, 
#                          function(x) !all(x==x[1])),
#                 ]
        concs <- concs[
        	     apply(concs[, -1], 1, function(x) any(x!=x[1])),
                 ]
    }
    write.table(concs, file=file, sep='\t')
}


saveFluxes <- function(dfba, file="dynFluxes.tab", all=FALSE)
{
    feats <- attributes(dfba)
    fluxes <- feats$all_fluxes

    if (all != TRUE) {
#        fluxes <- feats$all_fluxes[
#        	      apply(feats$all_fluxes[, -1], 1, function(x) !all(x==x[1])),
#                  ]
        fluxes <- fluxes[
        	      apply(fluxes[, -1], 1, function(x) any(x!=x[1])),
                  ]
    }

    write.table(fluxes, file=file, sep='\t')
                

}


# This function makes a plot like the one defined for class optsol_dynamicFBA
# but uses lines instead of splines.
#	Splines look good but, besides being unscientific, they can led
#	to unrealistic plots as they produce additional oscillations in
#	attemting to fit the data.
#	A line plot should be less misleading.
lplotOptSolDFBA <- function(x,y,
    ylim=50,
    xlab = "",
    ylab = "Value",
    type = "p",
    pch = 20,
    col = "black",             
    # collower, colupper, pchupper, pchlower,
    # dottedline = TRUE,
    plotRxns=NULL,
    baseline = 0,
       ...) {
    # we'll prefer line plots to spline plots as the latter may
    # introduce weird biases at inflection points
    if(missing(plotRxns) || is.null(plotRxns)) {
        plot(x@timeVec,x@biomassVec,
             main='Biomass',xlab='Time',ylab=ylab);
        }
    else {
	def.par <- par(no.readonly = TRUE);
	layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE))
	#layout.show(2);
	# first plot biomass
	#plot(spline(x@timeVec,x@biomassVec, n = 201, method = "natural"), col = 1
	#                    ,main='Cell density',xlab='Time(hrs)',ylab="X(g/l)",type="l",lwd=2);
	plot(x@timeVec, x@biomassVec, 
             col=1, type="l", lwd=2,
             main='Cell density',
             xlab='Time(hrs)',
             ylab="Biomass(g/l)");
	points(x@timeVec,x@biomassVec, col = "red",lwd=2);

	# plot concentrations
	##plot(x@timeVec,2*x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
	## define min/max ()plot(x@timeVec,
	ymin <- min(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
        ymax <- max(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
	for ( i in 1:length(plotRxns) ){
	    plotInd=(x@excRxnNames %in% plotRxns[i]);
	    #print( x@concentrationMatrix[plotInd]);
	    if (i==1) {
                #plot(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"),
                plot(x@timeVec, x@concentrationMatrix[plotInd],
                     type="l", col=i, 
                     main="Concentrations",
                     ylab="mmol",
                     xlab='Time(hrs)',
                     ylim=c(ymin,ymax));
            } else {
       	        lines(x@timeVec, x@concentrationMatrix[plotInd],
                      col=i);
       	      }
	}
	#legend(1.8, ymax, plotRxns, col=1:length(plotRxns), lty=1);
	legend("left", plotRxns, col=1:length(plotRxns), lty=1);
    }

    #if (!missing(plotRxns)){		   }
}



# Plot ALL extracellular metabolites in the currently open graphic device
#	NOTE: this code may require retOptSol=FALSE
#   NOTE: superseeded by plotExcMets() below.
plotAllExchangeMetabolites <- function(dfba) {
    # dfba is an optSol object

    feats <- attributes(dfba)
    cl <- heat.colors(dim(feats$concentrationMatrix)[1])
    cl <- terrain.colors(dim(feats$concentrationMatrix)[1])
    cl <- topo.colors(dim(feats$concentrationMatrix)[1])
    cl <- cm.colors(dim(feats$concentrationMatrix)[1])
    cl <- rainbow(dim(feats$concentrationMatrix)[1])
    col <- 0
    for (i in row.names(feats$concentrationMatrix)) {
        col <- col + 1
        lines(as.matrix(feats$concentrationMatrix)[i,], col=cl[col], 
        	type='l')
    }
    legend('topleft', legend=rownames(feats$concentrationMatrix), col=cl, lwd=2)
    #matplot(as.matrix(feats$concentrationMatrix), type=c("b"))
}


# Plot active exchange metabolites on the currently open graphic device.
#	An open graphic device must exist
plotExchMets <- function(df_sol, all=FALSE) {
    times <- df_sol@timeVec
    concs <- as.matrix(df_sol@concentrationMatrix)
    # if we want to plot only non-zero ExcRxns
    if (all == FALSE) {
        concs <- concs[apply(concs[, -1], 1, function(x) any(x != x[1])),]
    }
    
    # choose color palette
    nrxns <- dim(concs)[1]
    cl <- heat.colors(nrxns)
    cl <- terrain.colors(nrxns)
    cl <- topo.colors(nrxns)
    cl <- cm.colors(nrxns)
    cl <- rainbow(nrxns)
    col <- 0    

    # find plot limits
    xmin <- min(sapply(times, function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
    xmax <- max(sapply(times, function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
    ymin <- min(sapply(concs, function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
    ymax <- max(sapply(concs, function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
    # create inital empty plot
    plot(NULL, NULL, xlim=c(xmin, xmax), ylim=c(ymin,ymax), type="l", 
         main="Concentrations", ylab = "mmol",xlab='Time(hrs)')
    # plot the metabolites
    #lines(concs[20,], type='l', col='black')
    #lines(concs[21,], type='l', col='red')
    #cat('ROW NAMES')
    #print(row.names(concs))
    for (i in row.names(concs)) {
        col <- col + 1
        lines(concs[i, ], col=cl[col], type='l', lty=col)
        #points(concs[i, ], col=cl[col], pch=col, cex=0.5)
    }

    legend('topleft', legend=rownames(concs), col=cl, lwd=1, cex=0.5)
}


# This is the same but for all exchange metabolites.
#	But in this one we do not control colors or extra labels
matplotAllExchRxns <- function(dfs) {
	matplot(t(as.matrix(dfs@concentrationMatrix)), type=c("l"))
	legend('topleft', legend=rownames(dfs@concentrationMatrix), col=cl, lwd=1, cex=0.5)
}

# It is easy to adapt it to plot only non-zero as above:
#	In this one we do not control colors or extra labels
matplotExchRxns <- function(dfs, all=FALSE) {
	concs <- as.matrix(dfs@concentrationMatrix)
	# we want to plot only non-zero ExcRxns
	if (all == FALSE) {
	    concs <- concs[apply(concs[, -1], 1, function(x) any(x != x[1])),]
	}
        matplot(t(concs), type=c("l"))
}

# Plot exchange reactions, each one in a different PNG file
#	WARNING: We use filenames hard-coded in the function!!!
plotExchRxns1by1 <- function(dfs, outDir='.', prefix="TS_PLOT", all=FALSE) {
    times <- dfs@timeVec
    concs <- as.matrix(dfs@concentrationMatrix)
    # if we want to plot only non-zero ExcRxns
    # concs <- concs[apply(concs[, -1], 1, function(x) !all(x==0)),]

    for (i in row.names(concs)) {
        if (all == FALSE) 
	    # do not plot inactive exchange reactions (Rxn whose concentration
            # does not change during the simulation
            if (max(concs[i,]) == min(concs[i,])) 
        	next;

        cat(paste("\t", i, "\n"))
        png(paste(outDir, '/', prefix, '_', i, '.png', sep=""), width=1000, height=1000)
        plot(concs[i,], type='l', col="blue",
	     main="Concentration", 
	     ylab = "mmol", xlab='Time(hrs)')
        legend('topright', legend=i, col="blue", lwd=1)
	dev.off()
    }
}

# Plot active exchange reactions normalized to 1 (this is useful to
# better visualize relationships between metabolite consumption rate
# changes)
plotExchMetsNormalized <- function(dfs, all=FALSE) {
    times <- dfs@timeVec
    concs <- as.matrix(dfs@concentrationMatrix)
    # if we want to plot only non-zero ExcRxns
    if (all == FALSE) {
        concs <- concs[apply(concs[, -1], 1, function(x) any(x != x[1])),]
    }
    
    # normalize the row values
    #concs <- t(apply(concs, 1, function(x)(x-min(x))/(max(x)-min(x))))
    concs <- t(apply(concs, 1, function(x)(x)/(max(x))))


    # choose color palette
    nrxns <- dim(concs)[1]
    cl <- heat.colors(nrxns)
    cl <- terrain.colors(nrxns)
    cl <- topo.colors(nrxns)
    cl <- cm.colors(nrxns)
    cl <- rainbow(nrxns)
    col <- 0    

    # find plot limits
    xmin <- min(sapply(times, function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
    xmax <- max(sapply(times, function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
    ymin <- min(sapply(concs, function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
    ymax <- max(sapply(concs, function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
    # create inital empty plot
    plot(NULL, NULL, xlim=c(xmin, xmax), ylim=c(ymin,ymax), type="l", 
         main="Normalized concentrations", ylab = "mmol",xlab='Time(hrs)')
    # plot the metabolites
    for (i in row.names(concs)) {
        col <- col + 1
        lines(concs[i, ], col=cl[col], type='l', lty=col)
        #points(concs[i, ], col=cl[col], pch=col, cex=0.5)
    }

    legend('topleft', legend=rownames(concs), col=cl, lwd=1, cex=0.5)
}


#	WORKS IN PROGRESS
#	THIS IS STILL IN DEVELOPMENT, WORKS BUT GIVES LITTLE USEFUL INFORMATION
#
plotAllFluxes <- function(dfba) {
    cl <- heat.colors(dim(dfba@all_fluxes)[1])
    cl <- terrain.colors(dim(dfba@all_fluxes)[1])
    cl <- topo.colors(dim(dfba@all_fluxes)[1])
    cl <- cm.colors(dim(dfba@all_fluxes)[1])
    cl <- rainbow(dim(dfba@all_fluxes)[1])
    col <- 0
    # start with empty plot
    plot(NULL, NULL, xlim=c(0,dim(dfba@all_fluxes)[2]), ylim=c(-1000,1000), xlab="time", ylab="flux")
    for (i in 1:dim(dfba@all_fluxes)[1]) {
        col <- col + 1
        lines(dfba@all_fluxes[i,], col=cl[col], 
        	type='l', lwd=2)
    }
    #legend('topleft', legend=rownames(feats$all_fluxes), col=cl, lwd=2)
    col <- 0
    # start with empty plot
    plot(NULL, NULL, xlim=c(0,dim(dfba@all_fluxes)[1]), ylim=c(-1000,1000), xlab="reaction", ylab="flux")
    for (i in 1:dim(dfba@all_fluxes)[2]) {
        col <- col + 1
        lines(dfba@all_fluxes[,i], col=cl[col], 
        	type='l', lwd=2)
    }
    #legend('topleft', legend=rownames(feats$all_fluxes), col=cl, lwd=2)
}

plotFluxes <- function(dfs, all=FALSE) {
	flux <- as.matrix(dfs@all_fluxes)
        if (all == FALSE) {
            flux <- flux[apply(flux[, -1], 1, function(x) any(x != x[1])),]
        }
	# if we want to plot only non-zero fluxes or fluxes between
	# some limits (which is of very little utility)
	#flux <- flux[apply(flux[, -1], 1, function(x) !all(x==0)),]
        #flux <- flux[apply(flux[, -1], 1, function(x) all(x<500)),]
        #flux <- flux[apply(flux[, -1], 1, function(x) all(x>-500)),]
        matplot((flux), type=c("l"))
        matplot(t(flux), type=c("l"))
}

#
#
#	END WORKS IN PROGRESS



#################################
#                               #
#  THE ACTUAL WORK STARTS HERE  #
#                               #
#################################


adynfba <- function() {

    defobj <- 'Biomass_SLI'	# default objective function
    defbiomass <- 'Biomass_SLI'	# default biomass reaction name
    deftsmin <- 60.0		# default time step in minutes
    defst <- 60			# default number of time steps
    definoc <- 0.10		# default inoculum size in g/L (initBiomass)
    defrarxn <- "EX_o2(e)"	# default reaction to use for RA
    defpppax <- "EX_o2(e)"	# default for first reaction in PPPA
    defpppay <- "EX_co2(e)"	# default for second reaction in PPPA
    defexch <- ""		# default time-variation of exchange rates
    defdeltas <- ""		# default substrate deltas per time step
    verbose <- 0

    option_list <- list(
	#-h / --help is added by default

        # options for basic setup
        make_option(c("-m", "--model"), 
            type="character", 
            default="model", 
            help="metabolic model file name [default '%default']"),

        make_option(c("-c", "--const"), 
            type="character", 
            default="constraints", 
            help="constraints file name [default '%default']"),

        make_option(c("-g", "--goal"), 
            type="character", 
            default=defobj, 
            help="goal production (ONLY ONE objective function) 
            [default '%default']"),

        # options for analysis selection
        make_option(c("-f", "--fba"), 
            action="store_true", 
            default=FALSE,
            help="whether to do FBA [default %default]"),

        make_option(c("-o", "--mtf"), 
            action="store_true", 
            default=FALSE,
            help="whether to optimize FBA by MTF [default %default]
            implies FBA"),

        make_option(c("-v", "--fva"), 
            action="store_true", 
            default=FALSE,
            help="whether to do a variability analysis using FVA [default %default]
            implies MTF and FBA"),

        make_option(c("-r", "--ra"), action="store_true", default=FALSE,
            help="whether to do a Robustness Analysis [default %default]"),

        make_option(c("-p", "--pppa"), 
            action="store_true", 
            default=FALSE,
            help="whether to do a Phenotypic Phase Plane Analysis (PPPA)
            [default %default]"),

        make_option(c("-d", "--dfba"), 
            action="store_true", 
            default=FALSE,
            help="whether to do dynamic analysis with DFBA [default %default] 
            implies FVA, MTF and FBA"),

        make_option(c("-a", "--adfba"), 
            action="store_true", 
            default=FALSE,
            help="whether to do an adaptive dynamic analysis with ADFBA [default %default] 
            implies FVA, MTF and FBA"),

        # options to control specific analyses
        make_option(c("-k", "--ctrl-react"), 
            type="character", 
            default=defrarxn, 
            help="control reaction to variate for RA [default '%default'] 
            (only used for RA)"),

        make_option(c("-x", "--pppax"), 
            type="character", 
            default=defpppax, 
            help="first controlled reaction for PPPA [default '%default'] 
            (only used for PPPA)"),

        make_option(c("-y", "--pppay"), 
            type="character", 
            default=defpppay, 
            help="second controlled reaction for PPPA [default '%default'] 
            (only used for PPPA)"),

        make_option(c("-b", "--biomass"), 
            type="character", 
            default=defbiomass, 
            help="name of the biomass reaction [default '%default']
            (only used for DFBA and ADFBA)"),

        make_option(c("-n", "--media"), 
            type="character", 
            default="medium", 
            help="nutrient media file name [default '%default']
            (only used for DFBA and ADFBA)"),

        make_option(c("-t", "--timestep"), 
            type="double", 
            default=deftsmin,
            help=paste("length of the time step in minutes [default %default' (",
            deftsmin/60.0, "h)]
            (only used for DFBA and ADFBA)", sep="")),

        make_option(c("-s", "--steps"), 
            type="integer", 
            default=defst,
            help=paste("maximum number of time steps [default %default (=",
            defst*deftsmin, "' or ", defst*deftsmin/60.0, "h)]
            (only used for DFBA and ADFBA)", sep="")),

        make_option(c("-i", "--inoculum"), 
            type="double", 
            default=definoc,
            help=paste("initial biomass (inoculum size) in g/L [default %default g/L]
            (only used for DFBA and ADFBA)", sep="")),

        make_option(c("-e", "--exch"), 
            type="character", 
            default=defexch,
            help=paste("exchange rate values over time [default '%default']
            (only used for ADFBA)", sep="")),

        make_option(c("-l", "--load"), 
            type="character", 
            default=defdeltas,
            help=paste("substrate concentration deltas over time (mmol/L) 
            [default '%default'] (only used for ADFBA)", sep="")),

	# debugging options (print additional information)
    	#make_option(c("-q", "--quiet"), 
        #    action="store_false", 
        #    default=TRUE,
        #    dest="verbose", 
        #    help="Print little output [default FALSE]"),

        make_option(c("-V", "--verbose"), 
            type="integer", 
            default=0,
            help="Verbosity level from 0 to 5 [default %default]
            (use larger values to increase the amount of information produced)")

    )
    # get command line options, if help option encountered print help and exit,
    # otherwise if options not found on command line then set defaults, 
    opt <- parse_args(OptionParser(option_list=option_list))

    # Now assign values to variables
    modelName <- opt$model
    constraintsName <- opt$const	# additional constraints
    obj <- opt$goal			# obhective function (only one)
    rarxn <- opt$ctrlreact		# control variable for RA
    pppax <- opt$pppax			# x and y control variables for PPPA
    pppay <- opt$pppay
    substrateName <- opt$media		# substrate to use for DFBA
    ts <- opt$timestep / 60		# convert to hours which is the reference time unit
    ns <- opt$steps			# number of steps
    inoculum <- opt$inoculum		# initial Biomass
    biomass <- opt$biomass		# name of biomass function
    exchangeRates <- opt$exch		# "" if not specified
    substrateDeltas <- opt$load		# "" if not specified
    
    if ( opt$fba  ) { doFBA  <- TRUE } else { doFBA  <- FALSE }
    if ( opt$mtf  ) { doMTF  <- TRUE } else { doMTF  <- FALSE }
    if ( opt$fva  ) { doFVA  <- TRUE } else { doFVA  <- FALSE }
    if ( opt$ra   ) { doRA   <- TRUE } else { doRA   <- FALSE }
    if ( opt$pppa ) { doPPPA <- TRUE } else { doPPPA <- FALSE }
    if ( opt$dfba ) { doDFBA <- TRUE } else { doDFBA <- FALSE }
    if ( opt$adfba ) { doADFBA <- TRUE } else { doADFBA <- FALSE }

    if ( opt$verbose ) {
        verbose <- opt$verbose	#5 is the maximum value possible right now
    } else {
        verbose <- 0
    }
    if (verbose > 3) {
        cat("\nadynfba called with option arguments\n")
        print(opt)
    }
    #print(verbose)  
    # some forced calculations
    if (verbose > 1) { 
        if ( doADFBA ) { doFVA <- TRUE ; cat('\nNOTE: ADFBA selected: will also do FVA') }
        if ( doDFBA ) { doFVA <- TRUE ; cat('\nNOTE: DFBA selected: will also do FVA') }
        if ( doFVA  ) { doMTF <- TRUE ; cat('\nNOTE: FVA  selected: will also do MTF') }
        if ( doMTF  ) { doFBA <- TRUE ; cat('\nNOTE: MTF  selected: will also do FBA') }
        cat('\n')
    }
    
    
    # Convert to full path names
    #modelName <- normalizePath(modelName)
    #constraintsName <- normalizePath(constraintsName)
    #substrateName <- normalizePath(substrateName)


    # set up output directory and log file
    # ------------------------------------
    #
    #	Prepare output directory 'out/model_name-etc..../'
    # dir.create will create a directory if it doesn't already exist,
    # printing a warning if it does exist (and showWarnings=TRUE)
    outName <- modelName
    if (constraintsName != "")
        #outName <- paste(outName, "_", constraintsName, sep="")
        outName <- paste(outName, "+", constraintsName, sep="")
    if (substrateName != "")
        outName <- paste(outName, "@", substrateName, sep="")
    if (exchangeRates != "")
        outName <- paste(outName, "=", file_path_sans_ext(basename(exchangeRates)), sep="")
    if (substrateDeltas != "")
        outName <- paste(outName, "+", substrateDeltas, sep="")

    # collapse the biomass vector in case we want to use
    # more than one objective function
    baseout <- paste('obj_', obj, sep="", collapse='_')
    dir.create(baseout, showWarnings = FALSE)
    outDir <- paste(baseout, "/", outName, sep="")
    dir.create(outDir, showWarnings = FALSE)

    if (log == 1) {
        # save log as 'out/mod_name/mod_name+calculations.Rlog'
	ext <- ""
        #if (exchangeRates != "") { ext <- paste(ext, "+", exchangeRates, sep="") }
        if (doFBA == TRUE) { ext <- paste(ext, "+FBA", sep="") }
        if (doMTF == TRUE) { ext <- paste(ext, "+MTF", sep="") }
        if (doFVA == TRUE) { ext <- paste(ext, "+FVA", sep="") }
        if (doRA == TRUE) { ext <- paste(ext, "+RA", sep="") }
        if (doPPPA == TRUE) { ext <- paste(ext, "+PPPA", sep="") }
        if (doDFBA == TRUE) { ext <- paste(ext, "+DFBA", sep="") }
        if (doADFBA == TRUE) { ext <- paste(ext, "+ADFBA", sep="") }
        logName   <- paste(outDir, "/", outName, ext, ".Rlog", sep="");
        #print(logName)
        #print(ext)
	cat('\n\nSaving output as', logName, "\n\n")
        logf <- openLogFile(logName)
    }

    cat("\n============================================================\n");
    cat("\n  ", modelName, "  ", constraintsName, "  ", substrateName, "  ", "\n");
    cat("\n============================================================\n");


    # Load (and modify) Model
    # -----------------------
    model <- loadModel(modelName)

    # this modifies the model, so, it would be wise to save it
    model <- customizeModel(model, obj, ci, constraints=constraintsName)

    # this will save tsv and SBML versions in separate sub-directories
    if ( verbose > 0 ) {
        saveModel(model, overwrite=TRUE, tsv=TRUE, sbml=TRUE, outDir=outDir)
        #saveModel(model, overwrite=TRUE, tsv=FALSE, sbml=FALSE, outDir=outDir)
    }
    # Print out In/Out diagnostic info
    checkModel(model)

    # Prepare the substrate
    # ---------------------
    # only if it will be required by DFBA, otherwise it will not be used
    if ((doDFBA == TRUE) || (doADFBA == TRUE)) {
        if (substrateName == "") {
            # provide a default substrate
            # ### j ###
            # We sometimes use a function setUpSubstrate to set the substrate
            # programatically. This has been removed to reduce complexity.
	    # define setUpSubstrate() function in case we need it (we shouldn't)
	    if (file.exists('setUpSubstrate.R')) {
	        source('setUpSubstrate.R')
                medium <- setUpSubstrate()
                cat("\nSubstrate has been set up with the provided R code.\n\n")
            } else {
	        cat('\nERROR: no substrate specified\n\n')
            }
        } else {
            medium <- readSubstrate(substrateName)
            cat("\nSubstrate has been set up from ", substrateName, ".\n\n")
        }
    }


    # =============
    # CALCULATIONS:
    # =============

    tryCatch( 
    {
        # Standard FBA 
        #
        if (doFBA) {
            cat('\n\n'); banner('FBA')
	    fba_obj_value <- fba(model, outDir=outDir)
        }

        # Minimize total flux
        #
        if (doMTF) {
            cat('\n\n'); banner('MTF')
	    mtf_sol <- mtf(model, outDir=outDir)
        }

        # Flux Variability Analysis (FVA)
        #
        if (doFVA) {
            cat('\n\n'); banner('FVA')
	    fva(model, outDir=outDir, verbose=verbose)
        }

        # Robustness analysis
        #
        if (doRA) {
            cat('\n\n'); banner('RA')
	    reaction <- rarxn
	    ra(model, reaction, outDir=outDir)
        }

        # Phenotypic phase plane analysis
        # We should select a meaningful pair of reactions here (e.g. BIOMASS/AMLB)
        # the two reactions should have predefined, sensible, ul and ll limits
        #
        if (doPPPA) {
            cat('\n\n'); banner('PPPA')
	    pppa(model, c(pppax, pppay), outDir=outDir)
	    #pppa(model, c('EX_mnl(e)', 'Biomass_SLI'))
	    #pppa(model, c('EX_aml(e)', 'Biomass_SLI'))
        }

        # Dynamic Flux Balance Analysis (DFBA)
        #
        if (doDFBA) {
            cat('\n\n'); banner('DFBA')
	    cat("\n--------------------------------------------------------------------------\n")
            #cat(paste("\nDoing DFBA: ", substrate, inoculum, ts, ns, "\n", sep=" "))
            cat("\nDoing DFBA\n")
            cat("\n--------------------------------------------------------------------------\n\n")
	    initBiomass <- inoculum
	    timeStep <- ts;
	    nSteps <- ns;

	    substrateRxns <- medium$substrateRxns
	    initConcentrations <- medium$initConcentrations
	    exclUptakeRxns <- medium$exclUptakeRxns

	    plotRxns <- c('EX_aml(e)', medium$substrateRxns);
	    plotRxns <- c(medium$substrateRxns);

	    dfbaProblem <- list(
		    substrateRxns=substrateRxns, 
		    initConcentrations=initConcentrations,
		    initBiomass=initBiomass, 
		    timeStep=timeStep, 
		    nSteps=nSteps, 
		    plotRxns=plotRxns,
		    exclUptakeRxns=exclUptakeRxns,
		    biomassRxn=biomass,
                    dynamicConstraints=NULL
		    )

	    dfs <- dfba(model, dfbaProblem, outDir=outDir, verbose=verbose)
            # use ADFBA instead
	    #dfs <- adfba(model, dfbaProblem, outDir=outDir, verbose=verbose)


	    #print(dfs)
        }

        # Adaptive Dynamic Flux Balance Analysis (ADFBA)
        #
        if (doADFBA) {
            cat('\n\n'); banner('aDFBA')
	    cat("\n--------------------------------------------------------------------------\n")
            cat("\nDoing Adaptive DFBA\n")
            cat("\n--------------------------------------------------------------------------\n\n")
	    initBiomass <- inoculum;
	    timeStep <- ts;
	    nSteps <- ns;

	    dynamicConstraints <- NULL
	    if (exchangeRates != "") {
		# we could check if it is a .R file and then, if it is
                # source it and set dynamicConstrints to the filename
                # without .R (i.e. to a function defined in the file,
                # which should be named after the control function)
                if ((file_ext(exchangeRates) == 'R') ||
                    (file_ext(exchangeRates) == 'r') ||
                    (file_ext(exchangeRates) == 'Rscript') ||
                    (file_ext(exchangeRates) == 'rscript')) {
                    source(exchangeRates)
                    # set nutrientChanges to file name without '.R'
                    func.name <- file_path_sans_ext(basename(exchangeRates))
                    # get object named func.name, which should be the
                    # desired function.
                    dynamicConstraints <- get(func.name)
                    dcfile <- paste(outDir, '/', modelName, "_dc.tab", sep="")
                    write.table(dynamicConstraints, file=dcfile, sep='\t')
                } else {
	            dynamicConstraints <- getRateChanges(exchangeRates,
		        timeStep, nSteps)
                    dcfile <- paste(outDir, '/', modelName, "_dc.tab", sep="")
                    write.table(dynamicConstraints, file=dcfile, sep='\t')
	        }
            } else {
            	dynamicConstraints <- NULL
            }

            if (substrateDeltas != "") {
                if ((file_ext(substrateDeltas) == 'R') ||
                    (file_ext(substrateDeltas) == 'r') ||
                    (file_ext(substrateDeltas) == 'Rscript') || 
                    (file_ext(substrateDeltas) == 'rscript')) {
                    # assume it is an R file defining the corresponding
                    # control function to use during calculations
    		    source(substrateDeltas)
                    # set nutrientChanges to file name without '.R'
                    func.name <- file_path_sans_ext(basename(substrateDeltas))
                    nutrientChanges <- get(func.name)
                    ncfile <- paste(outDir, '/', modelName, "_nc.tab", sep="")
                    write.table(nutrientChanges, file=ncfile, sep='\t')
                } else {
                    nutrientChanges <- getNutrientChanges(substrateDeltas,
                        timeStep, nSteps)
                    ncfile <- paste(outDir, '/', modelName, "_nc.tab", sep="")
                    write.table(nutrientChanges, file=ncfile, sep='\t')
                }
            } else {
            	nutrientChanges <- NULL
            }
	    substrateRxns <- medium$substrateRxns
	    initConcentrations <- medium$initConcentrations
	    exclUptakeRxns <- medium$exclUptakeRxns

	    # EX_rec(e) should be the recombinant protein
	    #plotRxns <- c('EX_rec(e)', medium$substrateRxns);
	    plotRxns <- c(medium$substrateRxns);
            
            if (doMTF) { method <- "MTF" } else { method <- "FBA" }
            
            method <- "FBA"
            #method <- "lpFBA"
            #method <- "directFBA"
            #method <- "MTF"		# CHECK: when failed it returns incorrectly
            #method <- "directMTF"
            #method <- "lpMTF"

	    adfbaProblem <- list(
		    substrateRxns=substrateRxns, 
		    initConcentrations=initConcentrations,
		    initBiomass=initBiomass, 
		    timeStep=timeStep, 
		    nSteps=nSteps, 
		    plotRxns=plotRxns,
		    exclUptakeRxns=exclUptakeRxns,
		    biomassRxn=biomass,
		    dynamicConstraints=dynamicConstraints,
                    nutrientChanges=nutrientChanges,
                    method=method
		    )
	    dfs <- adfba(model, adfbaProblem, outDir=outDir, verbose=verbose)

	    #print(dfs)
        }

        # for quick interactive diagnosis
        #
        cat('\nModel: ', model@mod_name, '\n')
        cat('\nObjective function:', obj, '\n')
        if (doFBA) {
            print( fba_obj_value )
        }

        cat('\n---------------------------------------------------')
        cat('\n  ', model@mod_name, ' done')
        cat('\n---------------------------------------------------\n')

    }, 
    warning=function(cond) { 
        message('WARNING:'); 
        message(cond) ; cat('\n\n') 
    },
    error=function(cond) { 
       traceback();
        message('ERROR:'); 
        message(cond) ; 
        cat('\n\n') },
    finally={
        message('INFO: closing files'); 
        while (sink.number() > 0 ) sink()
    }  
    )  # end TryCatch

    if (log == 1) {
        closeLogFile(logf)
    }

}

# for debugging
options(error=function() {traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

adynfba()

cat('\n\n'); banner('  bye')

# done: quit
q(save="no")

#
# NOTE: to debug typos use R and issue the command
# eval(parse('filename_to_check'))
# e.g.
# R -e 'eval(parse("adynDFBA.R"))'
# before calling a function
# options(error=recover)
# options(warn=2)
# inside a function
# browser()
