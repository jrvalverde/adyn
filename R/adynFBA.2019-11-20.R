#!/usr/bin/env Rscript

#
#	This is an R script file
#
#	You can use it by making it executable and using as a program or
# by calling it with Rscript or using 'source()' from within R
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
# the default objective is EX_BIOMASS (exchange of Biomass)
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
log <- 1


#----------------------------------------------------------------------
#	FUNCTION DEFINITIONS
#----------------------------------------------------------------------

#' myname()
#'
#' myname() obtains the name of the calling function
#'
#' This function will return the name of the calling function.
#' The CALLING function !!!
#' THE CALLING FUNCTION !!!!!!
#' It is intended as a way to obtain a function's own name for error
#' reporting. 
#  We could do it directly, within a function using 
#  myname <- deparse(sys.call())
#  but it would be difficult to understand. Isolating this in a function
#  call makes it more readable.
#'
#' 
#' @return	the name of the calling function (one up in the call stack)
#'
#' @usage	me <- myname()
#'
#' @examples
#'	f <- function() { print(myname()) }
#'	f()
#'	##[1] "print(myname())"
#'	## myname() is being passed to print() as an argument and is thus
#'	## called by 'print', hence the output
#'	##
#'	f <- function() { me <- myname() ; print(me) }
#'	f()
#'	##[1] "f()"
#'	## myname() is called by f() to obtain its own name and then the name
#'	## is handled to print().
#'
#' @author	(C) José R. Valverde, CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
myname <- function() { 
    deparse(sys.calls()[[sys.nframe()-1]]) 
}


#' cat.err()
#'
#' err() prints a generic error message using cat()
#'
#' @param	abort	(boolean) whether the program should be stopped
#'			(defaults to FALSE)
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.err('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.err(x, "is too big\n")}}
#'		f(13)
#'		##ERROR in  f(13) : 13 is too big
#'
#' @author	(C) José R. Valverde, CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
cat.err <- function(abort=FALSE, ...) { 
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    caller <- deparse( sys.call(-1) )
    cat('ERROR in ', caller, ":", ...)
    if (abort) {
        quit(save="no", status=1, runLast=FALSE)
    }
}


#' cat.warn()
#'
#' cat.warn() prints a generic warning message using cat()
#'
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.warn('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.warn(x, "is too big\n")}}
#'		f(13)
#'		##WARNING in  f(13) : 13 is too big
#'
#' @author	(C) José R. Valverde, CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
cat.warn <- function(...) { 
     #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
     #me <- as.list(sys.call())[[1]]
     #parent <- as.list(sys.call(-1))[[1]]
     caller <- deparse( sys.call(-1) )
     cat('WARNING in ', caller, ":", ...)
}


#' cat.info()
#'
#' cat.info() prints a generic information message using cat()
#'
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.warn('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.info(x, "is too big\n")}}
#'		f(13)
#'		##INFO  f(13) : 13 is too big
#'
#' @author	(C) José R. Valverde, CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
cat.info <- function(...) { 
     #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
     #me <- as.list(sys.call())[[1]]
     #parent <- as.list(sys.call(-1))[[1]]
     caller <- deparse( sys.call(-1) )
     cat('INFO ', caller, ":", ...)
}


# to be used instead of library(): this function ensures
# that the package is installed if not present in the system
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dependencies = TRUE)
    require(p, character.only = TRUE)
}


#' sourceDir
#'
#' sourceDir() sources all Q/R/S files in a directory
#'
#'	This may be convenient to source more than one file at once
#'	i.e. when getting dynamic constraints as a function,
#'	although it would likely make more sense if the file defining
#'	the function did include (source) itself all other required
#'	files explicitly.
#'	It may also be helpful to load all library files from a directory.
#'
#' @param	path	the path to the directory containing the files to source
#' @param	trace	Whether we want tracing output printed (default=TRUE)
#' @param	\dots	Any other parameters to pass on to source()
#'
#' @return	On end, any file ending in .(RrSsQq) present in the
#'		directory will have been sourced and the corresponding
#'		operations applied
#'
#' @usage	sourceDir('some/source/directory/')
#'
#' @examples
#'		sourceDir('./library/.')
#' 
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
sourceDir <- function(path, trace = TRUE, ...) {
   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if (trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if (trace) cat("\n")
   }
}


#' openLogFile
#' 
#' openLogFile() opens a file to save all output, deleting it first if it exists
#'
#' @param	logName	the name of the log file
#'
#' @return	a reference to the log file that can be used to close it later
#'
#' @usage	log <- openLogFile('output.log')
#'
#' @examples
#'		log <- openLogFile('output.log')
#'		print(log)
#'		closeLogFile(log)
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#'
openLogFile <- function(logName) {

    if (file.exists(logName)) file.remove(logName)

    logFile <- file(logName, open="wt");
    sink(logFile, type=c("output", "message"), split=TRUE);
    return (logFile)
}


#' closeLogFile
#'
#' closeLogFile() - close a log file previously opened by openLogFile
#' 
#' @param	logFile	a handle to the log file returned by a previous call
#'			to openLogFile
#' 
#' @return	on end, the log file is guaranteed to be closed
#'
#' @usage	closeLogFile(log)
#' 
#' @examples
#'		log <- openLogFile('output.log')
#'		print(log)
#'		closeLogFile(log)
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
closeLogFile <- function(logFile)  {
    close(logFile)
    while (sink.number() > 0) { sink() }
}


#' file.extension
#'
#' file.extension() returns the file extension suffix part of a file name
#' 
#' @param	path	the file name path
#' 
#' @return	the extension of the file (i.e. anything following after the
#'		last '.'
#'
#' @usage	ext <- file.extension('file/path.ext'
#' 
#' @examples
#'		ext <- file.extension('/some/path/to/a/file.ext')
#'		print(ext)
#'		## [1] "ext"
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
file.extension <- function (path) 
{
    parts <- strsplit(path, ".", fixed=T)[[1]]
    last <- parts[length(parts)]
    return(last)
}


#' loadModel
# ==========
#'
#' loadModel() loads a model in any supported format
#'
#' 	First try SBML and if that fails, try TSV
#
# NOTE: we should also check if the full filename has been given by testing
# its extension(s)
#' 
#' @param	modelName	the name of the file containing the model.
#'				It may be given with or without the extension.
#'				The reason is that for TSV format models, we
#'				need to load a number of files all starting with
#'				the same prefix. Thus, we will try first with
#'				modelName.xml, then with 'modelName' as the
#'				common root of a TSV model, and finally, use
#'				'modelName' as the full file name assuming it
#'				must necessarily be in XML/SBML format
#' 
#' @return	on end, loadModel() will return the read model structure
#'
#' @usage	mod <- readModel('model.sbml')
#' 
#' @examples
#'		## Assuming a file 'myModel.xml' exists in this directory
#'		mod <- readModel('myModel')
#'		## Assuming a file 'myModel_react.tsv' (and all other TSV
#'		## files defining a model exists
#'		mod <- readModel('myModel')
#'		## Assuming a file 'myModel.sbml' exists
#'		mod <- readModel('myModel.sbml')
#'
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
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
 	    #cat.err('No suitable model (.xml or .tsv) found', abort=TRUE);
    }
    return( model )
}


#' saveModel
#'
#' saveModel() is a commodity routine to ease saving an in-memory model
#'
#'	saveModel() can save the model as either TSV, SBML or both, and
#' overwrite (or refuse to) any previous file(s) with the same name if
#' they exist. Since TSV models will be saved as several files with different
#' terminations, only the root name (without extension) of the file needs to
#' be provided. By default, files will be created in the current directory,
#' but an alternate absolute or relative path to a directory may be provided
#' to save the files instead.
#'	NOTE that models will be saved in specific subdirectories named after
#' the format. This is to ensure output is duly comparmentalized.
#'	We normally use this function to save modified models associated to
#' any specific simulation. The SBML file is portable, but the TSV is easier
#' to search and modify (e.g. in an spreadsheet like LibreOffice Calc).
#'
#' 
#' @param	name	the root file name for the saved model: if the model
#'			will be saved in SBML format, it will be saved to
#'			a sub-folder 's_sbml' and named as 'name' plus the
#'			SBML level version and '.xml' (i.e. 
#'			s_sbml/"name"_L2v1.xml). If a TSV format model
#'			is to be produced, it will be saved in a sub-folder
#'			named 's_tsv', as three files named 
#'			's_tsv/"name"_desc.tsv', 's_tsv/"name"_met.tsv" and
#'			's_tsv/"name"_react.tsv'.
#' @param	tsv	boolean flag indicating whether a TSV format model
#'			should be output (defult=TRUE)
#' @param	sbml	boolean flag indicating whether a SBML format model
#'			should be output (default=TRUE)
#' @param	overwrite	boolean flag indicating whether to overwrite
#'			existing files with the same name if present (default
#'			=TRUE)
#' @param	outDir	the path to an output directory (default='.', the 
#'			current working directory
#' 
#' @return	on exit, the model(s) should have been saved to corresponding
#'		subdirectory(ies) in the output directory outDir.
#'
#' @usage	saveModel(mod, tsv=TRUE, sbml=TRUE, overwrite=TRUE, outDir='.')
#' 
#' @examples
#'		## read a model from 'iJV1220.xml' and save it again as
#'		## both as SBML and TSV
#'		mod <- readModel('iJV1220')
#'		## make some modification(s)
#'		mod <- changeObjFunc(mod, 'BIOMASS_SLI', 1)
#'		saveModel('iJV1220-obj-Biomass', tsv=T, sbml=T, overwrite=F)
#'		## should create two subdirectories with
#'		##	./s_sbml/iJV1220-obj-Biomass_L2v1.xml
#'		##	./s_tsv/iJV1220-obj-Biomass_desc.tsv
#'		##	./s_tsv/iJV1220-obj-Biomass_met.tsv
#'		##	./s_tsv/iJV1220-obj-Biomass_react.tsv
#'		
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
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

#' loadTableFile
#'
#' loadTableFile() loads a file containing a table in any supported format
#'
#'	This is a general routine. It will first detect which formats are
#' supported in your R installation by trying to load various foreign format
#' access packages. If any of these packages is found, then all corresponding
#' file extensions are considered.
#'	Once the supported formats (and extensions) have been detected, the
#' name provided will be successively combined with each of the extensions and
#' if a file named 'file'+'.'+ ext is found, it will be read using the 
#' appropriate method. If no combination of 'file' + '.' + any of the supported
#' file extensions is found, then we will try to read the file using its
#' name as provided (considering its own extension if any). If no extension
#' match is found, then we will simply try to load it with 'read.table()'.
#'	A special provision is made for reading .R or .Rscript files: if the
#' successfully found file extension corresponds to an R script, then it will
#' be sourced, and the file will be searched for either, a table that has
#' the same name as the file (without extension), which will then be considered
#' the table to read, or a function with the same name as the file (without
#' extension), which will then be considered a setup function to create the
#' desired table, on demand. Note that this function, in turn, might as well 
#' read an external file to load the table internally, that would be hidden 
#' from us and would not make any difference.
#'	
#' XXX JR XXX NOTE: I need to run more exhaustive tests on this function
#
# This one could be simplified by ignoring the extension and smply
# trying all known formats in order
#	That'd be safer but "dirtier"
#' 
#' @param	file	the name (with ot without extension) of the file
#'			that contains the table we want to load
#'		stringsAsFactors	flag indicating whether strings should
#'			be read as factors or not (default FALSE)
#' @param	header	boolean flag indicating if the table contains a header
#'			(defaults to FALSE)
#' @param	check.names boolean flag to tell if table names should be 
#'			checked for R conformance (defaults to FALSE)
#'
#' 
#' @return	the table read
#'
#' @usage	table <- loadTableFile(file="myTable")
#' 
#' @examples
#'		## this should be enough
#'		table <- loadTableFile("myTable")	
#'		## but any of these should also work if you want to be specific
#'		table <- loadTableFile("myTable.tsv")
#'		table <- loadTableFile("myTable.csv")
#'		table <- loadTableFile("myTable.dat")
#'		table <- loadTableFile("myTable.txt")
#'		table <- loadTableFile("myTable.xlsx")
#'		table <- loadTableFile("myTable.arff")
#'		table <- loadTableFile("myTable.dbf")
#'		table <- loadTableFile("myTable.dta")
#'		table <- loadTableFile("myTable.rec")
#'		table <- loadTableFile("myTable.mtp")
#'		table <- loadTableFile("myTable.mat")
#'		table <- loadTableFile("myTable.sav")
#'		table <- loadTableFile("myTable.ssd")
#'		table <- loadTableFile("myTable.sas")
#'		table <- loadTableFile("myTable.sys")
#'		table <- loadTableFile("myTable.syd")
#'		table <- loadTableFile("myTable.xpt")
#'		table <- loadTableFile("myTable.dump")
#'		table <- loadTableFile("myTable.R")
#'		table <- loadTableFile("myTable.Rscript")
#'		table <- loadTableFile("myTable.Rdata")
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
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
    suff <- c('tsv', 'csv', 'tab', 'dat', 'txt', 'R', 'Rscript', 'Rdata')
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

    # if none is available, then ask for a file name
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
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'tab') {
            table <- read.table(fnames[i], sep="\t", 
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'csv') {
            table <- read.table(fnames[i], sep=",", 
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'txt') {
            table <- read.table(fnames[i], 
                                header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
            if (! empty(table)) { return(table) }
        } else if (ext == 'dat') {
            table <- read.table(fnames[i], 
            			header=header, check.names=check.names,
            			fill=T, strip.white=T, blank.lines.skip=T)
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
            } else if (ext == 'dump') {
                #table <- data.restore(fnames[i]) 
                table <- read.S(fnames[i], 
                                header=header, check.names=check.names)
                if (! empty(table)) { return(table) }
            }
        } else if ((ext == 'R') || (ext == 'Rscript') ||
                   (ext == 'r') || (ext == 'rscript')) {
            # it is an R file that should contain a function named like the file
            source(fnames[i])
            # we could have used return(dget(data'.R') but then the file should
            # contain only an anonymous function, this is more versatile.
            # Check if an object named as data exists and if it is a function
            data <- file_path_sans_ext(basename(fname[i]))
            
            # if the content was a function
            if (exists(data) && is.function(get(data))) {
                # in principle, this is the same as return(get(data))
                # this way, we have together both ways of doing it
                return(eval(parse(text=data)))
            } 
            # if the content was an array or a data.frame
            else if (exists(data) && 
                       (is.array(get(data)) || is.data.frame(get(data)))) {
                # in case that instead of a function we got a data definition
                # (matrix is a kind of array)
                table <- data
            }else {
                cat('loadTableFile: ERROR - ', data, '.R{script} should contain R code\n',
                    '    as a function ', data, '(model, concentrations, fluxes, step)\n',
                    '    or a "', data, '" array or data.frame', sep='')
                return(NULL)
            }
            ### JR ### NEEDS TO BE TESTED
        } else if (((ext == "Rdata") || (ext == "rdata")) && exists(data)) {
            return(load(data))
        } else {
            # At this point we have tried all known suffixes. As a last resort, 
            # we will try to read it as a text file
            table <- read.table(fnames[i], header=header, check.names=check.names)
	    # and leave it as default if nothing else works
            if (! empty(table)) { return(table) }
        }
    }
    # if we reach here it is because no approach has worked
    cat("\n\nloadTableFile: File", data, "was not valid!\n\n\n")

    return(table)    
}


#' loadConstraints
#'
#' loadConstraints() opens and load a file with a constraints table
#'
#	The file name needs not carry an extension. First, the name will
#' be concatenated with a number of extensions and each will be tried (.tab, 
#' .tsv, .txt, .dat, .csv, .xlsx, .Rdata, .R and .Rscript). If none of the
#' composite names works, then the file name provided will be tried as a 
#' tab-separated file. If that also fails, then we will try to let the user
#' choose a file.
#'	If the file is an R script, then it will be loaded and we'll check
#' if a table or data.frame with the same name as the file (without extension)
#' has been created or, alternatively, if a function has been defined, which
#' will be called to build the desired table.
#' 
#' @param	constraints	the name of the file with the constraints table
#' 
#' @return	a table with the loaded constraints
#'
#' @usage	constr <- loadConstraints('filename')
#' 
#' @examples	
#'		constr <- loadConstraints('exponential_growth')
#'		print(constr)
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
loadConstraints <- function(constraints="")
{
    if (constraints == "") {
        return (NULL)
    }
    # Read mods. file (if any) and apply it
    options(stringsAsFactors=FALSE)		# we want strings as strings
    # by default, # will be used as comment.char
    if (file.exists(paste(constraints, ".tsv", sep=""))) {
    	limits <- read.table(paste(constraints, ".tsv", sep=""), sep="\t", header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, '.dat', sep=""))) {
    	limits <- read.table( paste(constraints, '.dat', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, ".tab", sep=""))) {
     	limits <- read.table(paste(constraints, ".tab", sep=""), sep="\t", header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, ".csv", sep=""))) {
    	limits <- read.table(paste(constraints, ".csv", sep=""), sep=",", header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, ".txt", sep=""))) {
    	limits <- read.table(paste(constraints, ".txt", sep=""), header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if ((is.element('xlsx', installed.packages()[,1])) &&
             (file.exists(paste(constraints, ".xlsx", sep="")))) {
        limits <- read.xlsx(paste(constraints, ".xlsx", sep=""), sheet.index=1)
    } 
    else if (file.exists(paste(constraints, '.Rdata', sep=""))) {
        dat <- load (paste(constraints, '.Rdata', sep=""))
    } 
    else if (file.exists(paste(constraints, '.R', sep=""))) {
        # it is an R file that should contain a function or array named like the file
        source(paste(constraints, '.R', sep=""))
        # we could have used return(dget(data'.R') but then the file should
        # contain only an anonymous function; this is more versatile.
        # Check if an object named as data exists and if it is a function
        if (exists(constraints) && is.function(get(constraints))) {
            return(eval(parse(text=constraints)))
        } else if (exists(constraints) && 
                   (is.array(get(constraints)) || is.data.frame(get(constraints)))) {
            # in case that instead of a function we got a data definition
            # (matrix is a kind of array)
            limits <- constraints
        } else {
            cat('loadConstraints: ERROR - ', constraints, '.R should contain\n',
                '    a function ', constraints, '(model, concentrations, fluxes, step)', sep='')
            return(NULL)
        }
        ### JR ### NEEDS TO BE TESTED
    } 
    # last resort: try to read it directly (without extension) as a TAB file
    else if (file.exists(constraints)) {
        limits <- read.table(constraints, sep="\t", header=FALSE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else {
        cat("\n\nloadConstraints: File ", constraints, "{.txt, .tsv, .tab, .csv} not found\n")
    	limits <- read.table(file.choose(), header=FALSE)
    }
    return(limits)
}

#' applyConstraints
#'
#' applyConstraints() returns a new model based on the one provided with the
#'		modifications specified
#'
#'	Given a model, applyConstraints() can modify the objective function(s)
#' using a vector of reaction names and a vector of coefficients, and can modify
#' any reaction limits according to a table of reaction names followed by the
#' lower and upper limit values.
#'	Reaction limits are specified in a table with three columns: reaction
#' name, lower limit and upper limit. Lines starting with a hash mark ('#')
#' are considered comments and ignored. Changes are applied successively, i.e.
#' if a reaction appears more than once, the last changes will prevail. e.g.:
#'
#'	##Fluxes set in units of mmol/h*gDW 
#'	##               Lower Bound     Upper Bound
#'	##EX_glc(e)      -2.49397        0.00001
#'	EX_glc(e)       -1.40000        0.00001
#'	EX_mnl(e)       0.000000        0.00001
#'	...
#'
#'	Any reaction (including the Biomass reaction) can be modified.
#' 
#' @param	model	the base model to be modified
#' @param	obj	a vector containing the names of the objective functions
#' @param	ci	a vector containing the coeeficients for each objective
#' @param	limits	a table with three columns specifying the reactions to
#'			modify and the new lower and upper limits; '#' indicate
#'			comments. The columns need not be named but must be
#'			in order 'reaction-name', 'lower-limit', 'upper-limit'
#' @param	verboe	the level of verbosity desired (the amount of output to
#'			generate while doing the work)
#' 
#' @return	a new model with the changes applied
#' 
#' @usage	mod <- applyConstraints(mod, obj=obj_vector, ci=ci_vector, limits=limits_table)
#' 
#' @examples
#'		## Load S.lividans model
#'		data(Slividans)
#'		## change the objective function
#'		Slividans <- applyConstraints(Slividans, 
#'			c('Biomass_SLI', 'EX_aml(e)'),
#'			c(      1      ,     1))
#'		## change the reaction limits
#'		constr <- matrix(c('Ex_glc(e)', -2.493, 0.0001,
#'		                   'Ex_mnl(e)', -1.4, 0.0001),
#'		                   nrow=2, ncol=3)
#'		Slividans <- customizeModel(Slividans, limits=constr)
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
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

    cat("Changing model name from", model@mod_name, "to", 
    	paste(model@mod_name, "_", constraints, sep=""), "\n\n")
    model@mod_name <- paste(model@mod_name, "_", constraints, sep="")

    return (model)
}

#' customizeModel
#'
#' customizeModel() applies the constraints specified and returns a new model
#
# Unused, redefined below. Work in progress (pending verification of correctness).
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
NEWcustomizeModel <- function(model, obj, ci, constraints="", verbose=0)
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


#' customizeModel
#'
#' customizeModel(): return a new model including our customizations
#'
# Limit uptake rates
# ------------------
# Reactions with upper and lower bounds set to zero are not functional.
# Exchange reactions with lower bound = 0 and upper bound > 0 only
#	allow secretion.
# To enable metabolite uptake, set lower bound < 0
#' 
#' @param	model	The base model to customize
#' @param	obj	A vector with the name of the objective reactions
#' @param	ci	A vector with the corresponding coefficients
#' @param	constraints	A file name that contains a matrix with three 
#'			columns corresponding to the reaction to modify and 
#'			its new lower and upper limits, e.g.:
#'		##Fluxes set in units of mmol/h*gDW 
#'		##               Lower Bound     Upper Bound
#'		##EX_glc(e)      -2.49397        0.00001
#'		EX_glc(e)       -1.40000        0.00001
#'		EX_mnl(e)       0.000000        0.00001
#'		...
#'
#' @return	a modified model with the requested changes
#'
#' @usage	mod <- customizeModel(mod, obj, ci, constraints)
#' 
#' @examples
#'		## Load S.lividans model
#'		data(Slividans)
#'		## change the objective function
#'		Slividans <- applyConstraints(Slividans, 
#'			c('Biomass_SLI', 'EX_aml(e)'),
#'			c(      1      ,     1))
#'		## change the reaction limits
#'		constr <- matrix(c('Ex_glc(e)', -2.493, 0.0001,
#'		                   'Ex_mnl(e)', -1.400, 0.0001),
#'		                 nrow=2, ncol=3)
#'		Slividans <- customizeModel(Slividans, limits=constr)
#'
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
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
    if (file.exists(paste(constraints, ".tsv", sep=""))) {
    	limits <- read.table(paste(constraints, ".tsv", sep=""), sep="\t", header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, '.dat', sep=""))) {
    	limits <- read.table( paste(constraints, '.dat', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, ".tab", sep=""))) {
     	limits <- read.table(paste(constraints, ".tab", sep=""), sep="\t", header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, ".csv", sep=""))) {
    	limits <- read.table(paste(constraints, ".csv", sep=""), sep=",", header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(constraints, ".txt", sep=""))) {
    	limits <- read.table(paste(constraints, ".txt", sep=""), header=FALSE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if ((is.element('xlsx', installed.packages()[,1])) &&
             (file.exists(paste(constraints, ".xlsx", sep="")))) {
        limits <- read.xlsx(paste(constraints, ".xlsx", sep=""), sheet.index=1)
    } 
    else if (file.exists(paste(constraints, '.Rdata', sep=""))) {
        dat <- load (paste(constraints, '.Rdata', sep=""))
    } 
    else if (file.exists(paste(constraints, '.R', sep=""))) {
        # it is an R file that should contain a function or array named like the file
        source(paste(constraints, '.R', sep=""))
        # we could have used return(dget(data'.R') but then the file should
        # contain only an anonymous function; this is more versatile.
        # Check if an object named as data exists and if it is a function
        if (exists(constraints) && is.function(get(constraints))) {
            return(eval(parse(text=constraints)))
        } else if (exists(constraints) && 
                   (is.array(get(constraints)) || is.data.frame(get(constraints)))) {
            # in case that instead of a function we got a data definition
            # (matrix is a kind of array)
            limits <- constraints
        } else {
            cat('customizeModel: ERROR - ', constraints, '.R should contain\n',
                '    a function ', constraints, '(model, concentrations, fluxes, step)', sep='')
            return(NULL)
        }
        ### JR ### NEEDS TO BE TESTED
    } 
    # last resort: try to read it directly (without extension) as a TAB file
    else if (file.exists(constraints)) {
        limits <- read.table(constraints, sep="\t", header=FALSE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else {
        cat("\n\ncustomizModel: File ", constraints, "{.txt, .tsv, .tab, .csv} not found\n")
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


#' checkModel
#'
#' checkModel() inspects the model and reports exchange reactions
#' 
#' @param	model	The model to check
#' 
#' @return	this function prints a report of the exchange and uptake
#'		reactions on standard output
#'
#' @usage	checkModel(mod)
#' 
#' @examples
#'		## Load S.lividans model
#'		data(Slividans)
#'		checkModel(Slividans)
#'		## change some reaction limits
#'		constr <- matrix(c('Ex_glc(e)', -2.493, 0.0001,
#'		                   'Ex_mnl(e)', -1.4, 0.0001),
#'		                   nrow=2, ncol=3)
#'		Slividans <- customizeModel(Slividans, limits=constr)
#'		checkModel(Slividans)
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
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


#
# ----------------- Various sigmoid functions -------------------------
#

# Here are a bunch of sigmoid functions that may (or may not) be more
# or less comfortable to use as lueprints for printing/display depending
# on the parameters available.

# Sigmoid fitting is a very general approach that may not be the
# most sensible. It may be used (with due caveats) to approximate
# MM or growth curves, but it is far from being the ideal approach.
#
# We have also tried other fits, like splines, etc. (see below) to generate
# missing intermediate Biomass data. 
#
# Try to choose the most adequate by yourself based on scientific grounds
# or, if you are not sure, select the one with less assumptions (i.e. linear
# fitting) as the less 'inaccurate' (Science 101). This is specially relevant
# when you have scarce growth data, as the missing data may hide relevant 
# patterns that would be obscured by the assumptions of a given model. 
# On the other hand, when you have enough (as in many) data points and no
# parametric model gets tightly close to your data, you may want to use
# a spline to reproduce better the observed behaviour. But, generally, you
# should FIRST study your data, decide on the best parametric model and
# use a THEORETICALLY SOUND PARAMETRIC MODEL with preference.
#
# You may be interested in reading the (highly recommended) paper
# grofit: Fitting Biological Growth Curves in R (2010) Kahn, M., et al.
# J. Stat. Soft. vol. 33, issue 7.
#
# Another interesting reference is
# https://en.wikipedia.org/wiki/List_of_probability_distributions
#
# Although we try to avoid external packages in the code for simplicity, we
# should certainly use them in later versions. Once we publish this, and if, 
# and when, we can gather any support for doing it.
#

# generic sigmoid, logistic function or Richard's curve
#
#                 K - A
# Y(t) = A + ------------------
#                          1/v
#             /        -Bt\  
#             \C + Q e    /
# where
#	Y = growth, weight, size...
#	A = lower asymptote
#	K = upper asymptote when C = 1. If A=0 and C=1 then K is called
#	    the carrying capacity
#	B = growth rate
#	v > 0 affects near which asymptote maximal growth occurs
#	Q = related to Y(0)
#	C = typically 1. If not, then the upper assymptote instead of K
#	    will be A + ((K -A) / (C^(1/v)))
#
# alternatively, it may be written as
#
#                  K -A
# Y(t) = A + ---------------------
#                              1/v
#             /        -B(t-M)\
#             \ C + Q e       /
# where 
#	M = starting time, t0, at which Y(t0) = A + ((K-A) / ((C+1)^(1/v)))
#
# or as
#                   K - A
# Y(t) = A + ---------------------
#                              1/v
#             /        -B(t-M)\
#             \ C + Q e       /
#
# where it is easier to set the starting time and the value of Y at t0
#
#
#                                    -a
#                           /     -x\
# It can be standardised as \1 + e  /   , a > 0
#
#
#
# The logistic function, with maximum growth rate at time M, is the case
# of the general equation where Q = v = 1
#
#               L
# f(x) = -----------------
#              -k(x - x0)
#         1 + e
# where
#	x0 = x value of the sigmoid's midpoint
#	L  = the maximum value of the curve
#	K  = logistic growth rate or steepness of the curve
#
# This leads to the standard logistic function (k=1, x0=0 and L=1)
#                          x
#            1           e       1      1         x
# f(x) = --------- = -------- = --- +  --- tanh( --- )
#              -x      x         2      2         2
#         1 + e       e  + 1
#
# here, 1 - f(x) = f(-x), and corresponds to an ODE f(x)(1 - f(x)), the
# continuous version of the logistic map.
#
#
#
# Gompertz's model for growth curves is a sigmoid following the formula
#
# N(t) = b + N0 exp(-c (exp(a·t) - 1))
#
# where
#	N0 = initial number of cells/organisms
#	a  = asymptote
#	b  = positive number, displacement across the X axis
#	c  = positive number, growth rate
#
#
# Another interesting function is the Von Bertalanffy growth function (VGBF)
#
# L(a) = L∞ ( 1 - exp(-k(a-a0)) )
#
# where
#	a  = time
#	k  = growth coefficient
#	a0 = the value used to calculate size when time is zero
#	L∞ = asymptotic size
# try (e.g.) with L∞ = 275.2, K = 0.52, a0 = -0.47 (values for 
# size change of Girella nigricans with age)
# or
# y = sigmoid.vbgf(x, Linf=45.78855, K=0.23959, t0=-0.08253) ; plot(x,y)
#
sigmoid.vbgf <- function(x, K, t0, Linf) {
    y = Linf * (1 - exp(-K * (x - t0)))
    return(y)
}


# five parameter logistic function
# See Liao and Liu (2009) Re-parameterization of the five-parameter
# logistic function. J. of Chemometrics
# https://onlinelibrary.wiley.com/doi/abs/10.1002/cem.1218
#
# Starting from the four-parameter function:
#
#               A - D
# y = D + --------------------
#                    B
#              /  x \
#          1 + | --- |
#              \  C /
# 
# its advantage is that the parameters have useful interpretations: A and D
# are the two asymptotes, B is the shape parameter or slope, and C is the 
# ED50 value with is the concentration/dose/x corresponding to a response
# hafway between the asymptotes.
#
# We can move to a 5PL (five-parameter logistic function) 
#
#               A - D
# y = D + --------------------
#                           g
#          +-           B -+
#          |     /  x \    |
#          | 1 + | --- |   |
#          |     \  C /    |
#          +-             -+
#
# which introduces a fifth parameter 'g' for assymmetry in the function.
# The EDC50 value from this 5PL is C * (2^(1/g) - 1)^(1/B)
# Since the EDC50 is often a parameter of interest, a new re-parameterized
# function has been proposed:
#
#                     A - D
# y = D + -------------------------------
#                                      g
#         +-                      B -+ 
#         |       1/g       /  x \   |
#         | 1 + (2    -1) * | --- |  |
#         |                 \  C*/   |
#         +-                        -+
#
# where C* is ED50, as in the 4PL, and hence C* = C(2^(1/g) - 1)^(1/B) and
# C is the parameter of the standard 5PL.
#
# The advantage here is that when parameterizing the curve, it is possible
# to obtain directly C* (ED50) and its variability, instead of indirectly
# with a convolute operation. See the paper by Liao and Liu for details.
#

#' sigmoid.4pl()
#' 
#' sigmoid.4pl() computes the Y values for a four parameter logistic (4PL) sigmoid
#'
#
#               A - D
# y = D + --------------------
#                    B
#              /  x \
#          1 + | --- |
#              \  C /
# 
#' its advantage is that the parameters have useful interpretations: A and D
#' are the two asymptotes, B is the shape parameter or slope, and C is the 
#' ED50 value with is the concentration/dose/x corresponding to a response
#' hafway between the asymptotes.
#'
#' @param	x	The value(s) whose Y must be computed
#' @param	A	Lower (left) asymptote
#' @param	B	Shape parameter
#' @param	C	ED50 value, x value corresponding to a response alfway
#'			between the asymptotes
#' @param	D	Upper (right) asymptote
#'
#' @return	The Y value(s) corresponding to the specified sigmoid for X
#'
#' @usage	sigmoid.4pl(x, A, B, C, D)
#'
#' @examples
#'	## e.g. A=0, D=4, C=1.6, B=1.5 (g=1)
#'	x=seq(0, 5, 0.1)
#'	y = sigmoid.4pl(x, A=0, D=4, C=1.6, B=1.5); plot(log(x), y)
#'
#' @author	José R. Valverde, CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
sigmoid.4pl <- function(x, A, B, C, D) {
    y = D + ( (A - D) / (1 + ((x / C)^B)) )
    return(y)
}


#' sigmoid.5pl()
#'
#' sigmoid.5pl() computes the five-parameter logistic sigmoid
#'
#
#               A - D
# y = D + --------------------
#                           g
#          +-           B -+
#          |     /  x \    |
#          | 1 + | --- |   |
#          |     \  C /    |
#          +-             -+
#
#' As in the 4 parameter logistic formula for a sigmoid, here, A and D
#' are the two asymptotes, B is the shape parameter or slope, C has become
#' a new coefficient dependent on ED50, and g is a fifth parameter that 
#' introduces an asymmetry on the function shape.
#'
#' The ED50 value in this version is now  C * (2^(1/g) - 1)^(1/B) which makes
#' it more difficult to estimate from a curve fit as it no longer enters
#' directly the formula.
#'
#' @param	x	The value(s) whose Y must be computed
#' @param	A	Lower (left) asymptote
#' @param	B	Shape parameter
#' @param	C	ED50 value, x value corresponding to a response alfway
#'			between the asymptotes
#' @param	D	Upper (right) asymptote
#' @param	g	Asymmetry factor
#'
#' @return	The Y value(s) corresponding to the specified sigmoid for X
#'
#' @usage	sigmoid.5pl(x, A, B, C, D, g)
#'
#' @examples
#'	## e.g. A=0, B=2, C=2.942, D=4, g=3
#'	x=seq(0, 5, 0.1)
#'	y = sigmoid.5pl(x, A=0, B=2, C=2.942, D=4, g=3) ; plot(log(x), y)
#'
#'	## e.g. A=0, D=4, C=1.6, B=1.6, g=2
#'	x=seq(0, 5, 0.1)
#'	y = sigmoid.5pl(x, A=0, D=4, C=1.6, B=1.5, g=2); plot(log(x), y)
#'
#' @author	José R. Valverde, CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#'
#
sigmoid.5pl <- function(x, A, B, C, D, g) {
    y = D + ( (A - D) / ((1 + (x/C)^B)^g) )
    return(y)
}


#' sigm.5pl.LiaoLiu()
#'
#' sigmoid.5pl() computes the five-parameter logistic sigmoid of Liao & Liu
#'
#                     A - D
# y = D + -------------------------------
#                                      g
#         +-                      B -+ 
#         |       1/g       /  x \   |
#         | 1 + (2    -1) * | --- |  |
#         |                 \  C*/   |
#         +-                        -+
#
#' As in the 4 parameter logistic formula for a sigmoid, here, A and D
#' are the two asymptotes, B is the shape parameter or slope, C* is the 
#' ED50 value with is the concentration/dose/x corresponding to a response
#' hafway between the asymptotes, and g is a fifth parameter that introduces
#' an asymmetry on the function shape.
#'
#'
#' The ED50 value in this version is now C* =  C * (2^(1/g) - 1)^(1/B) which 
#' makes it now easy egain to estimate from a curve fit as it is now obtained
#' directly on reparametrization.
#'
#' @param	x	The value(s) whose Y must be computed
#' @param	A	Lower (left) asymptote
#' @param	B	Shape parameter
#' @param	C	ED50 value, x value corresponding to a response alfway
#'			between the asymptotes
#' @param	D	Upper (right) asymptote
#' @param	g	Asymmetry factor
#'
#' @return	The Y value(s) corresponding to the specified sigmoid for X
#'
#' @usage	sigmoid.5pl(x, A, B, C, D, g)
#'
#' @examples
#'	## e.g. A=0, B=2, C*=1.5, D=4, g=3
#'	x=seq(0, 5, 0.1)
#'	y = sigm5pl.LiaoLiu(x, A=0, B=2, Cstar=1.6, D=4, g=3) ; plot(log(x), y)
#'
#'	## e.g. A=0, D=4, C=1.6, B=1.6, g=2
#'	x=seq(0, 5, 0.1)
#'	y = sigmoid.5pl.LiaoLiu(x, A=0, D=4, Cstar=1.6, B=1.5, g=2); plot(log(x), y)
#'
#' @author	José R. Valverde, CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#'
#
sigmoid.5pl.LiaoLiu <- function(x, A, B, Cstar, D, g) {
    y = D + ( (A - D) / ((1 + (2^(1/g) - 1) * (x / Cstar)^B)^g) )
    return(y)
}



#' sigm
#'
#' sigm() General four-parameter sigmoid function
# This is a sigmoid function that we will model. We'll use this function
# to generate values using the predicted model.
#' 
#' @param	p	A four valued vector with the parameters for the sigmoid
#'			p[1] lower asymptote
#'			p[2] carrying capacity
#'			p[3] growth rate
#'			p[4] time_max
#' @param	x	The x value whose corresponding y value we want
#' 
#' @return	the y value for x of the sigmoid function defined by p
#'
#' @usage	sigm(p, x)
#' 
#' @examples
#'		sigm(1, 3)
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
sigm <- function(x, p) {
    p[1] + ((p[2] - p[1]) / (1 + exp(-p[3] * (x - p[4]))))
}

# the same function explained
### Sigmoid function ### create a function to generate a sigmoid pattern
#' sigmoid4
#'
#' sigmoid4() General four-parameter sigmoid function
# This is the sigmoid function that we will model. We'll use this function
# to generate values using the predicted model.
#' 
#' @param	x	The x value(s) whose corresponding y value(s) we want
#' @param	lower_asymptote
#' @param	carrying capacity
#' @param	growth rate
#' @param	time_max	time of maximal growth rate
#' 
#' @return	the y value for x of the sigmoid function defined by p
#'
#' @usage	sigm(p, x)
#' 
#' @examples
#'		sigm(1, 3)
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
sigmoid.4 <- function(x, lower_asymptote, carrying_capacity, growth_rate, time_max) {
    return(lower_asymptote + (
    		(carrying_capacity - lower_asymptote)
                /
		(1 + exp(-growth_rate * (x - time_max)))
		)
	   )
}


# another version which uses significant parameters to calculate the
# correct parameters for the curve
# https://stats.stackexchange.com/questions/265266/adjusting-s-curves-sigmoid-functions-with-hyperparameters
#
#' s_curve()
#'
#' s_curve() computes a sigmoid curve using specific convenience parameters
#'
#' s_curve takes convenient parameters and generates a sigmoid curve that
#' complies with them. It does it by using these parameters to derive the
#' actual parameters needed to fit the sigmoid curve.
#' The slope of the curve is defined by the difference between x50L and x50U
#' (the closer they are, the steeper the curve will be).
#'
#' @param	x	the values for which the function will be evaluated
#' @param	ymin	floor (or lower asymptote), actually left asymptote
#' @param	ymax	ceiling (or upper asymptote), rather right asymptote
#' @param	x50L	lower bound containing 50% of the values
#' @param	x50U	upper bound containing 50% of the values
#'
#' @return	the y value corresponding to x
#'
#' @usage	s_curve(x, ymin, ymax, x50L, x50U)
#'
#' @examples 
#'	library(tidyverse)
#'	df <- tibble(x = seq(0, 30, 0.5)) %>% 
#'	      mutate(y = s_curve(x, 
#'                   ymin=5000, ymax=20000, x50L=11, x50U=14))
#'	## Plot graph:
#'	ggplot(df, aes(x, y)) + geom_line() + xlim(0, 30) + ylim(0, 35000)
#'	
#' @author: José R. Valverde, CNB-CSIC
#'
#' @license: EU-GPL
#'
#' @export
#	
s_curve <- function(x, ymin, ymax, x50L, x50U) {
    # Where ymin is the floor, ymax is the ceiling, x50L & x50U are the 
    # lower and upper bounds for containing 50% of the values.
    #
    # Example:
    # library(tidyverse)
    # df <- tibble(x = seq(0, 30, 0.5)) %>% 
    #    mutate(y = s_curve(x, ymin = 5000, ymax = 20000, x50L = 11, x50U = 14))
    # Graph:
    # ggplot(df, aes(x, y)) + geom_line() + xlim(0, 30) + ylim(0, 35000)
    #

    a = (x50L + x50U) / 2	# mid point
    b = 2 / abs(x50L - x50U)	# slope/growth rate at mid point
    c = ymin			# lower asymptote
    d = ymax - c		# step to upper asymptote

    y = c + ( d / ( 1 + exp( b * (x - a) ) ) )

    return(y)
}


# for dose-response curves, we may use (assuming a constant slope)
#
#               (Top - Bottom)
# Y = Bottom + -----------------
#                    LogEC50 - X
#              1 + 10
#
# X is the logarithm of agonist concentration and Y the response, LogEC50 is
# the logarithm of the EC50 (effective concentration, 50%) --sometimes ED50
# (effective dose 50%) or IC50 (inhibitory concentration, 50%).
#
# A more generic curve would be
#
#                     (Top - Bottom)
# Y = Bottom + -----------------------------
#                    (LogEC50 - X)*HillSlope
#              1 + 10
#
# where Hill Slope defines the steepness of the curve, slope factor or Hill
# coefficient. If positive the curve increases, if negative, decreases. A 
# standard curve would have a Hill slope of 1. When > 1 the curve is steeper.
#


# general 3-parameter sigmoid function. Check the reference for interpretation.
#    a / (1 + exp(b - (x * c)))
# this one would be equivalent to the next function with different meaning
#    a / (1 + exp( -cx + b)) = 
#    K / (1 + exp(-B x + (BM))) = 
#    K / (1 + exp(-B * (x-M)))
# thus c=max growth rate, b=c * midpoint, a=asymptote
sigmoid.3 <- function(x, params) {
    a = params[1]
    b = params[2]
    c = params[3]
    
    a / (1 + exp( b - (x * c) ) )
}

# alternate formulation (note the change in meaning of the paramters)
# general 3-parameter sigmoid function, where
#	K = upper asymptote
#	B = growth rate at the midpoint, when growth is steepest
#	M = time of maximum growth
sigmoid.3 <- function(x, params) {
    K = params[1]
    B = params[2]
    M = params[3]
    
    K / (1 + exp( -B * (x - M) ) )
}
# or, alternatively (note the reversion of -B and (x-M) to B and (M-x)
#                 asymptote
# f(t) = ---------------------------
#         1 + exp(scale * (mid - t)
# where
#	mid is the x (or t) coordinate of the midpoint, where growth is
#	steepest
#	asymptote is the maximum value
#	scale indicates how steep is the curve
#
# since growth rate changes with time, we can compute it taking f'(x)
# which, at the middle point (where mid=t) becomes
#	maximal growth rate = slope(mid) = f'(mid) = asymptote * (scale / 4 )
#


sigmoid.normalized <- function(x) {
# example:
# 	x <- seq(-5, 5, 0.01)
#	plot(x, normalized_sigmoid(x), col='blue')

	1 / (1 + exp(-x))
}

#
# Note that 
# package pracma has a two parameter sigmoid function
# sigmoid(x, a, b) that implements the formula y = 1/(1+e^[-a(x-b)])
# whic is a solution to the ODE y' = y(1-y)
#
# tropFishR contains a VBGF Von Bertalanffy Growth Function.
#
# drc contains a general gompertz function and dose-response functions.
#
# e1071 contains a normalized sigmoid function for use in IA/ML applications,
# and package sigmoid contains functions for use in ML
#

#
#--------------------------- End of sigmoid functions ----------------------
#

# unused
# Alternate implementation (the one afterwards will preempt this one)
# A function to convert "observed" values (note that they might as well be
# the result of interpolation) to rates (Delta(x) / Delta(y))
#	Note that growth uses a different "rate" (mu, see below)
# This implementation offsets (shifts) all values once to the left and may make
# more sense under some circumstances. It all depends on how do we provide
# the values
#
obs2rates <- function(x, y) {
    n <- length(x)
    mu <- 1:n

    # There is no value previous to i=1
    for (i in 2:n) {
    	# this should allow us to specify jumps
        if (x[i] == x[i-1]) next
        mu[i] <- (y[i] - y[i-1]) / (x[i] - x[i-1]);
    }
    mu[1] <- mu[2] # start from a non-zero value

    return( data.frame(time=x, rate=mu) )
}

# unused
# A function to convert "observed" values (note that they might as well be
# the result of interpolation between a few key points). There must be
# one value per time point, and a time point may be repeated.
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
        # for a time point, then after the calculation, we would add/remove
        # nutrients: the next time point will start from the new concentration,
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
# DFBA. This assumes that we start at time zero
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


#' loadDynamicRates
#'
#' loadDynamicRates() reads from a file the changes to reaction rates to be
#'	applied at which simulation steps
#'
#' In an ADFBA simulation reaction limits need not neccessarily be fixed
#' beforehand and may be changed at any time. The dynamic evolution of
#' reaction rates can be specified in a file containing only the times at
#' which rates change and, for each time, the new upper and lower limits
#' of the reactions whose limits change.
#'
#' Since reactions have two limits (upper and lower), we need to specify
#' whenever there must be a change, the reaction name, which limit is
#' changed and its new value. The limit is specified by appending to the
#' reaction name [upp] or [low].
#'
#' A table may look like the following example
#'
#' time    Biomass[low]   EX_glc[upp]
#' 10      1.3	          0.5
#' 20      0.9            0.5
#'
#' This implies that, from time 10 onwards, the lower limit for the reaction
#' named 'Biomass' in the model will be 1.3, and the upper limit for the 
#' reaction named 'EX_glc' will be 0.5. At time 20, these limits will be
#' changed again, 'Biomass' lower limit will be set to 0.9 and 'EX_glc' upper
#' imit will be maintained at 0.5
#'
#' @param	data	the name of the file that contains the data
#' @param	verbose	the verbosity level when producing output
#' 
#' @return	the table read
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
loadDynamicRates <- function(data='', verbose=0) {
    if (data == '') {
    	return (NULL)
    }
    if (file.exists(data)) {
    	dat <- read.table(data, sep='\t', header=TRUE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.dat', sep=""))) {
    	dat <- read.table(paste(data, '.dat', sep=""), sep='\t', header=TRUE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.tsv', sep=""))) {
    	dat <- read.table(paste(data, '.tsv', sep=""), sep='\t', header=TRUE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.tab', sep=""))) {
    	dat <- read.table(paste(data, '.tab', sep=""), sep='\t', header=TRUE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.csv', sep=""))) {
        dat <- read.table(paste(data, '.csv', sep=""), sep=',', header=TRUE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if ((is.element('xlsx', installed.packages()[,1])) &&
             (file.exists(paste(data, ".xlsx", sep="")))) {
        dat <- read.xlsx(paste(data, ".xlsx", sep=""), sheet.index=1)
    } 
    else {
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

# Setup experimental data
# -----------------------
#

#' interpolateDynamicRates
#
#' interpolateDynamicRates() takes as input a short table of rates to be
#' set at specific time points and interpolates intermediate values.
#'
#' This function reads a file with experimental values measured at specific
#' times, and interpolates the values for the times that we are going to
#' simulate, returning a data frame that we can use to steer the dynamically
#' adjusted FBA calculations.
#'
#' Data is read from a file in TAB-separated format
#'	Data MUST be organized as follows
#'	First line contains 'Time' and reaction names
#'	Next values containe observed values for each reaction at each time
#'	Missing values are indicated by NA
#'
#' Values are interpolated from the data in the file to fill in
#' a table with nSteps values separated tStep hours each.
#'
#' The Biomass reaction may be specified so it can be approximated
#' to a Sigmoid function.
#'
#' Any non-Biomass reaction will be approximated by linear interpolation
#' by default, though other interpolation methods may be specified.
#'
#' If you also want Biomass interpolated using the same method as all
#' other reactions, simply do not specify it as the Biomass reaction.
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

    # create output data frame with target times (we start at tStep * 1 = tStep
    out = data.frame( Time = seq(tStep, nSteps * tStep, by=tStep))

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
                # Sigmoid interpolation only makes sense if Biomass
                # has been given as observed values
		#
                # Do sigmoid interpolation followed by conversion to rate.
                # We use here a formula with four parameters corresponding
                # to function 'sigm()' above
	        mdl <- nls(y ~ a + ((b - a) / (1 + exp(-c * (x - d)))), 
	            start=list(a = min(y), b=max(y), c= 1, d=median(x)), 
		    trace=FALSE, algorithm="port")
		p <- coef(mdl)
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
	    # In principle, linear is recommended as it makes no assumptions
            # and slines may lead to anomalies like negative concentrations.
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
    	dat <- read.table(data, sep='\t', header=TRUE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.dat', sep=""))) {
    	dat <- read.table( paste(data, '.dat', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.tsv', sep=""))) {
    	dat <- read.table( paste(data, '.tsv', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.tab', sep=""))) {
    	dat <- read.table( paste(data, '.tab', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.txt', sep=""))) {
    	dat <- read.table( paste(data, '.txt', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.csv', sep=""))) {
        dat <- read.table( paste(data, '.csv', sep=""), sep=',', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.Rdata', sep=""))) {
        dat <- load (paste(data, '.Rdata', sep=""))
    } 
    else if (file.exists(paste(data, '.R', sep=""))) {
        # it is an R file that should contain a function named like the file
        source(paste(data, '.R', sep=""))
        # we could have used return(dget(data'.R') but then the file should
        # contain only an anonymous function, this is more versatile.
        # Check if an object named as data exists and if it is a function
        if (exists(data) && is.function(get(data))) {
            return(eval(parse(text=data)))
        } else if (exists(data) && 
                   (is.array(get(data)) || is.data.frame(get(data)))) {
            # in case that instead of a function we got a data definition
            # (matrix is a kind of array)
            dat <- data
        }else {
            cat('getRateChanges: ERROR - ', data, '.R should contain\n',
                '    a function ', data, '(model, concentrations, fluxes, step)', sep='')
            return(NULL)
        }
        ### JR ### NEEDS TO BE TESTED
    } 
    else {
        cat("\n\ngetRateChanges: File", data, "could not be read!\n\n\n")
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
                # we use a four parameter formula here that corresponds
                # to function 'sigm()' above.
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
    	dat <- read.table(data, sep='\t', header=TRUE, check.names='FALSE',
            fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.dat', sep=""))) {
    	dat <- read.table( paste(data, '.dat', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.tsv', sep=""))) {
    	dat <- read.table( paste(data, '.tsv', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.tab', sep=""))) {
    	dat <- read.table( paste(data, '.tab', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.txt', sep=""))) {
    	dat <- read.table( paste(data, '.txt', sep=""), sep='\t', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.csv', sep=""))) {
        dat <- read.table( paste(data, '.csv', sep=""), sep=',', header=TRUE,
            check.names='FALSE', fill=T, strip.white=T, blank.lines.skip=T)
    } 
    else if (file.exists(paste(data, '.Rdata', sep=""))) {
        dat <- load (paste(data, '.Rdata', sep=""))
    } 
    else if (file.exists(paste(data, '.R', sep=""))) {
        # it is an R file that should contain a function named like the file
        source(paste(data, '.R', sep=""))
        # we could have used return(dget(data'.R') but then the file should
        # contain only an anonymous function, this is more versatile.
        # Check if an object named as data exists and if it is a function
        if (exists(data) && is.function(get(data))) {
            return(eval(parse(text=data)))
        } else if (exists(data) && 
                   (is.array(get(data)) || is.data.frame(get(data)))) {
            # in case that instead of a function we got a data definition
            # (matrix is a kind of array)
            dat <- data
        }else {
            cat('getNutrientChanges: ERROR - ', data, '.R should contain\n',
                '    a function ', data, '(model, concentrations, fluxes, step)', sep='')
            return(NULL)
        }
        ### JR ### NEEDS TO BE TESTED
    } 
    else {
        cat("\n\ngetNutrientChanges: File", data, "could not be read!\n\n\n")
	return(NULL)
    }
    
    # provide feedback
    cat("\ngetNutrientChanges: processing", data, "\n\n")
    print(dat)
    
    # get list of nutrient names
    # (the corresponding exchange reactions would be of the form
    # 'EX_'+nutrient+'(e)': e.g. ala_L -> EX_ala_L(e)
    nutrients <- names(dat)
    nutrients <- nutrients[ -which(nutrients %in% 'Time') ]

    tStep = as.numeric(tStep)	# in case it is passed as a string
    nSteps = as.numeric(nSteps)	# in case it is passed as a string

    # create output data frame with target times
    # starting at first time step
    out = data.frame( Time = seq(tStep, nSteps * tStep, by=tStep))
    #out = data.frame( Time = seq(1, nSteps * tStep, by=tStep))

    if ('yes' == 'yes') {
        # do it nutrient by nutrient
        times <- dat['Time'][,1]	# Prepare a vector with time values

        for (n in nutrients) {
            # create an empty vector of nSteps zeros
            changes=c(rep(0, nSteps))
            # set the values we got in 'dat'
            changes[times]=dat[,n]
            # join the vector with the output data frame
            # since we get nutrient names, we need to add the reation afixes
            out[ paste('EX_', n, '(e)', sep='') ] <- changes
        }
    } else  {
        # do it time by time
        times <- dat['Time'][,1]	# Prepare a vector with time values

        # initialize the dataframe to all zeros
        for (n in nutrients) {
	    out[n] <- rep(0.0, nSteps)
        }

        for (i in  1:length(dat$Time)) {
            t <- dat$Time[i]
            if ("mode" == "closest") {
                target.time <- which.min(abs(out$Time - t))
                # the data is in dat[t, reactions]
                # should go in out[target.time, reactions]
                # where reactions are of the form 'EX_'+nutrient-name+'(e)'
                out[target.time, paste('EX_', nutrients, '(e)', sep='')] <- dat[t, reactions]
            } else if ("mode" == "closest-after") {
                # set in the first time that is >= the specified time
                # that is: if the time is matched, then use it, if not at the
                # closest time afterwards
                #target.time <- Position(function(x) {x >= t}, out$Time)
                target.time <- which(out$Time >= t)[1]
                out[target.time, paste('EX_', nutrients, '(e)', sep='')] <- dat[x, nutrients]
            }
        }
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
# FBA only returns a single flux distribution that corresponds to maximal growt
# under given growth conditions. However, alternate optimal solutions may exist
# which correspond to maximal growth. FVA calculates the full range of numerical
# values for each reaction flux within the network.
# 
# The function fluxVar performs a flux variability analysis with a given model 
# [Mahadevan and Schilling, 2003].  The minimum and maximum flux values for 
# each reaction in the model are calculated, which still support a certain 
# percentage of a given optimal functional state Z_opt.
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
    	subs <- read.table(paste(substrate, ".tsv", sep=""), sep="\t", 
        	header=FALSE, 
            	check.names=FALSE, fill=TRUE, strip.white=T, blank.lines.skip=T)
    } else if (file.exists(paste(substrate, ".tab", sep=""))) {
     	subs <- read.table(paste(substrate, ".tab", sep=""), sep="\t", 
        	header=FALSE, 
            	check.names=FALSE, fill=T, strip.white=T, blank.lines.skip=T)
    } else if (file.exists(paste(substrate, ".csv", sep=""))) {
    	subs <- read.table(paste(substrate, ".csv", sep=""), sep=",", 
        	header=FALSE, 
            	check.names=FALSE, fill=T, strip.white=T, blank.lines.skip=T)
    } else if (file.exists(paste(substrate, ".txt", sep=""))) {
    	subs <- read.table(paste(substrate, ".txt", sep=""), header=FALSE)
    } else if (file.exists(paste(substrate, ".dat", sep=""))) {
    	subs <- read.table(paste(substrate, ".dat", sep=""), sep="\t", 
        	header=FALSE, 
            	check.names=FALSE, fill=T, strip.white=T, blank.lines.skip=T)
    } else if (file.exists(substrate)) {
        # last resort, try to read it as a TAB file
        subs <- read.table(substrate, sep="\t", header=FALSE, 
            	check.names=FALSE, fill=T, strip.white=T, blank.lines.skip=T)
    } else {
        cat("File ", substrate, "{.txt, .tsv, .tab, .csv} not found\n")
    	subs <- read.table(file.choose(), header=FALSE)
    }
    
    # loop over all the substrates and apply them
    # We'll read in
    #    Metabolite	Mol.Weight	mg/L
    # and compute the concentration in mmol/L
    # When g/L < 0 we'll take it as an excluded uptake reaction
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
    outConcentrations <- paste (dfbaDir, '/', modelName, "_concs.tab", sep="")
    outFluxes <- paste (dfbaDir, '/', modelName, "_fluxes.tab", sep="")
    outBiomass <- paste (dfbaDir, '/', modelName, "_biomass.tab", sep="")


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

    # write out the concentrations and fluxes obtained
#    write.table(df_sol@concentrationMatrix, file=outConcentrations, sep='\t', 
#    	quote=FALSE, row.names=TRUE, col.names=TRUE)
#    write.table(df_sol@all_fluxes, file=outFluxes, sep='\t', 
#    	quote=FALSE, row.names=TRUE, col.names=TRUE)
    saveConcentrations(af_sol, outConcentrations)
    saveFluxes(af_sol, outFluxes)
    write.table(df_sol@biomassVec, file=outBiomass, sep='\t', 
    	quote=FALSE, row.names=TRUE, col.names=TRUE)

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
    outConcentrations <- paste (adfbaDir, '/', modelName, "_concs.tab", sep="")
    outFluxes <- paste (adfbaDir, '/', modelName, "_fluxes.tab", sep="")
    outBiomass <- paste (adfbaDir, '/', modelName, "_biomass.tab", sep="")


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
    # write out the concentrations and fluxes obtained
#    write.table(af_sol@concentrationMatrix, file=outConcentrations, sep='\t', 
#    	quote=FALSE, row.names=TRUE, col.names=TRUE)
#    write.table(af_sol@all_fluxes, file=outFluxes, sep='\t', 
#    	quote=FALSE, row.names=TRUE, col.names=TRUE)
    saveConcentrations(af_sol, outConcentrations)
    saveFluxes(af_sol, outFluxes)
    write.table(af_sol@biomassVec, file=outBiomass, sep='\t', 
    	quote=FALSE, row.names=TRUE, col.names=TRUE)


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



lplot <- function(x,y,
                 ylim=50,
                   xlab = "",
                   ylab = "Value",
                   type = "p",
                   pch = 20,
                   col = "black",             
    #               collower, colupper, pchupper, pchlower,
    #               dottedline = TRUE,
                  plotRxns=NULL,
                   baseline = 0,
                   legend.pos="left",
                    ...) {
    if(missing(plotRxns)){
       plot(x@timeVec,x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
       }
    else {
	    def.par <- par(no.readonly = TRUE);
	    layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE))
	    #layout.show(2);
	    # first plot biomass
	     plot(x@timeVec,x@biomassVec, col = 1,
			        main='Cell density',
                                xlab='Time(hrs)',
                                ylab="X(g/l)",
                                type="l", lwd=2);
	     points(x@timeVec,x@biomassVec, col = "red",lwd=2);

	    # plot concentrations
	    ##plot(x@timeVec,2*x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
	    ## define min/max ()plot(x@timeVec,
	    ymin <- min(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
            ymax <- max(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
	    for ( i in 1:length(plotRxns) ) {
		    plotInd=(x@excRxnNames %in% plotRxns[i]);
		    #print( x@concentrationMatrix[plotInd]);
		    if(i==1) {
                       plot(x@timeVec, x@concentrationMatrix[plotInd],
                       type="l", col =i,main="Concentrations",ylab = "mmol/l",xlab='Time(hrs)'
                       ,ylim=c(ymin,ymax));
                       
                    } else {
       		      lines(x@timeVec, x@concentrationMatrix[plotInd], col=i);
       		    }
	    }
	    legend(legend.pos, plotRxns, col=1:length(plotRxns), lty=1);
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
    defst <- 90			# default number of time steps
    definoc <- 0.10		# default inoculum size in g/L (initBiomass)
    defrarxn <- "EX_o2(e)"	# default reaction to use for RA
    defpppax <- "EX_o2(e)"	# default for first reaction in PPPA
    defpppay <- "EX_co2(e)"	# default for second reaction in PPPA
    defexch <- ""		# default time-variation of exchange rates
    defdeltas <- ""		# default substrate deltas per time step
    defmedium <- "medium"	# default growth medium composition
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

        make_option(c("-n", "--medium"), 
            type="character", 
            default=defmedium, 
            help="nutrient medium file name [default '%default']
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

        make_option(c("-l", "--deltas"), 
            type="character", 
            default=defdeltas,
            help=paste("substrate concentration deltas at specific time points (mmol/L) 
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
    substrateName <- opt$medium		# substrate to use for DFBA
    ts <- opt$timestep / 60		# convert to hours which is the reference time unit
    ns <- opt$steps			# number of steps
    inoculum <- opt$inoculum		# initial Biomass
    biomass <- opt$biomass		# name of biomass function
    exchangeRates <- opt$exch		# "" if not specified
    substrateDeltas <- opt$deltas	# "" if not specified
    
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
    if (verbose > 3) { 
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
    if (exchangeRates != "")
	outName <- paste(outName, "+", exchangeRates, sep="")
        #outName <- paste(outName, "=", file_path_sans_ext(basename(exchangeRates)), sep="")
    if (substrateName != "")
        outName <- paste(outName, "@", substrateName, sep="")
    if (substrateDeltas != "")
        outName <- paste(outName, "+", substrateDeltas, sep="")

    outName <- paste(outName, ":", ts, "~", ns, sep='')
#    outName <- paste(outName, ":", ts, "~", ns, sep='')

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
                    dcfile <- paste(outDir, '/', modelName, "_dc.R", sep="")
                    # TEST THIS. This will ONLY copy the named file, but
                    # it won't copy any additional files sourced by it, hence,
                    # it might save an incomplete copy... corollary: in this
                    # case, and as an exception, you better keep all your eggs
                    # in the same basket :-{)}
                    file.copy(exchangeRates, dcfile)
                } else {
                    # try to read a data (non-R-source) file into a data.frame
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
                    # TEST THIS. This will ONLY copy the named file, but
                    # it won't copy any additional files sourced by it, hence,
                    # it might save an incomplete copy... corollary: in this
                    # case, and as an exception, you better keep all your eggs
                    # in the same basket :-{)}
                    file.copy(substrateDeltas, ncfile)
                } else {
                    nutrientChanges <- getNutrientChanges(substrateDeltas,
                        timeStep, nSteps)
                    ncfile <- paste(outDir, '/', modelName, "_nc.tab", sep="")
                    write.table(nutrientChanges, file=ncfile, sep='\t')
                    # Note that you can provide delayed concentration-dependent 
                    # changes by keeping a log of back concentrations in a
                    # global variable and using the <<- operator to update
                    # it inside this function and then reacting to a past value
                    # instead of the current one.
                    # e.g. delay three time steps in reacting to O2 depletion
                    # lastO2conc <- c(0.0, 0.0, 0.0 0.0)
                    # then inside the function
                    #	lastO2conc[ step %% 4 ] <- currconc
                    #	if (lastO2conc[ (step - 3) %% 4 ] < 1) 
                    #		currcon <- currconc + supplement
                    # a similar trick will allow detecting trends and changes
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
            #method <- "lpMTF"
            #method <- "directMTF"

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
