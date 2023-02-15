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
