library(adfba)
library(sybilSBML)

load('ref/SlividansWT.RData')
verbose <- 4
sink('ref/REF')
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
                    verboseMode=verbose);
plot(af_sol, plotRxns=plotRxns);
sink()
png('ref/plotADFBA.png', width=750, height=750)
plot(af_sol, plotRxns=plotRxns);
dev.off



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


write.table(dynamicConstraints, "rates_ref.dat", sep='\t')
names(initConcentrations) <- substrateRxns
write.table(initConcentrations, "medium.dat", sep='\t')
saveModel(model)
