################################################
# Function: dynamicFBA
#
# Performs a dynamic flux balance analysis
# 
# The function dynamicFBA() is inspired by the function
# dynamicFBA() contained in the COBRA Toolbox.
# The algorithm is the same.

dynamicFBAplus <- function (model,
			substrateRxns,
                        initConcentrations,
                        initBiomass,
                        timeStep,
                        nSteps,
                        exclUptakeRxns,
			retOptSol = TRUE,
                  	fld = FALSE,
                        verboseMode = 2, ...){
    #PARAMETERS:
    #===========
    # model                 Sybil model structure (class modelorg)
    # substrateRxns         List of exchange reaction names for substrates
    #                       initially in the media that may change (e.g. not
    #                       h2o or co2)
    # initConcentrations    Initial concentrations of substrates (in the same
    #                       structure as substrateRxns)
    # initBiomass           Initial biomass (must be non zero)
    # timeStep              Time step size
    # nSteps                Maximum number of time steps
    # fld                   indicates if all fluxes at all steps will be returned.
    # retOptSol             indicates if optsol calss will be returned or simple list
    #
    #OPTIONAL PARAMETERS
    #===================
    # exclUptakeRxns        List of uptake reactions whose substrate
    #                       concentrations do not change (Default =
    #                       {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'})
    # 
    #RETURN VALUES:
    #=============
    # concentrationMatrix   Matrix of extracellular metabolite concentrations
    # excRxnNames           Names of exchange reactions for the EC metabolites
    # timeVec               Vector of time points
    # biomassVec            Vector of biomass values
    # all_fluxes            Matrix containing the fluxes of all reactions at different steps

    # aids for debugging/verbose messages
    myName <- 'adaptiveDFBA '
    info <- paste(myName, 'INFO: ')
    warn <- paste(myName, 'WARNING: ')
    err  <- paste(myName, 'ERROR: ')

    optsol.OK = 5

    # Try to find the real Biomass function before accepting to use only
    # the first objective function (which might even be multiple objectives)
    #
    # THIS IS VERY, VERY RISKY
    #	It might be that the biomass reaction is called something else
    #	and the one named *biomass* is not the real Biomass reaction
    #	or that there is more than one reaction called '*biomass*'
    #   and the proper one is not the first. Hopefully we will get it
    #	right
    #
    biomassIdx <- grep('biomass', react_id(model), ignore.case=TRUE)[1]
    # if found use it
    if (biomassIdx != 0)
        biomassRxn <- react_id(model)[biomassIdx]
    else {
        # if not, revert to first objective function
        biomassIdx <- which(sybil::obj_coef(model) != 0)[1]
        biomassRxn <- sybil::react_id(model)[biomassIdx][1]
    }
    cat(info, 'Biomass reaction:', biomassRxn, "(",biomassIdx, ')\n')
    cat(warn, 'IF THIS IS WRONG, USE INSTEAD adaptiveDFBA() TO SET IT\n\n')

    #
    # If no initial concentration is given for a substrate that has an open
    # uptake in the model (i.e. model.lb < 0) the concentration is assumed to
    # be high enough to not be limiting. If the uptake rate for a nutrient is
    # calculated to exceed the maximum uptake rate for that nutrient specified
    # in the model and the max uptake rate specified is > 0, the maximum uptake 
    # rate specified in the model is used instead of the calculated uptake
    # rate.

    if (verboseMode > 4) {
        cat(info, 'Provided Concentrations\n')
        #print(substrateRxns)
        #print(initConcentrations)
        for (i in 1:length(substrateRxns))
            cat(i, substrateRxns[i], initConcentrations[i], '\n', sep="	")
    }

    ##--------------------------------------------------------------------------##
    # check prerequisites 
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }
    ##--------------------------------------------------------------------------##

    # Uptake reactions whose substrate concentrations do not change
    if (missing(exclUptakeRxns)){
        exclUptakeRxns = c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)');
    }
    if (verboseMode > 2){
       cat(info, 'Excluded uptake reactions\n')
       cat(info, '(substrate concentration does not change:\n')
       for (i in exclUptakeRxns) cat(i, '	')
       cat('\n\n')
    }

    ##--------------------------------------------------------------------------##
    # Find exchange reactions
    excReact = findExchReact(model);
    excReactInd=(react_id(model) %in% react_id(excReact));#excReact$exchange
    #represent extra cellular reaction with boolean vector.
    exclUptakeRxnsInd=is.element(react_id(model) ,exclUptakeRxns);
    #Exclude reactions with concentrations that will not be changed 
    excReactInd = excReactInd & !exclUptakeRxnsInd;   #excInd & ~ismember(model.rxns,exclUptakeRxns);
    #get reaction names
    excRxnNames =react_id(model)[excReactInd];                #excRxnNames = model.rxns(excInd);

    ##--------------------------------------------------------------------------##
    # Find substrate reactions
    substrateRxnsInd=(react_id(model) %in% substrateRxns)
    # Figure out if substrate reactions are correct: all substrate reactions should be exchange reactions.
    missingSub = substrateRxnsInd & !excReactInd;
    if (sum(missingSub) != 0){
        print(sum(missingSub));
        print(react_id(model)[missingSub]);
        cat(warn, 'Invalid substrate uptake reaction!\n\n');
    }

    ##--------------------------------------------------------------------------##
    # Initialize concentrations
    #substrateMatchInd = intersect(excRxnNames,substrateRxns);
    #table(excRxnNames);
    #vector(length=length(excRxnNames),mode="numeric");
    #concentrations[1:length(concentrations)]=0;
    #
    # we'll use one concentration per every possible reaction
    # and we'll start from a concentration of zero by default
    concentrations=rep(0,length(react_id(model)))
    # XXX JR XXX
    # THERE IS AN INCONSISTENCY HERE: WHEN SUBSTRATE RXNS ARE NOT IN THE SAME
    # ORDER AS THEY ARE IN THE MODEL, THIS ASSIGNMENT WILL ASSIGN THE
    # CONCENTRATIONS ERRONEOUSLY SINCE THE ORDER IN initConcentrations IS 
    # NOT SORTED BUT THE INDEX OF SUBSTRATE REACTIONS IS!!!
    #concentrations[substrateRxnsInd] = initConcentrations;
    for (i in 1:length(substrateRxns)) {
        # note that there MUST be a concentration for each substrate
        # and that concentrations and substrates MUST have been provided 
        # in the same order!!!
        concentrations[react_id(model) == substrateRxns[i]] <- initConcentrations[i]
    }

    # save original lower bounds for all reactions
    originalBound <- -lowbnd(model);# take all to be able to directly update

    # Deal with uptake reactions for which there are no initial concentrations
    #	Corollary: for every uptake reaction there MUST be an initial
    #	concentration that is different from zero, otherwise 1000 will be used 
    #
    #	This should likely be removed, requiring everything to be explicit
    #
    #(concentrations == 0 & originalBound > 0);
    # this will set the concentration for ALL metabolites to 1000 except
    # for those in the medium !!!
    # it will cancel also any substrate explicitly set to zero
    noInitConcentration <- (concentrations==0)&(lowbnd(model)<0)
    #
    # it would likely be better to use
    #	concentrations is for all reactions in the model
    #	select those that
    #		are zero
    #		had a lower bound < 0
    #		are exchange reactions
    #		and were NOT explicitly provided
    noInitConcentration <- (concentrations==0) &
    			   (lowbnd(model)<0) &
                           (excReactInd) &
                           (!(substrateRxnsInd))
                          
    concentrations[noInitConcentration] <- 1000;
    # still, this will set the concentration of any potential substrate
    # whose concentration has not been specified to 1000, providing potentially
    # unwanted substrates


    biomass = initBiomass;

    ##--------------------------------------------------------------------------##
    # Initialize bounds
    #	compute maximal uptake bound as the concentration available per
    #	biomass unit per unit of time lapse
    #	I.e.: this is the maximum amount of each metabolite available,
    #	hence cells cannot consume more than this even if the maximum
    #	uptake limit in the model permits higher potential consumption
    #
    #	My it be that this could cause a problem if the computed uptake
    #	where smaller than the upperbound?
    #	e.g. if EX_X is between (-)1.5 and (-)0.5 and the computed availability
    #   is < 0.5, then the new bounds would be (-)<0.5 and (-0.5) with
    #	lowbnd > upbnd
    #
    #	As long as the uptake limits are of the form <0 .. >=0 this shouldn't
    #   be a problem. If it were, the simulation should fail and then we
    #	can fix this below.
    # we can take at most biomass*timeStep per hour
    uptakeBound =  concentrations/(biomass*timeStep);

    # Make sure bounds are not higher than those that were specified in the model
#uptakeBound[uptakeBound > 1000] = 1000;
#uptakeBound <- ifelse(abs(uptakeBound) < 1e-9,0,uptakeBound);

    # if the computed maximal uptake bound is greater than the original
    # limit (and the original uptake bound was not zero)
    aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
    
    # then we will use the smaller original bound
    uptakeBound[aboveOriginal] = originalBound[aboveOriginal];
    
    # set the lower bound (which, since uptakes are neative, must be
    # the maximum uptake possible) to the corresponding value (i.e., either the
    # original bound (if availability is larger) or the computed bound
    # (if availability is smaller)
    lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];
    
    if (verboseMode > 2) {
        cat(info, 'Step', 0, 'uptake rate corrections:\n');
        cat('Name	PreviousBound	UptakeBound	LowBnd	UppBnd	aboveOrig\n')
        for (i in 1:length(excReactInd)) {
            if (excReactInd[i]) {
                cat(react_id(model)[i], 
                    originalBound[i], 
                    -uptakeBound[i], 
                    lowbnd(model)[i],
                    uppbnd(model)[i],
                    aboveOriginal[i],
                    '\n', sep="	")
            }
        }
        cat('\n')
    }

    ##--------------------------------------------------------------------------##
    # Initialize simulation
    concentrationMatrix = concentrations[excReactInd];
    biomassVec = biomass;
    timeVec = 0;
    
    ##------------------------------------- Prepare Problem object --------------------##
    # get OptObj instance
    #lpmod <- prepProbObj(model,
    #                         nCols      = react_num(model),
    #                         nRows      = met_num(model),
    #             #            alg        = "FBA",
    #                         solver     = solver,
    #                         method     = method,
    #                         lpdir      = lpdir
    #                         #solverParm = solverParm
    #             )
    lpmod <- sybil::sysBiolAlg(model, algorithm = "fba", ...)

    if (verboseMode > 2) {
        cat(info, 'Step number    Biomass\n');
        cat('0    ', biomass, '\n\n', sep='	')
    }
    if (verboseMode > 3) {
        cat(info, 'Concentration of substrates\n')
        for (i in substrateRxns)
            cat(i, concentrations[react_id(model) == i], '\n', sep="	")
        cat('\n')
    }
    if (verboseMode > 4) {
        cat(info, 'All Concentrations\n')
        for (i in 1:length(substrateRxnsInd))
            cat(i, react_id(model)[i], concentrations[i], '\n', sep="	" )
        cat('\n')
    }


    # Inititialize progress bar ...');	## FROM MATLAB OpenCOBRA
    #if (verboseMode == 2)  progr <- .progressBar();
    if (verboseMode == 1) { pb <- txtProgressBar(0, nSteps) }
    if (verboseMode > 1) {
        # requires package tcltk
        library(tcltk)
        pb <- tkProgressBar(title = "DFBA progress", min = 0,
                            max = nSteps, width = 300)
    }

    all_fluxes=NULL;
    all_stat=NULL;

    ##-----------------------------------------------------------------------------##
    # Do the simulation
    for (stepNo in 1:nSteps){

        #if (verboseMode == 2)  progr <- .progressBar(stepNo, nSteps, progr);
	if (verboseMode == 1) { setTxtProgressBar(pb, stepNo) }
        if (verboseMode > 1) { 
            setTkProgressBar(pb, stepNo, 
            	label=paste( round(stepNo/nSteps*100, 0),
                "% done"))
	}

        # Run FBA
        if (verboseMode > 2) {
            fba <- sybil::optimizeProb(model, algorithm="fba");
            cat(info, "Objective function =", mod_obj(fba), "\n")
            computed_biomass = biomass * exp(mod_obj(fba) * timeStep);
            cat(info, "Computed biomass =", computed_biomass, "\n\n")
            if (verboseMode > 3) {
                cat(info, "Computing MTF\n\n")
                mtf_sol <- sybil::optimizeProb(model, algorithm="mtf", wtobj = mod_obj(fba)); 
                cat('Net flux of the exchange reactions\n');
                cat('----------------------------------\n');
                #excReact = findExchReact(model);
	        #excReactInd=(react_id(model) %in% react_id(excReact));#excReact$exchange
                fd <- getFluxDist(mtf_sol, excReact);
                print( getNetFlux(fd) )
                cat('\n')
            }
        }

        sol = sybil::optimizeProb(lpmod);
        # if the objective function is not Biomass, this is plainly wrong!
        mu =  sol$obj;  ##objvalue sol.f
	# XXX JR XXX these alternatives should fix the issue
	mu_bmnam <- sol$fluxes[biomassIdx];
        mu_bmrxn <- sol$fluxes[react_id(model) == biomassRxn]
        mu <- mu_bmnam
        if (verboseMode > 2) {
	    cat(info, 'SOL$obj:', sol$obj, '\n')
            cat(info,
                'mu(obj)', mu_obj, 
            	'mu(biomassIdx)', mu_bmnam, 
                'mu(biomassRxn)', mu_bmrxn, '\n')
            cat(info, 'using MU =', mu, '\n\n')
        }
        ## checkSolStat
        if ( length(checkSolStat(sol$stat,solver(problem(lpmod))))!=0 ){
	    # checkSolStat tells us which problem(s) could not be solved
            cat(warn, 'No feasible solution - nutrients exhausted\n');
            break;
        }
        all_stat = c(all_stat,sol$stat)
	# compute new biomass
        biomass = biomass*exp(mu*timeStep);
        #biomass = biomass*(1+mu*timeStep);
        biomassVec = c(biomassVec,biomass);
	# add the fluxes to the table of flux changes over time
        if(fld){
            if (stepNo == 1) {
	            all_fluxes = sol$fluxes;
            }else{
	            all_fluxes = cbind(all_fluxes,sol$fluxes);
            }
        }
        
	# get uptake fluxes after FBA
        uptakeFlux = sol$fluxes[excReactInd];
        # Update concentrations
        concentrations[excReactInd]= concentrations[excReactInd] - uptakeFlux/mu*biomass*(1-exp(mu*timeStep));
        #concentrations = concentrations + uptakeFlux*biomass*timeStep;
        concentrations[concentrations <= 0] = 0;
        concentrationMatrix = cbind(concentrationMatrix,concentrations[excReactInd]);
        if (verboseMode > 2) {
            cat(info, 'Concentration of substrates after', stepNo, 'steps\n')
            for (i in substrateRxns)
                cat(i, concentrations[react_id(model) == i], '\n', sep="	")
            #for (i in react_id(excReact))
            #    print(concentrations[excReactInd])
            cat('\n')
        }

        # Update bounds for uptake reactions
        uptakeBound[excReactInd] =  concentrations[excReactInd]/(biomass*timeStep);
        # This is to avoid any numerical issues
     ###uptakeBound[uptakeBound > 1000] = 1000;
        # Figure out if the computed bounds were above the original bounds
        aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
        # Revert to original bounds if the rate was too high
        uptakeBound[aboveOriginal] = originalBound[aboveOriginal];# uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
     ###uptakeBound=ifelse(abs(uptakeBound) < 1e-9,0,uptakeBound);
        ## Change lower bounds according to the result of last step
        #lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];  
        uppb_tmp <- getColsUppBnds(problem(lpmod), which(excReactInd));
        changeColsBnds(problem(lpmod),which(excReactInd),lb=-uptakeBound[excReactInd],ub=uppb_tmp);

        if (verboseMode > 2) {
            cat(info, 'Step', stepNo, 'uptake rate corrections:\n');
            cat('Name	PreviousBound	UptakeBound	LowBnd	UppBnd	aboveOrig\n')
            for (i in 1:length(excReactInd)) {
                if (excReactInd[i]) {
                    cat(react_id(model)[i], 
                        originalBound[i], 
                        -uptakeBound[i], 
                        lowbnd(model)[i],
                        uppbnd(model)[i],
                        aboveOriginal[i],
                        '\n', sep="	")
                }
            }
            cat('\n')
        }

        #waitbar(stepNo/nSteps,h);
        timeVec = c(timeVec,stepNo*timeStep);

	if (verboseMode > 2) cat(info, 'Objective at step', stepNo, ':', sol$obj, '\n');
        if (verboseMode > 2) cat(info, "Biomass at t =", stepNo, ":", biomass, '\n', sep="	");
    }# end loop

    ##--------------------------------------------------------------------------##
    # Simulation completed... prepare RETURN output
    # browser();
    row.names(concentrationMatrix)=react_id(model)[excReactInd];
    ## Preparing OUTPUT
    #concentrationMatrix,excRxnNames,timeVec,biomassVec
    if (isTRUE(retOptSol)) {
        if(is.null(all_fluxes)) all_fluxes=as.matrix(NA);
	    return (optsol_dynamicFBA(solver = solver(problem(lpmod)),
				      method = method(problem(lpmod)),
				      nprob  = stepNo,
				      ncols  = react_num(model),
				      nrows  = met_num(model),
				      fld    = fld,
				      all_fluxes = all_fluxes,
				      concmat=concentrationMatrix,
				      exRxn=excRxnNames,
						      tmVec=timeVec,  
						      bmVec=biomassVec
						      )
	      )
    }else{
		    return(optsol <- list(		  nprob  = stepNo,
				      ncols  = react_num(model),
				      nrows  = met_num(model),
				      all_fluxes = all_fluxes,
				      all_stat=all_stat,
				      concentrationMatrix=concentrationMatrix,
				      excRxnNames=excRxnNames,
						      timeVec=timeVec,  
						      biomassVec=biomassVec
			                      ))
    }
}
