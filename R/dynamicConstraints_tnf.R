# THESE ARE BASED ON TNF USING A LINEAR MODEL

# Estimate new uptake rate for ALANINE 
# using final concentrations calculated at the previous step
ala_L.rate <- function(conc) {

    return (as.numeric((-1.599804 + 476.949126 * conc["EX_tnf(e)"] + 
                    	0.021925 * conc["EX_pyr(e)"] + 
			0.538590 * conc["EX_succ(e)"] + 
			0.026429 * conc["EX_glc(e)"] + 
			0.002828 * conc["EX_nh4(e)"])))

}


# Estimate new uptake rate for TNF
# using final concentrations calculated at the previous step
tnf.rate <- function(conc) {

    return (as.numeric((1.102038e-05 + -2.784397e-07 * conc["EX_glc(e)"] + 
			9.065259e-08 * conc["EX_nh4(e)"] + 
			8.217468e-07 * conc["EX_ala_L(e)"])))

}


# Estimate new uptake rate for NH4
# using final concentrations calculated at the previous step
nh4.rate <- function(conc) {

    return (as.numeric((-51.9854350 + 0.1823635 * conc["EX_ala_L(e)"] + 
                    	0.3387304 * conc["EX_glc(e)"] + 
			-4.3104901 * conc["EX_gly(e)"] + 
			-10.1093257 * conc["EX_lys_L(e)"] + 
			-217.0650661 * conc["EX_met_L(e)"] +
			210.5804690 * conc["EX_phe_L(e)"] +
			-2.3891587 * conc["EX_pyr(e)"] +
			-5.8809109 * conc["EX_succ(e)"] +
			11.1485758 * conc["EX_thr_L(e)"] +
			0.5827572 * conc["EX_urea(e)"])))


}


# Estimate new uptake rate for GLUCOSE
# using final concentrations calculated at the previous step
glc.rate <- function(conc) {

    return (as.numeric((-0.4952746758 + 0.0005174833* conc["EX_fe3(e)"] + 
                    	-0.0332264707 * conc["EX_ala_L(e)"] + 
			-0.0066566305 * conc["EX_pyr(e)"] + 
			33.3482587003 * conc["EX_tnf(e)"])))
				
}


# Estimate new uptake rate for PYRUVATE
# using final concentrations calculated at the previous step

pyr.rate <- function(conc) {

    return (as.numeric((1.946748126 + 0.016501685 * conc["EX_asp_L(e)"]+ 
                    	0.012096337 * conc["EX_glc(e)"] + 
			-0.006298109* conc["EX_nh4(e)"] + 
			-0.219101510 * conc["EX_succ(e)"] + 
			-3.758139542 * conc["EX_phe_L(e)"] +
			-0.101717502 * conc["EX_ser_L(e)"])))
				
}


# Estimate new uptake rate for SUCCINATE
# using final concentrations calculated at the previous step

succ.rate <- function(conc) {

    return (as.numeric((0.025225789 + -0.008852523* conc["EX_glc(e)"] + 
                    	0.859733864 * conc["EX_met_L(e)"] + 
			0.002508077 * conc["EX_nh4(e)"] + 
			-72.345352214 * conc["EX_tnf(e)"] + 
			-0.007430037 *  conc["EX_pyr(e)"]+
			0.022415539 * conc["EX_thr_L(e)"])))
				
}


# Estimate new uptake rate for ASPARTIC ACID
# using final concentrations calculated at the previous step

asp_L.rate <- function(conc) {

    return (as.numeric((-0.40607880 + 174.82484063 * conc["EX_tnf(e)"] + 
                    	 0.01428527 * conc["EX_glc(e)"] + 
			-0.03793353 * conc["EX_glu_L(e)"] + 
			-0.34711620 * conc["EX_lys_L(e)"] + 
			 0.08508223 * conc["EX_pyr(e)"] +
			 0.19665256 * conc["EX_ser_L(e)"] +
			 0.29939134 * conc["EX_succ(e)"] +
			-0.36956391 * conc["EX_thr_L(e)"])))


}


# Estimate new uptake rate for GLUTAMINE
# using final concentrations calculated at the previous step

glu_L.rate <- function(conc) {

    return (as.numeric((-0.90581610 + 354.94076176 * conc["EX_tnf(e)"] + 
                    	0.02295263 * conc["EX_glc(e)"] + 
			-0.51435023 * conc["EX_lys_L(e)"] + 
			0.14282992 * conc["EX_pyr(e)"] + 
			0.31658572 * conc["EX_ser_L(e)"] +
			0.81269716 * conc["EX_succ(e)"] +
			-0.58471202 * conc["EX_thr_L(e)"])))


}


# Estimate new uptake rate for SERINE
# using final concentrations calculated at the previous step

ser_L.rate <- function(conc) {

    return (as.numeric((-0.30816902 + 150.30296353 * conc["EX_tnf(e)"] + 
                    	0.01471203 * conc["EX_glc(e)"] + 
			-0.04295528 * conc["EX_glu_L(e)"] + 
			-0.42382146 * conc["EX_lys_L(e)"] + 
			0.08510758 * conc["EX_pyr(e)"] +
			0.12974364 * conc["EX_ser_L(e)"] +
			0.14826292 * conc["EX_succ(e)"] + 
			-0.33638735 * conc["EX_thr_L(e)"])))


}


# Estimate new uptake rate for HISTIDINE
# using final concentrations calculated at the previous step

his_L.rate <- function(conc) {

    return (as.numeric((-2.015036e-04 + 4.351762e-06* conc["EX_nh4(e)"] + 
                    	-4.950312e-03 * conc["EX_tnf(e)"])))


}


# Estimate new uptake rate for VALINE
# using final concentrations calculated at the previous step

val_L.rate <- function(conc) {

    return (as.numeric((-2.176089e-03 + 4.699638e-05 *  conc["EX_nh4(e)"] + 
                    	-5.350770e-02  * conc["EX_tnf(e)"])))


}


# Estimate new uptake rate for LEUCINE
# using final concentrations calculated at the previous step

leu_L.rate <- function(conc) {

    return (as.numeric((-0.0073745352  + 0.0001592659 * conc["EX_nh4(e)"]+ 
                    	-0.1813291797  * conc["EX_tnf(e)"])))


}


# Estimate new uptake rate for ISOLEUCINE
# using final concentrations calculated at the previous step

ile_L.rate <- function(conc) {

    return (as.numeric((-2.256686e-03 + 4.873704e-05*  conc["EX_nh4(e)"] + 
                    	-5.548912e-02  * conc["EX_tnf(e)"])))


}


# Estimate new uptake rate for PROLINE
# using final concentrations calculated at the previous step

pro_L.rate <- function(conc) {

    return (as.numeric((-4.485168e-03 + 9.686498e-05 * conc["EX_nh4(e)"] + 
                    	-1.102838e-01  * conc["EX_tnf(e)"])))


}


# Estimate new uptake rate for TYROSINE
# using final concentrations calculated at the previous step

tyr_L.rate <- function(conc) {

    return (as.numeric((-3.022410e-04  + 6.527503e-06 * conc["EX_nh4(e)"] + 
                    	-7.430670e-03   * conc["EX_tnf(e)"])))


}


# Estimate new uptake rate for PHENYLALANINE
# using final concentrations calculated at the previous step

phe_L.rate <- function(conc) {

    return (as.numeric((2.448324e-02  + -1.962992e-05 *  conc["EX_fe3(e)"] + 
                    	-1.049245e-04  * conc["EX_nh4(e)"] + 
			-9.434223e-05 * conc["EX_ser_L(e)"] + 
			-2.475855e+00  * conc["EX_tnf(e)"])))


}


# Estimate new uptake rate for METHIONINE
# using final concentrations calculated at the previous step

met_L.rate <- function(conc) {

    return (as.numeric((-7.680819e-03  + -3.175957e-05 * conc["EX_glc(e)"] + 
                    	1.176688e-02  * conc["EX_lys_L(e)"] + 
			-1.290836e-04  * conc["EX_pyr(e)"] + 
			-2.578283e-04 * conc["EX_ser_L(e)"] + 
			5.522181e-04 * conc["EX_thr_L(e)"])))


}


# Estimate new uptake rate for LYSINE
# using final concentrations calculated at the previous step

lys_L.rate <- function(conc) {

    return (as.numeric((0.0627605438  + -0.0013970106 * conc["EX_asp_L(e)"] + 
                    	0.0014362438 * conc["EX_glc(e)"] + 
			-0.3431860815  * conc["EX_met_L(e)"]+ 
			-0.0000766168 * conc["EX_nh4(e)"] + 
			0.0003009253  * conc["EX_ser_L(e)"] +
			-0.0035808775 * conc["EX_succ(e)"])))


}


# Estimate new uptake rate for GLYCINE
# using final concentrations calculated at the previous step

gly.rate <- function(conc) {

    return (as.numeric((0.971052132  + 0.006033858 * conc["EX_glu_L(e)"] + 
                    	1.878482743 * conc["EX_met_L(e)"] + 
			-2.838069276  * conc["EX_phe_L(e)"] + 
			0.040088058 * conc["EX_pyr(e)"] + 
			0.117735498 * conc["EX_succ(e)"] +
			-0.112714554 * conc["EX_thr_L(e)"])))


}


# Estimate new uptake rate for THREONINE
# using final concentrations calculated at the previous step

thr_L.rate <- function(conc) {

    return (as.numeric((1.485657670  + 0.001066872 * conc["EX_ala_L(e)"] + 
                    	0.005483605 * conc["EX_glc(e)"] + 
			0.030829724  * conc["EX_gly(e)"] + 
			-5.163714999  * conc["EX_his_L(e)"] + 
			-0.121899734 * conc["EX_lys_L(e)"] +
			-1.482683308  * conc["EX_met_L(e)"] +
			0.007359645 * conc["EX_pyr(e)"] +
			1.642848585 * conc["EX_phe_L(e)"] +
			0.008218513 * conc["EX_succ(e)"])))


}


dynamicConstraints_tnf <- function(model, 
                               concentrations=vector(), 
                               fluxes=vector(), 
                               stepNo=0) {
    
    # list of rate functions to evaluate 
    
    rxns <- c('ala_L', 'asp_L', 'glc', 'glu_L', 'gly', 'his_L', 'ile_L',
              'leu_L', 'lys_L', 'met_L', 'nh4', 'phe_L', 'pro_L', 'pyr', 
              'ser_L', 'succ', 'thr_L',  'tnf', 'tyr_L', 'val_L')

    for (i in rxns) {
        # compute new constraints based on previous concentrations
        # fn <- get(paste(i, "rate", sep='.')) ; rate <- fn(concentrations)
        rate <- get(paste(i, "rate", sep='.'))(concentrarions) )
        # set constraints to new values
        lo <- rate * 0.8
        hi <- rate * 1.2
        rxnName <- paste('EX_', i, '(e)', sep='')
        lowbnd(model)[react_id(model) == rxnName] <- lo
        uppbnd(model)[react_id(model) == rxnName] <- hi
    }
    
    # return new model
    return(model)
}
