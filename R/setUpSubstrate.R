### JR ### should be sourced
setUpSubstrate <- function(useGLC, useMNT, useEXP) {
    #
    # Initial concentrations
    # ----------------------
    # initial concentrations are in millimoles per litre
    #	i.e the substrate concentration is scaled to define the amount
    #	of substrate available per unit of biomass per unit of time.
    #	Varma and Palsson:
    #	For E.coli, the maximum O2 utilization rate was 15 mmol of O2 per
    #	g (dry weight) per h.
    #	For E.coli aerobic cultures, the ratio of the growth rate to the
    #	biomass yield was 10.5 mmol of Glc per g (dry weight) per h
    #	For E. coli anaerobic cultures the maximum glucose utilization
    #	rate was determined to be 18.5 mmol per g (gry weight) per h
    # 
    # Some examples are:
    # 
    Mnt =	     c('EX_mnl(e)');
    Mnt_mg_per_l = c( 10000 );
    Mnt_mwt = c( 182.172  );
    Mnt_mmol_per_l = Mnt_mg_per_l / Mnt_mwt;

    Glc =    c('EX_glc(e)');
    Glc_mg_per_l = c( 10000);
    Glc_mwt = c( 180.15588 );
    Glc_mmol_per_l = Glc_mg_per_l / Glc_mwt;

    #GlcMnt = 	c('EX_glc(e)', 'EX_mnl(e)');
    #GlcMnt_g_per_l = c( 5,     5 );
    #GlcMnt_mwt = c( 180.15588,	182.172	 );
    GlcMnt = 	c('EX_glc(e)', 'EX_mnl(e)');
    GlcMnt_mg_per_l = c( 10000,    10000 );
    GlcMnt_mwt = c( 180.15588,	182.172	 );
    GlcMnt_mmol_per_l = GlcMnt_mg_per_l / GlcMnt_mwt;

    # SIMULATE CULTURE MEDIUM
    # =======================
    #
    # The culture medium used contains
    # Major components
    #	D-Glucose	1%	10   g/l	 55.5075 mmol/l 
    #	D-Mannitol	1%	10   g/l	 54.8931 mmol/l 
    #	NH4 2(SO4)		 2   g/l	  9.5163 mmol/l
    #	Bacto casamino acids	 5   g/l	
    #	Mg SO4 7(H2O)		 0.6 g/l	  1.7516 mmol/l	
    #	NaH2PO4			 6.8 g/l	 56.6775 mmol/l
    #	K2HPO4			11.4 g/l	 65.4511 mmol/l
    #
    #
    # Chemical formulas and molecular weights are
    #
    # GLC	D-Glucose	C6H12O6		180.1559 g/mol
    # MNT	D-Mannitol	C6H14O6		182.172  g/mol
    #	NH4		NH4		 18.04   g/mol
    #	CO2		CO2		 44.01   g/mol
    #	O2		O2		 31.9988 g/mol
    #	SO4		SO4		 96.0626 g/mol 
    #	PO4		PO4		 94.9714 g/mol
    #	MG (atom)	MG		 24.305  g/mol
    #	Na (atom)	Na		 22.9898 g/mol
    #	K (atom)	K		 39.0983 g/mol
    #	Ca (atom)	Ca		 40.078  g/mol
    #	Cl (atom)	Cl		 35.453  g/mol
    #
    #	NH4(SO4)2			210.1637 g/mol
    #	Mg(SO4)2 (H2O)7			342.5372 g/mol
    #	Na H2 PO4			119.9770 g/mol
    #	K2 H PO4			174.1759 g/mol
    #
    # Aminoacid	Chemical formula	Molecular weight(g/mol)	#*5g/l->mmol/l
    # ALA	Alanine		C3H7NO2		89.0935			
    # ARG	Arginine	C6H14N4O2	174.2017
    # ASN	Asparagine	C4H8N2O3	132.1184
    # ASP	Aspartate	C4H7NO4		133.1032
    # CYS	Cysteine	C3H7NO2S	121.1590
    # GLU	Glutamate	C5H9NO4		147.1299
    # GLN	Glutamine	C5H10N2O3	146.1451
    # GLY	Glycine		C2H5NO2		75.0669
    # HIS	Histidine	C6H9N3O2	155.1552
    # ILE	Isoleucine	C6H13NO2	131.1736
    # LEU	Leucine		C6H13NO2	131.1736
    # LYS	Lysine		C6H14N2O2	146.1882
    # MET	Methionine	C5H11NO2S	149.2124
    # PHE	Phenylalanine	C9H11NO2	165.1900
    # PRO	Proline		C5H9NO2		115.1310
    # SER	Serine		C3H7NO3		105.0930
    # THR	Threonine	C4H9NO3		119.1197
    # TRP	Tryptophan	C11H12N2O2	204.2262
    # TYR	Tyrosine	C9H11NO3	181.1894
    # VAL	Valine		C5H11NO2	117.1469
    #	
    #	To convert, use rule (i.e. divide weight by molar mass [mol. weight])
    #
    #	    grams   molar mass in grams
    #       ----- = -------------------
    #	    moles       one mole
    #
    #    Hence
    #	[NH3] = 2 * [(NH4)2 SO4] = 2 * 15.1355 = 30.2710
    #	[Pi] = [NaH2PO4] + [K2HPO4] = 56.6775 + 65.4511 = 122.1286
    #	[SO4] = [(NH4)2 SO4] + [Mg SO4 7(H2O)] = 15.1335 + 1.7516 = 16.8851
    #	[K] = 2 * [K2HPO4] = 2 * 65.4511 = 130.9022
    #
    major_components_minus_aa = c(
	    'EX_nh4(e)', 'EX_na1(e)', 'EX_pi(e)', 
            'EX_so4(e)', 'EX_k(e)',   'EX_mg2(e)'
    );

    concentration_of_major_components_minus_aa = c(
	    30.2710,  56.6775,  122.1286,
            16.8851,    130.9022,  1.7516
    );

    #
    # K, MG are not found in the metabolites or reaction lists, so we'll remove them
    #
    major_components_minus_aa = c(
	    'EX_nh4(e)', 'EX_na1(e)', 'EX_pi(e)', 
            'EX_so4(e)'
    );

    concentration_of_major_components_minus_aa = c(
	    30.2710,  56.6775,  122.1286, 
            16.8851
    );

    if (useGLC == 1) {
	    major_components_minus_aa = c(
		    Glc, major_components_minus_aa)
	    concentration_of_major_components_minus_aa = c(
		    Glc_mmol_per_l, concentration_of_major_components_minus_aa)
    }

    if (useMNT == 1) {
	    major_components_minus_aa = c(
		    Mnt, major_components_minus_aa)
	    concentration_of_major_components_minus_aa = c(
		    Mnt_mmol_per_l, concentration_of_major_components_minus_aa)
    }

    # minor components
    #	ZnSO4 7(H2O)		 1  mg/l
    #	FeSO4 7(H2O)		 1  mg/l
    #	MnCl2 4(H2O)		 1  mg/l
    #	CaCl2			 1  mg/l
    #
    # These concentrations (but for the carbon source and, maybe some a.a. if
    # protein overproduction is very large) should not be rate-limiting (i.e.
    # there should be enough of each component for all the cultivation period
    # to occur without signs of depletion). We can leave, thus, most of them
    # unspecified and assume they never run out.
    #
    #	Note also that casamino acids may contain up to 50-60% salts whose 
    # contribution has not been accounted for. We'll generally consider
    # that these componentns are in excess and do not need to be specified.
    # When not specified, the simulation system will provide them in excess.
    #
    # They can be listed as excluded uptake reactions and we may ignore their
    # actual concentrations
    #
    minor_components_in_excess = c('EX_zn2(e)', 'EX_fe2(e)', 'EX_mn2(e)', 'EX_cl(e)', 'EX_ca2(e)');
    #
    minor_components_in_excess = c()

    # casamino acids
    # --------------
    #
    # 1. using composition from www.myopure.com
    # -----------------------------------------
    #
    # aa	name		formula		molar mass	#HC%	in 5g	g -> mmol
    # ALA	Alanine		C3H7NO2		89,0935		2,9	0,1450	1,6275
    # ARG	Arginine	C6H14N4O2	174,2017	3,4	0,1700	0,9759
    # ASN	Asparagine	C4H8N2O3	132,1184	0	0,0000	0,0000
    # ASP	Aspartate	C4H7NO4		133,1032	6,4	0,3200	2,4041
    # CYS	Cysteine	C3H7NO2S	121,1590	0,3	0,0150	0,1238
    # GLU	Glutamate	C5H9NO4		147,1299	21,1	1,0550	7,1705
    # GLN	Glutamine	C5H10N2O3	146,1451	0	0,0000	0,0000
    # GLY	Glycine		C2H5NO2		75,0669		1,7	0,0850	1,1323
    # HIS	Histidine	C6H9N3O2	155,1552	3	0,1500	0,9668
    # ILE	Isoleucine	C6H13NO2	131,1736	4,4	0,2200	1,6772
    # LEU	Leucine		C6H13NO2	131,1736	9	0,4500	3,4306
    # LYS	Lysine		C6H14N2O2	146,1882	7,5	0,3750	2,5652
    # MET	Methionine	C5H11NO2S	149,2124	2,7	0,1350	0,9048
    # PHE	Phenylalanine	C9H11NO2	165,1900	5,3	0,2650	1,6042
    # PRO	Proline		C5H9NO2		115,1310	10,4	0,5200	4,5166
    # SER	Serine		C3H7NO3		105,0930	5,1	0,2550	2,4264
    # THR	Threonine	C4H9NO3		119,1197	4	0,2000	1,6790
    # TRP	Tryptophan	C11H12N2O2	204,2262	1,4	0,0700	0,3428
    # TYR	Tyrosine	C9H11NO3	181,1894	5,2	0,9615	5,3068
    # VAL	Valine		C5H11NO2	117,1469	6,3	0,7937	6,7748
    # 
    # giving
    #	ALA	ARG	ASN	ASP	CYS	GLU	GLN	GLY	HIS	ILE	LEU	LYS	MET	PHE	PRO	SER	THR	TRP	TYR	VAL
    #	1,6275	0,9759	0,0000	2,4041	0,1238	7,1705	0,0000	1,1323	0,9668	1,6772	3,4306	2,5652	0,9048	1,6042	4,5166	2,4264	1,6790	0,3428	5,3068	6,7748
    #
    #pct =,	c( 2.9,	3.4,	0,	6.4,	0.3,	21.1,	0,	1.7,	3,	4.4,	9,	7.5,	2.7,	5.3,	10.4,	5.1,	4,	1.4,	5.2,	6.3 )
    #aa_mwt = c(89.0935,	174.2017,	132.1184,	133.1032,	121.1590,	147.1299,	146.1451,	75.0669,	155.1552,	131.1736,	131.1736,	146.1882,	149.2124,	165.1900,	115.1310,	105.0930,	119.1197,	204.2262,	181.1894,	117.1469);
    #aa_mmolg = pct / aa_mwt;
    #aa_mmol5g = c( 1.6275,	0.9759,	0.0000,	2.4041,	0.1238,	7.1705,	0.0000,	1.1323,	0.9668,	1.6772,	3.4306,	2.5652,	0.9048,	1.6042,	4.5166,	2.4264,	1.6790,	0.3428,	5.3068,	6.7748 );
    #
    #substrateRxns = c(
    #'EX_glc(e)',		'EX_mnl(e)',	'EX_pi(e)', 	'EX_AMLB', 
    #'EX_ala_L(e)',		'EX_arg_L(e)',	'EX_asn_L(e)',	'EX_asp_L(e)', 
    #'EX_cys_L(e)',		'EX_glu_L(e)',	'EX_gln_L(e)',	'EX_gly(e)', 
    #'EX_his_L(e)',		'EX_ile_L(e)',	'EX_leu_L(e)',	'EX_lys_L(e)', 
    #'EX_met_L(e)',		'EX_phe_L(e)',	'EX_pro_L(e)',	'EX_ser_L(e)',  
    #'EX_thr_L(e)',		'EX_trp_L(e)',	'EX_tyr_L(e)',	'EX_val_L(e)' 
    #);
    #initConcentrations = c(
    #30.,	  		12.,		1.,		1e-5,	    
    #1.6275,	       0.9759,		1e-5,		2.4041, 
    #0.1238,	       7.1705,		1e-5,		1.1323, 
    #0.9668,	       1.6772,		3.4306,		2.5652, 
    #0.9048,	       1.6042,		4.5166,		2.4264, 
    #1.6790,	       0.3428,		5.3068,		6.7748 
    #);
    #
    #
    #
    # 2. using composition from D'Huys 2011
    # -------------------------------------
    #
    #	a.a.	mg/g
    #	ALA	17.8
    #	ARG	n.d.
    #	ASN	0
    #	ASP	42.56
    # 	CYS	n.d.
    # 	GLU	102.9
    # 	GLN	0
    # 	GLY	12
    # 	HIS	12.4
    # 	ILE	14.41
    # 	LEU	39.3
    # 	LYS	23.36
    # 	MET	11.92
    # 	PHE	19.8
    # 	PRO	46
    # 	SER	20
    # 	THR	10.7
    # 	TRP	0
    # 	TYR	14.48
    # 	VAL	31.8
    #
    # giving
    #
    # Casamino acids using composition from D'Huys 2011
    # -------------------------------------------------
    #
    aa =  c('ALA',	'ARG',	'ASN',	'ASP',	
	    'CYS',	'GLU',	'GLN',	'GLY',	
            'HIS',	'ILE',	'LEU',	'LYS',	
            'MET',	'PHE',	'PRO',	'SER',	
            'THR',	'TRP',	'TYR',	'VAL');
    exchange_aa =  c('EX_ala_L(e)',	'EX_arg_L(e)',	'EX_asn_L(e)',	'EX_asp_L(e)',
		    'EX_cys_L(e)',	'EX_glu_L(e)',	'EX_gln_L(e)',	'EX_gly(e)',
                    'EX_his_L(e)',	'EX_ile_L(e)',	'EX_leu_L(e)',	'EX_lys_L(e)',	
                    'EX_met_L(e)',	'EX_phe_L(e)',	'EX_pro_L(e)',	'EX_ser_L(e)',	
                    'EX_thr_L(e)',	'EX_trp_L(e)',	'EX_tyr_L(e)',	'EX_val_L(e)');

    #mg_per_g = c(17.8,	-0,	0,	42.56,	-0,	102.9,	0,	12,	12.4,	14.41,	39.3,	23.36,	11.92,	19.8,	46,	20,	10.7,	0,	14.48,	31.8);
    #
    #	NOTE: we'll use 1e-04 for 0 values so as to be able to specify and
    # follow them
    mg_per_g = c(17.8,	15.18,	1e-4,	42.56,	
		 1.32,	102.9,	1e-4,	12,	
        	 12.4,	14.41,	39.3,	23.36,	
        	 11.92,	19.8,	46,	20,	
        	 10.7,	1e-4,	14.48,	31.8);

    aa_mwt =   c(89.0935,	174.2017,	132.1184,	133.1032,	
		121.1590,	147.1299,	146.1451,	75.0669,
        	155.1552,	131.1736,	131.1736,	146.1882,	
        	149.2124,	165.1900,	115.1310,	105.0930,	
        	119.1197,	204.2262,	181.1894,	117.1469)

    aa_mmol_per_g_l = mg_per_g / aa_mwt; 	# .* 1 (per 1 litre)
    aa_mmol_per_5g_l = aa_mmol_per_g_l * 5;
    aa_mmol_per_15g_l = aa_mmol_per_g_l * 15;


    # We'll set these exchange reactions to zero
    zeroRxns = c('EX_2hymeph(e)', 'EX_acald(e)', 'EX_nformanth(e)', 'EX_glyc(e)', 
		 'EX_gua(e)',  'EX_hxan(e)', 'EX_na1(e)C', 'EX_ade(e)', 'EX_urea(e)');
    zeroConc = c(1.e-5,	1.e-5,	  1.e-5,    1.e-5,   1.e-5,	
        	 1.e-5,     1.e-5,	  1.e-5,    1.e-5,   1.e-5);

    # Other exchange reactions
    # ------------------------
    # we will include them to see if we can plot them as well.
    otherExcRxns = c('EX_co2(e)', 'EX_nh4(e)', 'EX_o2(e)') 

    #
    # For S. coelicolor (Kim et al.)
    # - - - - - - - - - - - - - - - -
    #substrateRxns = c('EX_glc(e)', 'EX_mnl(e)');
    #initConcentrations = c(10, 1.e-4);		# mmol/L
    #initConcentrations = c(1.e-4, 10);		# mmol/L
    #initConcentrations = c(10, 10);		# mmol/L
    #initConcentrations = c(.01, .27);		# mmol/L
    #exclUptakeRxns = c('EX_o2(e)');
    #
    # Mannitol
    # - - - - -
    if (useMNT == 1) {
	# Only mannitol
	#substrateRxns = c( Mnt );	# 5g GLC 5g MNT
	#initConcentrations = Mnt_mmol_per_l;

	# Mannitol + casamino acids
	substrateRxns = c( Mnt, exchange_aa );
	initConcentrations = c( Mnt_mmol_per_l, aa_mmol_per_5g_l );
    }

    # Glucose
    # - - - - -
    if (useGLC == 1) {
	# Only glucose
	#substrateRxns = c( Glc );	# 5g GLC 5g MNT
	#initConcentrations = Glc_mmol_per_l;

	# Glucose + casamino acids
	substrateRxns = c( Glc, exchange_aa );
	initConcentrations = c( Glc_mmol_per_l, aa_mmol_per_5g_l );
    }

    # Glucose/Mannitol 
    # - - - - - - - - -
    if (useGLC == 1 && useMNT == 1) {
	# only glucose + mannitol
	#substrateRxns = c( GlcMnt );	# 5g GLC 5g MNT
	#initConcentrations = GlcMnt_mmol_per_l;

	#substrateRxns = c( GlcMnt, zeroRxns );	# 5g GLC 5g MNT
	#initConcentrations = ( 1.e-4, GlcMnt_mmol_per_l(2), zeroConc );

	# Glucose + Mannitol + casamino acids
	substrateRxns = c( GlcMnt, exchange_aa );
	initConcentrations = c( GlcMnt_mmol_per_l, aa_mmol_per_5g_l );

	#substrateRxns = c( GlcMnt, exchange_aa zeroRxns );	# 5g GLC 5g MNT
	#initConcentrations = c( 1.e-4, GlcMnt_mmol_per_l(2), aa_mmol_per_5g_l zeroConc);
    }


    if (useEXP == 1) {
	# For our environmental set up
	# Note that this relies on useGLC and useMNT having set the
	# composition of the medium earlier.
	substrateRxns = c( major_components_minus_aa, 
			   exchange_aa );
	initConcentrations = c( concentration_of_major_components_minus_aa,
				aa_mmol_per_5g_l);

    }


    # Excluded uptake reactions will have no concentration limits
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    #exclUptakeRxns=c('EX_o2(e)')
    #
    #
    #exclUptakeRxns = c(minor_components_in_excess, 
    #		    'EX_pi(e)','EX_na1(e)','EX_nh4(e)', 'EX_so4(e)');
    #
    #exclUptakeRxns = c(major_components_in_excess,
    #		minor_components_in_excess, 'EX_pi(e)','EX_na1(e)','EX_nh4(e)');
    #
    #exclUptakeRxns = c(minor_components_in_excess,
    #			'EX_co2(e)', 'EX_o2(e)', 'EX_h2o(e)', 'EX_h(e)');
    # WE CANNOT EXCLUDE THESE REACTIONS IF WE WANT TO TRACK THEM!!!
    #
    exclUptakeRxns=c('EX_nh4(e)', 'EX_so4(e)', 
    		'EX_co2(e)', 'EX_o2(e)', 'EX_h2o(e)', 'EX_h(e)');
    # hence, we'll use instead only
    exclUptakeRxns=c('EX_h2o(e)', 'EX_h(e)')

    medium <- list(
    	substrateRxns=substrateRxns, 
	initConcentrations=initConcentrations,
    	exclUptakeRxns=exclUptakeRxns)

    return (medium)
}

