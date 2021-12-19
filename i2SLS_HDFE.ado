program define i2SLS_HDFE, eclass
	syntax [anything] [if]  [in] [aweight pweight fweight iweight]  [, DELta(real 1) ABSorb(varlist) gmm2s Robust CLuster(varlist numeric)]
	marksample touse
	
	if "`gmm2s'" !="" {
		local opt0 = "`gmm2s' "
	}
	if  "`robust'" !="" {
		local opt1  = "`robust' "
	}
	if "`cluster'" !="" {
		local opt2 = "cluster(`cluster') "
	}
	local option = "`opt0'`opt1'`opt2'"
	
	local list_var `anything'
	* get depvar and indepvar
	gettoken depvar list_var : list_var
	gettoken indepvar list_var : list_var, p("(")
    * get endogenous variables and instruments
	gettoken endog list_var : list_var, bind
	gettoken endog endog : endog, p("(")
    gettoken endog instr_temp : endog , p("=")
    gettoken equalsign instr_temp : instr_temp , p("=")
	gettoken instr instr_temp : instr_temp, p(")")
	*di `"`indepvar'"'
	*di `"`indepvar'"'
	*di `"`endog'"'
	*di `"`instr'"'
	
		*** FWL Theorem Application
	cap drop Z0_*
	cap drop E0_*
	cap drop M0_*
	cap drop fe
	cap drop Y0_*
	cap drop xb_hat*
	quietly hdfe `indepvar'   [`weight'] , absorb(`absorb') generate(M0_)
	tempname DF_ADAPT
	scalar `DF_ADAPT' = e(df_a) //- e(N_hdfe)
	quietly hdfe `endog'  [`weight'] , absorb(`absorb') generate(E0_)
	quietly hdfe `instr'  [`weight'] , absorb(`absorb') generate(Z0_)
	*** Initialisation de la boucle
	tempvar fe
	quietly gen fe = 0
	tempvar y_tild  
	quietly gen `y_tild' = log(`depvar' + `delta'*exp(fe))
	quietly	hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_) 
	quietly ivreg2 Y0_  M0_*  (E0_* = Z0_*) [`weight'`exp'] if `touse' , `option' noconstant
	matrix beta_new = e(b)
	quietly predict xb_hat,xb
	local k = 0
	local eps = 1000	
	*** ItÃ©rations iOLS
	_dots 0
	while (`k' < 2000 & `eps' > 1e-5) {
		matrix beta_initial = beta_new
			* calcul de constante
		tempvar temp1
	quietly	gen `temp1' = `depvar'*exp(-xb_hat)
	quietly sum `temp1' if e(sample)
	
	tempname phi_hat
		scalar `phi_hat' = log(`r(mean)')
		cap drop `temp1'
		* Calcul de c_hat
		
		tempvar temp2
		quietly	gen `temp2' = log(`depvar' + `delta'*exp(`phi_hat' + xb_hat)) - xb_hat - `phi_hat'
		quietly sum `temp2' if e(sample)
		tempname c_hat
		scalar `c_hat' = `r(mean)'
		* Update d'un nouveau y_tild et regression avec le nouvel y_tild
		
		quietly replace `y_tild' = log(`depvar' + `delta'*exp(xb_hat)) - `c_hat'
		cap drop Y0_
	    quietly hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_) 
		* Update d'un nouveau y_tild et regression avec le nouvel y_tild
	quietly	 ivreg2 Y0_  M0_*  ( E0_* = Z0_*) [`weight'`exp'] if `touse', `option' noconstant 
			
	cap drop xb_hat_M xb_hat_N
	quietly predict xb_hat_M, xb		// predict with demeaned
	foreach var in `indepvar' { // rename variables for prediction
	quietly	rename M0_`var' TEMP`var'
	quietly	rename `var' M0_`var'
	}
	quietly predict xb_hat_N, xb // predict without demeaned
	foreach var in `indepvar' { // rename variables back to normal
	quietly	rename M0_`var' `var'
	quietly	rename TEMP`var' M0_`var'
	}
	quietly	replace fe = `y_tild'- Y0_ + xb_hat_M - xb_hat_N
	quietly	replace xb_hat = xb_hat_N + fe
		matrix beta_new = e(b)
		* DiffÃ©rence entre les anciens betas et les nouveaux betas
		matrix diff = beta_initial - beta_new
		mata : st_matrix("abs_diff", abs(st_matrix("diff")))
		mata : st_matrix("abs_diff2", st_matrix("abs_diff"):*st_matrix("abs_diff"))
		mata : st_matrix("criteria", rowsum(st_matrix("abs_diff2"))/cols(st_matrix("abs_diff2")))
		local eps = criteria[1,1]
		local k = `k'+1
		_dots `k' 0
	}

	*** Calcul de la bonne matrice de variance-covariance
	* Calcul du "bon" rÃ©sidu
	foreach var in `indepvar' {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename M0_`var' `var'
	}	
		foreach var in `instr'  {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename Z0_`var' `var'
	}
		foreach var in `endog'  {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename E0_`var' `var'
	}
quietly ivreg2 Y0_ `indepvar'  (`endog' = `instr') [`weight'`exp'] if `touse', `option' noconstant 
	foreach var in `indepvar' {      // rename variables back
	quietly	rename `var' M0_`var'
	quietly	rename TEMP_`var' `var'
	}
		foreach var in  `instr' {      // rename variables back
	quietly	rename `var' Z0_`var'
	quietly	rename TEMP_`var' `var'
	}
		foreach var in `endog'{      // rename variables back
	quietly	rename `var' E0_`var'
	quietly	rename TEMP_`var' `var'
	}
	preserve
	keep if e(sample)	
	local N_DF = e(Fdf2) -`DF_ADAPT' 
	tempvar ui
    gen `ui' = `depvar'*exp(- xb_hat - `phi_hat')
	matrix beta_final = e(b)
	quietly sum [`weight'`exp'] if e(sample)
	tempname nobs
	scalar `nobs' = r(N)
	* Calcul de Sigma_0, de I-W, et de Sigma_tild
	matrix Sigma = e(V)
	*tempname cste
	*gen `cste' = 1
	tempvar ui_bis
	gen `ui_bis' = 1 - `delta'/(`delta' + `ui')
	local var_list E0_* M0_*  // `cste'
	local instr_list  Z0_*  M0_* // `cste'
	mata : X=.
	mata : Z=.
	mata : IW=.
	mata : st_view(X,.,"`var_list'")
	mata : st_view(Z,.,"`instr_list'")
	mata : st_view(IW,.,"`ui_bis'")
	mata : W = (cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X))
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = W*Sigma_hat*W
	mata : invXpPzIWX = invsym(0.5*cross(X,Z)*invsym(cross(Z,Z))*(Z':*IW')*X+ 0.5*X'*(IW:*Z)*invsym(cross(Z,Z))*cross(Z,X))
	mata : Sigma_tild = invXpPzIWX*Sigma_0*invXpPzIWX
	
	/* Does not work with large sample 
	mata : IW = diag(IW) // preparation 
	mata : Pz = Z*invsym(Z'*Z)*Z'
	mata : Sigma_hat = st_matrix("Sigma") // retrouver la matrice de ivreg2
	mata : Sigma_0 = (X'*Pz*X)*Sigma_hat*(X'*Pz*X)
	mata : invXpPzIWX = invsym(0.5*X'*(Pz*IW+IW*Pz)*X)
	*/
	
	mata: st_matrix("Sigma_tild", Sigma_tild) // used in practice
   	matrix Sigma_tild = Sigma_tild*(`e(Fdf2))')/(`N_DF') // adjust DOF
	*** Stocker les rÃ©sultats dans une matrice
	local names : colnames e(b)
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
    ereturn post beta_final Sigma_tild , obs(`=r(N)') depname(`depvar') esample(`touse')  dof(`N_DF') 
	restore 
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn scalar widstat = e(widstat)
ereturn scalar arf = e(arf)
ereturn local cmd "i2SLS_HDFE"
ereturn local vcetype `option'
ereturn display

* drop 
	cap drop E0_*
	cap drop Z0_*
	cap drop M0_* 
	cap drop Y0_*
	cap drop xb_hat*
end
