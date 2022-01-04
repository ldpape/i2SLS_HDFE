* 15/12/2021 : corrected "cross" in S.E. which is complicated by the symmetrization
* 15/12/2021 : corrected iteration logical specification
* 16/12/2021 : corrected absorb from varlist to string, as in ppmlhdfe
* 22/12/2021 : coded with matrix multiplication instead of pre-canned program
* 22/12/2021 : added convergence control (limit and maximum)
* 04/01/2022 : added constant + checks for convergence + corrected problem with collinear variables affecting final 2SLS
cap program drop i2SLS_HDFE
program define i2SLS_HDFE, eclass
	syntax [anything] [if] [in] [aweight pweight fweight iweight]  [, DELta(real 1) ABSorb(varlist) LIMit(real 0.00001) MAXimum(real 1000) Robust CLuster(varlist numeric) ]
		marksample touse
	preserve 
	quietly keep if `touse'
	if  "`robust'" !="" {
		local opt1  = "`robust' "
	}
	if "`cluster'" !="" {
		local opt2 = "cluster(`cluster') "
	}
	local option = "`opt1'`opt2'"
	*** Obtain lists of variables 
	local list_var `anything'
	gettoken depvar list_var : list_var
	if (strpos("`list_var'","(")==2){  // anormal case : no X, only Z 
	local list_var: subinstr local list_var "(" "", all
	local indepvar 
	gettoken endog list_var : list_var, bind
	gettoken endog instr_temp : endog , p("=")
	local list_var: subinstr local list_var "=" "", all
	gettoken instr list_var : list_var, bind
	}
	else{  // normal case : X exists 
	gettoken indepvar list_var : list_var, p("(")
	gettoken endog list_var : list_var, bind
	gettoken endog endog : endog, p("(")
	gettoken endog instr_temp : endog , p("=")
	gettoken equalsign instr_temp : instr_temp , p("=")
	gettoken instr instr_temp : instr_temp, p(")")
	}		
	*** Initialisation de la boucle
	tempvar cste
	gen `cste' = 1
	** drop collinear variables
    _rmcoll `indepvar', forcedrop 
	local var_list `endog' `r(varlist)' 
	local instr_list `instr' `r(varlist)'
	local alt_varlist `r(varlist)'
	*** FWL Theorem Application
	cap drop Z0_*
	cap drop E0_*
	cap drop M0_*
	cap drop Y0_*
	cap drop xb_hat*
	if "`indepvar'"=="" { // case with no X , only FE 
	quietly hdfe `endog' [`weight'] , absorb(`absorb') generate(E0_)
	tempname DF_ADAPT
	scalar `DF_ADAPT' = e(df_a) //- e(N_hdfe)
	quietly hdfe `instr' [`weight'] , absorb(`absorb') generate(Z0_)
	tempvar y_tild  
	quietly gen `y_tild' = log(`depvar' + 1)
	quietly	hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_) 
	mata : X=.
	mata : PX=.
	mata : PZ=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`endog'")
	mata : st_view(PX,.,"E0_*")
	mata : st_view(PZ,.,"Z0_*")
	mata : st_view(y_tilde,.,"`y_tild'")
	mata : st_view(Py_tilde,.,"Y0_")
	mata : st_view(y,.,"`depvar'")	
	}
	else { // standard case with both X and FE
	quietly hdfe `r(varlist)' [`weight'] , absorb(`absorb') generate(M0_)
	tempname DF_ADAPT
	scalar `DF_ADAPT' = e(df_a) //- e(N_hdfe)
	quietly hdfe `endog'  [`weight'] , absorb(`absorb') generate(E0_)
	quietly hdfe `instr'  [`weight'] , absorb(`absorb') generate(Z0_)
	tempvar y_tild  
	quietly gen `y_tild' = log(`depvar' + 1)
	quietly	hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_) 
	mata : X=.
	mata : PX=.
	mata : PZ=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'")
	mata : st_view(PX,.,"E0_* M0_*")
	mata : st_view(PZ,.,"Z0_* M0_*")
	mata : st_view(y_tilde,.,"`y_tild'")
	mata : st_view(Py_tilde,.,"Y0_")
	mata : st_view(y,.,"`depvar'")	
	}
	* prepare  future inversions 
	mata : invPzX = invsym(cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX))*cross(PX,PZ)*invsym(cross(PZ,PZ))
	mata : beta_initial = invPzX*cross(PZ,Py_tilde)
	mata : beta_t_1 = beta_initial // needed to initialize
	mata : beta_t_2 = beta_initial // needed to initialize
	mata : q_hat_m0 = 0
	local k = 1
	local eps = 1000	
	mata: q_hat = J(`maximum', 1, .)
	*** Iterations iOLS
	_dots 0
	while ((`k' < `maximum') & (`eps' > `limit' )) {
	mata: xb_hat_M = PX*beta_initial 
	mata: xb_hat_N = X*beta_initial
	mata: fe = y_tilde - Py_tilde + xb_hat_M - xb_hat_N
	mata: xb_hat = xb_hat_N + fe		
		* update du alpha
	mata: alpha = log(mean(y:*exp(-xb_hat)))
	mata: y_tilde = log(y + `delta'*exp(xb_hat :+ alpha )) :-mean(log(y + `delta'*exp(xb_hat :+ alpha)) -xb_hat :- alpha  )
		* regression avec le nouvel y_tild
	cap drop `y_tild' 
	quietly mata: st_addvar("double", "`y_tild'")
	mata: st_store(.,"`y_tild'",y_tilde)
	cap drop Y0_
    quietly hdfe `y_tild' [`weight'] , absorb(`absorb') generate(Y0_)
	mata : st_view(Py_tilde,.,"Y0_")
		* 2SLS 
	mata : beta_new = invPzX*cross(PZ,Py_tilde)
	mata: criteria = mean(abs(beta_initial - beta_new):^(2))
	mata: st_numscalar("eps", criteria)
	mata: st_local("eps", strofreal(criteria))
		* safeguard for convergence.
	if `k'==`maximum'{
		  di "There has been no convergence so far: increase the number of iterations."  
	}
	if `k'>4{
	mata: q_hat[`k',1] = mean(log( abs(beta_new-beta_initial):/abs(beta_initial-beta_t_2)):/log(abs(beta_initial-beta_t_2):/abs(beta_t_2-beta_t_3)))	
	mata: check_3 = abs(mean(q_hat)-1)
		if mod(`k'-4,50)==0{
    mata: q_hat_m =  mm_median(q_hat[((`k'-49)..`k'),.] ,1)
	mata: check_1 = abs(q_hat_m - q_hat_m0)
	mata: check_2 = abs(q_hat_m-1)
	mata: st_numscalar("check_1", check_1)
	mata: st_local("check_1", strofreal(check_1))
	mata: st_numscalar("check_2", check_2)
	mata: st_local("check_2", strofreal(check_2))
	mata: st_numscalar("check_3", check_3)
	mata: st_local("check_3", strofreal(check_3))
	mata: q_hat_m0 = q_hat_m
		if ((`check_1'<1e-4)&(`check_2'>1e-2)) {
di "delta is too small to achieve convergence -- update to larger value"
	local k = `maximum'
		}
		if ((`check_3'>0.5) & (`k'>500)) {
	local k = `maximum'
di "q_hat too far from 1"
		}
					  }
	}
	if `k'>2 { // keep in memory the previous beta_hat for q_hat 
	mata:   beta_t_3 = beta_t_2
	mata:   beta_t_2 = beta_initial
	}
	mata: beta_initial = beta_new
	local k = `k'+1
	_dots `k' 0
	}
	*** Calcul de la bonne matrice de variance-covariance
	* Calcul du "bon" residu
	mata: ui = y:*exp(-xb_hat)
	mata: ui = ui:/(`delta' :+ ui)
	* Final 2SLS with ivreg2 
		foreach var in `alt_varlist' {     // rename variables for last ols
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
cap _crcslbl Y0_ `depvar' // label Y0 correctly
quietly ivreg2 Y0_ `alt_varlist' (`endog' = `instr') [`weight'`exp'] , `option' noconstant   // standard case with X and FE 
if "`indepvar'"=="" {
quietly ivreg2 Y0_ `alt_varlist' (`endog' = `instr') [`weight'`exp'] , `option' noconstant   // case with no X , only FE 
}
	foreach var in `alt_varlist' {      // rename variables back
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
	local N_DF = e(Fdf2) -`DF_ADAPT' 
	* Calcul de Sigma_0, de I-W, et de Sigma_tild
	matrix beta_final = e(b) // 	mata: st_matrix("beta_final", beta_new)
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX))*Sigma_hat*(cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX)) // recover original HAC 
	mata : invXpPzIWX = invsym(0.5*cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,ui,PX)+ 0.5*cross(PX,ui,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX))
	mata : Sigma_tild = invXpPzIWX*Sigma_0*invXpPzIWX
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
