{smcl}
{* *! version 1.0 22march2021}{...}
{vieweralsosee "[R] poisson" "help poisson"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "ivreghdfe" "help reghdfe"}{...}
{vieweralsosee "ppml" "help ppml"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{viewerjumpto "Syntax" "i2SLS_HDFE##syntax"}{...}
{viewerjumpto "Description" "i2SLS_HDFE##description"}{...}
{viewerjumpto "Citation" "i2SLS_HDFE##citation"}{...}
{viewerjumpto "Authors" "i2SLS_HDFE##contact"}{...}
{viewerjumpto "Examples" "i2SLS_HDFE##examples"}{...}
{viewerjumpto "Description" "i2SLS_HDFE##Testing"}{...}
{viewerjumpto "Stored results" "i2SLS_HDFE##results"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:i2SLS_HDFE} {hline 2}}Iterated Two Stage Least Squares (i2SLS) with High Dimensional Fixed Effect with delta {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:i2SLS_HDFE}
{depvar} [{indepvars}]
{ifin} {it:{weight}} {cmd:,} [{help i2SLS_HDFE##options:options}] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab: Standard Errors: Classical/Robust/Clustered}
{synopt:{opt vce}{cmd:(}{help i2SLS_HDFE##opt_vce:vcetype}{cmd:)}}{it:vcetype}
may be classical (assuming homoskedasticity), {opt r:obust}, or {opt cl:uster} (allowing two- and multi-way clustering){p_end}
{syntab: Delta}
{synopt:{opt delta}{cmd:(}{help i2SLS_HDFE##opt_vce:delta}{cmd:)}}{it:delta} is any strictly positive constant. {p_end}
{syntab: Absorb}
{synopt:{opt absorb}{cmd:(}{help i2SLS_HDFE##opt_vce:absorb}{cmd:)}}{it:absorb} takes any categorical or binary variable. You can generate a combination of fixed effects with egen var = group(var1 var2). {p_end}

{marker description}{...}
{title:Description}

{pstd}{cmd:i2SLS_HDFE} iterated Two Stage Least Squares (i2SLS) with High Dimensional Fixed effects and delta, as described by {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3444996":Bellego, Benatia, and Pape (2021)}.

{pstd}This package:

{pmore} 1. relies on Stata's ivreg2 routine for estimation. {p_end}

{pmore} 2. assumes the iOLS exogeneity condition with instrument Z, E(Z'log(delta+U))=constant.  {p_end}

{pmore} 3. allows for high dimensional fixed effects through the "absorb(.)" option.  {p_end}


{title:Background}

{pstd} {cmd: i2SLS_HDFE} estimates iOLS_delta, a solution to the problem of the log of zero, in the context of an endogenous regressor with high dimensional fixed effects.  This method relies on running the "ivreg2" function iteratively.
and requires the "HDFE" package to also be installed. This provides the reader with the final 2SLS estimates and allows the use the post-estimation commands available under ivreg2 (using Y_tilde = log(Y + delta*exp(xb+fe))) as a 
dependent variable.  

{synoptset 22}{...}
{synopthdr: variables}
{synoptline}
{synopt:{it:depvar}} Dependent variable{p_end}
{synopt:{it:indepvars}} List of explanatory variables {p_end}
{synoptline}
{p2colreset}{...}


{marker caveats}{...}
{title:Caveats}

{pstd} Convergence is decided based on coefficients (sum of squared coefficients < 1e-15) and not on the modulus of the contraction mapping. {p_end}

{pstd} Standard errors are modified by this procedure (and will differ from those in i2SLS_ivreg2) but are nonetheless still asymptotically consistent. {p_end}

{pstd} The {help reg postestimation##predict:predict}, {help test}, and {help margins} postestimation commands are available after {cmd:i2SLS_HDFE}. {p_end}

{marker contact}{...}
{title:Authors}

{pstd} Christophe Bellego, David Benatia, Louis Pape {break}
CREST - ENSAE - HEC Montréal - Ecole Polytechnique {break}
Contact: {browse "mailto:louis.pape@polytechnique.edu":louis.pape@polytechnique.edu}
{p_end}

{marker citation}{...}
{title:Citation}

{pstd}
Bellégo Christophe, Benatia David, and Pape Louis-Daniel, Dealing with Logs and Zeros in Regression Models (2019).
Série des Documents de Travail n° 2019-13.
Available at SSRN: https://ssrn.com/abstract=3444996

{marker examples}{...}
{title:Examples}

{pstd}First, we compare i2SLS_HDFE with IV-Poisson (using ivpois)
{browse "http://www.haghish.com/statistics/stata-blog/stata-programming/download/ivpois.html":ivpois help file}.
{p_end}
{hline}
{phang2}{cmd:. use "http://www.stata-press.com/data/r14/airline"}{p_end}
{phang2}{cmd:. i2SLS_HDFE injuries (XYZowned=n) , delta(1) vce(robust)}{p_end}
{phang2}{cmd:. ivpois injuries (XYZowned=n)}{p_end}
{hline}

{pstd} Second, we show how to test for the pattern of zeros with i2SLS. We rely on data on households' trips away from home, as used in {browse "https://www.stata.com/manuals/rivpoisson.pdf":ivpoisson manual}.
to study the effect of cost of transportation (tcost). 
{p_end}
{hline}
{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. webuse trip }{p_end}
{phang2}{cmd:. gen outside = trips>0 }{p_end}

{phang2}{cmd:. i2SLS_HDFE trips cbd ptn worker weekend (tcost=pt) , delta(1) robust }{p_end}
{phang2}{cmd:. i2SLS_HDFE trips cbd ptn worker weekend (tcost=pt) , delta(1000) robust  }{p_end}
{phang2}{cmd:. ivpois trips cbd ptn worker weekend endog(tcost) exog(pt) }{p_end}
{phang2}{cmd:. ivpoisson gmm trips cbd ptn worker weekend (tcost=pt), multiplicative}{p_end}

{phang2}{cmd:. cap program drop i2SLS_bootstrap  }{p_end}
{phang2}{cmd:. program i2SLS_bootstrap, rclass  }{p_end}
{phang2}{cmd:. i2SLS_HDFE trips cbd ptn worker weekend (tcost=pt) , delta(1) robust  }{p_end}
{phang2}{cmd:. scalar delta = 1  }{p_end}
{phang2}{cmd:. *lhs of test  }{p_end}
{phang2}{cmd:. predict xb_temp, xb  }{p_end}
{phang2}{cmd:. gen u_hat_temp = trips*exp(-xb_temp)  }{p_end}
{phang2}{cmd:. gen lhs_temp = log(delta+u_hat_temp) - log(delta)  }{p_end}
{phang2}{cmd:. * rhs of test  }{p_end}
{phang2}{cmd:. gen temp = log(trips + delta*exp(xb_temp)) - xb_temp  }{p_end}
{phang2}{cmd:. egen c_hat_temp = mean(temp)   }{p_end}
{phang2}{cmd:. logit outside cbd ptn worker weekend (tcost=pt) }{p_end}
{phang2}{cmd:. predict p_hat_temp, pr  }{p_end}
{phang2}{cmd:. gen rhs_temp = (c_hat_temp-log(delta))/p_hat_temp  }{p_end}
{phang2}{cmd:. * run the test  }{p_end}
{phang2}{cmd:. reg lhs_temp rhs_temp if outside, nocons   }{p_end}
{phang2}{cmd:. matrix b = e(b)  }{p_end}
{phang2}{cmd:. ereturn post b  }{p_end}
{phang2}{cmd:. * drop created variables  }{p_end}
{phang2}{cmd:. cap drop *temp  }{p_end}
{phang2}{cmd:. end  }{p_end}

{phang2}{cmd:. bootstrap lambda = _b[rhs_temp] , reps(50): i2SLS_bootstrap  }{p_end}
{phang2}{cmd:. test lambda==1  }{p_end}
{hline}

{pstd} Third, we show how to test for the pattern of zeros with IV-Poisson.
{p_end}
{hline}

{phang2}{cmd:.  cap program drop IVPoisson_boostrap  }{p_end}
{phang2}{cmd:.  program IVPoisson_boostrap, rclass  }{p_end}
{phang2}{cmd:.  * estimate the model  }{p_end}
{phang2}{cmd:.  ivpois trips cbd ptn worker weekend endog(tcost) exog(pt)  }{p_end}
{phang2}{cmd:.  * lhs of test   }{p_end}
{phang2}{cmd:.  predict xb_temp, xb  }{p_end}
{phang2}{cmd:.  gen u_hat_temp = trips*exp(-xb_temp)  }{p_end}
{phang2}{cmd:.  egen mean_u_temp = mean(u_hat_temp)  }{p_end}
{phang2}{cmd:. gen lhs_temp = u_hat_temp*exp(-xb_temp)  // obtain the correct residuals  }{p_end}
{phang2}{cmd:.  * rhs of test  }{p_end}
{phang2}{cmd:.  logit outside ivpois trips cbd ptn worker weekend pt  // use the instrument instead of the endogenous variable  }{p_end}
{phang2}{cmd:.  predict p_hat_temp, pr  }{p_end}
{phang2}{cmd:.  gen rhs_temp = (mean_u_temp)/p_hat_temp  }{p_end}
{phang2}{cmd:.  * run the test  }{p_end}
{phang2}{cmd:.  reg lhs_temp rhs_temp if outside, nocons   }{p_end}
{phang2}{cmd:.  matrix b = e(b)  }{p_end}
{phang2}{cmd:.  ereturn post b  }{p_end}
{phang2}{cmd:.  * drop created variables  }{p_end}
{phang2}{cmd:.  cap drop *temp  }{p_end}
{phang2}{cmd:.  end  }{p_end}
 
{phang2}{cmd:. bootstrap lambda = _b[rhs_temp] , reps(50): IVPoisson_boostrap  }{p_end}
{phang2}{cmd:. test lambda==1  }{p_end}


{pstd} Fourth, you can convert your results into latex using esttab where "eps" provides the convergence criteria:
{p_end}
{hline}
{phang2}{cmd:. eststo clear}{p_end}
{phang2}{cmd:. eststo: i2SLS_HDFE trips cbd ptn worker weekend (tcost=pt) , delta(1) robust}{p_end}
{phang2}{cmd:. eststo: i2SLS_HDFE trips cbd ptn worker weekend (tcost=pt) , delta(10) robust}{p_end}
{phang2}{cmd:. eststo: esttab * using table_iv.tex,  scalars(delta eps) }{p_end}
{hline}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:i2SLS_HDFE} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}} marks the sample used for estimation {p_end}
{synopt:{cmd:e(eps)}} sum of the absolute differences between the parameters from the last two iterations of i2SLS {p_end}
{synopt:{cmd:e(k)}} number of iterations of iOLS{p_end}
