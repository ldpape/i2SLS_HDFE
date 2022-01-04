# i2SLS_HDFE : Iterated Two Stage Least Squares (i2SLS) with High Dimensional Fixed Effects (HDFE) using the "ivreg2" function for any delta>0
# /!\ Still preliminary : use at your own risk.
To install this code into Stata, run the following (requires at least Stata 14) :

> cap ado uninstall i2SLS_HDFE

> net install i2SLS_HDFE, from("https://raw.githubusercontent.com/ldpape/i2SLS_HDFE/master/")

You will need to have the following packages installed using:
>ssc install reghdfe

>ssc install hdfe

>ssc install ivreg2
