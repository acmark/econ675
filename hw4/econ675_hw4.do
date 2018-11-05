// Erin Markiewitz
// ECON 675 Assignment 4
********************************************************************************
clear all
set more off, perm
set seed 12345
global dir "/Users/erinmarkiewitz/Dropbox/Phd_Coursework/Econ675/hw4"
global datadir $dir\data
global resdir $dir\results

cap log close
log using $resdir\pset4_stata.smcl, replace

****
* Question 1
****
****
* Question 2
****

****
* Question 3
****
clear all
set seed 12345 

*construct dgp variance covariance matrix 
matrix P = (1,.85 \.85, 1)
mat A = cholesky(P)

program modelsim, rclass
	args A
	drop _all 
	set obs 50

	*generate component normal variables 
	gen c1= invnorm(uniform())
	gen c2= invnorm(uniform())

	*use cholesky decomp to back out x,z 
	gen x = `A'[1,1] * c1 + `A'[1,2] * c2
	gen z = `A'[2,1] * c1 + `A'[2,2] * c2

	*general epsilon and outcome variable 
	gen e= invnorm(uniform())
	gen y = 0.5*x + z + e 


	*simulate model selection process (flip order of regs for speed) 
	reg y x 
	scalar beta_tilde = _b[x]
	
	reg y x z 
	scalar beta_hat = _b[x]


	if abs(_b[z]/_se[z])>=1.96 {
		scalar beta_check = beta_hat
	}
	else{
		scalar beta_check = beta_tilde
	}
	
end

simulate beta_hat = beta_hat beta_check=beta_check beta_tilde=beta_tilde, ///
seed(1234) reps(1000): modelsim A

**TODO: empircal coverage rate 

**

sum * 
estpost summarize *
estout using hw4_q3_1_stata.tex, cells("mean sd min max") style(tex)  replace

kdensity beta_hat, normal name(beta_hat,replace)
gr export hw4_q3_bhat_stata.png ,replace
kdensity beta_check, normal name(beta_check,replace)
gr export hw4_q3_bcheck_stata.png ,replace
kdensity beta_tilde,normal name(beta_tilde,replace)
gr export hw4_q3_btilde_stata.png ,replace

sum beta_hat
gen cov_upper = r(mean) + 1.96*r(sd)/sqrt(50)
gen cov_lower = r(mean) - 1.96*r(sd)/sqrt(50)
gen covrate_hat =  cond(beta_hat<= cov_upper & beta_hat>= cov_lower,1,0)
gen covrate_tilde =  cond(beta_tilde<= cov_upper & beta_tilde>= cov_lower,1,0)
gen covrate_check  = cond(beta_check<= cov_upper & beta_check>= cov_lower,1,0)
sum covrate*
estout using hw4_q3_2_stata.tex, cells("mean ") style(tex)  replace










cap log close
