// Erin Markiewitz
// ECON 675 Assignment 4
********************************************************************************
clear all
set more off, perm
set seed 12345
global dir "/Users/erinmarkiewitz/Dropbox/Phd_Coursework/Econ675/hw4"
global datadir $dir\data
global resdir $dir\results
cd $dir
cap log close
log using $resdir\pset4_stata.smcl, replace

****
* Question 1
****


****
* Question 2
****
/*
clear all
set seed 12345
scalar n1 = 185
scalar n0 = 260
scalar n2 = 2490
import delimited using LaLonde_all.csv, clear 
gen log_re74 = log(re74 + 1)
gen log_re75 = log(re75 + 1)
gen age2 = age^2
gen age3 = age^3
gen educ2 = educ^2
gen black_u74 = black*u74
gen educ_log_re74 = educ*log_re74

local covars_a = "age educ black hisp  married nodegr log_re74 log_re75"
local covars_b = "age educ black hisp  married nodegr log_re74 log_re75 age2 educ2 u74 u75"
local covars_c = "age educ black hisp  married nodegr log_re74 log_re75 age2 educ2 u74 u75 age3 black_u74 educ_log_re74"


*initialize matricies (0) Experimental data, (2) PSID Control

mat ate0 = J(19,4,.) 
mat att0 = J(19,4,.) 
mat ate2 = J(19,4,.) 
mat att2 = J(19,4,.) 

*diff in means 
reg re78 treat if treat ==1 | treat == 0, hc2
mat ate0[1,1] = _b[treat]
mat ate0[1,2] = _se[treat]
mat ate0[1,3] = ate0[1,1] - _se[treat] * 1.96
mat ate0[1,4] = ate0[1,1] + _se[treat] * 1.96

mat att0[1,1] = _b[treat]
mat att0[1,2] = _se[treat]
mat att0[1,3] = att0[1,1] - _se[treat] * 1.96
mat att0[1,4] = att0[1,1] + _se[treat] * 1.96

reg re78 treat if treat ==1 | treat == 2, hc2

mat ate2[1,1] = - _b[treat]
mat ate2[1,2] = _se[treat]
mat ate2[1,3] = ate2[1,1] - _se[treat] * 1.96
mat ate2[1,4] = ate2[1,1] + _se[treat] * 1.96

mat att2[1,1] = - _b[treat]
mat att2[1,2] = _se[treat]
mat att2[1,3] = att2[1,1] - _se[treat] * 1.96
mat att2[1,4] = att2[1,1] + _se[treat] * 1.96




*OLS
local base_count = 2
foreach num of numlist 0 2 { 
local count = `base_count'

foreach cv in a b c { 
di "`covars_`cv''"
di "`num'"
reg re78 treat `covars_`cv'' if treat ==1 | treat == `num', hc2
mat ate`num'[`count',1] = _b[treat] * (1 - `num') // stata thinks 2 is treatment
mat ate`num'[`count',2] = _se[treat]
mat ate`num'[`count',3] = ate`num'[`count',1] - _se[treat] * 1.96
mat ate`num'[`count',4] = ate`num'[`count',1] + _se[treat] * 1.96

mat att`num'[`count',1] = _b[treat] * (1 - `num') // stata thinks 2 is treatment
mat att`num'[`count',2] = _se[treat]
mat att`num'[`count',3] = att`num'[`count',1] - _se[treat] * 1.96
mat att`num'[`count',4] = att`num'[`count',1] + _se[treat] * 1.96

local ++count 
}
}


*Reg. Impute
local base_count = `count' 
foreach num of numlist 0 2 { 
local count = `base_count'

foreach cv in a b c { 
di "`covars_`cv''"
di "`num'"
teffects ra (re78 `covars_`cv'') (treat) if treat ==1 | treat == `num', ate
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat ate`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat ate`num'[`count',2] = _se[`colv']
mat ate`num'[`count',3] = ate`num'[`count',1] - ate`num'[`count',2] * 1.96
mat ate`num'[`count',4] = ate`num'[`count',1] + ate`num'[`count',2] * 1.96

teffects ra (re78 `covars_`cv'') (treat) if treat ==1 | treat == `num', atet
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat att`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat att`num'[`count',2] = _se[`colv']
mat att`num'[`count',3] = att`num'[`count',1] - att`num'[`count',2] * 1.96
mat att`num'[`count',4] = att`num'[`count',1] + att`num'[`count',2] * 1.96

local ++count 
}
}


*IPW 
local base_count = `count' 
foreach num of numlist 0 2 { 
local count = `base_count'

foreach cv in a b c { 
di "`covars_`cv''"
di "`num'"
capture teffects ipw (re78) (treat `covars_`cv'', probit) if treat ==1 | treat == `num', ate osample(otest) iter(50)
 if _rc==498 {
               display "Overlap Assumption Violated"
			   teffects ipw (re78) (treat `covars_`cv'', probit) if (treat ==1 | treat == `num') & otest ==0, ate iter(50)
}
drop otest
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat ate`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat ate`num'[`count',2] = _se[`colv']
mat ate`num'[`count',3] = ate`num'[`count',1] - ate`num'[`count',2] * 1.96
mat ate`num'[`count',4] = ate`num'[`count',1] + ate`num'[`count',2] * 1.96

capture teffects ipw (re78) (treat `covars_`cv'', probit) if treat ==1 | treat == `num', atet osample(otest) iter(50)
 if _rc==498 {
               display "Overlap Assumption Violated"
			    teffects ipw (re78) (treat `covars_`cv'', probit) if (treat ==1 | treat == `num') & otest ==0, atet iter(50)
}
drop otest
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat att`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat att`num'[`count',2] = _se[`colv']
mat att`num'[`count',3] = att`num'[`count',1] - att`num'[`count',2] * 1.96
mat att`num'[`count',4] = att`num'[`count',1] + att`num'[`count',2] * 1.96


local ++count 
}
}


*DR 
local base_count = `count' 
foreach num of numlist 0 2 { 
local count = `base_count'

foreach cv in a b c { 
di "`covars_`cv''"
di "`num'"
capture teffects ipwra (re78) (treat `covars_`cv'', probit) if treat ==1 | treat == `num', ate osample(otest) iter(50)
 if _rc==498 {
               display "Overlap Assumption Violated"
			   teffects ipw (re78) (treat `covars_`cv'', probit) if (treat ==1 | treat == `num') & otest ==0, ate iter(50)
}
drop otest
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat ate`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat ate`num'[`count',2] = _se[`colv']
mat ate`num'[`count',3] = ate`num'[`count',1] - ate`num'[`count',2] * 1.96
mat ate`num'[`count',4] = ate`num'[`count',1] + ate`num'[`count',2] * 1.96

capture teffects ipwra (re78) (treat `covars_`cv'', probit) if treat ==1 | treat == `num', atet osample(otest) iter(50)
 if _rc==498 {
               display "Overlap Assumption Violated"
			    teffects ipw (re78) (treat `covars_`cv'', probit) if (treat ==1 | treat == `num') & otest ==0, atet iter(50)
}
drop otest
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat att`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat att`num'[`count',2] = _se[`colv']
mat att`num'[`count',3] = att`num'[`count',1] - att`num'[`count',2] * 1.96
mat att`num'[`count',4] = att`num'[`count',1] + att`num'[`count',2] * 1.96


local ++count 
}
}


*Reg. Impute
local base_count = `count' 
foreach num of numlist 0 2 { 
local count = `base_count'

foreach cv in a b c { 
di "`covars_`cv''"
di "`num'"
teffects nnmatch (re78 `covars_`cv'') (treat) if treat ==1 | treat == `num', ate nneighbor(1) metric(maha)
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat ate`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat ate`num'[`count',2] = _se[`colv']
mat ate`num'[`count',3] = ate`num'[`count',1] - ate`num'[`count',2] * 1.96
mat ate`num'[`count',4] = ate`num'[`count',1] + ate`num'[`count',2] * 1.96

teffects nnmatch (re78 `covars_`cv'') (treat) if treat ==1 | treat == `num', atet nneighbor(1) metric(maha)
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat att`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat att`num'[`count',2] = _se[`colv']
mat att`num'[`count',3] = att`num'[`count',1] - att`num'[`count',2] * 1.96
mat att`num'[`count',4] = att`num'[`count',1] + att`num'[`count',2] * 1.96

local ++count 
}
}

*PS Matching 
local base_count = `count' 
foreach num of numlist 0 2 { 
local count = `base_count'

foreach cv in a b c { 
di "`covars_`cv''"
di "`num'"
capture teffects psmatch (re78) (treat `covars_`cv'', probit) if treat ==1 | treat == `num', ate osample(otest) iter(50)
 if _rc==498 {
               display "Overlap Assumption Violated"
			   teffects psmatch (re78) (treat `covars_`cv'', probit) if (treat ==1 | treat == `num') & otest ==0, ate iter(50)
}
cap drop otest
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat ate`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat ate`num'[`count',2] = _se[`colv']
mat ate`num'[`count',3] = ate`num'[`count',1] - ate`num'[`count',2] * 1.96
mat ate`num'[`count',4] = ate`num'[`count',1] + ate`num'[`count',2] * 1.96

capture teffects psmatch (re78) (treat `covars_`cv'', probit) if treat ==1 | treat == `num', atet osample(otest) iter(50)
 if _rc==498 {
               display "Overlap Assumption Violated"
			    teffects psmatch (re78) (treat `covars_`cv'', probit) if (treat ==1 | treat == `num') & otest ==0, atet iter(50)
}
cap drop otest
local colnmsb: coln e(b)
local colnmsv: coln e(V)
local colb: word 1 of `colnmsb'
local colv: word 1 of `colnmsv'
mat att`num'[`count',1] = _b[`colb'] * (1 - `num') // stata thinks 2 is treatment
mat att`num'[`count',2] = _se[`colv']
mat att`num'[`count',3] = att`num'[`count',1] - att`num'[`count',2] * 1.96
mat att`num'[`count',4] = att`num'[`count',1] + att`num'[`count',2] * 1.96


local ++count 
}
}



mat li ate0
mat li ate2

mat li att0
mat li att2
**TODO: put into charts 

*/


****
* Question 3
****
/*
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
