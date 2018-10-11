clear all
set more off
global dir "/Users/aaronmarkiewitz/Dropbox/Phd_Coursework/Econ 675"
cd "$dir"
set seed 675
set scheme s2mono

global dir $dir\data
global dir $dir\data
import delimited LaLonde_1986.csv


*Question 2.5

gen cons		= 1	
gen educ2 		= educ^2 
gen blacke74 	= black*earn74
local y earn78
local xlist treat black age educ educ2 earn74 blacke74 u74 u75 cons


local alpha .05
local n = _N
local d = wordcount("`xlist'")
*Procedure writen for Question 2.4
mata:
W = I(`d')
st_view(y=., ., "`y'")
st_view(X=., ., tokens("`xlist'"))

XW  	= cross(X', W)
XWX  	= cross(XW, X)
XWXinv  = invsym(XWX)
XWXinvc = cholinv(XWX)
b  = XWXinv*cross(X, y)
bc = XWXinvc*cross(X, y)

e  = y - X*b
v_e = diag(e:*e)
n  = rows(X)

dof = `n' - `d'
V  = XWXinv*(X'*v_e*X)*XWXinv*(n/dof)
se = sqrt(diagonal(V))

t_stat = (b):/se
p = 2*ttail(rows(X)-cols(X),abs(t_stat))
cil = b+invt(dof, (`alpha'/2))*se
ciu = b+invt(dof, (1-`alpha'/2))*se
st_matrix("b", b')
st_matrix("V", V')

end

mata: output = (b,se,t_stat,p,cil,ciu)
mmat2tex output using .\mata_out.tex, replace

reg `y' `xlist', nocons robust

esttab using .\q2_5outreg.tex, replace


** Question 3 ** 
clear all
import delimited LaLonde_1986.csv



*Question 3.1

sum earn78 if treat==0
local N0 = r(N)
local mu0 = r(mean)
local sd0 = r(sd)
local V0 = r(Var)/r(N)
local sig_sq0 = r(Var)

sum earn78 if treat==1
local N1 = r(N)
local mu1 = r(mean)
local sd1 = r(sd)
local V1 = r(Var)/r(N)
local sig_sq1 = r(Var)

local tau = `mu1'-`mu0'
local v = sqrt(`V1'+`V0')
local t_stat = `tau'/`v'
local p = 2*normal(-abs(`t_stat'))

local mu0 = round(`mu0', .01)
local mu1 = round(`mu1', .0001)
local sd0 = round(`sd0', .01)
local sd1 = round(`sd1', .0001)

di "`tau'"


local cil = `tau' - invnormal(0.975)*`v'
local ciu = `tau' + invnormal(0.975)*`v'

di "`CIlower'"
di "`CIupper'"



mata: output = (`tau',`v',`t_stat',`p',`cil',`ciu')
mmat2tex output using .\31_out_mata.tex, replace



*Question 3.2a

* difference in means estimator
permute treat diffmean=(r(mu_2)-r(mu_1)), reps(1999) nowarn: ttest earn78, by(treat) 
matrix pval_dm = r(p)
local p_dm = pval_dm[1,1]
di "DM p-value"
di "DM p-value= `p_dm'"


* KS statistic
permute treat ks=r(D), reps(1999) nowarn: ksmirnov earn78, by(treat) 
matrix pval_ks = r(p)
local p_ks = pval_ks[1,1]
di "KS p-value= `p_ks'"



*Question 3.2b

* Infer missing values under the null of constant treatment effect
gen     Y1_imputed = earn78
replace Y1_imputed = earn78 + `tau' if treat==0

gen     Y0_imputed = earn78
replace Y0_imputed = earn78 - `tau' if treat==1


* Write program to put into bootstrap function
program define meandiff, rclass
	summarize   Y1_imputed if treat==1
	local 		tau1 = r(mean)
	sum 		Y0_imputed if treat==0
	local 		tau0 = r(mean)
	return      scalar meandiff = `tau1' - `tau0'
end

* Run bootstrap function using meandiff program
bootstrap diff = r(meandiff), reps(1999): meandiff

/*


*Question 3.3
twoway function y= 1 - normal(invnormal(0.975)-x/`v') + normal(-invnormal(0.975)-x/`v'), range(-5000 5000)
graph export stata_power.png,replace

mata: mata clear
mata:
 function myfunc(N, s0, s1, p, tau){
 
   return(1 - normal(invnormal(0.975)-tau/sqrt(1/N*s1*(1/p)+1/N*s0*(1/(1-p)))) +
       normal(-invnormal(0.975)-tau/sqrt(1/N*s1*(1/p)+1/N*s0*(1/(1-p)))) -0.8)
 
 }
 s0 =  30072466.58373794
 s1 =  61896056.06715253
    p     = 2/3
   tau   = 1000
   
  mm_root(x=., &myfunc(), 1000, 1500, 0, 10000, s0,s1, p ,tau)
      
	x
	  
end 


