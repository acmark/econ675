// Erin Markiewitz
// ECON 675 Assignment 3
********************************************************************************
clear all
set more off, perm
set seed 12345
global dir "/Users/erinmarkiewitz/Dropbox/Phd_Coursework/Econ675/hw3"
global datadir $dir\data
global resdir $dir\results

cap log close
log using $resdir\pset2_stata.smcl, replace



*******
*** Problem 1 
*******
use pisofirme, clear
gen s = 1 - cond(danemia==.,1,0)
gen ls_incomepc =log(S_incomepc+1)
glm s S_age S_HHpeople ls_incomepc, family(binomial) link(logit) r
estout using hw3_q1_9a_stata.tex, cells("b var t p ci") style(tex)  replace



glm s S_age S_HHpeople ls_incomepc, family(binomial) link(logit) vce(bs, r(99)  seed(123) nodots)
estout using hw3_q1_9b_stata.tex, cells("b p ci") style(tex)  replace
predict prop_score 
kdensity prop_score
gr export  hw3_q1_9a_stata.png, replace



*******
*** Problem 2 a
*******
use pisofirme, clear
gen constant = 1 
gen s = 1 - cond(danemia==.,1,0)
gen ls_incomepc =log(S_incomepc+1)
gmm (danemia - logistic({xb: dpisofirm S_age S_HHpeople ls_incomepc})) , instruments(dpisofirm S_age S_HHpeople ls_incomepc, noconstant) winitial(identity) vce(bs, r(49)  seed(123) nodots)
estout using hw3_q2_2a_stata.tex, cells("b ci") style(tex)  replace


*******
*** Problem 2 3c
*******
glm s S_age S_HHpeople ls_incomepc dpisofirme, family(binomial) link(logit)
predict prop_score 
gen w_s_age = S_age/prop_score 
gen w_s_hhpeople = S_HHpeople/prop_score 
gen w_ls_incomepc = ls_incomepc/prop_score 
gen w_dpisofirme = dpisofirme/prop_score 
gmm (danemia - logistic({xb: S_age S_HHpeople ls_incomepc dpisofirme})), instruments(w_*, noconstant)
estout using hw3_q2_3c_stata.tex, cells("b ci") style(tex)  replace


*******
*** Problem 2 3d
*******
drop if prop_score < 0.1
gmm (danemia - logistic({xb: S_age S_HHpeople ls_incomepc dpisofirme}) ), instruments(w_*, noconstant) vce(bs, r(49)  seed(123) nodots)
estout using hw3_q2_3d_stata.tex, cells("b ci") style(tex)  replace



*******
*** Problem 3 a
*******
clear all
set obs 1000
set seed 123 
gen x = runiform()

sum x
local max_x = r(max)
bs max_x_star = r(max) , reps(599) saving(mbs,replace): sum x
use mbs, clear
gen stat = 1000*(`max_x' - max_x_star)
twoway (histogram stat, bin(16) ) (function exp(-x),range(0 8)) ,title("Distribution of Bootstrap Statistic")
gr export hw3_Q3_1_stata.png ,replace



*******
*** Problem 3 b
*******
clear all
set obs 1000
set seed 123 
gen x = runiform()
sum x
local max_x = r(max)
di `max_x'


program pbs, rclass
	args max_x
	drop _all
	set obs 1000
	gen x_pbs = runiform(0,`max_x')
	egen max_x_pbs = max(x_pbs)
	drop if _n>1
end
pbs `max_x'

simulate max_pbs= max_x_pbs, reps(599): pbs `max_x'


gen stat = 1000*(`max_x' - max_pbs)
twoway (histogram stat, bin(16) ) (function exp(-x),range(0 8)), title("Distribution of Parametric Bootstrap Statistic")
gr export hw3_Q3_2_stata.png ,replace






