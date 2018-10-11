clear all
set more off, perm
set seed 12345

global dir "/Users/erinmarkiewitz/Dropbox/Phd_Coursework/Econ675/hw2"
global datadir $dir\data
global resdir $dir\results

cap log close
log using $resdir\pset2_stata.smcl, replace


*******
*** Problem 1 
*******
/*
set obs 10000
timer on 1
program IMSEsim, rclass
drop _all
set obs 1000
gen x =  rnormal(-1/4, 5/8)
gen fx = normalden(-1/4, 5/8)
_kdens x, at(x) generate(fxh) bw(.5) kernel(epan2)
gen diffLI = (fx - fxh)^2 
gen diffL0 = 0


forvalues i = 1/1000 {
	_kdens x if _n != `i', at(x) generate(fxh`i') bw(.5) kernel(epan2) 
	replace diffL0 = (fx - fxh`i')^2 if _n == `i'
}

qui summ diffLI
return scalar data1 = r(mean)
qui summ diffL0
return scalar data2 = r(mean)
end



simulate IMSE_LI=r(data1) IMSE_L0 = r(data2), reps(1) nodots: IMSEsim
timer off 1
timer list
*/


*******
*** Problem 3
*******
set obs 500
drop _all
local theta = 1 
local d = 5
local n = 500

mata:
X 	= uniform(`n',`d'):*2 :-1
ep	= invnormal(uniform(`n',1)):*0.3637899:*(1 :+ rowsum(X:^2)) 
gx	= exp(rowsum(X:^2))
T	= invnormal(uniform(`n',1)) + rowsum(X:^2):^.5 :>= 0
Y   = T + gx + ep 


A = asarray_create("real",1)
cons= J(500,1,1)
X2 	= X:^2
X3 	= X:^3
X4 	= X:^4
X5 	= X:^5
X6 	= X:^6
X7 	= X:^7
X8 	= X:^8
X9 	= X:^9
X10 = X:^10
X1k = X#X
X2k = X2#X2
X3k = X3#X3
X4k = X4#X4
X1k = X1k[1::`n',2::5], X1k[1::`n', 8::10], X1k[1::`n',14::15], X1k[1::`n', 20]
X2k = X2k[1::`n',2::5], X2k[1::`n', 8::10], X2k[1::`n',14::15], X2k[1::`n', 20]
X3k = X3k[1::`n',2::5], X3k[1::`n', 8::10], X3k[1::`n',14::15], X3k[1::`n', 20]
X4k = X4k[1::`n',2::5], X4k[1::`n', 8::10], X4k[1::`n',14::15], X4k[1::`n', 20]



asarray(A,1,X)
asarray(A,2,(asarray(A,1),X2))
asarray(A,3,(asarray(A,2),X1k))
asarray(A,4,(asarray(A,3),X3))
asarray(A,5,(asarray(A,4),X2k))
asarray(A,6,(asarray(A,5),X4))
asarray(A,7,(asarray(A,6),X3k))
asarray(A,8,(asarray(A,7),X5))
asarray(A,9,(asarray(A,8),X4k))
asarray(A,10,(asarray(A,9),X6))
asarray(A,11,(asarray(A,10),X7))
asarray(A,12,(asarray(A,11),X8))
asarray(A,13,(asarray(A,12),X9))
asarray(A,14,(asarray(A,13),X10))


Z = qrsolve(cons,(T,asarray(A,1)))
ZZ  = Z*Z'
Yhat = ZZ*Y
W = diag(ZZ)
lolk = mean((Y-Yhat)/(1-W)^2)


ZQ = (cons,asarray(A,1))*invsym((cons,asarray(A,1))'*(cons,asarray(A,1)))*(cons,asarray(A,1))'

M = I(`n') - ZQ
YM = M*Y
TM = M*T
theta_hat = (TM'*YM) / (TM'*TM)
sigma = diag((YM - TM*theta_hat)'*(YM - TM*theta_hat))
se_hat = sqrt(TM'*sigma*TM / (TM'*TM))


end

/*
forv iter = 0/9{
svmat X`iter', n(x`iter'_)
}


forv iter = 1/4{
svmat X`iter'k, n(xk`iter')
}

svmat ep, n(ep)
svmat gx, n(gx)
svmat T, n(T)
svmat Y, n(Y)
g cons = 1 //only for npseries

local poly = "cons x1_* x2_* xk1* x3_* xk2* x4_* x3k* x5_* xk4* x6_* x7_* x8_* x9_* x10_*"
forv iter = 1/14{
local wc = 5 + 5*`iter'
local varlist = substr("`poly'",1,`wc')
di "`varlist'"
semipar Y1 T1, nonpar(`varlist') 


}

di "`poly'"


