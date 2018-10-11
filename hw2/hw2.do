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





**************
**** Problem 2 
**************

**************
**** Problem 2a 
***************
set obs 1000


* Define cross validation function: CV(list, i): vars=variable list, i = max polynomial
mata
	void CV(vars, i) {
		st_view(y=., ., "y")
		st_view(X=., ., tokens(vars))
		XpX  = cross(X, X)
		XpXinv  = invsym(XpX)
		b  = XpXinv*cross(X, y)
		w = diagonal(X*XpXinv*X')
		muhat = X*b
		num = (y - muhat):*(y - muhat)
		den= (J(1000,1,1) - w):*(J(1000,1,1) - w)
		div = num:/den
		CV = mean(div)
		CV
		st_numscalar("mCV"+strofreal(i), CV)
	}
end


* Program which runs the monte-carlo experiment
program CVsim, rclass
	drop _all
	set obs 1000
	forvalues i = 0/20 { 
		gen CV`i' = 0
	}
	gen x = runiform(-1,1)
	gen e = x^2*(rchi2(5)-5)
	gen y = exp(-0.1*(4*x-1)^2)*sin(5*x)+e
	forvalues i = 0/20 { 
		gen x`i' = x^`i'
	}
	forvalues i = 0/20 {
		global xlist = "x0-x`i'"
		di "$xlist"
		mata CV("$xlist", `i')
		replace CV`i' = mCV`i'
	}
end 



* Run the experiment
set seed 12345
simulate CV0=CV0 CV1=CV1 CV2=CV2 CV3=CV3 CV4=CV4 CV5=CV5 CV6=CV6 CV7=CV7 CV8=CV8 /// 
	CV9=CV9 CV10=CV10 CV11=CV11 CV12=CV12 CV13=CV13 CV14=CV14 CV15=CV15 ///
	CV16=CV16 CV17=CV17 CV18=CV18 CV19=CV19 CV20=CV20, reps(100) nodots: CVsim
collapse *
gen i = 1
reshape long CV, i(i) j(k)
sort CV
local min = k[1]
twoway scatter CV k, ytitle("Mean CV") xtitle("K") xlabel(0(2)20) xmtick(0(1)20) xline(`min')




***************
***Problem 2b 
***************

* Program which runs the monte-carlo experiment for mu_0
program muhatsim, rclass
	drop _all
	set obs 1000
	gen x = runiform(-1,1)
	gen e = x^2*(rchi2(5)-5)
	gen y = exp(-0.1*(4*x-1)^2)*sin(5*x)+e
	forvalues p = 0/7 { 
		gen x`p' = x^`p'
	}
	reg y x0-x7, nocons
	clear
	set obs 11
	gen n = _n
	gen foo = 1
	gen x = -1+(_n-1)/5
	forvalues p = 0/7 { 
		gen x`p' = x^`p'
	}
	predict muhat
	predict se, stdp
	generate lb = muhat - invnormal(0.975)*se
	generate ub = muhat + invnormal(0.975)*se

	
	
	keep n muhat foo lb ub 
	reshape wide muhat lb ub, i(foo) j(n)
end



set seed 12345
simulate muhat1=muhat1 muhat2=muhat2 muhat3=muhat3 muhat4=muhat4 muhat5=muhat5 ///
	muhat6=muhat6 muhat7=muhat7 muhat8=muhat8 muhat9=muhat9 muhat10=muhat10 muhat11=muhat11 ///
	ub1=ub1 ub2=ub2 ub3=ub3 ub4=ub4 ub5=ub5 ub6=ub6 ub7=ub7 ub8=ub8 ub9=ub9 ub10=ub10 ub11=ub11 ///
	lb1=lb1 lb2=lb2 lb3=lb3 lb4=lb4 lb5=lb5 lb6=lb6 lb7=lb7 lb8=lb8 lb9=lb9 lb10=lb10 lb11=lb11, reps(1000) nodots: muhatsim
gen i = _n
reshape long muhat ub lb, i(i) j(grid)
collapse muhat ub lb, by(grid)
gen x = -1+ (grid-1)/5
twoway (function y = exp(-0.1*(4*x-1)^2)*sin(5*x), range(-1 1) lcolor(red)) ///
	(line muhat x, lcolor(gs6)) (line lb x, lcolor(gs6) lpattern(dash)) (line ub x, lcolor(gs6) lpattern(dash)), ///
	legend(order(1 "DGP" 2 "Prediction" 3 "Confidence Interval") rows(1)) ytitle(Y) xtitle(X) title("Q2b")
graph export $resdir\pset2q2b.png, replace



* Program which runs the monte-carlo experiment for mu_1
program dmuhatsim, rclass
	drop _all
	set obs 1000
	gen x = runiform(-1,1)
	gen e = x^2*(rchi2(5)-5)
	gen y = exp(-0.1*(4*x-1)^2)*((0.8-3.2*x)*sin(5*x)+5*cos(5*x)) + e
	forvalues p = 0/7 { 
		gen x`p' = x^`p'
	}
	reg y x0-x7, nocons
	clear
	set obs 11
	gen n = _n
	gen foo = 1
	gen x = -1+(_n-1)/5
	forvalues p = 0/7 { 
		gen x`p' = x^`p'
	}
	predict dmuhat
	predict se, stdp
	generate lb = dmuhat - invnormal(0.975)*se
	generate ub = dmuhat + invnormal(0.975)*se

	
	
	keep n dmuhat foo lb ub 
	reshape wide dmuhat lb ub, i(foo) j(n)
end


set seed 12345
simulate dmuhat1=dmuhat1 dmuhat2=dmuhat2 dmuhat3=dmuhat3 dmuhat4=dmuhat4 dmuhat5=dmuhat5 ///
	dmuhat6=dmuhat6 dmuhat7=dmuhat7 dmuhat8=dmuhat8 dmuhat9=dmuhat9 dmuhat10=dmuhat10 dmuhat11=dmuhat11 ///
	ub1=ub1 ub2=ub2 ub3=ub3 ub4=ub4 ub5=ub5 ub6=ub6 ub7=ub7 ub8=ub8 ub9=ub9 ub10=ub10 ub11=ub11 ///
	lb1=lb1 lb2=lb2 lb3=lb3 lb4=lb4 lb5=lb5 lb6=lb6 lb7=lb7 lb8=lb8 lb9=lb9 lb10=lb10 lb11=lb11, reps(1000) nodots: dmuhatsim
gen i = _n
reshape long dmuhat ub lb, i(i) j(grid)
collapse dmuhat ub lb, by(grid)
gen x = -1+ (grid-1)/5
twoway (function y = exp(-0.1*(4*x-1)^2)*((0.8-3.2*x)*sin(5*x)+5*cos(5*x)), range(-1 1) lcolor(red)) ///
	(line dmuhat x, lcolor(gs6)) (line lb x, lcolor(gs6) lpattern(dash)) (line ub x, lcolor(gs6) lpattern(dash)), ///
	legend(order(1 "DGP" 2 "Prediction" 3 "Confidence Interval") rows(1)) ytitle(Y) xtitle(X) title("Q2d")
graph export $resdir\pset2q2d.png, replace




/*

*******
*** Problem 3
*******

drop _all
set obs 1000
local theta = 1 
local d = 5
local n = 500

forvalues p = 1/14 { 
	gen v_hat`p' = .
	gen theta_hat`p' = .

}


mata:
	void polyloop(i) {
		real matrix v_hat
		real matrix theta_hat

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


		theta_hat = I(1,14):*0
		v_hat = I(1,14):*0
		k_hat = I(1,14):*0

		for (j=1; j<=14; j++) {

			Z = qrsolve(cons,(T,asarray(A,j)))
			ZZ  = Z*Z'
			Yhat = ZZ*Y
			W = diag(ZZ)
			ZQ = (cons,asarray(A,j))*invsym((cons,asarray(A,j))'*(cons,asarray(A,j)))*(cons,asarray(A,j))'
			M = I(`n') - ZQ
			YM = M*Y
			TM = M*T
			theta_hat[1,j] = (TM'*YM) / (TM'*TM)
			sigma = diag(ZQ*(Y-T*theta_hat[1,j]))
			v_hat[1,j] = invsym(T'*ZQ*T)*(T'*ZQ*sigma*ZQ*T)*invsym(T'*ZQ*T)
			st_store(i, "v_hat"+strofreal(j), v_hat[1,j])
			st_store(i, "theta_hat"+strofreal(j), theta_hat[1,j])
		}

	}
end


forvalues i = 1/10 {
	mata polyloop(`i')
}
save output_q3.dta, replace

use output_q3,clear
gen obs = _n
reshape long v_hat theta_hat, i(obs) j(k)
collapse (mean) mean_v_hat= v_hat  mean_theta_hat=theta_hat (sd) sd_theta_hat = theta_hat, by(k) 
gen mean_bias = mean_theta_hat  - 1
gen v_theta_hat = sd_theta_hat^2

*coverage rate test
//mean_theta_hat + 1.96 * sd_theta_hat/sqrt(`n')

log close
translate $resdir\pset2_stata.smcl $resdir\pset2_stata.pdf, replace


