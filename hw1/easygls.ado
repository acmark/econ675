program drop _all

program define easyfgls, eclass
	syntax varlist [if] [in]
	marksample touse
	*gettoken gets the first token from a list store first token in depvar remaining tokens in indepvars
	gettoken depvar indepvars :varlist 
	tempname b W a
	mata: easyfgls_mata("`depvar'", "`indepvars'", "`W'", "`bn'", "`a'" , "`N'")
	matrix colnames `W' = `indepvars' _cons
	ereturn local cmd "easyfgls"
	display "measyfgls"
	ereturn display
end


********************************************************************************
** MATA FUNCTIONS USED
********************************************************************************

mata: mata clear
mata:
	void easyfgls_mata(string scalar depvar, string scalar indepvars, ///
					    string matrix ws, string scalar bns, string scalar as, string scalar ns) {
		real vector y, bh
		real matrix X, vh
		real matrix W, vh
		real vector bns, bh
		real scalar n
		real scalar a
		y = st_data(., depvar)
		X = st_data(., tokens(indepvars))
		W = st_data(.,ws)
		n = rows(y)
		X = (X,J(n,1,1))
		bh = cholinv(cross(X,X))*cross(X,y);
		
		S = optimize_init()
		optimize_init_evaluator(S, &MyPoisson2())
		optimize_init_evaluatortype(S, "gf0")
		optimize_init_argument(S, 1, y)
		optimize_init_argument(S, 2, X)
		optimize_init_params(S, bh')
		bh = optimize(S)
		vh = optimize_result_V_oim(S)
		st_matrix(bs, bh)
		st_matrix(vs, vh)
		st_numscalar(ns, n)
	}

	void MyPoisson2(real scalar todo, real vector beta, real vector y, real matrix X, lnf, grad, Hess) {
		real vector xb
		xb = X*beta'
		lnf = -exp(xb) + y :*xb - lnfactorial(y)
	}
end
*/
