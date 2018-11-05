


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

sum * 
kdensity beta_hat, normal name(beta_hat,replace)
kdensity beta_check, normal name(beta_check,replace)
kdensity beta_tilde,normal name(beta_tilde,replace)


