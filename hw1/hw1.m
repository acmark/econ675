clear all; close all; clc;
cd('/Users/aaronmarkiewitz/Dropbox/Phd_Coursework/Econ 675')
rng(675)
n           = 10;
d           = 2;
alpha       = .05;
beta_null   = [0 0]';
y           = randn(n,1);
ep          = randn(n,1);
X           = randn(n,d);
W           = eye(n);

%% 2.4.a: point-estimate
L = chol((X'*W*X),'lower');
chol = inv(L)'*inv(L);
beta_chol_hat = chol * (X'*W*y)
beta_hat = inv((X'*W*X)) * (X'*W*y) 

%% 2.4.b: t-test
V_hat = inv((X'*W*X)/n) * ((X'*W*ep*ep'*W*X)/sqrt(n)) * inv((X'*W*X)/n);
t_test = zeros(1,d);
ej     = zeros(d,1);
for iter = 1:d
   ej(iter) = 1;
   t_test(iter) = (ej'*beta_hat - ej'*beta_null) ./ (sqrt(ej'*V_hat*ej)/n);
   ej(iter) = 0;
end

%% 2.4.c: p_values
p_values = normcdf(t_test);
%% 2.4.d: Confidence Intervals
p_values = normcdf(t_test);
ej     = zeros(d,1);
ci     = zeros(2,d);
for iter = 1:d
    ej(iter) = 1;
    ci(1,iter) = ej'*beta_hat - norminv(alpha/2,0,sqrt(ej'*V_hat*ej)/n);
    ci(2,iter) = ej'*beta_hat + norminv(alpha/2,0,sqrt(ej'*V_hat*ej)/n);
    ej(iter) = 0;
end

%% output
output  = {'beta_hat', 'beta_null' ,'t_test', 'p_values', 'CI ('+string(1-alpha)+')'};
for iter = 1:d
output(end+1,:) = {beta_hat(iter), beta_null(iter),t_test(iter), p_values(iter), ci(iter,:)'};
end
output_mat = [beta_hat, beta_null,t_test', p_values', ci]
