%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework #3 ECON 512                                    %
% Written by Joonkyo (Jay) Hong, 4 Oct 2018               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% Load the data
load hw3.mat

% Set the options
options_nm = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000,'TolFun',10e-6,'TolX',10e-6);
options_qn = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
options_nls = optimoptions(@lsqnonlin,'Display','off','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',10000,'MaxIterations',5000,'TolX',10e-8);

% Initial Value: OLS estimate for log(y+1) = X*beta + u
% Of course this is bad specifiction when the DGP is true because
% E(log(y)|x) =/= log E(y|x)

% But it provides us with somewhat "reasonable" initial value for
% MLE,NLS,GMM...

init_theta = (X'*X)\X'*log(1+y);

%% Questions 1 and 2 (MLE)

% Negative Likelihood Function
lik_fcn = @(theta) -sum(-exp(X*theta)+y.*(X*theta));

% Nelder-Mead
tic
[thetamle_nm,~,~,~] = fminsearch(lik_fcn,init_theta,options_nm);
t_mlenm = toc;

% Quasi-Newton
tic
[thetamle_qn,~,~,~] = fminunc(lik_fcn,init_theta,options_qn);
t_mleqn = toc;

%% Questions 3 and 4 (NLS)

% Nonlinear Least Square objective function
obj_fcn_lsq = @(theta) (y-exp(X*theta));

% Least Square Non-linear 
tic
thetanls_lq = lsqnonlin(obj_fcn_lsq,init_theta,[],[],options_nls);
t_nlslq = toc;

% Nelder-Mead
obj_fcn = @(theta) (y-exp(X*theta))'*(y-exp(X*theta));
tic
[thetanls_nm,~,~,~] = fminsearch(obj_fcn,init_theta,options_nm);
t_nlsnm = toc;

%% Question 5 (Check the Robustness to Initial Value)

initbeta0 = 0:0.25:4;
initbeta6 = -0.5:0.05:0;

M = length(initbeta0);
N = length(initbeta6); 
mle_nm = zeros(M,N);
mle_qn = zeros(M,N);
nls_lq = zeros(M,N);
nls_nm = zeros(M,N);

    for m=1:M
         for n=1:N
         
              init_theta = [initbeta0(m);init_theta(2:5,1);initbeta6(n)];
              
              x = fminsearch(lik_fcn,init_theta,options_nm);
              mle_nm(m,n) = (x-thetamle_nm)'*(x-thetamle_nm)/6;
              
              x = fminunc(lik_fcn,init_theta,options_qn);
              mle_qn(m,n) = (x-thetamle_qn)'*(x-thetamle_qn)/6;
              
              x = lsqnonlin(obj_fcn_lsq,init_theta,[],[],options_nls);
              nls_lq(m,n) = (x-thetanls_lq)'*(x-thetanls_lq)/6;
              
              x = fminsearch(obj_fcn,init_theta,options_nm);
              nls_nm(m,n) = (x-thetanls_nm)'*(x-thetanls_nm)/6;
                  
         end         
     end
 

robustness_measure = 10e-10*(1./[mean(mean(mle_nm)) mean(mean(mle_qn)) mean(mean(nls_lq)) mean(mean(nls_nm))]);
convergence_time = [t_mlenm t_mleqn t_nlslq t_nlsnm];

 disp("              ");
 disp("RESULT");
 disp("Parameter estimates of each method");
 disp("MLE_NM        MLE_QN        NLS_LQ      NLS_NM");
 disp(num2str([thetamle_nm thetamle_qn thetanls_lq thetanls_nm]));
  disp("              ");
 disp("Robustness Measure of each method");
 disp("MLE_NM        MLE_QN        NLS_LQ      NLS_NM");
 disp(num2str(robustness_measure));
 disp("              ");
 disp("Time to Convergence of each method");
 disp("MLE_NM        MLE_QN        NLS_LQ      NLS_NM");
 disp(num2str(convergence_time));


