%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework #1 ECON 512                                    %
% Written by Joonkyo (Jay) Hong, 22 Aug 2018              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1.

X = [1 1.5 3 4 5 6 7 9 10];
Y1 = -2 + .5*X;
Y2 = -2 + .5*X.^2;

figure(1)
plot(X,[Y1; Y2]);
legend("Y1","Y2");
xlabel("X");

%% Problem 2

vec_problem2 = linspace(-10,20,200);
vec_problem2 = vec_problem2';
ans_problem2 = sum(vec_problem2);

%% Problem 3

A = [2 4 6;
     1 7 5;
     3 12 4];
b = [-2;3;10];

C = A'*b;
D = (A'*A)\b;
E = b'*A*[1;1;1];
F = A;
F(:,3)=[];
F(2,:)=[];
x = A\b;

%% Problem 4

B = kron(eye(5),A);

%% Problem 5

matrix_problem5 = normrnd(10,5,5,3);
ans_problem5 = matrix_problem5;
ans_problem5(ans_problem5<10)=0;
ans_problem5(ans_problem5>=10)=1;

%% Problem 6

dataset = csvread('datahw1.csv',0,0);

ymat = (dataset(:,5));
xmat = [ones(length(ymat),1) dataset(:,3:4) (dataset(:,6))];

ols_est = xmat'*xmat\xmat'*ymat;
k = length(ols_est);
emat = ymat-xmat*ols_est;
wave = 4;
obs= length(ymat)/4;
center_sand = zeros(k,k);

for i=1:obs
    ei = emat((i-1)*wave+1:i*wave,1); xi = xmat((i-1)*wave+1:i*wave,:);
    center_sand = center_sand + xi'*ei*ei'*xi;
end

ols_cov = (xmat'*xmat)\center_sand/(xmat'*xmat);
ols_se = sqrt(diag(ols_cov));

 disp("              ");
 disp("OLS estimates in Problem 6");
 disp("Parameter Estimates and Standard Errors");
 disp(" beta0        beta1       beta2      beta3");
 disp(num2str(ols_est'));
 disp(num2str(ols_se'));
 
 