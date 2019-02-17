function [D0, D1, D2] = computeD(P1,P2)
global L nu 

% This function returns back the market shares of each product
% given the strategy profile P1, P2

% P1 and P2: policy functions, L by L matrices
% Written by Joonkyo Hong
% 16 Feb 2019

 temp = 1 + exp(nu - P1) + exp(nu - P2);
 D0   = ones(L,L)./temp;
 D1   = exp(nu-P1)./temp;
 D2   = exp(nu-P2)./temp;
 
end

