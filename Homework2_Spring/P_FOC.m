function z = P_FOC(Vin, Pin)
global L betta cost 

% This function returns back the First Order Condition from Bellmann Eq 

% Written by Joonkyo Hong
% 16 Feb 2019


P1 = Pin;
P2 = Pin';

[D0, D1, D2] = computeD(P1,P2);
[W0, W1, W2] = computeW(Vin);

z = betta.*(D0.*W0 + D1.*W1 + D2.*W2) - betta*W1 - (1-D1).*(P1 - repmat(cost,1,L)) + 1;

end

