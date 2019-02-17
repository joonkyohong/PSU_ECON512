function Vout = updateV(Vin, Pnew, Pold)
global L betta cost 

% This function returns back the updated value function through the Bellmann
% Equation
% Vout: newly updated value function

% Vin: value function, L by L matrix
% Pnew: updated policy function, L by L matrix
% Pold: previous policy function, L by L matrix

% Written by Joonkyo Hong
% 16 Feb 2019

 P1 = Pnew;
 P2 = Pold';
 [D0, D1, D2] = computeD(P1,P2);
 [W0, W1, W2] = computeW(Vin);
 
    Vout = D1.*(P1-repmat(cost,1,L)) + betta.*(D0.*W0 + D1.*W1 + D2.*W2);
 
end

