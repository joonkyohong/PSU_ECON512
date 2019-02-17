function [W0, W1, W2] = computeW(Vin)
global  trans0 trans1

% This function returns back the continuation value given the ohter
% player's strategy and the market outcome in the current period
% W0: continuation when no one made a sale
% W1: when the player made a sale
% W2: when the opponent made a sale

% Vin: value function, L by L matrix

% Written by Joonkyo Hong
% 16 Feb 2019

 W0 = trans0*Vin*trans0';          % no one made a sale
 W1 = trans1*Vin*trans0';          % the player made a sale
 W2 = trans0*Vin*trans1';          % the opponent made a sale
 
end

