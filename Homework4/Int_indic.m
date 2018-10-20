function [ returnVal ] = Int_indic(a, b, Nx, Ny)
% This function integrates 1(x^2 + y^2 <= 1) over [0,1] * [0,1]
% Build upon int_simp.m
% built by Joonkyo Hong 20 Oct. 2018

%Int_indic: Calculate a definite integral of 1(x^2 + y^2 <= 1) over [0,1] * [0,1]
%by using the midpoint rule for N
%evenly spaces segments
%   a - lower bound (must be 2 by 1)
%   b - upper bound (must be 2 by 1)
%   N - Number of segments (must be 2 by 1)

    assert(a(1) < b(1));
    assert(a(2) < b(2));    
    assert(floor(Nx) == Nx);
    assert(floor(Ny) == Ny);    
    assert(Nx > 0);
    assert(Ny > 0);
    
    
    %Make sure Nx and Ny are even:
    if (mod(Nx,2) == 1)||(mod(Ny,2) == 1)
        Nx = Nx+1;
        Ny = Ny+1;
    end

    %Setup Grid: 
    hx = (b(1) - a(1))/Nx;
    X = a(1):hx:b(1);
    
    hy = (b(2) - a(2))/Ny;
    Y = a(2):hy:b(2);    

    %Setup weights
    wx = 2*ones(length(X),1);
    evens = 2:2:Nx;
    wx(evens) = 4;
    wx(1) = 1;
    wx(end) = 1;
    wx=wx*(hx/3);
    
    wy = 2*ones(length(Y),1);
    evens = 2:2:Ny;
    wy(evens) = 4;
    wy(1) = 1;
    wy(end) = 1;
    wy=wy*(hy/3);

    %Evaluate function on grid:
    fv = indic_fcn(X,Y);

    %Sum...
    returnVal = wx'*fv*wy;

end

