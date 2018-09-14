function [ X ] = chebspace( Xmin,Xmax,N )
%CHEBSPACE Summary of this function goes here
% Returns a set of N+1 points that is distributed linearly on the chebyshev
% nodes.

[~,X] = cheb(N);
X = (Xmin + Xmax) / 2 + (Xmax - Xmin )* X / 2;

end

