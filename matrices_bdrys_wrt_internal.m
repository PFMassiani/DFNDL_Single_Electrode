function [ Dmod,Emod,Umod,Vmod ] = matrices_bdrys_wrt_internal( D,E,U,V,bdry_coef )
%MATRICES_BDRYS_WRT_INTERNAL Summary of this function goes here
%   Assume the boundary conditions are of the form (the state is F(1:N+1), and U and V are state-independant scalar forcing terms) :
%           D * F = U + bdry_coef(1) * F(1)
%           E * F = V + bdry_coef(2) * F(N+1)
%   This function returns the line-vectors Dmod and Emod and the modified
%   scalar controls Umod and Vmod such that:
%       F(1) = Umod + Dmod * F(2:N)
%       F(N+1) = Vmod + Emod * F(2:N)

Dmod(:) = ( (D(end)/(E(end) - bdry_coef(2)))*E(2:end - 1) - D(2:end-1)) / (D(1) - bdry_coef(1) - D(end)*E(1)/(E(end) - bdry_coef(2) ) );
Emod(:) = ( (E(1)/(D(1) - bdry_coef(1)) ) * D(2:end - 1) - E(2:end-1) ) / (E(end) - bdry_coef(2) - (E(1)*D(end))/(D(1) - bdry_coef(1) ) );

Umod_left = U / ( D(1) - bdry_coef(1) - D(end)*E(1)/( E(end) - bdry_coef(1) ) );
Vmod_left = - V * (D(end)/(E(end) - bdry_coef(2) ) ) / ( D(1) - bdry_coef(1) - D(end)*E(1)/( E(end) - bdry_coef(1) ) );

Umod_right = - U * (E(1) / ( D(1) - bdry_coef(2) ) ) / ( E(end) - bdry_coef(2) - E(1)*D(end)/(D(1) - bdry_coef(1) ) );
Vmod_right = V / ( E(end) - bdry_coef(2) - E(1)*D(end)/(D(1) - bdry_coef(1) ) );

Umod = Umod_left + Vmod_left;
Vmod = Umod_right + Vmod_right;
end

