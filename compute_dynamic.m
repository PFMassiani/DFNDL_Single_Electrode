function [ dynamic ] = compute_dynamic( t,x,matrices,params )
%COMPUTE_DYNAMIC Summary of this function goes here
%   Detailed explanation goes here
matrices = matrices_linearized_t_dep(t,matrices,params);
[D,exo_bdry,nonlinear] = flatten_matrices(matrices,params);
dynamic = D * x  + exo_bdry + nonlinear(x);

end

