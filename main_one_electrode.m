clc;
clear;
close all;

params = params();

%% Getting the basic operators
matrices = matrices_linearized_t_indep(params);
matrices = @(t) matrices_linearized_t_dep(t,matrices,params);

%% Initial conditions
init = initial_conditions(matrices.adim.ephase.x,matrices.adim.sphase.r,params);
y0 = init.state;

%% Computing the discretized system in state-space form
D = @(t) flatten_matrices(matrices(t),params);
dynamic = @(t,x) (D(t) * x )
