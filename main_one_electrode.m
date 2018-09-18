clc;
clear;
close all;

params = params();

%% Getting the basic operators
matrices = matrices_linearized_t_indep(params);
matrices = matrices_linearized_t_dep(matrices,params);

%% Initial conditions
init = initial_conditions(matrices.adim.ephase.x,matrices.adim.sphase.r,params);
y0 = init.state;

%% Computing the discretized system in state-space form
[D,exo_bdry,nonlinear] = flatten_matrices(matrices,params);
dynamic = @(t,x) (D(t) * x  + exo_bdry(t) + nonlinear(t,x));

%% Running the simulation
tf = 60;
tspan = [0 tf];
[t,state] = ode15s(dynamic,tspan,y0);

%% Plotting the final concentration
NTOT = (matrices.N_s + 1)* (matrices.N_elyte - 1);
Nem2 = matrices.N_elyte - 1;
Nsm2 = matrices.N_s - 1;
mask_us = logical([zeros(1,2*Nem2), ones(1,Nsm2),zeros(1,NTOT - 2*Nem2 - Nsm2)]);
mask_dl = logical([ones(1,Nem2), zeros(1, NTOT - Nem2)]);
mask_elyte = logical([zeros(1,Nem2), ones(1,Nem2), zeros(1,NTOT - 2*Nem2)]);

x = matrices.adim.ephase.x(2:end-1);
r = matrices.adim.sphase.r(2:end-1);
tplot = size(state,1);
figure
plot(matrices.adim.ephase.x(2:end-1),state(end,mask_elyte));