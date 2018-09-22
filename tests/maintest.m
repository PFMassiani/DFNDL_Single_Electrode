clc;
clear;
close all;

params = params2();

%% Getting the basic operators
matrices = matrices_linearized_t_indep(params);
matrices = matrices_linearized_t_dep(matrices,params);

NTOT = (matrices.N_s + 1)* (matrices.N_elyte - 1) + matrices.N_elyte - 1;
Nem2 = matrices.N_elyte - 1;
Nsm2 = matrices.N_s - 1;
Ns = matrices.N_s + 1;

%% Initial conditions
init = initial_conditions(matrices,params);
y0 = init.state;

%% Computing the discretized system in state-space form
[D,exo_bdry,nonlinear] = flatten_matrices(matrices,params);
dynamic = @(t,x) (full(D(t)) * x  + exo_bdry(t) + nonlinear(t,x));

%% Running the simulation
tf = 1;
tspan = [0 tf];
Mass = blkdiag(speye(NTOT - Nem2), speye(Nem2,Nem2));
options = odeset('NormControl','on','MassSingular','yes','Mass',Mass,'AbsTol',1e-9,'RelTol',1e-7);
[t_K,state_K] = ode15s(dynamic,tspan,y0,options);

[t,x,r,dl,elyte,us,cs] = postproc(t_K,state_K,params,matrices);
% params.K = 0;
% matrices = matrices_linearized_t_indep(params);
% matrices = matrices_linearized_t_dep(matrices,params);
% init = initial_conditions(matrices.adim.ephase.x,matrices.adim.sphase.r,params);
% y0 = init.state;
% [D,exo_bdry,nonlinear] = flatten_matrices(matrices,params);
% dynamic = @(t,x) (full(D(t)) * x  + exo_bdry(t) + nonlinear(t,x));
% [t_noK,state_noK] = ode15s(dynamic,t_K,y0,options);

%% Plots


figure
for k = 1:length(t)
    plot(x,elyte(k,:),'o-'); grid on;
    title(sprintf('time = %0.7f',t(k)));
    ylabel('CE')
    ylim([0.7e3 1.2e3])
    pause(0.001);
end
figure
for k = 1:length(t)
    plot(x,dl(k,:),'o-'); grid on;
    title(sprintf('time = %0.7f',t(k)));
    ylabel('DL')
      ylim([-1 1])
    pause(0.001);
end
figure
for k = 1:length(t)
    plot(reshape(cs(k,:,:),Nem2*Ns,1),'o-'); grid on;
    title(sprintf('time = %0.7f',t(k)));
    ylabel('CS')
     ylim([-1000 params.csmax])
    pause(0.001);
end

figure
for k = 1:length(t)
    plot(reshape(us(k,:,:),Nem2*Ns,1),'o-'); grid on;
    title(sprintf('time = %0.7f',t(k)));
    ylabel('US')
%     ylim([0.7e3 1.2e3])
    pause(0.001);
end