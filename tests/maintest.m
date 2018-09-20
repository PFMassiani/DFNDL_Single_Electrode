clc;
clear;
close all;

params = params();

%% Getting the basic operators
matrices = matrices_linearized_t_indep(params);
matrices = matrices_linearized_t_dep(matrices,params);

NTOT = (matrices.N_s + 1)* (matrices.N_elyte - 1) + matrices.N_elyte - 1;
Nem2 = matrices.N_elyte - 1;
Nsm2 = matrices.N_s - 1;

%% Initial conditions
init = initial_conditions(matrices.adim.ephase.x,matrices.adim.sphase.r,params);
y0 = init.state;

%% Computing the discretized system in state-space form
[D,exo_bdry,nonlinear] = flatten_matrices(matrices,params);
dynamic = @(t,x) (full(D(t)) * x  + exo_bdry(t) + nonlinear(t,x));

%% Running the simulation
tf = 0.01;
tspan = [0 tf];
Mass = blkdiag(speye(NTOT - Nem2), speye(Nem2,Nem2));
options = odeset('NormControl','on','MassSingular','yes','Mass',Mass,'AbsTol',1e-9,'RelTol',1e-7);
[t_K,state_K] = ode15s(dynamic,tspan,y0,options);

params.adim.K_n = 0;
matrices = matrices_linearized_t_indep(params);
matrices = matrices_linearized_t_dep(matrices,params);
init = initial_conditions(matrices.adim.ephase.x,matrices.adim.sphase.r,params);
y0 = init.state;
[D,exo_bdry,nonlinear] = flatten_matrices(matrices,params);
dynamic = @(t,x) (full(D(t)) * x  + exo_bdry(t) + nonlinear(t,x));
[t_noK,state_noK] = ode15s(dynamic,t_K,y0,options);
%% Plotting the final concentration

mask_us = logical([zeros(1,2*Nem2), ones(1,Nsm2),zeros(1,NTOT - 2*Nem2 - Nsm2)]);
mask_dl = logical([ones(1,Nem2), zeros(1, NTOT - Nem2)]);
mask_elyte = logical([zeros(1,Nem2), ones(1,Nem2), zeros(1,NTOT - 2*Nem2)]);
mask_ussurf = logical([zeros(NTOT - Nem2,1); ones(Nem2,1)]);
mask_us1 = logical([zeros(2*Nem2+Nsm2,1);ones(Nsm2,1);zeros(NTOT - 2*Nem2 - Nsm2,1)]);
mask_uss = logical([zeros(1,2*Nem2), ones(1,Nsm2*Nem2), zeros(1,NTOT - 2*Nem2-Nsm2*Nem2)]);
%full_state = reconstruct_state(state,matrices);
x = matrices.adim.ephase.x(2:end-1);
r = matrices.adim.sphase.r(2:end-1);

%figure
%plot(matrices.adim.ephase.x(2:end-1),state(end,mask_elyte));
%ylabel('ce')
%figure
%plot(matrices.adim.ephase.x(2:end-1),state(end,mask_dl))
%ylabel('phi')
%figure
%plot(matrices.adim.sphase.r(2:end-1),state(end,mask_us1))
%ylabel('us')
%usquare = zeros(size(x,1),Nsm2,Nem2);
%usline = zeros(size(x,1),Nsm2*Nem2);
%for i = 1:size(x,1)
%    usquare(i,:,:) = x(i,mask_uss);
    %usline(i,:)= reshape(usquare(i,:,:),1,length(usquare(i,:,:)));
%end
figure
for k = 1:length(t_K)
    plot(state_K(k,mask_elyte),'o-'); grid on;
    title(sprintf('time = %0.7f',t_K(k)));
ylabel('With K')
    ylim([0.6 1.2])
    pause(0.001);
end
figure
for k = 1:length(t_noK)
    plot(state_noK(k,mask_elyte) - state_K(k,mask_elyte),'o-'); grid on;
    title(sprintf('time = %0.7f',t_K(k)));
ylabel('Error')
    ylim([-0.3 0.3])
    pause(0.001);
end
all(abs(state_noK(:,mask_elyte) - state_K(:,mask_elyte))<1e-6)

%{
figure
ylim([-0.002 0])
for k = 1:length(t)
    plot(state(k,mask_dl),'o-'); grid on;

    title(sprintf('time = %0.7f',t(k)));
 %   ylabel('dl time')
    ylim([0 5])
    pause(0.001);
end

figure
ylabel('us time')
%ylim([-0.002 0])
for k = 1:length(t)
    plot(state(k,mask_uss),'o-'); grid on;
    title(sprintf('time = %0.7f',t(k)));
ylabel('us time')
%    ylim([0.8 1])
    pause(0.001);
end
%}