function [ p ] = params2
%PARAMS2 Summary of this function goes here
%   Detailed explanation goes here

% Values taken from params_LCO.m on Github - FastDFN
p.L = 100e-6;
p.R_s = 10e-6;
p.epsilon_e = 0.3;
p.epsilon_s = 0.6;
p.a_s = 3*p.epsilon_s / p.R_s;
p.D_s = 3.9e-14;
p.sigma = 100;
p.t_0 = 0.4;
p.Faraday = 96487;
p.R = 8.314472;
p.alpha = 0.5; 
p.T = 298.15;
p.csmax =  3.6e3 * 372 * 1800 / p.Faraday;
p.volt_max = 4.7; %4.1113; %4.7;
p.volt_min = 3.105; %2.6;
p.k = 1e-5;


% Values taken from the simulations sent by Saehong
p.kappa = (1.5 + 0.1) / 2; % Mean value of the conductivity of electrolyte [S/m]
p.D_e = (3.5 + 2.4)*10^(-10) / 2; % Mean value of the diffusivity [m^2/s]

% Value taken from "State-of-Charge Estimation with a
% Doyle-Fuller-Newman Li-ion Battery Model", R. Drummond
p.aC = 22*1e3; % [F/m^2]


%Values to find
p.dlnfdlnce = 0;
p.i_0_ref = 1;%((1-p.t_0)/p.Faraday)^(p.alpha/(1 - p.alpha)) * p.csmax^(2*p.alpha/(1 - p.alpha)) * (p.a_s*p.R_s^2/p.D_s)^(1/(1 - p.alpha));

%Computing OCP characteristics
% p.OCP_0 = 4;
% p.beta = 0.1; % "Identifiability and Parameter Estimation of the Single Particle Lithium-Ion Battery Model", David Howey
% p.csmax = 30000; % [mol/m^3], "Identifiability and Parameter Estimation of the Single Particle Lithium-Ion Battery Model", David Howey
% p.A = p.epsilon_s * p.L / (4*pi*p.R_s^3 / 3);
% p.Qth = p.epsilon_s * p.L * p.csmax * p.Faraday * p.A;
% p.OCP_slope_wrt_stoechio = -p.epsilon_s * p.L * p.csmax * p.Faraday * p.A * p.beta;% "Identifiability and Parameter Estimation of the Single Particle Lithium-Ion Battery Model", David Howey
% p.OCP_slope = p.OCP_slope_wrt_stoechio / (p.R_s * p.csmax);
%p.OCP_slope = -1e4;

% Parameters appearing in the system
p.V_0 = ((1 - p.t_0)/p.Faraday)^(p.alpha/(1 - p.alpha)) * p.csmax^(2*p.alpha/(1 - p.alpha)) * (p.a_s * p.R_s^2/p.D_s)^(1/(1-p.alpha));
p.V_0 = (p.Faraday / p.aC) * (p.k * (p.R_s ^2 / p.D_s) * (p.a_s*(p.csmax^2)*(1 - p.t_0))^p.alpha)^(1/(1 - p.alpha));
p.theta_c = (p.R_s^2 * p.sigma * p.kappa ) / (p.L^2*p.aC * (p.sigma + p.kappa) * p.D_s);
p.theta_d = (p.D_e * p.R_s^2) / (p.D_s * p.L^2 * p.epsilon_e);
p.E = p.Faraday * p.V_0 / (p.R * p.T);
%p.E = p.E *10^(-floor(log10(p.E)));
p.K = 1e9*2 * (1 - p.t_0) * (1 + p.dlnfdlnce) / p.E;
p.mu = p.R_s^2*p.i_0_ref/(p.D_s*p.Faraday);
% p.i_0_ref = p.k * p.Faraday * (p.Faraday / ((1 - p.t_0) * p.V_0 * p.aC)) / p.csmax; % Expression to check
p.i_0_ref = p.aC * p.D_s * p.V_0/ (p.a_s * p.R_s^2);
%Current profile
p.i_dim = @(t) (10); % [A]
p.i = @(t) (p.i_dim(t * p.R_s^2 / p.D_s) * p.L / (p.sigma * p.V_0));

%% Simulation parameters
p.dscrtzn.N_e = 10; % Number of nodes in the electrolyte
p.dscrtzn.N_s = 5;% Number of nodes in the solid phase
end

