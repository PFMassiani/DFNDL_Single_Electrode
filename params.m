function [ params ] = params
%Returns the parameters for the DFN model with Double Layer capacitance

%% Copied value from params_LCO.m
% Geometric Params
% Thickness of each layer
p.L_n = 100e-6;     % Thickness of negative electrode [m]
p.L_s = 25e-6;     % Thickness of separator [m]
p.L_p = 100e-6;     % Thickness of positive electrode [m]

L_ccn = 25e-6;    % Thickness of negative current collector [m]
L_ccp = 25e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 10e-6;   % Radius of solid particles in negative electrode [m]
p.R_s_p = 10e-6;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.6;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.5;      % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 0.3;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 1.0;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.3;   % Volume fraction in electrolyte for pos. electrode

% make element to caclulate phi_{s} by Saehong Park 
p.epsilon_f_n = 0.1;  % Volume fraction of filler in neg. electrode
p.epsilon_f_p = 0.2;  % Volume fraction of filler in pos. electrode

epsilon_f_n = p.epsilon_f_n;  % Volume fraction of filler in neg. electrode
epsilon_f_p = p.epsilon_f_p;  % Volume fraction of filler in pos. electrode


% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]

% Mass densities
rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]
rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]
rho_e =  1324;    % Electrolyte [kg/m^3]
rho_f = 1800;     % Filler [kg/m^3]
rho_ccn = 8954;   % Current collector in negative electrode
rho_ccp = 2707;   % Current collector in positive electrode

% Compute cell mass [kg/m^2]
m_n = p.L_n * (rho_e*p.epsilon_e_n + rho_sn*p.epsilon_s_n + rho_f*epsilon_f_n);
m_s = p.L_s * (rho_e*p.epsilon_e_n);
m_p = p.L_p * (rho_e*p.epsilon_e_p + rho_sp*p.epsilon_s_p + rho_f*epsilon_f_p);
m_cc = rho_ccn*L_ccn + rho_ccp*L_ccp;

% Lumped density [kg/m^2]
p.rho_avg = m_n + m_s + m_p + m_cc;

% Transport Params
% Diffusion coefficient in solid
p.D_s_n0 = 3.9e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p0 = 1e-13;  % Diffusion coeff for solid in pos. electrode, [m^2/s]

% Diffusional conductivity in electrolyte
p.dactivity = 0;

p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 10;    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% Miscellaneous
p.t_plus = 0.4;       % Transference number
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 1;           % Electrode current collector area [m^2]

% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]

p.alph = 0.5;         % Charge transfer coefficients

p.R_f_n = 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_c = 0;         % Contact Resistance/Current Collector Resistance, [Ohms-m^2]

% Nominal Reaction rates
p.k_n0 = 1e-5;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 3e-7; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

% Thermodynamic Params

% Thermal dynamics
p.C_p = 2000;   % Heat capacity, [J/kg-K]
p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0

% Ambient Temperature
p.T_amb = 298.15; % [K]

% Activation Energies
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% All units are [J/mol]
p.E.kn = 37.48e3;
p.E.kp = 39.57e3;
p.E.Dsn = 42.77e3;
p.E.Dsp = 18.55e3;
p.E.De = 37.04e3;
p.E.kappa_e = 34.70e3;

% Reference temperature
p.T_ref = 298.15; %[K]

% Concentrations
% Maxima based on DUALFOIL 
% line 588 in DUALFOIL Fortran code

p.c_s_n_max = 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]
%p.c_s_n_max = 3.6e3 * 372 * 2260 / p.Faraday;   % Max concentration in anode, [mol/m^3]

%p.c_s_p_max = 3.6e3 * 247 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]
p.c_s_p_max = 3.6e3 * 274 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.n_Li_s = 2.5; %2.781;        % Total moles of lithium in solid phase [mol]
p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]

% Cutoff voltages
p.volt_max = 4.7; %4.1113; %4.7;
p.volt_min = 3.105; %2.6;

% Discretization parameters
% Discrete time step
p.delta_t = 1;

% Pade Order
p.PadeOrder = 3;

% Finite difference points along x-coordinate
p.Nxn = 70;
p.Nxs = 35;
p.Nxp = 70;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

% Modifications
p.R_f_n = 1000 * p.R_f_n;

params = p;
%% Custom parameters
%Linearization point
params.adim.c_0 = 1; % TO CHOOSE

% Conductivity
params.kappa_n = (1.5 + 0.1) / 2; % Mean value of the conductivity of electrolyte [S/m]
params.kappa_p = (1.5 + 0.1) / 2; % Mean value of the conductivity of electrolyte [S/m]

% Diffusivity
params.D_e_n0 = (3.5 + 2.4)*10^(-10) / 2; % Mean value of the diffusivity, negative electrode [m^2/s]
params.D_e_p0 = (3.5 + 2.4)*10^(-10) / 2; % Mean value of the diffusivity, positive electrode [m^2/s]

% Double layer volumic capacitance
% Value is taken from Drummond & al., "State-of-Charge Estimation with a
% Doyle-Fuller-Newman Li-ion Battery Model"
params.C_n = 22e3 / params.a_s_n; % Negative electrode [F/m^2]
params.C_p = 22e3 / params.a_s_p; % Positive electrode [F/m^2]

% Linearized open circuit potential constant
params.eta = 1; % TODO : needs to be specified

% Exchange current density
params.i_0 = 1; % TODO : needs to be specified [A/m^2]

% Typical voltage
params.V_0_n = params.R_s_n^2 * params.i_0 / (params.D_s_n0 * params.C_n); % Negative electrode [V]
params.V_0_p = params.R_s_p^2 * params.i_0 / (params.D_s_p0 * params.C_p); % Positive electrode [V]

%Adimensional coefficients
params.adim.mu_n = params.eta * params.C_n / (params.R_s_n * params.Faraday); % Adimensional linearized open circuit potential, Negative electrode
params.adim.mu_p = params.eta * params.C_p / (params.R_s_p * params.Faraday); % Adimensional linearized open circuit potential, Positive electrode
params.adim.E_n = params.Faraday * params.V_0_n / (params.R * params.T_amb); % Adimensional energy, Negative electrode
params.adim.E_p = params.Faraday * params.V_0_p / (params.R * params.T_amb); % Adimensional energy, Positive electrode
params.adim.theta_c_n = (params.R_s_n^2 / params.D_s_n0) * (params.sig_n + params.kappa_n) / (params.sig_n * params.kappa_n) * 1 / (params.a_s_n * params.C_n); % Ratio comparing solid diffusion characteristic time and double layer dynamics characteristic time, Negative electrode
params.adim.theta_c_p = (params.R_s_p^2 / params.D_s_p0) * (params.sig_p + params.kappa_p) / (params.sig_p* params.kappa_p) * 1 / (params.a_s_p * params.C_p); % Ratio comparing solid diffusion characteristic time and double layer dynamics characteristic time, Positive electrode
params.adim.theta_d_n = (params.D_e_n0 / params.D_s_n0) * params.R_s_n^2 / params.L_n^2; % Ratio comparing solid diffusion and electrolyte diffusion characteristic times, Negative electrode
params.adim.theta_d_p = (params.D_e_p0 / params.D_s_p0) * params.R_s_p^2 / params.L_p^2; % Ratio comparing solid diffusion and electrolyte diffusion characteristic times, Positive electrode
params.adim.theta_f_n = (params.D_s_n0 / params.R_s_n^2)*params.R_f_n * params.C_n; % Ratio comparing solid diffusion and filme resistance characteristic times, negative electrode
params.adim.theta_f_p = (params.D_s_p0 / params.R_s_p^2)*params.R_f_p * params.C_p; % Ratio comparing solid diffusion and filme resistance characteristic times, positive electrode
params.adim.rho_n = params.adim.theta_d_n / params.adim.theta_c_n; % Rho coefficient, negative electrode
params.adim.rho_n = params.adim.theta_d_p / params.adim.theta_c_p; % Rho coefficient, positive electrode
%% Modifications
params.adim.E_n = params.adim.E_n * 10^(-5);
%% End of modifs
params.adim.K_n = 2*(1-params.t_plus)*(1+0)/params.adim.E_n; %TODO

%Charging profile
params.adim.i = @(t) (1e13* sin(1000*2*pi*t) * params.L_n/(params.V_0_n * params.sig_n)); % TO CHOOSE

%% Simulation parameters
params.dscrtzn.N_e_n = 10; % Number of nodes in the electrolyte (negative electrode)
params.dscrtzn.N_s_n = 5;% Number of nodes in the solid phase (negative electrode)
params.dscrtzn.N_e_p = 10;% Number of nodes in the electrolyte (positive electrode)
params.dcrtzn.N_s_p = 10;% Number of nodes in the solid phase (positive electrode)

params.misc.xmin = -0.0015;
params.misc.xmax = 0.0015;


function y =door(t)
if t<10
    y=1e11;
else
    y=0;
end
end
end