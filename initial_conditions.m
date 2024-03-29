function [ initial ] = initial_conditions( matrices,params )
%INITIAL_CONDITIONS Returns initial conditions that are compatible
% with the boundary conditions.
% IMPORTANT : since the boundary conditions are included in the
% differentiation matrices and the dynamics, it is absolutely necessary to
% initialize the model with an initial condition that satisfies the
% boundary conditions !

%% Dimensionalized Initial conditions
% Part supposed to be modified to fit the desired initial conditions

x = matrices.ephase.x; % Abscissa across the electrode; x(1) = L, x(end) = 0
r = matrices.sphase.r; % Radius across the particle; r(1) = R_s, r(end) = 0

%Initial distributions of the states, in SI units
% Remark : these distributions HAVE TO satisfy the boundary conditions of
%                   the model.
%                   Initializing dl and cs is tricky because of this
%                   constraint. The following script initializes the
%                   electrode in the following state:
%                           - same surface concentration for all particles
%                           - us is exponential along the radius
%                           of a particle
%                           - dl affine in x, with slope current_profile(0)
%                           and origin value given by the OCP taken at the
%                           stoechiometry of the initial surface concentration.
% elyte_init = 1000 * ones(length(x),1); 
elyte_init = 1000 * ones(length(x),1);
cssurf_init = 0.6*params.csmax;
dl_init = ocp_dualfoil(cssurf_init/params.csmax)*ones(length(x),1) + x*current_profile(0);
cs_init = init_cs(dl_init,params);
us_init = zeros(size(cs_init));
for i = 1:size(us_init,1)
    us_init(i,:) = r(i)*cs_init(i,:);
end
% syms coef;
% warning('off','symbolic:solve:FallbackToNumerical');
% for j = 1:length(x)
%     bv = - params.mu*matrices.adim.ephase.dl.butler_volmer(dl_init(j)/params.V_0,ussurf_init);
%     c = solve(coef + coef/(exp(coef) - 1) == ussurf_init + bv,coef);
%     b = bv / (exp(c) - 1);
%     us_init(:,j) = b*(exp(c*r) - 1);
% end

%% MODIFICATION 
% To check positivity of cs depending on IC

for j = 1:length(x)
    us_init(:,j) = r(:)*cssurf_init;
end

%% Adimensionalized and flattened initial conditions
% This should not be modified
% Now that the dimensionalized IC have been specified, we adimensionalize
% them and put them in the right format so the rest of the program can
% handle them.

N_s = params.dscrtzn.N_s;

initial.elyte = elyte_init / (params.aC * params.V_0 * (1 - params.t_0) / params.Faraday);
initial.dl = dl_init / params.V_0;
initial.us = us_init / (params.csmax*params.R_s);
initial.ussurf = initial.us(1,:)';

initial.elyte2N = initial.elyte(2:end-1);
initial.dl2N = initial.dl(2:end-1);
initial.us2N2R = zeros((length(r)-2) * (length(x)-2),1);
for i = 1:length(r)-2
    for j = 1:length(x)-2
        index = (j-1)*(N_s-1) + i;
        initial.us2N2R(index) = initial.us(i+1,j+1);
    end
end
initial.ussurf2N = initial.ussurf(2:end-1);

initial.state = [initial.dl2N;initial.elyte2N;initial.us2N2R;initial.ussurf2N];
end

