function [ initial ] = initial_conditions( matrices,params )
%INITIAL_CONDITIONS Returns initial conditions that are compatible
% with the boundary conditions.
% IMPORTANT : since the boundary conditions are included in the
% differentiation matrices and the dynamics, it is absolutely necessary to
% initialize the model with an initial condition that satisfies the
% boundary conditions !
x = matrices.adim.ephase.x;
r = matrices.adim.sphase.r;

N_elyte = params.dscrtzn.N_e;
N_s = params.dscrtzn.N_s;

elyte0 = 1000; % [mol/m^3]
ussurf0 = 10000*params.R_s;
i_initial = params.i_dim(0);

elyte0_adim = elyte0 / (params.aC * params.V_0 * (1 - params.t_0) / params.Faraday);
ussurf0_adim = ussurf0 / (params.OCP_slope / params.V_0);
dl0_adim = ussurf0_adim;
i_initial_adim = params.i(0);

% initial.elyte = elyte0 * ones(N_elyte+1,1);
initial.elyte = elyte0_adim*(1+ (x.^3)/3 - (x.^2)/2);
initial.dl = x *i_initial_adim + dl0_adim;
%initial.dl = cos(2*pi*x);
initial.elyte2N = initial.elyte(2:end-1);
initial.dl2N = initial.dl(2:end-1);

initial.us2N2R = zeros((length(r)-2) * (length(initial.dl)-2),1);
for i = 1:length(r)-2
    for j = 1:length(initial.dl)-2
        index = (j-1)*(N_s-1) + i;
        initial.us2N2R(index) = initial.dl(j+1) * r(i+1);
%         initial.us2N2R(index) = (initial.dl(j+1)-3)*r(i+1) - 2*r(i+1)^2 + r(1+i)^3;
    end
end
%for j = 1:length(initial.dl)-2
%    range = (j-1)*length(r(2:end-1)) + 1:(length(r)-2);
%    initial.us2N2R(range) = r(2:end-1)*initial.dl(j+1);
%end
initial.ussurf = initial.dl2N;
initial.state = [initial.dl2N;initial.elyte2N;initial.us2N2R;initial.ussurf];
end

