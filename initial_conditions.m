function [ initial ] = initial_conditions( x,r,params )
%INITIAL_CONDITIONS Returns initial conditions that are compatible
% with the boundary conditions.
% IMPORTANT : since the boundary conditions are included in the
% differentiation matrices and the dynamics, it is absolutely necessary to
% initialize the model with an initial condition that satisfies the
% boundary conditions !

N_elyte = params.dscrtzn.N_e_n;
N_s = params.dscrtzn.N_s_n;

elyte0 = 1;
i0 = params.adim.i(0);

initial.elyte = elyte0 * ones(N_elyte+1,1);
initial.dl = squeeze(x .* (x - 1) * i0);

initial.elyte2N = initial.elyte(2:end-1);
initial.dl2N = initial.dl(2:end-1);

initial.us2N2R = zeros((length(r)-2) * (length(initial.dl)-2),1);
for i = 1:length(r)-2
    for j = 1:length(initial.dl)-2
        index = (j-1)*(N_s-1) + i;
        initial.us2N2R(index) = r(i) * initial.dl(j);
    end
end
initial.state = [initial.dl2N;initial.elyte2N;initial.us2N2R];
end

