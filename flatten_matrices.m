function [ D,exo_bdry,nonlinear ] = flatten_matrices( matrices,params )
%FLATTEN_MATRICES Summary of this function goes here
%   The natural structure of the problem (2 spatial dimensions) is not
%   permitted by the MATLAB solver. For instance, the state variable us
%   cannot be a matrix : it needs to be a column vector. This function
%   takes in parameter the operators following the natural structure of the
%   problem, and returns a matrix that can be used in a solver.

%% Constants definition
NTOT = (matrices.N_s + 1)* (matrices.N_elyte - 1);
Nem2 = matrices.N_elyte - 1;
Nsm2 = matrices.N_s - 1;

%% State masks

mask.elyte = logical([zeros(Nem2,1);
                                        ones(Nem2,1);
                                        zeros(NTOT - 2*Nem2,1)]);
mask.dl(:,:) = logical([eye(Nem2);
                                          zeros(NTOT-Nem2,Nem2)]);
mask.all_dl = logical([ones(Nem2,1);
                                          zeros(NTOT-Nem2,1)]);
mask.us2R(:,:) = logical([zeros(2*Nem2,Nem2);
                                            kron(eye(Nem2),ones(Nsm2,1))]);

%% Linear dyanmics matrix
base.dl =@(t) ([params.adim.theta_c_n * sparse(matrices.adim.ephase.dl.D2_BC(t)) sparse(Nem2,NTOT - Nem2)]);
base.elyte = @(t) ([params.adim.theta_c_n * sparse(matrices.adim.ephase.dl.D2_BC(t))  params.adim.theta_d_n * sparse(matrices.adim.ephase.elyte.D2_BC) sparse(Nem2, NTOT - 2*Nem2)]);
base.sphase = [sparse(NTOT - 2*Nem2,2*Nem2), kron(speye(Nem2), matrices.adim.sphase.us.D2_BC)];

indexlist_linear_dl_us = zeros(Nem2,2);

for i = 1:size(indexlist_linear_dl_us,1)
    indexlist_linear_dl_us(i,1) = (Nsm2 -1)*(i - 1) + i;
    indexlist_linear_dl_us(i,2) = i;
end
for i = 1:size(indexlist_linear_dl_us,1)
    i_real = indexlist_linear_dl_us(i,1);
    j_real = indexlist_linear_dl_us(i,2);
    base.sphase(i_real:i_real+Nsm2-1,j_real) = matrices.adim.sphase.us.exogeneous_linear_dl(:);
end

D = @(t) ([base.dl(t);
        base.elyte(t);
        base.sphase]);
        
%% Nonlinear dynamics terms
% Bulter-Volmer term in dl's dynamics
% Computed in one line in the synthesis

% Logarithmic term in dl's and elyte's dynamics
% Computed in one line in the synthesis


%% Source terms because of boundaries

% Linear boundary conditions
exo_bdry = @(t) ([params.adim.theta_c_n * matrices.adim.ephase.dl.exogeneous_bdry(t);
                                params.adim.theta_d_n * matrices.adim.ephase.elyte.exogeneous_bdry(t);
                                sparse(NTOT - 2*Nem2,1)]);

% Nonlinear boundary condition : Butler-Volmer term in us dynamics
nonlinear_us = @(t,x) (- apply_elementwise_masked( matrices.adim.sphase.us.exogeneous_nonlinear, Nsm2, x, mask.dl,mask.us2R ));

%% Nonlinear terms synthesis

nonlinear = @(t,x) ([params.adim.K_n * params.adim.theta_c_n * matrices.adim.ephase.elyte.logterm_BC(t,x(mask.elyte));
                                      - matrices.adim.ephase.dl.butler_volmer(x(mask.all_dl),compute_ussurf2N(x,matrices,mask)) + params.adim.K_n * params.adim.theta_c_n * matrices.adim.ephase.elyte.logterm_BC(t,x(mask.elyte));
                                      nonlinear_us(t,x)]);
        
end

function y = apply_elementwise_masked(f,dimOut_f,x,mask1,mask2)
    y = zeros(length(f)*dimOut_f,1);
    size(mask1);
    for i = 1:size(mask1,2)
        p = (i-1)*(dimOut_f) + (1:dimOut_f);
        y(p) = feval(f,x(mask1(:,i)),x(mask2(:,i)));
    end
end