function [ D,exo_bdry,nonlinear ] = flatten_matrices( matrices,params )
%FLATTEN_MATRICES Summary of this function goes here
%   The natural structure of the problem (2 spatial dimensions) is not
%   permitted by the MATLAB solver. For instance, the state variable us
%   cannot be a matrix : it needs to be a column vector. This function
%   takes in parameter the operators following the natural structure of the
%   problem, and returns a matrix that can be used in a solver.

%% Constants definition
NTOT = (matrices.N_s + 1)* (matrices.N_elyte - 1) + matrices.N_elyte - 1;
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
                                            kron(eye(Nem2),ones(Nsm2,1));
                                            zeros(NTOT - 2*Nem2 - Nsm2*Nem2)]);
mask.ussurf = logical([zeros(NTOT - Nem2,1);
                                            ones(Nem2,1)]);
%% Linear dynamics matrix
base.dl = [params.theta_c * sparse(matrices.adim.ephase.dl.D2_BC) sparse(Nem2,NTOT - Nem2)];
base.elyte = [params.theta_c * sparse(matrices.adim.ephase.dl.D2_BC)  params.theta_d * sparse(matrices.adim.ephase.elyte.D2_BC) sparse(Nem2, NTOT - 2*Nem2)];
base.us = [sparse(NTOT - 3*Nem2,2*Nem2), kron(speye(Nem2), matrices.adim.sphase.us.D2_BC), sparse(NTOT - 3* Nem2, NTOT - 2*Nem2 - Nsm2*Nem2);];
base.ussurf = [sparse(Nem2,2*Nem2), sparse(matrices.adim.sphase.ussurf.linear_term_us), sparse(matrices.adim.sphase.ussurf.linear_term_ussurf)];

%indexlist_linear_dl_us = zeros(Nem2,2);

%for i = 1:size(indexlist_linear_dl_us,1)
%    indexlist_linear_dl_us(i,1) = (Nsm2 -1)*(i - 1) + i;
%    indexlist_linear_dl_us(i,2) = i;
%end
%for i = 1:size(indexlist_linear_dl_us,1)
%    i_real = indexlist_linear_dl_us(i,1);
%    j_real = indexlist_linear_dl_us(i,2);
%    base.sphase(i_real:i_real+Nsm2-1,j_real) = matrices.adim.sphase.us.exogeneous_linear_dl(:);
%end

D = [base.dl;
        base.elyte;
        base.us;
        base.ussurf];
        
%% Nonlinear dynamics terms
% Bulter-Volmer term in dl's dynamics
% Computed in one line in the synthesis

% Logarithmic term in dl's and elyte's dynamics
% Computed in one line in the synthesis


%% Source terms because of boundaries

% Exogeneous boundary conditions
exo_bdry = [params.theta_c * matrices.adim.ephase.dl.exogeneous_bdry;
                                params.theta_c * matrices.adim.ephase.dl.exogeneous_bdry+params.theta_d * matrices.adim.ephase.elyte.exogeneous_bdry;
                                sparse(NTOT - 2*Nem2,1)];

% Nonlinear boundary condition : Butler-Volmer term in us dynamics
% Computed in one line in the synthesis

%% Nonlinear terms synthesis

nonlinear = @(x) ([-matrices.adim.ephase.dl.butler_volmer(x(mask.elyte),x(mask.all_dl),x(mask.ussurf)) + params.K * params.theta_c * matrices.adim.ephase.elyte.logterm_BC(x(mask.elyte));
                                       params.K * params.theta_c * matrices.adim.ephase.elyte.logterm_BC(x(mask.elyte));
                                      reshape(matrices.adim.sphase.us.exogeneous_nonlinear(x(mask.ussurf)), Nsm2 * Nem2, 1);%flatten(matrices.adim.sphase.us.exogeneous_nonlinear(x(mask.ussurf)));
                                      params.mu / (1 - matrices.adim.sphase.D1_noBC(1)) * matrices.adim.ephase.dl.butler_volmer(x(mask.elyte),x(mask.all_dl),x(mask.ussurf))]);
        
end

function flat = flatten(M)
    flat = zeros(size(M,1)*size(M,2),1);
    for i = 1:size(M,2)
        p = (i-1)*size(M,1) + (1:size(M,1));
        flat(p) = M(:,i);
    end
end
%function y = apply_elementwise_masked(f,dimOut_f,x,mask1,mask2)
%    y = zeros(length(f)*dimOut_f,1);
%    size(mask1);
%    for i = 1:size(mask1,2)
%        p = (i-1)*(dimOut_f) + (1:dimOut_f);
%        y(p) = feval(f,x(mask1(:,i)),x(mask2(:,i)));
%    end
%end