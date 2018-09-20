function [ matrices ] = matrices_linearized_t_indep( params )
%MATRICES Summary of this function goes here
%   Detailed explanation goes here

matrices.N_s = params.dscrtzn.N_s_n;
matrices.N_elyte = params.dscrtzn.N_e_n;
L = params.L_n;
Rs = params.R_s_n;

%% Scaled differentiation operators (no Boundary Conditions)

[ephase.D1_cheb,ephase.x_cheb] = cheb(matrices.N_elyte);
[sphase.D1_cheb,sphase.r_cheb]  = cheb(matrices.N_s);

%       For the dimensional system
matrices.ephase.D1_noBC = (2/L)*ephase.D1_cheb;
matrices.sphase.D1_noBC = (2/Rs)*sphase.D1_cheb;
matrices.ephase.x = (1 + ephase.x_cheb)*L/2;
matrices.sphase.r = (1 + sphase.r_cheb) * Rs/2;

%       For the adimensional system
matrices.adim.ephase.D1_noBC = 2 * ephase.D1_cheb;
matrices.adim.sphase.D1_noBC = 2* sphase.D1_cheb;
matrices.adim.ephase.D2_noBC = matrices.adim.ephase.D1_noBC ^ 2;
matrices.adim.sphase.D2_noBC = matrices.adim.sphase.D1_noBC ^ 2;

matrices.adim.ephase.x = (1 + ephase.x_cheb)/2;
matrices.adim.sphase.r = (1 + sphase.r_cheb)/2;

%% Boundary conditions : electrolyte phase

% Getting the modifications that should be made to the _noBC matrices to get
% the _BC matrices :
[matrices.adim.ephase.elyte.Dmod,matrices.adim.ephase.elyte.Emod,matrices.adim.ephase.elyte.Umod,matrices.adim.ephase.elyte.Vmod] = ...
                    matrices_bdrys_wrt_internal(...
                                                                            matrices.adim.ephase.D1_noBC(1,:),...
                                                                            matrices.adim.ephase.D1_noBC(end,:),...
                                                                            @(t)(0),...
                                                                            @(t)(0),...
                                                                            [0,0]...
                                                                            );
% Remark : since there is no exogeneous boundary condition, matrices.adim.ephase.elyte.Umod and matrices.adim.ephase.elyte.Vmod are 0.

% Actually computing the _BC matrix, and the all the terms on the right side of dc_e/dt = ... and ddl/dt = ... :
% These terms are 
%           dc_e/dt = params.adim.theta_d * (elyte.D2_BC * ce + elyte.exogeneous_bdry ) ...
%                                  +  params.adim.theta_c * (dl.D2_BC * dl + dl.exogeneous_bdry) ...
%                                  +  params.adim.K * params.adim.theta_c *elyte.logterm_BC(ce)
%           ddl/dt =   params.adim.theta_c * (dl.D2_BC * dl + dl.exogeneous_bdry) ...
%                                  - lambda(dl,ussurf) ...
%                                  +  params.adim.K * params.adim.theta_c *elyte.logterm_BC(ce)
% The matrices depending on dl depend on time (because they depend on i(t)),
%       and then are defined in the file matrices_linearized_t_dep.m

matrices.adim.ephase.elyte.D2_BC = matrices.adim.ephase.D2_noBC(2:end-1,2:end-1)...
                                                                           + matrices.adim.ephase.D2_noBC(2:end-1,1) *  matrices.adim.ephase.elyte.Dmod...
                                                                           + matrices.adim.ephase.D2_noBC(2:end-1,end) *  matrices.adim.ephase.elyte.Emod;
matrices.adim.ephase.elyte.exogeneous_bdry_left = @(t) (matrices.adim.ephase.D2_noBC(2:end-1,1) * matrices.adim.ephase.elyte.Umod(t)); % This is 0 bcz we do not have an exogeneous term at the bdry for ce, but we still define it to have consistent notations
matrices.adim.ephase.elyte.exogeneous_bdry_right = @(t)(matrices.adim.ephase.D2_noBC(2:end-1,end) * matrices.adim.ephase.elyte.Vmod(t)); % Same
matrices.adim.ephase.elyte.exogeneous_bdry = @(t)(matrices.adim.ephase.elyte.exogeneous_bdry_left(t) + matrices.adim.ephase.elyte.exogeneous_bdry_right(t)); % Same
matrices.adim.ephase.elyte.logterm_BC = ...
            @(t,ELYTE2N) (matrices.adim.ephase.D2_noBC(2:end-1,2:end-1) * log(1 + ELYTE2N / params.adim.c_0)...
                                            + matrices.adim.ephase.D2_noBC(2:end-1,1) * log( 1 + (matrices.adim.ephase.elyte.Umod(t) + matrices.adim.ephase.elyte.Dmod * ELYTE2N)/params.adim.c_0)...
                                            + matrices.adim.ephase.D2_noBC(2:end-1,end) * log( 1 + (matrices.adim.ephase.elyte.Vmod(t) + matrices.adim.ephase.elyte.Emod * ELYTE2N)/params.adim.c_0)...
                                        );

%% Boundary conditions : solid phase
%The operators defined below give the dynamics for u_s, with automatically
%satisfied boundary conditions. The equation is (k is a parameter that designs the x-absciss, k is in [2, N_elyte]:
%   du_s/dt = D2_BC us + exogeneous_linear_dl * dl(k)...
%                       + exogeneous_nonlinear(dl(k),us)
matrices.adim.sphase.us.D2_BC = matrices.adim.sphase.D2_noBC(2:end-1,2:end-1);
matrices.adim.sphase.us.exogeneous_nonlinear = ...
            @(USSURF) (matrices.adim.sphase.D2_noBC(2:end-1,1) * USSURF');
%matrices.adim.sphase.us.compute_ussurf = ...
%           @(DLK,US2R) (DLK - numeric_inverse_x_min_sinh(params.misc.xmin,params.misc.xmax,coef,...
%                                                        DLK...
%                                                        + matrices.adim.sphase.D1_noBC(1,2:end-1) * US2R / (matrices.adim.sphase.D1_noBC(1) - 1)...
%                                                          )...
%                                                   );
                                               
%% Algebraic equation
matrices.adim.sphase.ussurf.linear_term_ussurf = - eye(matrices.N_elyte - 1);
matrices.adim.sphase.ussurf.linear_term_us = kron(eye(matrices.N_elyte - 1),matrices.adim.sphase.D1_noBC(1,2:end-1) / (1 - matrices.adim.sphase.D1_noBC(1)));

end
