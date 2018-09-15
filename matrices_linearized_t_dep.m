function [ t_dep_matrices ] = matrices_linearized_t_dep( t,t_indep_matrices,params )
%MATRICES_LINEARIZED_T_DEP Summary of this function goes here
%   Detailed explanation goes here
[t_dep_matrices.Dmod,t_dep_matrices.Emod,t_dep_matrices.Umod,t_dep_matrices.Vmod] = ...
                   matrices_bdrys_wrt_internal(...
                                                                            t_indep_matrices.adim.ephase.D1_noBC(1,:),...
                                                                            t_indep_matrices.adim.ephase.D1_noBC(end,:),...
                                                                            params.adim.i(t),...
                                                                            -params.adim.i(t),...
                                                                            [0,0]...
                                                                            );
t_dep_matrices.D2_BC = t_indep_matrices.adim.ephase.D2_noBC(2:end-1,2:end-1)...
                                                                           + t_indep_matrices.adim.ephase.D2_noBC(2:end-1,1) *  t_dep_matrices.Dmod...
                                                                           + t_indep_matrices.adim.ephase.D2_noBC(2:end-1,end) *  t_dep_matrices.Emod;
t_dep_matrices.exogeneous_bdry = t_indep_matrices.adim.ephase.D2_noBC(2:end-1,end) * t_dep_matrices.Umod...
                                                                           + t_indep_matrices.adim.ephase.D2_noBC(2:end-1,1) * t_dep_matrices.Vmod;

t_dep_matrices.butler_volmer = ...
                @(DL2N, USSURF2N) (...
                                   sinh(params.adim.alpha * params.adim.E * (DL2N - USSURF2N) )...
                                   );
end

