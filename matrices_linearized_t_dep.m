function [ matrices ] = matrices_linearized_t_dep( matrices,params )
%MATRICES_LINEARIZED_T_DEP Summary of this function goes here
%   Detailed explanation goes here
[matrices.adim.ephase.dl.Dmod,matrices.adim.ephase.dl.Emod,matrices.adim.ephase.dl.Umod,matrices.adim.ephase.dl.Vmod] = ...
                   (matrices_bdrys_wrt_internal(...
                                                                            matrices.adim.ephase.D1_noBC(1,:),...
                                                                            matrices.adim.ephase.D1_noBC(end,:),...
                                                                            params.i,...
                                                                            @(t)(params.i(t)),...
                                                                            [0,0]...
                                                                            ));
matrices.adim.ephase.dl.D2_BC = ...
                                    @(t) (matrices.adim.ephase.D2_noBC(2:end-1,2:end-1)...
                                                       + matrices.adim.ephase.D2_noBC(2:end-1,1) *  matrices.adim.ephase.dl.Dmod...
                                                       + matrices.adim.ephase.D2_noBC(2:end-1,end) *  matrices.adim.ephase.dl.Emod);
matrices.adim.ephase.dl.exogeneous_bdry = ...
                                    @(t) (matrices.adim.ephase.D2_noBC(2:end-1,1) * matrices.adim.ephase.dl.Umod(t)...
                                                                           + matrices.adim.ephase.D2_noBC(2:end-1,end) * matrices.adim.ephase.dl.Vmod(t));

matrices.adim.ephase.dl.butler_volmer = ...
                @(DL2N, USSURF2N) (...
                                  (params.alpha * params.E * (DL2N - USSURF2N) )...
                                  );
end

