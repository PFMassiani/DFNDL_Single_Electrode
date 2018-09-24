function [ u,x,r,dl,elyte,us,cs ] = postproc( t,state,params,matrices )
%POSTPROC Summary of this function goes here
%   Detailed explanation goes here
NTOT = (params.dscrtzn.N_s + 1)* (params.dscrtzn.N_e - 1) + params.dscrtzn.N_e - 1;
Nem2 = params.dscrtzn.N_e - 1;
Nsm2 = params.dscrtzn.N_s - 1;
Ne = Nem2 + 2;
Ns = Nsm2 + 2;

mask_us = logical([zeros(2*Nem2,1);ones(NTOT - 3*Nem2,1);zeros(Nem2,1)]);
mask_dl = logical([ones(1,Nem2), zeros(1, NTOT - Nem2)]);
mask_elyte = logical([zeros(1,Nem2), ones(1,Nem2), zeros(1,NTOT - 2*Nem2)]);
mask_ussurf = logical([zeros(NTOT - Nem2,1); ones(Nem2,1)]);

u = t * params.R_s^2/params.D_s;
x =matrices.adim.ephase.x * params.L;
r = matrices.adim.sphase.r * params.R_s;

dl = zeros(length(t),Ne);
elyte = zeros(length(t),Ne);
us = zeros(length(t),Ns,Nem2);
for i = 1:length(t)
    matrices = matrices_linearized_t_dep(t(i),matrices,params);
    dl(i,2:end-1) = state(i,mask_dl);
    elyte(i,2:end-1) = state(i,mask_elyte);
    us(i,2:end-1,:) = reshape(state(i,mask_us),Nsm2,Nem2);
    
    dl(i,1) = matrices.adim.ephase.dl.Umod + matrices.adim.ephase.dl.Dmod * dl(i,2:end-1)';
    dl(i,end) = matrices.adim.ephase.dl.Vmod + matrices.adim.ephase.dl.Emod * dl(i,2:end-1)';
    
    elyte(i,1) = matrices.adim.ephase.elyte.Umod + matrices.adim.ephase.elyte.Dmod * elyte(i,2:end-1)';
    elyte(i,end) = matrices.adim.ephase.elyte.Vmod + matrices.adim.ephase.elyte.Emod * elyte(i,2:end-1)'; 
    
    us(i,end,:) = 0;
    us(i,1,:) = state(i,mask_ussurf);
end
dl = dl * params.V_0;
elyte =  elyte * params.aC * params.V_0 * (1 - params.t_0) / params.Faraday;
us = us * (params.csmax*params.R_s);

cs = zeros(size(us));
for v = 1:length(u)
    for j = 1:size(us,3)
        cs(v,1:end-1,j) = us(v,1:end-1,j) ./ r(1:end-1)';
        cs(v,end,j) = (matrices.sphase.D1_noBC(end,:)*us(v,:,j)')';
    end
end
end

