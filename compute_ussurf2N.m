function [ ussurf ] = compute_ussurf2N( x,matrices,mask )
%COPMUTE_USSURF Summary of this function goes here
%   Detailed explanation goes here
Nem2 = matrices.N_elyte - 1;

warning('The function compute_ussurf2N is used, and it shouldn''t');

ussurf = zeros(Nem2,1);
for k=1:Nem2
    ussurf(k) = matrices.adim.sphase.us.compute_ussurf(x(mask.dl(:,k)),x(mask.us2R(:,k)));
end

end

