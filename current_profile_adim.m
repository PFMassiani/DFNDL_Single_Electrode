function [ i_adim ] = current_profile_adim( t_adim,params )
%CURRENT_PROFILE_ADIM Summary of this function goes here
%   Detailed explanation goes here

i_adim = current_profile(t_adim * params.R_s^2 / params.D_s) * params.L / (params.sigma * params.V_0);
end

