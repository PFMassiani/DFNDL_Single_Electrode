function [ cs ] = init_cs( dl,params )
%INIT_CS Summary of this function goes here
%   Detailed explanation goes here

sto_guess = 0.7*ones(length(dl),1);

TOLERANCE = 1e-8;
for k = 1:length(dl)
    sto_high = 1;
    sto_low = 0.2;
    if(dl(k) < ocp_dualfoil(1))
        warning('DFNDL:InitialConditions',['No compatible initial conditions found for node #',num2str(k),' : double layer potential is too low. Try putting a lower input current, raising the maximal concentration in solid phase, or changing the OCP function.']);
        sto_guess(k) = 1;
    elseif (dl(k) > ocp_dualfoil(0.2))
        warning('DFNDL:InitialConditions',['No compatible initial conditions found for node #',num2str(k),' : double layer potential is too high. Try putting a lower input current or changing the OCP function.']);
        sto_guess(k) = 0.2;
    else
        while(abs(dl(k) - ocp_dualfoil(sto_guess(k))) > TOLERANCE)
            if dl(k) - ocp_dualfoil(sto_guess(k)) > 0
                sto_high = sto_guess(k);
            else
                sto_low = sto_guess(k);
            end
            sto_guess(k) = (sto_high + sto_low) / 2;
        end
    end
end
cs = params.csmax * ones(params.dscrtzn.N_s+1,1) * sto_guess';

end

