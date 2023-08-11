function [pout, emission] = microturbine(p,status,prated)
    emissionSO2 = 0.206*10^(-6); % kg/Wh
    emissionCO2 = 649*10^(-6);
    emissionNOx = 9.89*10^(-6);
    
    if(status==1)
        pout(p>prated) = prated;
        pout(p<=prated && p>=0)=p;
    else
        pout = 0;
    end

    emission = [emissionSO2*pout emissionCO2*pout emissionNOx*pout];
end