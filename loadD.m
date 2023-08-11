function [pout, emission] = loadD(P)
    emissionSO2 = 0.42/10^6; % kg/Wh
    emissionCO2 = 1230/10^6;
    emissionNOx = 2.35/10^6;
    
    pout = P;
    
    emission = [emissionSO2*pout emissionCO2*pout emissionNOx*pout];
end