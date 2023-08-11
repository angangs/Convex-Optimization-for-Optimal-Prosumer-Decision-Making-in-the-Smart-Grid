function [ObjFun,Futility,Fgen,suc] = objFun(putil, pmt , status_mt, toff,cutil,index)

    if(index==1)
        statusMT = status_mt(1);
        statusMT_prev = 0;
    else
        statusMT = status_mt(index);
        statusMT_prev = status_mt(index-1);
    end
    
    %% UTILITY
    Futility = cutil' .* putil;
       
    %% GENERATOR MAINTAINANCE AND START-UP
    stp = 3600; %
    hot_startup = 30; %sec
    cold_startup = 200; %sec
    
    cooling_time = 520; %sec
    
    if statusMT_prev==0 && statusMT==1
        suc = ( ( hot_startup / stp ) + ( cold_startup / stp ) * (1 - exp ( - toff / (cooling_time/stp) )) ); % start-up cost
    else
        suc=0;
    end
   
    alpha=0.0074; beta=0.2333; cu=0.4333;    
    Fgen = ( alpha * (pmt*10^(-3)).^2 + beta * pmt*10^(-3) + cu + suc) .* statusMT;
    
    %% OBJECTIVE FUNCTION
     ObjFun =  sum( Futility +  Fgen); % $/Wh
end
