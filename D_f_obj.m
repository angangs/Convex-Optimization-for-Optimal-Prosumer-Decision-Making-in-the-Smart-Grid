function [gradient] = D_f_obj(pmt, statusMT,cutil)       
    alpha=0.0074*10^(-3); beta=0.2333*10^(-3);
    
    if statusMT==1
        gradient =  [cutil ( 2*alpha*pmt + beta ) 0]';
    else
        gradient =  [cutil 0]';
    end
end
