function [Hess] = D2_f_obj(statusMT)
    alpha=0.0074*10^(-3); 
    if statusMT==1
        Hess =  [0 0 0; 0 2*alpha 0; 0 0 0];
    else
        Hess =  [0 0;0 0];
    end
end
