function r_t = func_res(x,lamda,v,A,b,m,tau,statusMT,Vb,MTpmax,cutil,BATTpmax,BATTpmin,soc)
    
    if statusMT==1
        f = [ x(2)-MTpmax;
            -x(2);
            -x(3)*Vb-abs(BATTpmin);
            -x(3)*Vb-(1-soc)*BATTpmax;
            x(3)*Vb-(soc-0.2)*BATTpmax;
            x(3)*Vb-soc*BATTpmax;
            x(3)*Vb-abs(BATTpmin)];

        Df = [0 1 0;
            0 -1 0;
            0 0 -Vb;
            0 0 -Vb;
            0 0 Vb;
            0 0 Vb;
            0 0 Vb];

        r_t_d = D_f_obj(x(2), statusMT, cutil) + Df'*lamda + A'*v;      

    else
       f = [ -x(2)*Vb-abs(BATTpmin);
            -x(2)*Vb-(1-soc)*BATTpmax;
            x(2)*Vb-(soc-0.2)*BATTpmax;
            x(2)*Vb-soc*BATTpmax;
            x(2)*Vb-abs(BATTpmin)];

        Df = [0 -Vb;
            0 -Vb;
            0 Vb;
            0 Vb;
            0 Vb];
          
        r_t_d = D_f_obj(0, statusMT, cutil) + Df'*lamda + A'*v;      
        
    end  
    
    r_t_c = -diag(lamda)*f - (1/tau)*ones(m,1);
    r_t_p = A*x-b;
    r_t = [r_t_d; r_t_c; r_t_p];
end