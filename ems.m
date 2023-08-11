clear all; close all; clc;

%% READ CSV FILE
row=1;

numm=[  120000	-100000	120000	120000	40000
        80000	-75000	80000	80000	750000
        50000	-45000	120000	120000	80000  ];

for row=1:3
%% READ DATA
fprintf('Simulation Started\n\n') ;

num = csvread('Data.csv');

BATTpmax = num(row,1);
BATTpmin = num(row,2);
PVpmax = num(row,3);
WTpmax = num(row,4);
MTpmax = num(row,5);

%% TIME STEPS, AC BUS
time_step = 24;
totalUnits = 5;
p_bus = zeros(totalUnits,time_step);
index_GenMT = 1; index_GenPV = 2; index_GenWT = 3; index_GenBTR = 4; index_Load = 5; index_Util = 6;

% MICROURBINE 
mut = 6000;% minimum up time (sec)
mdt = 200;% minimum down time (sec)
ton=0;% counting the times turned on
toff=0;% counting the times turned off
statusMT=ones(1,time_step);

%% WIND TURBINE
statusWT=ones(1,time_step);

%% BATTERY
Vb=12;% initial voltage
pbold = 0;
[soc_cvx, qq1, qq2,num_of_batt] = battery( 0, pbold, BATTpmax);

%% UTILITY
vwind = [10.12 11.654 12.43 12.7654 12.123 12.5643 12.626 13.845 15.113 16.322 17.11 10.12 11.654 12.43 12.7654 10.12 11.654 12.43 12.7654 24.6363 25.346 26.7562 27.23 28.88];
irrad = [0.0075 0.3948 0.674 0.789 0.519 0.555  0.493  0.935  0.891 0.583 0.514 0.35 0.169 0.6472  0.1508 0.1324 0.10385 0.8206 0.1069 0.3603 0.7578 0.8679 0.0767 0.5666];
cutil = 10^(-3)*[0.0787 0.0785 0.0785 0.0789 0.088 0.088 0.089 0.09 0.0885 0.0805 0.0815 0.0825 0.079 0.078 0.089 0.0866 0.0799 0.082 0.08 0.085 0.09 0.0788 0.0788 0.085];
tempr = [33.57  66.38  54.61  65.75 42.07 70.09 30.81  54.57 57.26  28.02 65.66  72.92 62.93  49.59  54.64 43.63 32.40  39.59  56.84  47.46 68.34 59.68 55.18 44.52];
stp=3600;
%% EVALUATION
for i=1:time_step 
    %% CONSTRAINTS MT
   if( ton >= mut/stp )
        statusMT(i) = 0;
    end

    if( toff >= mdt/stp )
        statusMT(i) = 1;
    end

    if( statusMT(i) == 0 )
        ton=0;
        toff=toff+1;
        if(i<length(statusMT))
            statusMT(i+1)=0;
        end
    else
        toff = 0;
        ton = ton+1;
        if(i<length(statusMT))
            statusMT(i+1)=1;
        end
    end

    %% AC BUS
    p_bus(index_GenWT,i) = windturbine(vwind(i),statusWT(i),WTpmax);
    p_bus(index_Load,i) = power_bus_per_dayhour(time_step,i);
    if(mod(i,time_step)<=19 && mod(i,time_step)>=7)
        p_bus(index_GenPV,i) = pv_array(irrad(i),tempr(i),PVpmax);
    end
    
    if soc_cvx<=0.201
      soc_cvx=0.201;
    end
    %% CVX
%         n=3;
%         cvx_begin
%             variable x(n,1)
%             minimize objFun(x(1), x(2), statusMT, toff, cutil(i), i);
%             subject to 
%             x(1)+x(2)+x(3)*Vb==p_bus(index_Load,i)-(p_bus(index_GenWT,i)+p_bus(index_GenPV,i));
%             %-x(3)<=0;
%             %##########################
%             -x(3)*Vb-abs(BATTpmin)<=0;
%             -x(3)*Vb-(1-soc_cvx)*BATTpmax<=0;
%             %##########################
%             x(3)*Vb-(soc_cvx-0.2)*BATTpmax<=0;
%             x(3)*Vb-soc_cvx*BATTpmax<=0;
%             x(3)*Vb-abs(BATTpmin)<=0;
%             x(2)-MTpmax*statusMT(i)<=0;
%             -x(2)<=0;
%         cvx_end
%         x_cvx_x = cvx_optpnt.x;
%         fprintf('\nCVX finished. ');
            
    %% PRIMAL DUAL INTERIOR POINT METHOD
    clear r_t_d r_t_c r_t_p r_t 
    clear array_mtr 
    clear Dy Dx_pd Dlamda_pd Dv_pd 
    clear x lamda  v 
    clear eta tau f Df

    p=1; 
    mu=10;
    alpha = 1e-2; 
    beta = 5e-1;
    epsilon_feas = 1e-6;
    epsilon = 1e-6;

    b = p_bus(index_Load,i)-p_bus(index_GenWT,i)-p_bus(index_GenPV,i);

    v = rand(p,1);

    if statusMT(i)==1
        m=7;
        n=3;
        A = [1 1 Vb];
        x = rand(n,1);

        f = [ x(2)-MTpmax;
            -x(2);
            -x(3)*Vb-abs(BATTpmin);
            -x(3)*Vb-(1-soc_cvx)*BATTpmax;
            x(3)*Vb-(soc_cvx-0.2)*BATTpmax;
            x(3)*Vb-soc_cvx*BATTpmax;
            x(3)*Vb-abs(BATTpmin)];

        Df = [0 1 0;
            0 -1 0;
            0 0 -Vb;
            0 0 -Vb;
            0 0 Vb;
            0 0 Vb;
            0 0 Vb];
    else
        m=5;
        n=2;
        A = [1 Vb];
        x = rand(n,1);

        f = [ -x(2)*Vb-abs(BATTpmin);
            -x(2)*Vb-(1-soc_cvx)*BATTpmax;
            x(2)*Vb-(soc_cvx-0.2)*BATTpmax;
            x(2)*Vb-soc_cvx*BATTpmax;
            x(2)*Vb-abs(BATTpmin)];

        Df = [0 -Vb;
            0 -Vb;
            0 Vb;
            0 Vb;
            0 Vb];
    end

    Hessian = zeros(m,n);
    Hessian_f0 = D2_f_obj(statusMT(i));

    lamda = -1./f;

    k = 1;

    %% Repeat
    while (1)

        %% DETERMINE tau
        if statusMT(i)==1
            f = [ x(2,k)-MTpmax;
            -x(2,k);
            -x(3,k)*Vb-abs(BATTpmin);
            -x(3,k)*Vb-(1-soc_cvx)*BATTpmax;
            x(3,k)*Vb-(soc_cvx-0.2)*BATTpmax;
            x(3,k)*Vb-soc_cvx*BATTpmax;
            x(3,k)*Vb-abs(BATTpmin)];

        else
             f = [ -x(2,k)*Vb-abs(BATTpmin);
            -x(2,k)*Vb-(1-soc_cvx)*BATTpmax;
            x(2,k)*Vb-(soc_cvx-0.2)*BATTpmax;
            x(2,k)*Vb-soc_cvx*BATTpmax;
            x(2,k)*Vb-abs(BATTpmin)];
        end

        eta(k) = -(f)' * lamda(:,k);
        tau = (mu * n)/eta(k);
        t_max = 1;
        t = 0.99*t_max;
        %% UPDATE residuals
        if statusMT(i)
            r_t_d(:,k) = D_f_obj(x(2,k), statusMT(i), cutil(i)) + Df'*lamda(:,k) + A'*v(:,k);
            r_t_c(:,k) = -diag(lamda(:,k))*f - (1/tau)*ones(m,1);
            r_t_p(:,k) = A*x(:,k)-b;
            r_t(:,k) = [r_t_d(:,k); r_t_c(:,k); r_t_p(:,k)];
        else
            r_t_d(:,k) = D_f_obj(0, statusMT(i), cutil(i)) + Df'*lamda(:,k) + A'*v(:,k);
            r_t_c(:,k) = -diag(lamda(:,k))*f - (1/tau)*ones(m,1);
            r_t_p(:,k) = A*x(:,k)-b;
            r_t(:,k) = [r_t_d(:,k); r_t_c(:,k); r_t_p(:,k)];
        end

        %% UPDATE the invertible matrix for Δy
        array_mtr = [Hessian_f0 Df' A';
            -diag(lamda(:,k))*Df -diag(f) zeros(m,p);
            A zeros(p,m) zeros(p,p)];

        %% UPDATE Δy
        Dy = -inv(array_mtr)*r_t(:,k);
        Dx_pd = Dy(1:n);
        Dlamda_pd = Dy(n+1:n+m);
        Dv_pd = Dy(n+m+1:size(Dy,1));

        %% λ+Dλ>0
        while ~( sum( (lamda(:,k)+t*Dlamda_pd) > 0 ) == m )
            t=beta*t;
        end

        %% f(x+)<0
        if statusMT(i)==1
            f = [ (x(2,k)+t*Dx_pd(2))-MTpmax;
            -(x(2,k)+t*Dx_pd(2));
            -(x(3,k)+t*Dx_pd(3))*Vb-abs(BATTpmin);
            -(x(3,k)+t*Dx_pd(3))*Vb-(1-soc_cvx)*BATTpmax;
            (x(3,k)+t*Dx_pd(3))*Vb-(soc_cvx-0.2)*BATTpmax;
            (x(3,k)+t*Dx_pd(3))*Vb-soc_cvx*BATTpmax;
            (x(3,k)+t*Dx_pd(3))*Vb-abs(BATTpmin)];
        else
           f = [  -(x(2,k)+t*Dx_pd(2))*Vb-abs(BATTpmin);
                -(x(2,k)+t*Dx_pd(2))*Vb-(1-soc_cvx)*BATTpmax;
               (x(2,k)+t*Dx_pd(2))*Vb-(soc_cvx-0.2)*BATTpmax;
            (x(2,k)+t*Dx_pd(2))*Vb-soc_cvx*BATTpmax;
            (x(2,k)+t*Dx_pd(2))*Vb-abs(BATTpmin)];  
        end

        while ~( sum ( f < 0 ) == m )
            t=beta*t;

            if statusMT(i)==1
                f = [ (x(2,k)+t*Dx_pd(2))-MTpmax;
                -(x(2,k)+t*Dx_pd(2));
                -(x(3,k)+t*Dx_pd(3))*Vb-abs(BATTpmin);
                -(x(3,k)+t*Dx_pd(3))*Vb-(1-soc_cvx)*BATTpmax;
                (x(3,k)+t*Dx_pd(3))*Vb-(soc_cvx-0.2)*BATTpmax;
                (x(3,k)+t*Dx_pd(3))*Vb-soc_cvx*BATTpmax;
                (x(3,k)+t*Dx_pd(3))*Vb-abs(BATTpmin)];
            else
               f = [ -(x(2,k)+t*Dx_pd(2))*Vb-abs(BATTpmin);
                -(x(2,k)+t*Dx_pd(2))*Vb-(1-soc_cvx)*BATTpmax;
                    (x(2,k)+t*Dx_pd(2))*Vb-(soc_cvx-0.2)*BATTpmax;
                    (x(2,k)+t*Dx_pd(2))*Vb-soc_cvx*BATTpmax;
                    (x(2,k)+t*Dx_pd(2))*Vb-abs(BATTpmin)];  
            end

        end   

        %% BACKTRACKING
        while norm( func_res( x(:,k) + t*Dx_pd, lamda(:,k) + t*Dlamda_pd, v(:,k) + t*Dv_pd, A, b, m, tau, statusMT(i), Vb, MTpmax, cutil(i), BATTpmax,BATTpmin,soc_cvx ) ) >...
                (1-alpha*t) * norm(r_t(:,k))
            t = beta*t;
        end

        %% UPDATE x,λ,v
        x(:,k+1) = x(:,k) + t*Dx_pd ;
        lamda(:,k+1) = lamda(:,k) + t*Dlamda_pd;
        v(:,k+1) = v(:,k) + t*Dv_pd;

        %% TERMINATION CONDITION
        if norm(r_t_p(:,k))<=epsilon_feas && norm(r_t_d(:,k))<=epsilon_feas && eta(k)<=epsilon
            break; 
        end
        k = k+1;
    end

    if statusMT(i)==1
        [optimum_primal_dual_interior_point , Futility_prim, Fgen_prim, suc_prim] = objFun(x(1,k), x(2,k), statusMT, toff,cutil(i),i);
    else
        [optimum_primal_dual_interior_point , Futility_prim, Fgen_prim, suc_prim] = objFun(x(1,k), 0, statusMT, toff,cutil(i),i);
    end

    fprintf('\n\nPrimal Dual interior point method finished. ');

    fprintf('\n\nOptimum Primal Dual Interior Point: %d',optimum_primal_dual_interior_point);

%     fprintf('\n\nOptimal Points from CVX: [ %d %d %d ]',x_cvx_x(1),x_cvx_x(2),x_cvx_x(3));

    if statusMT(i)==1
        fprintf('\n\nOptimal Points from Primal Dual Interior Point: [ %d %d %d ]',x(1,k),x(2,k),x(3,k));
    else
        fprintf('\n\nOptimal Points from Primal Dual Interior Point: [ %d 0 %d ]',x(1,k),x(2,k));
    end
    %% WRITE TO BUS
%         p_bus(index_Util,i) = x_cvx_x(1);
%         p_bus(index_GenMT,i) = x_cvx_x(2);
%         p_bus(index_GenBTR,i) = x_cvx_x(3)*Vb;
%         val_obj(i) = cvx_optval;

        p_bus(index_Util,i) = x(1,k);
        if statusMT(i)==1
            p_bus(index_GenMT,i) = x(2,k);
            p_bus(index_GenBTR,i) = x(3,k)*Vb; 
        else
            p_bus(index_GenMT,i) = 0;
            p_bus(index_GenBTR,i) = x(2,k)*Vb;
        end
        val_obj(i) = optimum_primal_dual_interior_point;
    %% UPDATE BATTERY
    pbold = [pbold  p_bus(index_GenBTR,i)];
    [soc_cvx, pout_btr, qq3, num_of_batt] = battery( p_bus(index_GenBTR,i), pbold(1:end-1), BATTpmax);
end

loss_24_hour = sum(val_obj)
num(row,5+i+1) = loss_24_hour;
num(row,6:5+i) = val_obj;
csvwrite('Data.csv',num)
FileName =  sprintf('%d_%d_%d_%d_DATA.csv', BATTpmax, PVpmax, WTpmax, MTpmax);
csvwrite(FileName,p_bus)
plot_from_excel(FileName);
fprintf('Simulation Ended\n\n') ;
end
