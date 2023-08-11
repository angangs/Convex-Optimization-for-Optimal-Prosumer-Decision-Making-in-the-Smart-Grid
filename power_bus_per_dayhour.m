function [p_bus] = power_bus_per_dayhour(pertime,stp_time)

    if(stp_time<=pertime)
        if (stp_time==1 || stp_time==8 || stp_time==9 || stp_time==16 || stp_time==17 || stp_time==24)
            p_bus = 30000;
        elseif (stp_time==11 || stp_time==12 || stp_time==15 || stp_time==18 || stp_time==19 || stp_time==23)
            p_bus = 60000;
        elseif (stp_time==10)
            p_bus = 90000;
        elseif (stp_time==21 || stp_time==22)
            p_bus = 150000;
        elseif (stp_time==14 || stp_time==20)
            p_bus = 210000;
        elseif (stp_time==13)
            p_bus = 300000;
        else
            p_bus = 3000;
        end
    else
        if (mod(stp_time,pertime)==1 || mod(stp_time,pertime)==8 || mod(stp_time,pertime)==9 || mod(stp_time,pertime)==16 || mod(stp_time,pertime)==17 || mod(stp_time,pertime)==24)
            p_bus = 30000;
        elseif (mod(stp_time,pertime)==11 || mod(stp_time,pertime)==12 || mod(stp_time,pertime)==15 || mod(stp_time,pertime)==18 || mod(stp_time,pertime)==19 || mod(stp_time,pertime)==23)
            p_bus = 60000;
        elseif (mod(stp_time,pertime)==10)
            p_bus = 90000;
        elseif (mod(stp_time,pertime)==21 || mod(stp_time,pertime)==22)
            p_bus = 150000;
        elseif (mod(stp_time,pertime)==14 || mod(stp_time,pertime)==20)
            p_bus = 210000;
        elseif (mod(stp_time,pertime)==13)
            p_bus = 300000;
        else
            p_bus = 3000;
        end
    end
    
end