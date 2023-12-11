function [c , ceq]= top_ejem_pnl_res3(l)
    x=[1.3028,1.6972];
    d=[0.4343,1];
    sust=x+l*d;
    ceq=[];
    c=[-sust(1).^2 + sust(2);(sust(1) - 1).^2 + sust(2) - 5;-sust(2)];
end