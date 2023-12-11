function f= top_ejem_pnl_obj3(l)
%% Objective function for the second iteration of the Zoutendijk algorithm
    x=[1.3028, 1.6972];
    d=[0.4343, 1];
    sust=x + l*d;
    f = 2*sust(1) -sust(2); 
end

