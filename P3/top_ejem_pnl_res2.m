function [c , ceq]= top_ejem_pnl_res2(l)
%% Constraints associated to the objective function in the first iteration of the Zoutendijk algorithm
    ceq=[ ];
    c=[-(2-l)^2+1+l;(1-l)^2+l-4;-l-1];
end