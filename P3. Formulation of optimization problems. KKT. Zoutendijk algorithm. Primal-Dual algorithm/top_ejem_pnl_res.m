 function [c , ceq]= top_ejem_pnl_res(x)
 %% Constraints associated to the objective function of the exercise 1, practice 3
  ceq=[ ];
  c=[-x(1).^2+x(2);(x(1)-1).^2+x(2)-5;-x(2)];
 end