function [performance_rate, points_bet_hyp] = svm_soft(x1, x2, y, H, f, Aeq, beq, C, plt)
%% Description:
% Own-made MATLAB function implements a Support Vector Machine (SVM) with a 
% soft margin for binary classification.
%
% Syntax:
% [performance_rate, points_bet_hyp] = svm_soft(x1, x2, y, H, f, Aeq, beq, C, plt)
%
% Inputs:
% - x1: Feature values for the first dimension.
% - x2: Feature values for the second dimension.
% - y: Class labels (+1 or -1).
% - H: Hessian matrix in the quadratic programming problem.
% - f: Linear coefficient vector in the quadratic programming problem.
% - Aeq: Coefficients for linear equality constraints.
% - beq: Values for linear equality constraints.
% - C: Regularization parameter for the soft margin.
% - plt: Boolean indicating whether to plot the decision boundaries and support vectors (true or false).
%
% Outputs:
% 
% - performance_rate: The rate of correctly classified data points.
% - points_bet_hyp: The proportion of data points located between the decision boundaries.


% Lower and upper bound for the alpha value.
m= size(x1,1);
lb = zeros(m,1);
ub = C * ones(m,1);

[alpha, fun, exitflag] = quadprog(H, f, [], [], Aeq, beq, lb, ub);

if exitflag == -2
    return
else
    % normal vector to the hyperplanes
    w = [x1,x2]'*(alpha.*y);

    % identify support vectors for the class +1:
    for i=1:m
        if y(i) > 0
            indk = i;
        else
            indh = i;
        end
    end

    b = ([x1(indk),x2(indk)]*w + [x1(indh),x2(indh)]*w)/2;

    % Plotting
    % indneg = y<0; indpos = y>0;
    % H^+
    xx1 = min(x1):0.1:max(x1);
    xx2p = (b + 1 - w(1) * xx1)./ w(2);
    % H^-
    xx2n = (b - 1 - w(1) * xx1)./ w(2);
    % H
    xx2 = (b - w(1) * xx1)./ w(2);

    if plt == true
        gscatter(x1, x2, y)
        hold on
        plot(xx1,xx2p,'-m',xx1,xx2n,'-k', xx1,xx2,'-g','LineWidth',2);
        title("Binary classification with C="+num2str(C));
        xlabel('X_1');
        ylabel('X_2');
        ylim([min(x2) - abs(min(x2)), max(x2) + abs(max(x2))]);
        xlim([min(x1) - abs(min(x1)), max(x1) + abs(max(x1))])
        legend('class +1','class -1','H^+','H^-','H')
        grid on
        hold off
    end

    % parameters of each hyperplane
    % H1
    Xx1 = [ones(length(xx1),1), xx1'];
    bn = Xx1\xx2n';
    bp = Xx1\xx2p';
    bh = Xx1\xx2';
    hn = @(x1) bn(1) + bn(2).*x1;
    hp = @(x1) bp(1) + bp(2).*x1;
    h = @(x1) bh(1) + bh(2).*x1;

    % points in a wrong class
    missclassified = 0;
    % points between hyperplanes
    points_hyperplane = 0;
    % validation criterion

    % checking if wether the points fits in the hyperplane or are
    % misclassified
    for k = 1:length(x1)
        if bn(2) < 0
            if (x2(k) > hn(x1(k))) || (x2(k) < hp(x1(k)))
                points_hyperplane = points_hyperplane + 1;
            end
            if (x2(k) > h(x1(k))) && (y(k) ~= -1) || (x2(k) < h(x1(k))) && (y(k) ~= 1)
                missclassified = missclassified + 1;
            end
        else
            if (x2(k) < hn(x1(k))) || (x2(k) > hp(x1(k)))
                points_hyperplane = points_hyperplane + 1;
            end
            if (x2(k) < h(x1(k))) && (y(k) ~= -1) || (x2(k) > h(x1(k))) && (y(k) ~= 1)
                missclassified = missclassified + 1;
            end
        end
    end
    performance_rate = (length(x1)-missclassified)/length(x1);
    points_bet_hyp = points_hyperplane/length(x1);
end
end