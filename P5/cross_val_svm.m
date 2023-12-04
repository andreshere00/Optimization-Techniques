function [] = cross_val_svm(y,X,n,svm_method)
%% Cross-validation of a svm model.

% Create a cross-validation partition
cvp = cvpartition(y, 'KFold', n); % 5-fold cross-validation

% Initialize an array to store the trained SVM models
models = cell(cvp.NumTestSets, 1);

% Create a new figure and specify its size
fig = figure('Position', [100, 100, 1000, 600]);

for i = 1:cvp.NumTestSets
    % Split the data into training and testing sets
    trainIdx = cvp.training(i);
    testIdx = cvp.test(i);
    X_train = X(trainIdx, :);
    y_train = y(trainIdx);
    X_test = X(testIdx, :);
    y_test = y(testIdx);

    % Train the SVM model with an RBF kernel
    SVMModel = fitcsvm(X_train, y_train, 'KernelFunction', svm_method, 'BoxConstraint', 1, 'KernelScale', 'auto');

    % Store the trained model
    models{i} = SVMModel;
    
    if mod(n,2)==1
        m = floor(n/2) + 1;
    else
        m = n/2;
    end
       
    % Create subplots with larger size
    subplot(2, m, i, 'Parent', fig); % Adjust the subplot arrangement as needed

    % Create a grid for plotting
    x1 = linspace(min(X(:, 1)), max(X(:, 1)), 100);
    x2 = linspace(min(X(:, 2)), max(X(:, 2)), 100);
    [X1, X2] = meshgrid(x1, x2);
    XGrid = [X1(:), X2(:)];

    % Calculate scores for the grid points
    scores = predict(SVMModel, XGrid);
    scores = reshape(scores, size(X1));

    % Plot the decision boundary
    contour(x1, x2, scores, [0, 0], 'k', 'LineWidth', 2);
    hold on;

    % Plot the data points
    gscatter(X(:, 1), X(:, 2), y, 'rb');
    
    title(['Fold ', num2str(i)]);
    xlabel('Feature 1');
    ylabel('Feature 2');
end
end