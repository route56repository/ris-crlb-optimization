function [x_estimate] = orthogonalMatchingPursuit(A, y, k)
    % A: Measurement matrix
    % y: Measurement vector
    % k: Sparsity level

    [m, n] = size(A);
    x_estimate = zeros(n, 1); % Initialize the estimated sparse vector
    residual = y; % Initialize the residual

    indices = []; % Keep track of selected indices

    for iter = 1:k
        % Calculate correlations
        correlations = abs(A' * residual);

        % Find the index of the maximum correlation
        [~, index] = max(correlations);

        % Add the index to the selected indices
        indices = [indices, index];

        % Update the selected columns of A
        selectedA = A(:, indices);

        % Solve least squares problem to update x_estimate
        x_estimate(indices) = pinv(selectedA) * y;

        % Update the residual
        residual = y - selectedA * x_estimate(indices);
    end
end