function main(parameters)

EZ=parameters.EZ

tau0 = 6667; 
patient = 'P1';
connectivity_path = strcat('../data/connectivity_', patient, '/');
K = load(strcat(connectivity_path, 'weights.txt')); % Load connectivity matrix
K = normal(K); % Normalize
[N, ~] = size(K); % Get the number of nodes

%% Generate K0 and extended range of K matrices
results = cell(1, N);  % Adjusted to store more matrices
results{1} = K;  % K0, no modification

for i = 1:N
    K_modified = K;  % Copy the original matrix
    K_modified(i, :) = [];  % Delete the i-th row
    K_modified(:, i) = [];  % Delete the i-th column
    results{i+1} = K_modified;  % Store modified matrix
end

%% Initialize
format long
x0_range1 = -2.20; % Fixed value, not used in loops
x0_EZ_range1 = -2.2:0.1:-2.1
disp(x0_EZ_range1);

% Initialize a 3D cell array to collect all StabilityMatrix1 results for each K
AllStabilityMatrix1 = cell(1, N+1); % Expanded to match the new K ranges

%% Loop over all K matrices (K0 to KN)
for k_idx = 1:length(results)  % Iterate over all K matrices
    K_current = results{k_idx}; % Current matrix K0 or modified versions
    [N_current, ~] = size(K_current); % Update the size for the current matrix
    StabilityMatrix1 = NaN(1, length(x0_EZ_range1)); % Initialize for each K matrix
    
    %% Single loop over x0(EZ)
    for ix0_EZ = 1:length(x0_EZ_range1)
        x0 = x0_range1 + zeros(N_current, 1); % Set x0 as a vector with the fixed value
        EZ_index = EZ - 61; % Adjust EZ index to match current matrix size
        if EZ_index <= N_current % Ensure EZ index is within the range of the current matrix
            x0(EZ_index) = x0_EZ_range1(ix0_EZ); % Adjust the EZ node
        end

        % Solve fixed-point
        Z0 = 3 + zeros(N_current, 1); 
        one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K_current, tau0); 
        opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
        [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);

        % Calculate Coupling Matrix for numerical solution
        C1 = CouplingMatrix(z_fixed_numerical, K_current, tau0);

        % Calculate Stability
        [covariance1, err] = con2cov(C1, false, 10000000, 10, 1);
        stability1 = Stability(covariance1); 

        % Save Stability into matrix
        StabilityMatrix1(ix0_EZ) = stability1;
        fprintf('EZ = %d, K%d Stability1: %.3f\n', EZ, k_idx-1, stability1);
    end
    
    % Save StabilityMatrix1 for current K
    AllStabilityMatrix1{k_idx} = StabilityMatrix1; % Store results in temporary variable
end

% Nested Functions

function w = normal(W)
    [N, ~] = size(W);
    u = zeros(N * (N - 1) / 2, 1);
    k = 0;
    for i = 1:N
        for j = i + 1:N
            k = k + 1;
            u(k, 1) = W(i, j);
        end
    end
    u = sort(u);
    uu = u(floor(length(u) * 0.95));
    w = W - diag(diag(W));
    idx = w > uu;
    w(idx) = uu;
    w = w / (max(max(w)));
end
