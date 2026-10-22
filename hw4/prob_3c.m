%% prob 3c: Hebbian network
clear;
clc;
close all;

N = 500;                     
alpha_values = 0.02:0.04:0.62;
                           
num_cues = 10;          
flip_fraction = 0.10;       
max_sweeps = 1000;           
rng(1);

num_alpha = numel(alpha_values);

unstable_fraction = zeros(num_alpha, 1);
theory_fraction = zeros(num_alpha, 1);
mean_overlap = zeros(num_alpha, 1);
mean_sweeps = zeros(num_alpha, 1);


for a = 1:num_alpha

    
    alpha = alpha_values(a);  % current loading ratio
    P = round(alpha * N);
    P_values(a) = P;

    fprintf('\nalpha = %.2f, P = %d patterns\n', alpha, P);
    X = 2 * (rand(N, P) < 0.5) - 1;
    W = (X * X') / N; % hebbian weight matrix
    W(1:N+1:end) = 0; % set diagonals to zero
    H = W * X; % calculate local fields
    is_unstable = (X .* H) < 0;
    unstable_fraction(a) = nnz(is_unstable) / (N * P);
    theory_fraction(a) = normcdf(-1 / sqrt(alpha));

    fprintf('  Measured unstable fraction: %.5f\n', unstable_fraction(a));
    fprintf('  CLT prediction:            %.5f\n', ...
        theory_fraction(a));

    cue_overlaps = zeros(num_cues, 1);
    cue_sweeps = zeros(num_cues, 1);
    

    for trial = 1:num_cues % test recall

        mu = randi(P);
        target = X(:, mu);
        cue = target;
        num_flips = round(flip_fraction * N);
        flip_indices = randperm(N, num_flips); % use to avoid accidentally flipping the same bit twice
        cue(flip_indices) = -cue(flip_indices); % reverse signs

        [recalled, sweeps] = asynchronous_recall(W, cue, max_sweeps);
        cue_overlaps(trial) = (target' * recalled) / N; % m = (1/N)*target'*recalled

        cue_sweeps(trial) = sweeps;

    end
    mean_overlap(a) = mean(cue_overlaps);
    mean_sweeps(a) = mean(cue_sweeps);

    fprintf('  Mean final overlap:         %.4f\n', ...
        mean_overlap(a));
    fprintf('  Mean number of sweeps:      %.2f\n', ...
        mean_sweeps(a));

end

figure('Color', 'w');
plot(alpha_values, unstable_fraction, 'ko-', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'simulation');
hold on;
plot(alpha_values, theory_fraction, 'r--', ...
    'LineWidth', 1.8, ...
    'DisplayName', '\Phi(-1/\surd\alpha)');
xlabel('\alpha');
ylabel('fraction of unstable bits');
title('stability of stored patterns');
legend('Location', 'northwest');
grid on;
box on;


figure('Color', 'w');
plot(alpha_values, mean_overlap, 'bo-', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
yline(1 - 2*flip_fraction, 'k--', ...
    'initial cue overlap', 'LineWidth', 1.2);
xlabel('\alpha');
ylabel('\langle m \rangle');
title('testing recall');
grid on;
box on;

function [s, sweeps] = asynchronous_recall(W, s, max_sweeps)
% inputs: 
% W: symmetric weight matrix with zero diagonal
% s: initial binary state (column vector)
% max_sweeps: maximum number of complete sweeps

    N = length(s);

    for sweeps = 1:max_sweeps
        update_order = randperm(N);
        changed = false;

        for k = 1:N
            i = update_order(k);
            h_i = W(i,:) * s;
            % sign(0) = +1
            if h_i >= 0
                new_state = 1;
            else
                new_state = -1;
            end
            if new_state ~= s(i)
                s(i) = new_state;
                changed = true;
            end

        end
        if ~changed
            return;
        end

    end
end
