% 5c: drift + white noise with euler-maruyama with reset
% du/dt = mu + sqrt(2*sigma^2) * xi(t)
% threshold = theta, reset u -> 0

clear; close all; clc;

%% parameters
mu = 1;
theta = 1;

sigma2_values = [0.25, 1, 4];

dt = 1e-3;
Tburn = 100;   
T = 1000;          

nu_theory = mu/theta;

fprintf('Theoretical firing rate: %.4f\n\n', nu_theory);

% noise levels
for k = 1:length(sigma2_values)

    sigma2 = sigma2_values(k);
    sigma = sqrt(sigma2);

    % simulation time
    Ttotal = Tburn + T;
    nt = round(Ttotal/dt) + 1;

    t = (0:nt-1)' * dt;
    u = zeros(nt,1);

    spike_times = [];

    %% Euler-Maruyama simulation
    for n = 1:nt-1

        % Drift + noise
        u(n+1) = u(n) + mu*dt + sqrt(2*sigma2*dt)*randn;

        % threshold and reset
        if u(n+1) >= theta
            spike_times(end+1,1) = t(n+1); %#ok<SAGROW> from zack's code to ignore preallocation warning
            u(n+1) = 0;
        end
    end

    % discard transient
    keep = t >= Tburn;
    t_stat = t(keep);
    u_stat = u(keep);
    % spike times occurring after transient
    spikes_stat = spike_times(spike_times >= Tburn);
    % measured firing rate
    nu_measured = length(spikes_stat) / T;
    % measured fraction below reset
    frac_below = mean(u_stat < 0);
    % theoretical fraction below reset
    frac_theory = (sigma2/(mu*theta)) * (1 - exp(-mu*theta/sigma2));

    fprintf('sigma^2 = %.2f\n', sigma2);
    fprintf('  firing rate:       measured = %.4f, theory = %.4f\n', nu_measured, nu_theory);
    fprintf('  fraction u < 0:    measured = %.4f, theory = %.4f\n\n', frac_below, frac_theory);
    % analytical stationary density
    u_theory = linspace(-5, theta, 1000);
    p_theory = zeros(size(u_theory));
    positive = (u_theory >= 0);
    negative = (u_theory < 0);
    p_theory(positive) = (nu_theory/mu).*(1 - exp(mu*(u_theory(positive)-theta)/sigma2));

    p_theory(negative) = (nu_theory/mu).*(1-exp(-mu*theta/sigma2)).* exp(mu*u_theory(negative)/sigma2);

    figure;
    histogram(u_stat, 100, 'Normalization', 'pdf', 'DisplayStyle', 'bar');
    hold on;

    plot(u_theory, p_theory, 'LineWidth', 2);

    xline(0, '--k');
    xline(theta, '--r');

    xlabel('u');
    ylabel('density');
    title(sprintf('simulated model and theoretical p(u) with \\sigma^2 = %.2f', sigma2));
    legend('simulation', 'theory', 'reset', 'threshold', 'Location', 'northwest');
    box on;

end
