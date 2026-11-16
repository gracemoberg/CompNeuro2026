% prob 1
clear; close all; clc;
rng(1);
if ~exist('images','dir'), mkdir('images'); end

% parameters
rho = 0.5; pH = 0.5; tau = 5;

% long-run intake rate R(T) and belief-weighted rate lambda~(T)
Rfun = @(T) pH*(1-exp(-rho*T)) ./ (pH*(1-exp(-rho*T))/rho + (1-pH)*T + tau);
lam = @(T) rho*pH*exp(-rho*T) ./ (pH*exp(-rho*T) + 1 - pH);

% rate if the patch type were known on arrival
Rknown = pH / (pH/rho + tau);

% analytic values and optimal wait
fprintf('R(2) = %.5f, R(5) = %.5f, Rknown = %.5f\n', Rfun(2), Rfun(5), Rknown);

Tstar = fzero(@(T) Rfun(T) - lam(T), [0.1 20]);   % R(T) = lambda~(T)
Rstar = Rfun(Tstar);
fprintf('T* = %.4f, R* = %.5f, lambda~(T*) = %.5f\n', ...
        Tstar, Rstar, lam(Tstar));
fprintf('R*/Rknown = %.4f\n', Rstar/Rknown);

% simulate N patch visits at T = 2, T*, 10
N = 4e5;
Tlist = [2, Tstar, 10];
Rsim = zeros(1,3);
for k = 1:3
    T = Tlist(k);
    high = rand(N,1) < pH; % patch type
    X = -log(rand(N,1)) / rho; % time to find chunk (high patches only)
    found = high & (X < T); % chunk found before giving up
    stay = T*ones(N,1); % time spent in patch
    stay(found) = X(found); % leave at once after eating
    Rsim(k) = sum(found) / (sum(stay) + N*tau);
    fprintf('T = %7.4f: simulated R = %.5f, R(T) = %.5f\n', T, Rsim(k), Rfun(T));
end

% R(T), Rknown, and lambda~(T), marking T* and R*
Tgrid = linspace(0, 25, 1000);
figure('Position',[100 100 700 450]); hold on; box on;
plot(Tgrid, Rfun(Tgrid), 'b-', 'LineWidth', 1.8);
plot(Tgrid, lam(Tgrid), 'r-', 'LineWidth', 1.8);
yline(Rknown, 'k--', 'LineWidth', 1.5);
plot(Tstar, Rstar, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
xline(Tstar, ':', 'Color', 'black');
yline(Rstar, ':', 'Color', 'black');
plot(Tlist, Rsim, 'pentagram', 'MarkerSize', 9, 'MarkerFaceColor', 'black');
xlabel('wait time T');
ylabel('rate');
legend({'R(T)', '\tilde\lambda(T)', 'R_{\rm known}', ...
        '(T^*,R^*)', '', '', 'simulated'}, 'Location', 'northeast');
title(sprintf('T^*=%.3f, R^*=%.4f', Tstar, Rstar));
ylim([0 0.26]);
