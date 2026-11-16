% prob 2
clear; close all; clc;
if ~exist('images','dir'), mkdir('images'); end

% P(i,j): probability the shot goes in
% rows: shoot 2, shoot 3; cols: D2, D3
P = [0.40 0.60; 0.40 0.25];
A = diag([2 3]) * P; % expected points
xs = 8/17;  ys = 9/17; % equilibrium: x = P(shoot 3), y = P(guard paint)
V = 84/85;
fprintf('x* = %.4f, y* = %.4f, value = %.4f\n', xs, ys, V);

sig = @(t) 1 ./ (1 + exp(-t));
logit = @(p) log(p ./ (1-p));

% mean-flow orbit from (0.2, 0.3)  
meanflow = @(t, th) [ sig(th(1))*(1-sig(th(1))) * 0.85 * (sig(th(2)) - ys); -sig(th(2))*(1-sig(th(2))) * 0.85 * (sig(th(1)) - xs)];
th0 = [logit(0.2); logit(0.3)];
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);
tt = linspace(0, 3000, 6001);
[~, TH] = ode45(meanflow, tt, th0, opts);
Xo = sig(TH(:,1));  Yo = sig(TH(:,2));

Kint = @(c, z) integral(@(s) (s - c) ./ (s .* (1 - s)).^2, c, z, 'RelTol', 1e-12, 'AbsTol', 1e-14);
idx = 1:10:numel(tt);
K = zeros(size(idx));
for m = 1:numel(idx)
    K(m) = Kint(xs, Xo(idx(m))) + Kint(ys, Yo(idx(m)));
end
fprintf('mean flow: x in [%.4f, %.4f], y in [%.4f, %.4f]\n', ...
        min(Xo), max(Xo), min(Yo), max(Yo));
fprintf('K: mean = %.6f, max-min = %.3e, relative variation = %.3e\n', ...
        mean(K), max(K)-min(K), (max(K)-min(K))/mean(K));

figure('Position',[100 100 520 480]); hold on; box on;
plot(Xo, Yo, 'b-', 'LineWidth', 1.5);
plot(xs, ys, 'ro', 'MarkerFaceColor', 'r');
plot(0.2, 0.3, 'ks', 'MarkerFaceColor', 'k');
xlabel('x'); 
ylabel('y');
xlim([0 1]); ylim([0 1]);
title('mean-flow orbit from (0.2,0.3)');


% stochastic simulation, tau_theta = 0 and 0.05
eta = 0.02; 
N = 4e5;
tauvals = [0, 0.05];
th_init = [0.3; -0.3];
XH = cell(1,2); 
YH = cell(1,2);
for k = 1:2
    rng(1);
    [XH{k}, YH{k}] = simrun(tauvals(k), eta, N, th_init, P);
end

q = N/4;
for k = 1:2
    x = XH{k}; y = YH{k};
    d = sqrt((x - xs).^2 + (y - ys).^2);
    fprintf('tau = %.2f: last-quarter avg x = %.4f, y = %.4f\n', tauvals(k), mean(x(end-q+1:end)), mean(y(end-q+1:end)));
    fprintf('mean distance early = %.4f, late = %.4f\n', mean(d(1:q)), mean(d(end-q+1:end)));
end

figure('Position',[100 100 1000 450]);
for k = 1:2
    subplot(1,2,k); hold on; box on;
    plot(XH{k}(1:50:end), YH{k}(1:50:end), '-', 'Color', [0.2 0.4 0.9 0.5]);
    plot(Xo, Yo, 'k--', 'LineWidth', 1);
    plot(xs, ys, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    plot(XH{k}(1), YH{k}(1), 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    xlabel('x'); 
    ylabel('y');
    xlim([0 1]); 
    ylim([0 1]);
    title(sprintf('tau theta=%.2f', tauvals(k)));
end

function [xh, yh] = simrun(tauth, eta, N, th, P)
% policy-gradient learners with running-average baselines.
% a1 = 1: offense shoots a 3, a2 = 1: defense guards the paint (D2).
% update: dtheta = eta * [(r - b)(a - p) - tauth*theta*p*(1-p)]
    b1 = 0; b2 = 0;
    xh = zeros(N,1); yh = zeros(N,1);
    for i = 1:N
        x = 1/(1 + exp(-th(1)));
        y = 1/(1 + exp(-th(2)));
        xh(i) = x; 
        yh(i) = y;
        a1 = rand < x;
        a2 = rand < y;
        row = 1 + a1; % 1: shoot 2, 2: shoot 3
        col = 2 - a2; % 1: D2 (paint), 2: D3 (arc)
        r1 = (2 + a1) * (rand < P(row, col));   % points scored
        r2 = -r1;
        th(1) = th(1) + eta*((r1 - b1)*(a1 - x) - tauth*th(1)*x*(1-x));
        th(2) = th(2) + eta*((r2 - b2)*(a2 - y) - tauth*th(2)*y*(1-y));
        b1 = b1 + (r1 - b1)/i;% running averages
        b2 = b2 + (r2 - b2)/i;
    end
end
