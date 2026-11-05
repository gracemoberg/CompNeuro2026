% prob 3: bumps in a Mexican-hat neural field
clear; close all; clc;

G.L = 20;  G.dx = 0.02;  G.dt = 0.02;
G.x = (-G.L:G.dx:G.L)';
G.Nx = numel(G.x);

% kernel sampled on all lags and transformed once
wfun = @(x) 2*exp(-abs(x)) - exp(-abs(x)/2);
lags = (-(G.Nx-1):(G.Nx-1))'*G.dx;
G.Nfft = 2^nextpow2(3*G.Nx - 2);
G.Wf = fft(wfun(lags), G.Nfft);

Wfun = @(x) sign(x).*(2*exp(-abs(x)/2) - 2*exp(-abs(x))); % W(x) = int_0^x w

% kappa = 0.3, four initial half-widths
kappa = 0.3;
[am, ap, lam_m, lam_p] = bump_theory(kappa);
fprintf('kappa = %.2f: a_- = %.5f, a_+ = %.5f, lambda_e(a_-) = %.4f, lambda_e(a_+) = %.4f\n', ...
        kappa, am, ap, lam_m, lam_p);

a0list = [am - 0.02, am + 0.02, ap - 0.05, ap + 0.05];
labels = {'a_- - 0.02', 'a_- + 0.02', 'a_+ - 0.05', 'a_+ + 0.05'};
Tend = 60;
tgrid = (0:round(Tend/G.dt))'*G.dt;

figure(19); hold on;
for k = 1:numel(a0list)
    u0 = Wfun(G.x + a0list(k)) - Wfun(G.x - a0list(k));
    a = run_field(u0, kappa, Tend, G);
    plot(tgrid, a, 'LineWidth', 1.5);
    fprintf('a0 = %s: a(0) = %.4f, a(T_end) = %.4f\n', labels{k}, a(1), a(end));
end
plot([0 Tend], [am am], 'k:');
plot([0 Tend], [ap ap], 'k:');
xlabel('t'); ylabel('half-width of active region');
legend(labels, 'Location', 'east'); grid on;

% relaxation rate of the wide bump, started 0.05 too wide
kappaList = [0.1 0.3 0.45];
fitFloor = 5e-3; % only fit while |a - a_+| is above this level
tFit0 = 2; % skip the initial transient
fprintf('%6s %10s %12s %12s %14s\n', 'kappa', 'a_+', 'fitted rate', 'lambda_e', '|a_end - a_+|');
figure(27); hold on;
for k = 1:numel(kappaList)
    kp = kappaList(k);
    [~, ap, ~, lam_p] = bump_theory(kp);
    Tend  = min(400, 12/abs(lam_p));
    tgrid = (0:round(Tend/G.dt))'*G.dt;
    a0 = ap + 0.05;
    u0 = Wfun(G.x + a0) - Wfun(G.x - a0);
    a  = run_field(u0, kp, Tend, G);

    delta = abs(a - ap);
    sel   = delta > fitFloor & tgrid >= tFit0;
    pf    = polyfit(tgrid(sel), log(delta(sel)), 1);
    fprintf('%6.2f %10.4f %12.4f %12.4f %14.2e\n', kp, ap, pf(1), lam_p, delta(end));

    semilogy(tgrid, delta, 'LineWidth', 1.5);
    set(gca, 'YScale', 'log');
end
xlabel('t'); ylabel('|a(t) - a_+|');
legend(arrayfun(@(q) sprintf('\\kappa = %.2f', q), kappaList, 'UniformOutput', false));
grid on;

function [am, ap, lam_m, lam_p] = bump_theory(kappa)
% half-widths a_-, a_+ and even-mode eigenvalues for the bumps.
S     = sqrt(1 - 2*kappa);
am    = -log((1 + S)/2);
ap    = -log((1 - S)/2);
lam_m =  S*(1 + S)/(1 - S*(1 + S)/2);
lam_p = -S*(1 - S)/(1 + S*(1 - S)/2);
end

function a = run_field(u0, kappa, Tend, G)
% FE for u_t = -u + w * H_steep(u - kappa) returns half width of active region
nSteps = round(Tend/G.dt);
a = zeros(nSteps + 1, 1);
u = u0;
a(1) = half_width(u, kappa, G);
for n = 1:nSteps
    H = 0.5*(1 + tanh(50*(u - kappa))); % = 1/(1+exp(-100(u-kappa)))
    r = G.dx*real(ifft(G.Wf .* fft(H, G.Nfft)));
    r = r(G.Nx : 2*G.Nx - 1); % sum_j w(x_i - x_j) H_j dx
    u = u + G.dt*(-u + r);
    a(n+1) = half_width(u, kappa, G);
end
end

function a = half_width(u, kappa, G)
% half the distance between the two edge crossings u = kappa
idx = find(u > kappa);
if isempty(idx)
    a = 0; % bump has died
    return
end
i1 = idx(1);  i2 = idx(end);
if i1 == 1 || i2 == G.Nx
    a = G.L; % reached the domain boundary
    return
end
xl = G.x(i1-1) + (kappa - u(i1-1))*G.dx/(u(i1) - u(i1-1));
xr = G.x(i2)   + (u(i2) - kappa)*G.dx/(u(i2) - u(i2+1));
a  = (xr - xl)/2;
end
