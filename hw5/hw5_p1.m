% prob 1: E-I oscillation near a Hopf bifurcation
% Wilson-Cowan with r_E = r_I = 0, kappa = 0, w_II = 0.

clear; close all; clc;

P.eta = 4;  
P.tauE = 5;  
P.wEI = 10;  
P.wIE = 10;
P.Eb  = 0.1; 
P.Ib  = 0.1;

% tau_I = 10, vary epsilon
tauI = 10;
epsList = [0.005 0.01 0.02 0.05 0.1 0.2];
res1 = zeros(numel(epsList), 4);   % period, amplitude, lag (deg), settled flag
for k = 1:numel(epsList)
    r = simulate_wc(tauI, epsList(k), P);
    res1(k,:) = [r.period, r.amp, r.lagDeg, r.settled];
end
T1 = table(epsList(:), res1(:,1), res1(:,2), res1(:,3), res1(:,4), ...
    'VariableNames', {'eps','period_ms','amplitude','lag_deg','settled'});
disp('tau_I = 10');
disp(T1);

pfit = polyfit(log(epsList(:)), log(res1(:,2)), 1);
fprintf('fitted log-log slope of amplitude vs eps: %.4f\n', pfit(1));

figure;
loglog(epsList, res1(:,2), 'o-', 'LineWidth', 1.5); hold on;
loglog(epsList, exp(pfit(2))*epsList.^pfit(1), '--', 'LineWidth', 1.5);
xlabel('\epsilon'); ylabel('amplitude  max E - min E');
legend('simulation', sprintf('fit, slope = %.3f', pfit(1)), 'Location', 'northwest');
grid on;

% epsilon = 0.02, vary tau_I
tauList = [2 3 5 7.5 10 15 20];
eps2    = 0.02;
gE = P.eta*P.Eb*(1-P.Eb);  gI = P.eta*P.Ib*(1-P.Ib);
Bc = gE*gI*P.wEI*P.wIE/P.tauE;
res2 = zeros(numel(tauList), 3);
for k = 1:numel(tauList)
    r  = simulate_wc(tauList(k), eps2, P);
    om = sqrt(Bc/tauList(k) - 1/tauList(k)^2);
    res2(k,:) = [1000/r.period, 1000*om/(2*pi), r.settled];
end
T2 = table(tauList(:), res2(:,1), res2(:,2), res2(:,3), ...
    'VariableNames', {'tauI','f_sim_Hz','f_lin_Hz','settled'});
disp('eps = 0.02');
disp(T2);

figure;
plot(tauList, res2(:,1), 'o-', 'LineWidth', 1.5); hold on;
plot(tauList, res2(:,2), 's--', 'LineWidth', 1.5);
xlabel('\tau_I (ms)'); ylabel('frequency (Hz)');
legend('simulated', '\omega/2\pi from (b)'); grid on;

% time series, tau_I = 10, eps = 0.1
r = simulate_wc(10, 0.1, P);
sel = r.tw >= r.tw(end) - 5*r.period;
figure;
plot(r.tw(sel), r.E(sel), 'LineWidth', 1.5); hold on;
plot(r.tw(sel), r.I(sel), 'LineWidth', 1.5);
xlabel('t (ms)'); ylabel('activity'); legend('E(t)', 'I(t)'); grid on;
fprintf('period = %.4f ms, lag = %.3f deg\n', r.period, r.lagDeg);

function r = simulate_wc(tauI, eps, P)
% w_EE = (1+eps)*w_EE^Hopf and measure period, amplitude and lag of the I peaks behind the E peaks once settled.
eta = P.eta; tauE = P.tauE; wEI = P.wEI; wIE = P.wIE; Eb = P.Eb; Ib = P.Ib;
gE  = eta*Eb*(1-Eb);  gI = eta*Ib*(1-Ib);

wEE = (1+eps)*(1 + tauE/tauI)/gE;
pE  = log(Eb/(1-Eb))/eta - wEE*Eb + wEI*Ib;
pI  = log(Ib/(1-Ib))/eta - wIE*Eb;

f   = @(x) 1./(1 + exp(-eta*x));
rhs = @(t,y) [ (-y(1) + f(wEE*y(1) - wEI*y(2) + pE))/tauE;
               (-y(2) + f(wIE*y(1) + pI))/tauI ];

% linear growth rate of the amplitude, and linear period
mu   = 0.5*eps*(1 + tauE/tauI)/tauE;
Bc   = gE*gI*wEI*wIE/tauE;
om   = sqrt(Bc/tauI - 1/tauI^2);
Tlin = 2*pi/om;

tmax = 25/mu; % about 50 relaxation times of the amplitude
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
sol  = ode45(rhs, [0 tmax], [Eb + 1e-3; Ib], opts);

% fine uniform grid on the last 40 periods
dtg = 0.005;
tw  = (tmax - 40*Tlin : dtg : tmax)';
Y   = deval(sol, tw)';
E   = Y(:,1);  I = Y(:,2);

tpE = peak_times(tw, E);
tpI = peak_times(tw, I);
period = mean(diff(tpE));

lags = nan(size(tpE));
for k = 1:numel(tpE)
    j = find(tpI > tpE(k), 1, 'first');
    if ~isempty(j), lags(k) = tpI(j) - tpE(k); end
end
lagDeg = 360*mean(lags, 'omitnan')/period;

% check the amplitude is no longer changing (first vs last 10 periods)
nA = round(10*Tlin/dtg);
ampA = max(E(1:nA))       - min(E(1:nA));
ampB = max(E(end-nA:end)) - min(E(end-nA:end));

r.period  = period;
r.amp     = max(E) - min(E);
r.lagDeg  = lagDeg;
r.settled = abs(ampA - ampB)/ampB < 1e-3;
r.tw = tw; r.E = E; r.I = I;
end

function tp = peak_times(t, y)
% times of local maxima of y, refined by parabolic interpolation.
k  = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;
dt = t(2) - t(1);
ym = y(k-1); y0 = y(k); yp = y(k+1);
d  = 0.5*(ym - yp)./(ym - 2*y0 + yp);
tp = t(k) + d*dt;
end
