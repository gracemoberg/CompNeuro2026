% 2c: numerically compute the voltage where the slope vanishes and plot
% I_N(V) over the specified interval

V0 = 16.13; % mV
EN = 0; % mV
Mg = 1; % mM

B = @(V) 1 ./ (1 + (Mg/3.57).*exp(-V./V0));

% brackets in dI_N/dV
dIdV_factor = @(V) 1 + ((V - EN) .* (1 - B(V))) ./ V0;

% compute slope zero
V_zero = fzero(dIdV_factor, [-90 0]);

fprintf('Slope vanishes at V = %.3f mV\n', V_zero);


clear;
clc;
close all;

V0 = 16.13; % mV
EN = 0; % mV
gbar = 1; % nS

Mg_values = [0, 0.1, 1]; % mM
V = linspace(-90, 0, 1000); % mV

figure;
hold on;

for Mg = Mg_values

    B = 1 ./ (1 + (Mg/3.57).*exp(-V./V0));

    % NMDA current
    % nS * mV = pA
    IN = gbar .* B .* (V - EN);

    plot(V, IN, 'LineWidth', 2, 'DisplayName', sprintf('[Mg^{2+}] = %.1f mM', Mg));

end

% Find slope-zero voltage for Mg = 1 mM
Mg = 1;

Bfun = @(V) 1 ./ (1 + (Mg/3.57).*exp(-V./V0));
dIdV_factor = @(V) 1 + ((V - EN).*(1 - Bfun(V)))./V0;

V_zero = fzero(dIdV_factor, [-90 0]);

% current at the slope-zero point
B_zero = Bfun(V_zero);
I_zero = gbar .* B_zero .* (V_zero - EN);
plot(V_zero, I_zero, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k',  'DisplayName', sprintf('zero slope at %.3f mV', V_zero));

% Plot formatting
xlabel('membrane voltage V (mV)');
ylabel('NMDA current I_N (pA)');
title('I_N(V) vs. Voltage');

legend('Location', 'best');
grid on;
box on;

xlim([-90 0]);

hold off;

fprintf('For [Mg2+] = 1 mM:\n');
fprintf('dI_N/dV = 0 at V = %.3f mV\n', V_zero);