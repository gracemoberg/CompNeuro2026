%% hw1 problem 5
% find rheobase with bisection method 

%% part a

clear;
close all;
clc;

% parameters
tmax = 300;
gK = 36;
gNa = 120;
skipfrac = 0.4;

% repetitive firing criterion: at least 2 spikes after skipfrac

% initial bounds
Id_low = 5;      % no sustained firing
Id_high = 7;     % sustained repetitive firing

tol = 1e-3; % bisection tolerance

% bisect
while (Id_high - Id_low) > tol

    Id_mid = (Id_low + Id_high)/2;

    % run hh sim
    [~, ~, spk, ~] = hh_sim(Id_mid, tmax, gK, gNa, skipfrac);

    % do we meet the repetitive firign criterion?
    if length(spk) >= 2
        % repetitive
        Id_high = Id_mid;
    else
        % not repetitive
        Id_low = Id_mid;
    end
end

% estimate rheobase - bisection midpoint
Id_rheo = (Id_low + Id_high)/2;

fprintf('Rheobase = %.6f \n', Id_rheo);
fprintf('Lower bound = %.6f \n', Id_low);
fprintf('Upper bound = %.6f \n', Id_high);


%% plot voltage traces just below and just above rheobase

% choose currents immediately below and above the threshold
Id_below = Id_low;
Id_above = Id_high;

[t_below, V_below, spk_below, ~] = hh_sim(Id_below, tmax, gK, gNa, skipfrac);

[t_above, V_above, spk_above, ~] = hh_sim(Id_above, tmax, gK, gNa, skipfrac);


figure;

plot(t_below, V_below, 'LineWidth', 1.2);
hold on;
plot(t_above, V_above, 'LineWidth', 1.2);

xlabel('time (ms)');
ylabel('voltage (mV)');
title('Voltage traces near rheobase');

legend(sprintf('Below: I_d = %.4f', Id_below), sprintf('Above: I_d = %.4f', Id_above));

% explicitly give both traces identical limits
xlim([0 tmax]);

ymin = min([V_below; V_above]);
ymax = max([V_below; V_above]);
ylim([ymin ymax]);

grid on;

fprintf('\nBelow rheobase: %d spikes after transient\n', length(spk_below));
fprintf('Above rheobase: %d spikes after transient\n', length(spk_above));

%% part b
% compute and plot the f-I curve

% currents from just above rheobase to twice rheobase
Id_values = linspace(Id_rheo + 0.01, 2*Id_rheo, 20);

% firing rates
rates = zeros(size(Id_values));

% compute firing rate at each current
for k = 1:length(Id_values)

    [~, ~, ~, rates(k)] = hh_sim(Id_values(k), tmax, gK, gNa, skipfrac);

end

% onset rate
onset_rate = rates(1);

fprintf('\nOnset firing rate = %.3f Hz\n', onset_rate);

% plot f-I curve
figure;

plot(Id_values, rates, 'o-', 'LineWidth', 1.2);

xlabel('injected current I_d (\muA/cm^2)');
ylabel('firing rate (Hz)');
title('f-I curve');

grid on;

%% part c
% change gK and recompute rheobase

gK_old = gK;
gK = 30;

% initial bounds
Id_low = 1;
Id_high = 7;

% bisect
while (Id_high - Id_low) > tol

    Id_mid = (Id_low + Id_high)/2;

    % run hh sim
    [~, ~, spk, ~] = hh_sim(Id_mid, tmax, gK, gNa, skipfrac);

    % repetitive firing criterion
    if length(spk) >= 2
        Id_high = Id_mid;
    else
        Id_low = Id_mid;
    end

end

% rheobase for gK = 30
Id_rheo_new = (Id_low + Id_high)/2;


%% f-I curve for gK = 30

Id_values_new = linspace(Id_rheo_new + 0.01, 2*Id_rheo_new, 20);

rates_new = zeros(size(Id_values_new));

for k = 1:length(Id_values_new)

    [~, ~, ~, rates_new(k)] = hh_sim(Id_values_new(k), tmax, gK, gNa, skipfrac);

end


%% rheobase comparison table

fprintf('\nRheobase comparison:\n');
fprintf('gK (mS/cm^2)    Rheobase\n');
fprintf('36              %.3f\n', Id_rheo);
fprintf('30              %.3f\n', Id_rheo_new);


%% rheobase comparison figure
figure;

plot(Id_values, rates, 'o-', 'LineWidth', 1.2);
hold on;

plot(Id_values_new, rates_new, 's-', 'LineWidth', 1.2);

% mark rheobases
xline(Id_rheo, '--');
xline(Id_rheo_new, '--');
xlim([2,13]);
xlabel('injected current I_d (\muA/cm^2)');
ylabel('firing rate (Hz)');
title('f-I curves for g_K = 36 and g_K = 30');

legend('g_K = 36', 'g_K = 30', 'Rheobase, g_K = 36', 'Rheobase, g_K = 30', 'Location', 'northwest');

grid on;