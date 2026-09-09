%% HW1 problem 4 part c 
% PLIAF model w synaptic transient input

clear
clc
close all

%% set parameter values

R = 1;              % resistance
tau_m = 10;         % membrane time constant (ms)
tau_s = 2;          % synaptic time constant (ms)
u_rest = 0;         % resting voltage (mV)
u_th = 1;           % threshold voltage (mV)

%% critical current and theoretical peak time

r = tau_s/tau_m;

Icrit = (u_th-u_rest)/R * r^(1/(r-1));

tpeak_theory = (tau_m*tau_s)/(tau_s-tau_m) * log(tau_s/tau_m);

fprintf('Theoretical t_peak = %.4f ms\n', tpeak_theory);
fprintf('Critical current I0_crit = %.4f\n', Icrit);


I0_values = [6 1.0*Icrit 9]; % I0 amplitudes in vector


h = 0.0001;           % time step (ms)
tfinal = 50;        % total time (ms)
t = 0:h:tfinal;

%% voltage plot
figure
hold on

peak_voltages = zeros(1,3);
peak_times = zeros(1,3);

for j = 1:3

    I0 = I0_values(j);

    % initialize voltage
    u = zeros(size(t));
    u(1) = u_rest;

    I = @(t) I0*exp(-t/tau_s); % input current

    % LIAF right-hand side
    f = @(t,u) (-(u-u_rest) + R*I(t))/tau_m;

    %% RK4 integration
    for n = 1:length(t)-1

        k1 = f(t(n),u(n));

        k2 = f(t(n)+h/2, u(n)+h*k1/2);

        k3 = f(t(n)+h/2, u(n)+h*k2/2);

        k4 = f(t(n)+h, u(n)+h*k3);

        u(n+1) = u(n) + h*(k1 + 2*k2 + 2*k3 + k4)/6;

    end

    % find numerical peak
    [peak_voltages(j),peak_index] = max(u);
    peak_times(j) = t(peak_index);

    plot(t,u,'LineWidth',1.5)

end

%% threshold line
yline(u_th,'--','Threshold','LineWidth',1.5)

xlabel('t (ms)')
ylabel('u (mV)')
title('PLIAF single-neuron voltage')

legend('I0 < I_0^{crit}', 'I_0 = I_0^{crit}', 'I_0 > I_0^{crit}', 'Threshold', 'Location','best')

grid on
hold off

%% display numerical results

fprintf('\nNumerical results:\n')
fprintf('I0              t_peak (ms)       u_max (mV)\n')

for j = 1:3
    fprintf('%.4f          %.4f           %.4f\n', I0_values(j),peak_times(j),peak_voltages(j));
end


%% critical current as a function of tau_s
tau_s_values = logspace(log10(0.5),log10(100),500);

Icrit_values = zeros(size(tau_s_values));

for j = 1:length(tau_s_values)

    ts = tau_s_values(j);
    r = ts/tau_m;

    % at r = 1, use the limiting value e
    if abs(r-1) < 1e-8
        factor = exp(1);
    else
        factor = r^(1/(r-1));
    end

    Icrit_values(j) = (u_th-u_rest)/R * factor;

end

%% plot Icrit versus tau_s
figure

loglog(tau_s_values,Icrit_values,'LineWidth',1.5)

xlabel('log(\tau_s) (ms)')
ylabel('log(I_0^{crit})')
title('I_0^{crit} vs. \tau_s')

grid on