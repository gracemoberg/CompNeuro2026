% 3c - qif neurons with gap junction coupling
%
% dv1/dt = 1 + v1^2 + epsilon*(v2-v1)
% dv2/dt = 1 + v2^2 + epsilon*(v1-v2)


clear; close all; clc;


L = 200; % voltage cutoff
phi0 = 0.4; % phase lag

dt = 1e-4; % RK4 time step               
Tmax = 200;              

eps_plot = 0.02;
eps_fit = [0.005, 0.02, 0.05];
v1 = -L; % reset for neuron 1
theta2 = 1 - phi0; % phase lag for neuron 2
v2 = -cot(pi*theta2);

% simulate with eps = 0.02
[t,v1_hist,v2_hist,spike1,spike2] = simulate_QIF_RK4(eps_plot,L,dt,Tmax,v1,v2);

% measured phase lag: phi_k = (t_k^(2) - t_k^(1))/(t_{k+1}^(1) - t_k^(1))
ncycles = length(spike1)-1;
phi_k = nan(ncycles,1);
t_k = nan(ncycles,1);

for k = 1:ncycles

    t1 = spike1(k);
    t1_next = spike1(k+1);
    idx = find(spike2 > t1 & spike2 < t1_next,1,'first'); % find neuron 2 spike

    if ~isempty(idx)
        phi_k(k) = (spike2(idx)-t1)/(t1_next-t1);
        t_k(k) = t1;
    end
end

valid = ~isnan(phi_k); % if no neuron 2 spike found remove cycle

phi_k = phi_k(valid);
t_k = t_k(valid);
phi_theory = (1/pi)*atan(tan(pi*phi0).*exp(-2*eps_plot*t_k)); % theoretical prediction

% voltage trajectories for eps = 0.02 - does it work?
figure;
plot(t,v1_hist,'LineWidth',1);
hold on;
plot(t,v2_hist,'LineWidth',1);
xlabel('Time');
ylabel('Voltage');
legend('v_1','v_2','Location','best');
title('visualize sim for sanity check');
ylim([-L L]);
grid on;

% measured phase lag vs. prediction
figure;
plot(t_k,phi_k,'o','MarkerSize',4);
hold on;
plot(t_k,phi_theory,'-','LineWidth',2);
xlabel('time');
ylabel('phase lag');
legend('voltage','phase prediction', 'Location','best');
title(sprintf('measured phase lag \\epsilon = %.3f',eps_plot));
grid on;


% fit decay rates for eps = 0.005, 0.02, 0.05
decay_rates = zeros(size(eps_fit));

figure;
hold on;

for j = 1:length(eps_fit)

    epsilon = eps_fit(j);

    [t,v1_hist,v2_hist,spike1,spike2] = simulate_QIF_RK4(epsilon,L,dt,Tmax,v1,v2);

    % measure phase lag
    ncycles = length(spike1)-1;

    phi_k = nan(ncycles,1);
    t_k = nan(ncycles,1);

    for k = 1:ncycles

        t1 = spike1(k);
        t1_next = spike1(k+1);

        idx = find(spike2 > t1 & spike2 < t1_next,1,'first');

        if ~isempty(idx)

            phi_k(k) = (spike2(idx)-t1)/(t1_next-t1);

            t_k(k) = t1;

        end
    end

    valid = ~isnan(phi_k);

    phi_k = phi_k(valid);
    t_k = t_k(valid);

    % linearized theoretical prediction: log(tan(pi*phi_k)) = log(tan(pi*phi_0)) - 2*epsilon*t_k
    % theoretical slope = -2*epsilon

    y = log(tan(pi*phi_k));

    % remove noise and points too close to zero to prevent instability
    valid_fit = isfinite(y) & (phi_k > 1e-6) & (phi_k < 0.49);

    t_fit = t_k(valid_fit);
    y_fit = y(valid_fit);


    % linear regression to fit
    p = polyfit(t_fit,y_fit,1);
    slope = p(1);
    decay_rates(j) = -slope;
    plot(t_fit,y_fit,'o', 'MarkerSize',4, 'DisplayName', sprintf('simulated decay rate, \\epsilon = %.3f',epsilon));
    y_fit_line = polyval(p,t_fit);
    plot(t_fit,y_fit_line,'--', 'LineWidth',1.2, 'HandleVisibility','off');
    y_theory = log(tan(pi*phi0)) - 2*epsilon*t_fit;
    plot(t_fit,y_theory,'-', 'LineWidth',2, 'DisplayName', sprintf('theoretical: -2\\epsilon = %.3f',-2*epsilon));

end

xlabel('time');
ylabel('log(tan(\pi\phi_k))');
title('simulated and theoretical decay in phase diff');
legend('Location','best');
grid on;


% fitted decay rates
fprintf('\n');
fprintf('Electrical coupling decay rates\n');
fprintf(' epsilon       fitted rate       2*epsilon\n');

for j = 1:length(eps_fit)

    fprintf(' %7.3f       %10.6f       %10.6f\n', ...
        eps_fit(j), ...
        decay_rates(j), ...
        2*eps_fit(j));

end


% function: evolve system with rk4
function [t_hist,v1_hist,v2_hist,spike1,spike2] = simulate_QIF_RK4(epsilon,L,dt,Tmax,v1,v2)

    Nt = ceil(Tmax/dt); % timesteps

    t_hist = zeros(Nt+1,1);
    v1_hist = zeros(Nt+1,1);
    v2_hist = zeros(Nt+1,1);

    t_hist(1) = 0;
    v1_hist(1) = v1;
    v2_hist(1) = v2;

    % spike-time arrays
    spike1 = [];
    spike2 = [];

    for n = 1:Nt

        t = (n-1)*dt;
        y = [v1;v2];

        k1 = qif_rhs(y,epsilon);
        k2 = qif_rhs(y + 0.5*dt*k1, epsilon);
        k3 = qif_rhs(y + 0.5*dt*k2, epsilon);
        k4 = qif_rhs(y + dt*k3, epsilon);
        ynew = y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

        v1new = ynew(1);
        v2new = ynew(2);


        % detect threshold for neuron 1
        if v1new >= L

            frac = (L-v1)/(v1new-v1);  % linear interp for spike time
            tspike = t + frac*dt;
            spike1(end+1) = tspike;
            v1new = -L; % reset if spike

        end

        % detect threshold for neuron 2
        if v2new >= L

 
            frac = (L-v2)/(v2new-v2); % linear interp for spike time
            tspike = t + frac*dt;
            spike2(end+1) = tspike; 
            v2new = -L; % reset if spike

        end

        v1 = v1new;
        v2 = v2new;
        t_hist(n+1) = t + dt;
        v1_hist(n+1) = v1;
        v2_hist(n+1) = v2;

    end

end

function dydt = qif_rhs(y,epsilon)

    v1 = y(1);
    v2 = y(2);

    dv1 = 1 + v1^2 + epsilon*(v2-v1);
    dv2 = 1 + v2^2 + epsilon*(v1-v2);

    dydt = [dv1;dv2];

end