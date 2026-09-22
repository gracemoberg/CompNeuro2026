% 4d: simulate poisson shot noise directly
% V(t) = sum_k w exp(-(t-t_k)/tau)

clear; close all; clc;

% parameters
tau = 20e-3; % 20 ms
meanV = 10; % desired mean voltage (mV)

nu_tau = [2 20 200]; % nu_tau vals to simualte

dt = 0.05e-3;  % 0.05 ms
Tburn = 1; % burn-in time (s)
Tstat = 100; % statistics time (s)

results = zeros(length(nu_tau), 6);

figure;
for k = 1:length(nu_tau)
    
    % simulation-specific parameters
    ntau = nu_tau(k);
    nu = ntau/tau;
    w = meanV/ntau;
    Nburn = round(Tburn/dt);

    V = 0;

    decay = exp(-dt/tau);

    for j = 1:Nburn

        % Poisson events during this timestep
        N = poissrnd(nu*dt);

        % existing voltage decay + new voltage injections
        V = V*decay + w*N;

    end

    Nstat = round(Tstat/dt);

    Vtrace = zeros(Nstat,1);

    for j = 1:Nstat

  
        N = poissrnd(nu*dt);
        % shot-noise update
        V = V*decay + w*N;

        Vtrace(j) = V;

    end


    % computed stats
    Vmean_measured = mean(Vtrace);
    Vvar_measured = var(Vtrace,1);
    centered = Vtrace - Vmean_measured;
    Vskew_measured = mean(centered.^3) / Vvar_measured^(3/2);

    % expected stats
    Vmean_theory = nu*w*tau;
    Vvar_theory = nu*w^2*tau/2;
    Vskew_theory = (2^(3/2)) / (3*sqrt(ntau));

    results(k,:) = [ntau, nu, w,  Vmean_measured, Vvar_measured, Vskew_measured];

    fprintf('\nnu*tau = %g\n',ntau);
    fprintf('nu = %.2f Hz\n',nu);
    fprintf('w = %.4f mV\n',w);
    fprintf('Mean:     measured = %.4f, theory = %.4f\n', Vmean_measured, Vmean_theory);
    fprintf('Variance: measured = %.4f, theory = %.4f\n', Vvar_measured, Vvar_theory);
    fprintf('Skewness: measured = %.4f, theory = %.4f\n', Vskew_measured, Vskew_theory);


    % histogram and gaussian for \nu \tau = 2, 200
    if k == 1 || k == 3

        if k == 1
            subplot(1,2,1);
        else
            subplot(1,2,2);
        end

        histogram(Vtrace,  'Normalization','pdf', 'NumBins',60);

        hold on;


        % gaussian with same theoretical mean and variance
        sigmaV = sqrt(Vvar_theory);
        x = linspace(min(Vtrace),max(Vtrace),500);
        gaussian =  1/(sqrt(2*pi)*sigmaV) .* exp(-(x-Vmean_theory).^2/(2*sigmaV^2));
        plot(x,gaussian,'LineWidth',2);
        xlabel('V (mV)');
        ylabel('density');
        title(sprintf('\\nu\\tau = %g',ntau));
        legend('simulation','matching gaussian');
        box on;

    end

end


%% Results table
fprintf('\n\n');
fprintf('============================================================\n');
fprintf('                 SIMULATION RESULTS\n');
fprintf('============================================================\n');

fprintf(' nu*tau       w (mV)       Mean       Variance     Skewness\n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(nu_tau)

    fprintf(' %5.0f      %8.4f      %8.4f     %8.4f     %8.4f\n', ...
        results(k,1), ...
        results(k,3), ...
        results(k,4), ...
        results(k,5), ...
        results(k,6));

end

fprintf('============================================================\n');