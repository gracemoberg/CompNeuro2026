%% prob 2d: two-population decision model
clear; close all; clc;
rng(1);

% Parameters
mu = 0.2;
sigma = 0.5;
I0 = 5;
dt = 0.002;
nTrials = 10000;

theta = 0.5:0.5:4;
lambda = [0, 0.2];

% Store simulation results
accuracy = zeros(2,length(theta));
meanTime = zeros(2,length(theta));

for j = 1:length(lambda)
    for k = 1:length(theta)

        [correct, times] = simulateRates(lambda(j), theta(k), ...
            mu, sigma, I0, dt, nTrials);

        accuracy(j,k) = mean(correct);
        meanTime(j,k) = mean(times);

        fprintf('lambda = %.1f, theta = %.1f: accuracy = %.3f, time = %.3f\n', ...
            lambda(j), theta(k), accuracy(j,k), meanTime(j,k));
    end
end

%% Theoretical predictions
thetaFine = linspace(0.01,4,150);
pTheory = zeros(2,length(thetaFine));
tTheory = zeros(2,length(thetaFine));

for j = 1:2
    for k = 1:length(thetaFine)
        [pTheory(j,k), tTheory(j,k)] = ...
            predictions(lambda(j), thetaFine(k), mu, sigma);
    end
end

% Accuracy limit for the unstable circuit
pLimit = 0.5*erfc(-mu/(sigma*sqrt(2*lambda(2))));

%% Plot accuracy and decision time
colors = [0 0.35 0.7; 0.85 0.25 0.1];


% Accuracy
figure('Color', 'w');
hold on;

for j = 1:2
    plot(thetaFine,pTheory(j,:),'-','Color',colors(j,:), ...
        'LineWidth',1.5);

    plot(theta,accuracy(j,:),'o','Color',colors(j,:), ...
        'MarkerFaceColor',colors(j,:));
end

plot([0 4],[pLimit pLimit],'k:');
xlabel('threshold \theta');
ylabel('P(correct)');
title('accuracy');
legend('theory, \lambda=0','simulation, \lambda=0', ...
    'theory, \lambda=0.2','simulation, \lambda=0.2', ...
    'accuracy limit for \lambda = 0.2','Location','southeast');
xlim([0 4.1]);
ylim([0.48 1]);
grid on;
box on;

% Decision time
figure('Color', 'w');
hold on;

for j = 1:2
    plot(thetaFine,tTheory(j,:),'-','Color',colors(j,:), ...
        'LineWidth',1.5);

    plot(theta,meanTime(j,:),'o','Color',colors(j,:), ...
        'MarkerFaceColor',colors(j,:));
end

xlabel('threshold \theta');
ylabel('mean decision time');
title('decision time');
legend('theory, \lambda=0','simulation, \lambda=0', ...
    'theory, \lambda=0.2','simulation, \lambda=0.2', ...
    'Location','northwest');
xlim([0 4.1]);
grid on;
box on;

%% Compare the two circuits at matched accuracy
% First find the unstable circuit's theoretical performance at theta = 4.
[pUnstable,tUnstable] = predictions(0.2,4,mu,sigma);

% Find the perfect-integrator threshold with the same accuracy.
thetaMatch = (sigma^2/mu)*log(pUnstable/(1-pUnstable));
[pPerfect,tPerfect] = predictions(0,thetaMatch,mu,sigma);

fprintf('\nmatched-accuracy comparison:\n');
fprintf('unstable circuit: theta = 4, P = %.4f, T = %.4f\n', ...
    pUnstable,tUnstable);
fprintf('perfect integrator: theta = %.4f, P = %.4f, T = %.4f\n', ...
    thetaMatch,pPerfect,tPerfect);
fprintf('time saved by perfect integrator: %.1f%%\n', ...
    100*(tUnstable-tPerfect)/tUnstable);

%% Functions

function [correct,times] = simulateRates(lambda,theta,mu,sigma,I0,dt,N)

    ws = (1+lambda)/2;
    wc = ws;

    r1 = I0*ones(N,1);
    r2 = I0*ones(N,1);

    correct = false(N,1);
    times = zeros(N,1);
    active = (1:N)';
    step = 0;

    while ~isempty(active)
        step = step+1;

        % Update both rates using their values at the start of the step.
        a = r1(active);
        b = r2(active);

        r1new = a+dt*(-a+ws*a-wc*b+I0+mu/2) ...
            +sigma*sqrt(dt)*randn(length(active),1);
        r2new = b+dt*(-b+ws*b-wc*a+I0-mu/2) ...
            +sigma*sqrt(dt)*randn(length(active),1);

        r1(active) = r1new;
        r2(active) = r2new;

        % Stop trials whose rate difference reaches either threshold.
        d = r1new-r2new;
        hit = abs(d)>=theta;
        finished = active(hit);

        correct(finished) = d(hit)>0;
        times(finished) = step*dt;

        active = active(~hit);
    end
end

function [p,t] = predictions(lambda,theta,mu,sigma)

    if lambda == 0
        % Perfect integrator: analytical formulas.
        p = 1/(1+exp(-mu*theta/sigma^2));
        t = (theta/mu)*tanh(mu*theta/(2*sigma^2));
        return
    end

    % Unstable integrator: evaluate the scale-function integrals.
    x = linspace(-theta,theta,20001);
    g = exp(-(lambda*x.^2+2*mu*x)/(2*sigma^2));

    S = cumtrapz(x,g);
    p = S((length(x)+1)/2)/S(end);

    q = 1./g;
    left = trapz(x(1:(length(x)+1)/2), ...
        S(1:(length(x)+1)/2).*q(1:(length(x)+1)/2));

    R = fliplr(cumtrapz(fliplr(x),fliplr(g)));
    right = trapz(x((length(x)+1)/2:end), ...
        R((length(x)+1)/2:end).*q((length(x)+1)/2:end));

    t = ((1-p)*left+p*right)/sigma^2;
end

