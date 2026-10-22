%% prob 1c: ocular dominance stripes on a ring
clear; clc; close all;
rng(1);

N = 200;
dt = 0.02;
T = 40;
nSteps = round(T/dt);

sigmaE = 1;
sigmaI = [2 3];
nTrials = 20;

x = 0:N-1;
d = min(x,N-x);     % periodic distance on the ring

meanPeriod = zeros(size(sigmaI));
stdPeriod = zeros(size(sigmaI));
continuumPeriod = zeros(size(sigmaI));
ringPeriod = zeros(size(sigmaI));
modeNumber = zeros(size(sigmaI));

figure('Color','w');

for s = 1:length(sigmaI)
    si = sigmaI(s);

    % Mexican hat kernel and its Fourier transform
    K = exp(-d.^2/(2*sigmaE^2)) ...
        - (sigmaE/si)*exp(-d.^2/(2*si^2));
    Khat = fft(K);

    % Predicted period from the continuous model
    kStar = sqrt(4*log(si/sigmaE)/(si^2-sigmaE^2));
    continuumPeriod(s) = 2*pi/kStar;

    % Find the fastest-growing mode available on the ring
    growth = real(Khat(2:N/2+1));
    [~,m] = max(growth);
    modeNumber(s) = m;
    ringPeriod(s) = N/m;

    % Simulate from 20 small random initial conditions
    periods = NaN(nTrials,1);

    for trial = 1:nTrials
        w = 0.01*randn(1,N);

        for it = 1:nSteps
            w = w + dt*real(ifft(Khat.*fft(w)));
            w = max(-1,min(1,w));
        end

        % Count transitions between right- and left-eye preference
        eye = sign(w);
        eye(eye == 0) = 1;
        nBoundaries = sum(eye ~= circshift(eye,1));
        nPairs = nBoundaries/2;

        if nPairs > 0
            periods(trial) = N/nPairs;
        end

        if trial == 1
            wFinal = w;
        end
    end

    meanPeriod(s) = mean(periods,'omitnan');
    stdPeriod(s) = std(periods,'omitnan');

    % Plot the final pattern from the first trial
    subplot(2,1,s);
    hold on;
    area(x,max(wFinal,0),'FaceColor',[0.85 0.4 0.4], ...
        'EdgeColor','none');
    area(x,min(wFinal,0),'FaceColor',[0.4 0.5 0.85], ...
        'EdgeColor','none');
    plot(x,wFinal,'k','LineWidth',1);

    xlabel('Cortical position a');
    ylabel('w_-(a)');
    title(sprintf('\\sigma_i = %g',si));
    legend('Right-eye preference','Left-eye preference', ...
        'Location','best');
    xlim([0 N-1]);
    ylim([-1.1 1.1]);
    box on;
end

% Compare the measured and predicted stripe periods
fprintf('\nStripe period summary\n');
fprintf('sigma_i   Continuum   Ring mode   Measured mean +/- SD\n');

for s = 1:length(sigmaI)
    fprintf('%6.1f   %9.3f   %9.3f   %7.3f +/- %.3f\n', ...
        sigmaI(s),continuumPeriod(s),ringPeriod(s), ...
        meanPeriod(s),stdPeriod(s));
end
