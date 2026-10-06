%% 1c: firing number for periodically forced LIF 
% evolve subthreshold and evaluate first threshold crossing with fzero

clear; close all; clc;
epsilon = 0.5;
Ivalues = linspace(0.9,3,421);
transient = 500; % discard forcing period / burn in   
observation = 1000; % forcing periods for spike counts
firingNumber = zeros(size(Ivalues));

for j = 1:numel(Ivalues)
    spikes = lifSpikes(Ivalues(j),epsilon,transient+observation);
    firingNumber(j) = sum(spikes >= transient & spikes < transient+observation)/observation;
end

unforced = zeros(size(Ivalues));
mask = Ivalues > 1;
unforced(mask) = 1./log(Ivalues(mask)./(Ivalues(mask)-1));
Istar = 1/(1-exp(-1));
boundaries = Istar + [-1,1]*epsilon/sqrt(1+4*pi^2);

figure('Color','w');
plot(Ivalues,firingNumber,'b-','LineWidth',1.6); hold on;
plot(Ivalues,unforced,'--','Color',[.45 .45 .45],'LineWidth',1.3);
ylimits = [0,2.6];
plot([boundaries(1),boundaries(1)],ylimits,'r:','LineWidth',1.3);
plot([boundaries(2),boundaries(2)],ylimits,'r:','LineWidth',1.3);
xlabel('I'); 
ylabel('spikes/forcing period');
xlim([0.9,3]); ylim(ylimits); grid on;
legend('forced: epsilon = 0.5','unforced','1:1 boundaries',  'Location','northwest');
title('firing number (plateaus are mode-locked)');

% locate five other plateaus on the sampled grid, a p/q plateau means p spikes during q forcing periods
% resolution is 1/observation
ratios = [1/3,1/2,2/3,3/4,2];
labels = {'1/3','1/2','2/3','3/4','2/1'};
fprintf('candidate 1:1 boundaries: %.8f, %.8f\n',boundaries);
for j = 1:numel(ratios)
    ix = find(abs(firingNumber-ratios(j)) <= 1.01/observation);
    if ~isempty(ix)
        % longest contiguous set of sampled currents
        starts = [1,find(diff(ix)>1)+1];
        ends = [starts(2:end)-1,numel(ix)];
        [~,k] = max(ends-starts+1);
        group = ix(starts(k):ends(k));
        mid = group(ceil(numel(group)/2));
        text(Ivalues(mid),ratios(j)+0.08,labels{j}, 'HorizontalAlignment','center');
        fprintf('%s plateau: sampled I from %.4f to %.4f\n', labels{j},Ivalues(group(1)),Ivalues(group(end)));
    end
end
ix = find(abs(firingNumber-1) < 0.5/observation);
if ~isempty(ix)
    fprintf('Measured 1:1 plateau: sampled I from %.4f to %.4f\n', Ivalues(ix(1)),Ivalues(ix(end)));
end

% measure phase after the transient at I = 1.6, starting from v(0) = 0
spikes = lifSpikes(1.6,epsilon,transient+20);
phases = mod(spikes(spikes >= transient),1);
fprintf('Last five spike phases at I = 1.6:\n');
disp(phases(max(1,end-4):end).');

% function to compute lif spikes
function spikes = lifSpikes(I,epsilon,tEnd)
    G = @(t) I + epsilon/(1+4*pi^2)* (sin(2*pi*t)-2*pi*cos(2*pi*t));
    options = optimset('TolX',1e-11,'Display','off'); % high precision and suppress comments in fzero
    phaseBreaks = 1;
    q = (1-I)/epsilon;
    if abs(q) < 1
        a = asin(q);
        phaseBreaks = unique(sort([mod(a/(2*pi),1), mod((pi-a)/(2*pi),1),1]));
        phaseBreaks = phaseBreaks(phaseBreaks > 1e-12);
    end

    t0 = 0; v0 = 0;              
    spikes = zeros(ceil(5*tEnd),1);
    n = 0;
    for k = 0:ceil(tEnd)-1
        for phase = phaseBreaks
            tRight = min(k+phase,tEnd);
            if tRight <= t0
                continue;
            end
            voltage = @(t) G(t)+(v0-G(t0))*exp(-(t-t0)); % exact voltage simulation between spikes
            while voltage(tRight) >= 1
                tSpike = fzero(@(t) voltage(t)-1,[t0,tRight],options); % rootfinding
                n = n+1;
                spikes(n) = tSpike;
                t0 = tSpike; v0 = 0;     
                voltage = @(t) G(t)+(v0-G(t0))*exp(-(t-t0));
            end
            v0 = voltage(tRight);
            t0 = tRight;
        end
    end
    spikes = spikes(1:n);
end
