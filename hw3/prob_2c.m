% 2c: direct perturbation PRCs for LIF and FitzHugh--Nagumo

clear; close all; clc;
delta = 1e-3;
phi = (0:400)/401;% exclude phase 1 (the next spike/reset)

%% LIF: tau = 1, reset = 0, threshold = 1, I = 1.5
tau = 1; vr = 0; vth = 1; I = 1.5;
DeltaL = tau*log((I*tau-vr)/(I*tau-vth));
RL = zeros(size(phi));
for j = 1:numel(phi)
    tKick = phi(j)*DeltaL;
    vKick = I*tau+(vr-I*tau)*exp(-tKick/tau)+delta;
    if vKick >= vth
        Tprime = tKick;         
    else
        [~,~,te] = ode45(@(t,v) -v/tau+I, [tKick,tKick+2*DeltaL],vKick, odeset('RelTol',2e-10,'AbsTol',2e-12, 'Events',@lifThreshold,'MaxStep',DeltaL/50));
        Tprime = te(1); % absolute time since phase zero
    end
    RL(j) = (DeltaL-Tprime)/(DeltaL*delta);
end
Rexact = tau*exp(phi*DeltaL/tau)/(DeltaL*(I*tau-vr));
fprintf('LIF period: %.9f\n',DeltaL);
fprintf('LIF numerical min: %.6f at phase %.6f\n',min(RL),phi(find(RL==min(RL),1)));
fprintf('LIF numerical max on grid: %.6f at phase %.6f\n',max(RL),phi(find(RL==max(RL),1)));
fprintf('LIF analytic min and sup: %.6f, %.6f\n', tau/(DeltaL*(I*tau-vr)),tau/(DeltaL*(I*tau-vth)));

%% FHN: settle onto the cycle and define phase zero by upward v=0 crossing
fhn = @(t,x) [x(1)-x(1)^3/3-x(2)+0.5; 0.08*(x(1)+0.7-0.8*x(2))];
[~,~,te,xe] = ode45(fhn,[0,1200],[-1;1], odeset('RelTol',2e-10,'AbsTol',2e-12, 'Events',@upwardCrossing,'MaxStep',0.2));
lastPeriods = diff(te(end-5:end));
DeltaF = mean(lastPeriods);
x0 = xe(end,:).';               % state on the upward crossing
x0(1) = 0;
cycle = ode45(fhn,[0,DeltaF],x0, odeset('RelTol',2e-10,'AbsTol',2e-12, 'Events',@upwardCrossing,'MaxStep',0.2));

RF = zeros(size(phi)); RF8 = RF;
for j = 1:numel(phi)
    tKick = phi(j)*DeltaF;
    xKick = deval(cycle,tKick);
    xKick(1) = xKick(1)+delta;   % kick voltage; leave recovery variable fixed
    [~,~,tCross] = ode45(fhn,[tKick,8.5*DeltaF],xKick,odeset('RelTol',2e-10,'AbsTol',2e-12, 'Events',@upwardCrossing,'MaxStep',0.2));
    [~,k6] = min(abs(tCross-6*DeltaF));
    [~,k8] = min(abs(tCross-8*DeltaF));
    RF(j)  = (6*DeltaF-tCross(k6))/(DeltaF*delta);
    RF8(j) = (8*DeltaF-tCross(k8))/(DeltaF*delta);
end
[Rmin,jmin] = min(RF); [Rmax,jmax] = max(RF);
fprintf('\nFHN period: %.9f\n',DeltaF);
fprintf('FHN min: %.6f at phase %.6f\n',Rmin,phi(jmin));
fprintf('FHN max: %.6f at phase %.6f\n',Rmax,phi(jmax));

% estimate zero crossings by linear interp
z = [];
for j = 1:numel(phi)-1
    if RF(j)*RF(j+1) < 0
        z(end+1) = phi(j)-RF(j)*(phi(j+1)-phi(j))/(RF(j+1)-RF(j)); %#ok<SAGROW>
    end
end
fprintf('FHN zero-crossing phases:\n'); disp(z);

figure('Color', 'w');
plot(phi,RL,'b-','LineWidth',1.5); hold on;
plot(phi,Rexact,'r--','LineWidth',1.3);
xlabel('\phi'); ylabel('R(\phi)'); xlim([0,1]); grid on;
title('LIF');
legend('direct perturbation','analytic iPRC','Location','northwest');

figure('Color', 'w');
plot([phi,1],[RF,RF(1)],'b-','LineWidth',1.5); hold on;
plot([0,1],[0,0],'k:');
xlabel('\phi'); ylabel('R(\phi)'); xlim([0,1]); grid on;
title('FHN');

function [value,isterminal,direction] = lifThreshold(~,v)
    value = v-1; isterminal = 1; direction = 1;
end
function [value,isterminal,direction] = upwardCrossing(~,x)
    value = x(1); isterminal = 0; direction = 1;
end
