% prob 2: binocular rivalry with adaptation

clear; close all; clc;

P.w = 0.8;  P.beta = 1;  P.g = 2;  P.tau_a = 120;  P.eta = 1000;
y0   = [1; 0; 0; 0];  % (uL, uR, aL, aR)
tEnd = 30*P.tau_a;
nDiscard = 3; % discard first three dominance durations

% symmetric case I_L = I_R = 1.3
P.IL = 1.3;  P.IR = 1.3;
[t, Y] = run_rivalry(P, y0, tEnd);
B = (P.IL - P.beta)/P.g;

figure;
subplot(2,1,1);
plot(t, Y(:,1), t, Y(:,2), 'LineWidth', 1.2);
ylabel('rates'); legend('u_L', 'u_R'); ylim([-0.1 1.1]);
subplot(2,1,2);
plot(t, Y(:,3), t, Y(:,4), 'LineWidth', 1.2); hold on;
plot([0 tEnd], [B B],     'k--');
plot([0 tEnd], [1-B 1-B], 'k--');
text(0.02*tEnd, B + 0.03,     'B');
text(0.02*tEnd, 1 - B + 0.03, '1-B');
xlabel('t'); ylabel('adaptation');
legend('a_L', 'a_R', 'Location', 'east');

% I_L = 1.3, 1.45, 1.6 at I_R = 1.3
ILlist = [1.3 1.45 1.6];
P.IR = 1.3;
fprintf('\n(d)(ii): w = %.2f\n', P.w);
fprintf('%6s %10s %10s %10s %10s\n', 'I_L', 'T_L sim', 'T_R sim', 'T_L (c)', 'T_R (c)');
for k = 1:numel(ILlist)
    P.IL = ILlist(k);
    [t, Y] = run_rivalry(P, y0, tEnd);
    [TL, TR] = dominance_times(t, Y(:,1), Y(:,2), nDiscard);
    BL = (P.IL - P.beta)/P.g;  BR = (P.IR - P.beta)/P.g;
    TLth = P.tau_a*log((1 - BL)/BR);
    TRth = P.tau_a*log((1 - BR)/BL);
    fprintf('%6.3f %10.3f %10.3f %10.3f %10.3f\n', P.IL, TL, TR, TLth, TRth);
end

% w = 0.3 (below the window), I_L = I_R = 1.3
P2 = P;  P2.w = 0.3;  P2.IL = 1.3;  P2.IR = 1.3;
[t, Y] = run_rivalry(P2, y0, tEnd);
[TL, TR] = dominance_times(t, Y(:,1), Y(:,2), nDiscard);
Tth = P2.tau_a*log((P2.w + P2.IL)/(P2.g - P2.w - P2.IL));
fprintf('\n(d)(iii): w = %.2f\n', P2.w);
fprintf('T_L sim = %.3f, T_R sim = %.3f, formula = %.3f\n', TL, TR, Tth);
fprintf('relative differences (sim - formula)/formula: L %.3f, R %.3f\n', ...
        (TL - Tth)/Tth, (TR - Tth)/Tth);

figure;
subplot(2,1,1); plot(t, Y(:,1), t, Y(:,2)); ylabel('rates'); legend('u_L','u_R');
subplot(2,1,2); plot(t, Y(:,3), t, Y(:,4)); ylabel('adaptation');
xlabel('t'); legend('a_L','a_R');


function [t, Y] = run_rivalry(P, y0, tEnd)
% integrate the four-variable rivalry model 
f   = @(x) 0.5*(1 + tanh(0.5*P.eta*x)); 
rhs = @(t,y) [ -y(1) + f(P.w*y(1) - P.beta*y(2) - P.g*y(3) + P.IL);
               -y(2) + f(P.w*y(2) - P.beta*y(1) - P.g*y(4) + P.IR);
               (-y(3) + y(1))/P.tau_a;
               (-y(4) + y(2))/P.tau_a ];
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
sol  = ode15s(rhs, [0 tEnd], y0, opts); % use stiff solver
t    = (0:0.05:tEnd)';
Y    = deval(sol, t)';
end

function [TL, TR] = dominance_times(t, uL, uR, nDiscard)
% percept is dominant when its rate exceeds the other's by more than 1/2;
% in between, the previous dominant percept is carried forward
D = uL - uR;
state = zeros(size(D));
state(D >  0.5) =  1;
state(D < -0.5) = -1;
if state(1) == 0, state(1) = 1; end
cur = state(1);
for k = 2:numel(state)
    if state(k) == 0
        state(k) = cur;
    else
        cur = state(k);
    end
end
idx  = find(diff(state) ~= 0) + 1; % hand-over indices
tH   = t(idx);
sNew = state(idx); % percept that takes over at each hand-over
dur  = diff(tH); % duration of the percept taking over
sDur = sNew(1:end-1);
dur  = dur(nDiscard+1:end);
sDur = sDur(nDiscard+1:end);
TL = mean(dur(sDur ==  1));
TR = mean(dur(sDur == -1));
end
