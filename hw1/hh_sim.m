function [t, V, spk, rate] = hh_sim(Id, tmax, gK, gNa, skipfrac)
% HH_SIM  Classical Hodgkin-Huxley model under constant current injection.
%
%   [t,V,spk,rate] = HH_SIM(Id, tmax, gK, gNa, skipfrac)
%
%   Id        injected current density (uA/cm^2); defaults to 7, which
%             is above rheobase, so HH_SIM with no arguments gives
%             repetitive firing
%
%   Call it with a semicolon and capture what you need:
%       [t, V, spk, rate] = hh_sim;
%   A bare HH_SIM at the prompt displays ans, which is the first
%   output: 30001 elements of 0:0.01:300 scrolling past.
%   tmax      simulation duration (ms), default 300
%   gK        max potassium conductance (mS/cm^2), default 36
%   gNa       max sodium conductance (mS/cm^2), default 120
%   skipfrac  fraction of the trace discarded before counting spikes,
%             default 0.4, so that the onset transient is not mistaken
%             for repetitive firing
%
%   Returns the time vector, the voltage trace, the spike times (upward
%   crossings of 40 mV), and the mean firing rate in Hz over the counted
%   window (0 if fewer than two spikes). Potentials use the shifted
%   convention of the 1952 papers, so rest is near V = 0 and the spike
%   peak near V = 100.
%
%   The reported rheobase depends on how you define repetitive firing:
%   both skipfrac and the spike count you demand of it will move the
%   answer, because the underlying bifurcation is subcritical. State the
%   criterion you used rather than treating the number as absolute.
%
%   Regression checks with default parameters:
%     Id = 7 gives sustained firing at 58.3 Hz
%     Id = 5 gives no sustained firing
%   The 58.3 is solver independent: ode15s, ode45 and several stiff
%   solvers all agree to three digits, so a port that lands more than a
%   Hz away has a real bug. If Id = 5 fires, the gating rate functions
%   are wrong, most likely a sign in an exponent. The removable
%   singularities in am and an at V = 25 and V = 10 are patched below;
%   if you rewrite them, check that am(25) and an(10) come out finite.

if nargin < 1 || isempty(Id),       Id       = 7;    end
if nargin < 2 || isempty(tmax),     tmax     = 300;  end
if nargin < 3 || isempty(gK),       gK       = 36;   end
if nargin < 4 || isempty(gNa),      gNa      = 120;  end
if nargin < 5 || isempty(skipfrac), skipfrac = 0.4;  end

C = 1; gL = 0.3; ENa = 115; EK = -12; EL = 10.6;

V0 = 0;                                    % steady state gating at rest
y0 = [V0; ...
      am(V0)/(am(V0)+bm(V0)); ...
      ah(V0)/(ah(V0)+bh(V0)); ...
      an(V0)/(an(V0)+bn(V0))];

opts   = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);
[t, Y] = ode15s(@rhs, 0:0.01:tmax, y0, opts);
V      = Y(:,1);

i0  = max(2, round(skipfrac*numel(t)));
up  = find(V(i0:end-1) < 40 & V(i0+1:end) >= 40) + i0 - 1;
spk = t(up);

if numel(spk) >= 2
    rate = 1000*(numel(spk)-1)/(spk(end)-spk(1));   % Hz
else
    rate = 0;
end

    function dy = rhs(~, y)
        V_ = y(1); m = y(2); h = y(3); n = y(4);
        INa = gNa*m^3*h*(V_-ENa);
        IK  = gK *n^4  *(V_-EK);
        IL  = gL       *(V_-EL);
        dy = [ (Id - INa - IK - IL)/C; ...
               am(V_)*(1-m) - bm(V_)*m; ...
               ah(V_)*(1-h) - bh(V_)*h; ...
               an(V_)*(1-n) - bn(V_)*n ];
    end
end

% ---- gating rate functions (ms^-1), with removable singularities patched ----
function a = am(V)
if abs(V-25) < 1e-6, a = 1.0; else, a = 0.1*(25-V)/(exp((25-V)/10)-1); end
end
function b = bm(V), b = 4*exp(-V/18); end
function a = ah(V), a = 0.07*exp(-V/20); end
function b = bh(V), b = 1/(exp((30-V)/10)+1); end
function a = an(V)
if abs(V-10) < 1e-6, a = 0.1; else, a = 0.01*(10-V)/(exp((10-V)/10)-1); end
end
function b = bn(V), b = 0.125*exp(-V/80); end