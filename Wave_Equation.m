%% Simulation Class
% 
%   Wave equations
%   Standing Wave 
%   Voltage and Current time evolution
%   Reflection Coefficient 
%   Impedance and Admittance

%% Init
clear;
close all;
clc;

%% System Parameters
% Generator voltage [V]
Vg = 2;     

% Generator impedance [ohm]
Zg = 50;    

% Line length [m]
l = 1.25;

% Line impedance [ohm]
Z0 = 50;

% Frequency [MHZ]
f = 300;

% Attenation constant
alpha = 0;

% Load Impedance [ohm]
ZL = 12.5;

%% 1. Simulation Variables

% wave length
lambda = 300 / f;

% Period
T = 1/f;

% reflection coefficient in load
if ZL == Inf
    % in case of open circuit
    pl = 1;
else
    pl = (ZL - Z0)/(ZL + Z0);
end;

% points in the line
x = linspace(0, l, 1000);

% phase constant
beta = 2*pi/lambda;

% propagation constant
gamma = alpha + 1j*beta;

% Distance from load
d = l - x;

%% 2. General phasor equations

% Incident & reflection voltage waves amplitude
V1 = Z0 / (Z0 + Zg) * Vg;
V2 = V1 * pl *  exp(-2*gamma*l);

% Incident & reflected voltage waves
Vi = V1 .* exp(-gamma .* x);
Vr = V2 .* exp(+gamma.*x);

% Incident & reflected current waves
Ii = Vi / Z0;
Ir = -Vr / Z0;

% Voltage & Current Waves
V = Vi + Vr;
I = Ii + Ir;


%% 3-4. Standing waves
% Amplitude of voltage & current standing wave
Vsw = abs(Vi + Vr);
Isw = abs(Ii + Ir);

% Voltage and Current Standing waves
figure(1)
plot(x/lambda, Vsw, 'LineWidth', 2)
hold on
plot(x/lambda, Isw*Z0, 'LineWidth', 2)
hold off
title(['Standing wave amplitude' sprintf('\n') 'ZL = ' num2str(ZL) ...
       '\Omega | \alpha = ' num2str(alpha) ' Np/m'])
xlabel('x/\lambda [m]');
ylabel('Amplitude');
legend('Voltage standing wave [V]', 'Current standing wave [A*50]')
xlim([0 l])
grid on

% Voltage Standing Wave Ratio
if ZL == 0 || ZL == Inf
    VSWR = Inf
else
    VSWR = (1 + abs(pl)*exp(-2*alpha*d)) / (1 - abs(pl)*exp(-2*alpha*d))
end;

%% 5-6. Wave propagation animation

% Animation parameters
% Number of periods to simulate
numT = 5;

% Number of points per period
pointsPerT = 32;

% Total number of simulation points
numSimPoints = numT * pointsPerT;

% Time step and time init
Tstep = T/pointsPerT;
t = 0;

% Envelope plotting
figure(2)
hVS = plot(x/lambda, Vsw, 'b--');
hold on
hIS = plot(x/lambda, Isw*Z0, 'r--');
plot(x/lambda, -Vsw, 'b--')
plot(x/lambda, -Isw*Z0, 'r--')
title(['Propagation of Voltage & Current Waves' sprintf('\n') ...
      'ZL = ' num2str(ZL) '\Omega | \alpha = ' num2str(alpha) ' Np/m'])
xlabel('x/\lambda [m]');
ylabel('Amplitude');
xlim([0 l])
grid on

% Voltage & Current waves @t=0s
hV = plot(x/lambda, real(V), 'LineWidth', 2);
hI = plot(x/lambda, real(I), 'LineWidth', 2); 
legend([hVS hIS hV hI], 'Standing voltage wave amplitude [V]', ...
                        'Standing current wave amplitude [A*50]', ...
                        'Voltage Wave [V]', '50 \times Current Wave [A]')

% Animation
for k = 1 : numSimPoints
    % increment time
    t = t + Tstep;
    
    % Set new data for voltage & current using the new time
    set(hV, 'Ydata', real(V.* exp(1j.*2*pi*f*t)));
    set(hI, 'Ydata', real(I*Z0.* exp(1j.*2*pi*f*t)));
    
    % delay 100 ms
    pause(0.1);
end

% For each voltage maximum there is a minimum of current.
% For a open circuit (ZL = Inf), the voltage is at is maximum at the load 
% but is zero at the generator. The current is zero at the load but is
% maximum in the generator.  The VSWR = Inf, and reflection coefficient = 1.
% For a short circuit, the voltage is at is zero at the load but is maximum
% at the generator. The current is  maximum at the load but is zero in the
% generator. The VSWR = Inf, and reflection coefficient = -1
% For ZL = Z0, the current and voltage waves are in phase, and there is no
% standing wave. The VSWR = 1, and reflection coefficient = 0.

%% 8. Reflection coefficient
p = Vr./Vi;

% Indexes for unique points 
idx = 1:lambda/2 * length(x);

figure(3)
plot(p(idx), 'r')
title(['Reflection coefficient along the line' sprintf('\n') ...
      'ZL = ' num2str(ZL) '\Omega | \alpha = ' num2str(alpha) ' Np/m'])
xlabel('Re(\rho)')
ylabel('Im(\rho)')




