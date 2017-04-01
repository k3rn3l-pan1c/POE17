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
alpha = 0.0;

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
legend('Voltage standing wave [V]', 'Current standing wave [A\times50]')
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
numT = 3;

% Number of points per period
pointsPerT = 64;

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
                        'Standing current wave amplitude [A\times50]', ...
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

%% 1. Line Impedance
Zin = V ./ I;

figure(4)
plot(x/lambda, real(Zin), x/lambda, imag(Zin))
xlim([0 l])
title(['Input Impedance along the line' sprintf('\n') ...
      'ZL = ' num2str(ZL) '\Omega | \alpha = ' num2str(alpha) ' Np/m'])
xlabel('x/\lambda')
ylabel('Impedance (\Omega)')
legend('Resistance', 'Reactance')
grid on

%% 2. Match the Line Impedance
% To match the line using a serial reactance, the reflection coeficient
% must be zero. The distances x/\lambda adequated to match the line by
% using a serial reactance must have a resistance matching Z0, the
% caracterisitic impedance

% Input resistance equals characteristic resistance. Round Zin to units to
% accuracy in results

% Postive reactance implies a inductive reactance (inductor)
Lidx = find(round(real(Zin), 0) == real(Z0) & imag(Zin) > 0);

% Negative reactance implies a capacitive reactance (capacitor)
Cidx = find(round(real(Zin), 0) == real(Z0) & imag(Zin) < 0);


figure(5)
plot(x/lambda, real(Zin), x/lambda, imag(Zin))
hold on
plot(x(Lidx)/lambda, zeros(1, length(Lidx)), 'kx')
plot(x(Cidx)/lambda,  zeros(1, length(Cidx)), 'c*')
hold off
xlim([0 l])
title(['Impedance matching using a serial reactance' sprintf('\n') ...
      'ZL = ' num2str(ZL) '\Omega | \alpha = ' num2str(alpha) ' Np/m'])
xlabel('x/\lambda')
ylabel('Impedance (\Omega)')
legend('Resistance', 'Reactance', 'inductive reactance - Inductor', ...
       'Capacitive reactance -  Capacitor')
grid on



%% 3. Line Admittance
Yin = 1./Zin;

figure(6)
plot(x/lambda, real(Yin), x/lambda, imag(Yin))
xlim([0 l])
title(['Input Admittance along the line' sprintf('\n') ...
      'ZL = ' num2str(ZL) '\Omega | \alpha = ' num2str(alpha) ' Np/m'])
xlabel('x/\lambda')
ylabel('Admittance (S)')
legend('Conductance', 'Susceptance')
grid on

%% 4. Match the Line Admittance
% To match the line using a parallel susceptance, the reflection coeficient
% must be zero. The distances x/\lambda adequated to match the line by
% using a parallel susceptance must have a conductance matching Y0, the
% caracterisitic admittance

% Characteristic Admittance
Y0 = 1/Z0;

% Input admittance equals characteristic admittance. Round Yin to units to
% accuracy in results

% Positive susceptance implies a capacitive susceptance (capacitor)
Cidx = find(round(real(Yin), 3) == real(Y0) & imag(Yin) > 0);

% Negative susceptance implies a inductive susceptance (inductor)
Lidx = find(round(real(Yin), 3) == real(Y0) & imag(Yin) < 0);


figure(7)
plot(x/lambda, real(Yin), x/lambda, imag(Yin))
hold on
plot(x(Lidx)/lambda, zeros(1, length(Lidx)), 'kx')
plot(x(Cidx)/lambda,  zeros(1, length(Cidx)), 'c*')
hold off
xlim([0 l])
title(['Admittance matching using a parallel susceptance' sprintf('\n') ...
      'ZL = ' num2str(ZL) '\Omega | \alpha = ' num2str(alpha) ' Np/m'])
xlabel('x/\lambda')
ylabel('Admittance (S)');
legend('Conductance', 'Susceptance', 'inductive reactance - Inductor', ...
       'Capacitive reactance -  Capacitor')
grid on

