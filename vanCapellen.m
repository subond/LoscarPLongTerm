%%% Attempt at coding van Capellen C-Iron-Phosphorous model 1996
clear all

global koa vmix k54 CPred k36 kup  CPoxic CPanoxic k1011...
     k58 k59 k25 kt Fcbur
addpath('/home/komar/Dropbox/geocarb/Loscar d44Ca/myode');

kt =1;

% initial massed
M1i = 65e20;
M2i = 14e18;
M3i = 35e14;
M4i = 18e12;
M5i = 2e15;
M6i = 13e20;
M7i = 52e17;
M8i = 52e17;
M9i = 22e17;
M10i = 13e20;
%M11 = 
M12i = 1.2e20;
M13i = 1.529572300117367e+20;
M14i = 38.28e18;

% constants
CPred = 106;
% Cappellen scenarions:
% C:P oxic vs. C:P anoxic 
% without redox dependent P recycling 250 250
% minimum estimate of redox fractionation 250 500
% With redox dependent P cycling and 200 4000
CPoxic = 200; % 250
CPanoxic = 4000; %250 %500
k25 = 2.57*1e-9; %year^-1
k36 = 1.199999999999918e-26; % (mol/year)^-1.5
k54 = 3.930817610063*1e-3;
k58 = 5.56*1e-24;
k59 = 1.71*1e-2;
k1011 = 6.08*1e-10;
kup = 2.884615384615385e-09;
koa = 1.73*1e-5;

vmix = 3.000; %not sure what value they used!!!

F54 = k54*vmix*M5i;
F13 = CPred*F54;
F36 = k36*F13^2.5;



DOA = 1-koa*(vmix*M14i/F13);
x = 0.44+0.56*(DOA);
CPbur = (CPoxic*CPanoxic)/((1-DOA)*CPanoxic+DOA*CPoxic);

%Iron
F1011 = k1011*M10i;
F1113 = x*F1011;
F1112 = (1-x)*F1011;
F1210 = kup*M12i;
F1310 = kup*M13i;


%Phosphorous

F47 = F36/CPbur;
F45 = F54 - F47;

F58 = k58*F45^2.5;
F59 = k59*F1112;
F72=kup*M7i;
F82 = kup*M8i;
F92 = kup*M9i;
F25 = k25*M2i;


% Carbon

F31 = F13-F36;

F61 = kup*M6i;


t0 =0;
tfinal = 250e6;
y0 = [M14i M1i M2i M3i M4i M5i M6i M7i M8i M9i M10i M12i M13i];
tspan = tfinal -t0;

format compact
options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
    'InitialStep',1.,'Maxstep',tspan/1000.);
tic
[T,Y] = myode15s(@CappellenDiff,[t0 tfinal],y0,options);
toc
Fcbur(1) = F36;
Ox = Y(:,1);

figure(1)
plot(T/1e6, Ox/1e18)
ylabel('Oxygen')

figure(2)
plot(T,Fcbur)
ylabel('Org Carb Burial')
