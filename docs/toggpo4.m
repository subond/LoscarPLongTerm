function toggpo4
%--------------------------------------
%
% file: toggpo4.m
%
% OCN310L, R. E. Zeebe
% 
% cd L:\Teach\OCN310L
%
% 3-box ocean model (+1 for atmosphere)
%
% updates: 
%
% 10/24/06 REDTCP
% 09/20/05 revised
% 10/13/04 new file
%
% 
%--------------------------------------

clear all; % clear memory (all variables)
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%			3 BOX: PO4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global variables, known to all functions
global fhd Aoc fH PPH REDTCP TH V PPL PLPr;

format compact

Aoc    = 3.49e14; % 3.49e14 (m2) area ocean
REDTCP =   162.5; % 162.5 Redfield  P: Total C

% BOXES
%        Low  High Deep
%
%         1    2     3
%         L    H     D

H     = [100.;250.;3579.]; % 100 3579 250 (m) height of boxes
lH    = length(H);
fH    = .15;               % high-latitude fraction
A     = Aoc*[(1-fH);fH;1]; % m2
V     = A.*H;              % m3
Voc   = sum(V);            % m3

% Overturning
TH  = 20.e6*3600*24*365;   % (m3/y) 20 Sv conveyor transport 

% mixing term fhd
fhd = 60.e6*3600*24*365;   % (m3/y) 60 Sv (60e6 m3/s)

po4d0 = 2.144;  % (mumol/kg) initial PO4 deep 2.144
CPH   = 1.*A(2);% (mol/y) 1 high-lat C Export (CPH*13, Gt C/y)
                % 1 mol C/m2/y * A = mol/y
                     
% (mol/y) high-lat P Export                  
PPH = CPH/REDTCP; 

rho = 1.025e3;  % kg/m3 (1.025: Toggweiler 1999)
                % ftp://kosmos.agu.org/apend/pa/1999PA900033

% Define initial conditions.

po40 = [0.5;0.5;po4d0]*1e-6; % mol/kg (PO4 at t=0)
po40 = po40*rho;             % mol/m3
Y0   = po40;                 % initial value Y0


% initial inventory (mol/m3*m3)/(m3)*1e3 = mmol/m3
Mpi  = sum(Y0(1:3).*V)/Voc*1e3;	 % PO4 inventory/Voc


% set integration time                      

t0     =    0;  % (y) time (start)
tfinal = 1500;  % (y) time (end)

% solve differential equations.
% use matlab routine ode (Runge-Kutta)
% with function 'fdif' containing the difeq.
%
% set options, 'RelTol' (1e-3 by default) 
%              'AbsTol' (all components 1e-6 by default).
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
   
[tv,Y] = ode23t(@fdif,[t0 tfinal],Y0,options);

lt = length(tv);

% final inventories/average concentrations 
% (mol/m3*m3)/(m3)*1e3 = mmol/m3

Mp = sum(Y(lt,1:3).*V')/Voc*1e3;

% check inventories: difference to initial values
DPINV = [Mp-Mpi];
DPINV

% output LowLat P export (m3/y * mol/m3 = mol/y)
PPL
% output LowLat C export (m3/y * mol/m3 = mol/y)
%PLPr/REDTCP/(0.85*Aoc)
%PLPr/REDTCP/(3600*24*365)
PLPr
PLPr*12/1e15 % Gt C y

% rename results 'po4' and convert to mol/kg
po4(:,1) = Y(:,1)/rho;
po4(:,2) = Y(:,2)/rho;
po4(:,3) = Y(:,3)/rho;

% PO4 final in the three boxes, convert to mumol/kg
po4v  = [po4(lt,1) po4(lt,2) po4(lt,3)]*1e6;

% display result in command window
po4v

% plot results
figure(1)
plot(tv  ,po4  (:,1)*1e6,'g-');
hold on;
plot(tv  ,po4  (:,2)*1e6,'b-');
plot(tv  ,po4  (:,3)*1e6,'r-');
hold off;
grid;
xlabel('Time (y)');
ylabel('PO_4 (\mumol kg^{-1})');


%-------------------------------------
% function dydt = fdif(t,y)
% provides difequations 
%-------------------------------------

function dydt = fdif(t,y)

% global variables, known to all functions
global fhd Aoc fH PPH REDTCP TH V PPL PLPr;

% Low/high latitude areas
AL = Aoc*(1-fH);
AH = Aoc*   fH;

%    PO4 
% L   1 
% H   2 
% D   3 

% rename variables for differential equations

            % mol/m3
pl  = y(1); % PO4 low  lat surf
ph  = y(2); % PO4 high lat surf
pd  = y(3); % PO4 deep

% Volume, m3

V1 = V(1);
V2 = V(2);
V3 = V(3);

% Low Lat Carbon Export
EPL   = TH*pd*REDTCP; % (m3/y * mol/m3 = mol/y)
PLPr  = EPL;          % LowLat C-Export for output

% Low Lat PO4    Export
PPL   = TH*pd;        % (m3/y * mol/m3 = mol/y)
PP    = PPL+PPH;

% differential equations

dydt=[                          % PO4
     ( TH*(pd-pl)-PPL )/V1      % L  
     ( TH*(pl-ph)-PPH       ... % H 
     +fhd*(pd-ph)     )/V2      % H
     ( TH*(ph-pd)+PP        ... % D
     +fhd*(ph-pd)     )/V3      % D 
     ];                                 

return;
                              
