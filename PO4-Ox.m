
Yf=123; % +-24
A=-112; % +-24
r=32;   % +-19

DO=[150:0]; % dissolved oxygen in muM
Cex = 450e12; %moles of C per year
Crain = Cex*0.1; %lets assume 2% of Cex reaches the sediments
Prain = Crain/148; %148 is redfiled ratio in marine plankton at 1000ppm;

rREG = Yf+A*exp(-DO/r);

Cbur = Crain*0.02; %assume 2% of rain that reaches the floor gets buried
Cben = Crain-Cbur;