%%% Floegel and Wallmann
Yf=123; % +-24
A=-112; % +-24
r=32;   % +-19

DO=[170:-1:0]; % dissolved oxygen in muM
Cex = 450e12; %moles of C per year
Crain = Cex*0.2; %lets assume 20% of Cex reaches the sediments
Prain = Crain/148; %148 is redfiled ratio in marine plankton at 1000ppm;

rREG = Yf+A*exp(-DO./r);

Cbur = Crain*0.08; %assume 6% of rain that reaches the floor gets buried
Cben = Crain-Cbur;

Pben = min(Cben./rREG*106/148, Prain);
Pbur = Prain-Pben;

%%% Van Cappellen, Slomp, Tsandev

popb0 = 0.25*Pbur(1);
capb0 = 0.5*Pbur(1);
fepb0 = 0.25*Pbur(1);
pbtot0 = popb0+capb0+fepb0;


fepb = fepb0*DO/DO(1);
popb = popb0*(0.25+0.75*DO/DO(1));
k=capb0/(Prain-popb0);
capb = k*(Prain-popb).*(0.1+0.9*DO/DO(1));

pbtot = popb+capb+fepb;
figure(1)
% oxygen vs. Pbur
plot(DO, Pbur)
xlabel('[O_2] \mu M')
ylabel('P burial')
figure(2)
plot(DO, Pben)
xlabel('[O_2] \mu M')
ylabel('P benthic flux')
figure(3)
% oxygen vs. Pbur
plot(DO, pbtot,'k',DO,popb,'g',DO,fepb,'r',DO,capb,'b' )
xlabel('[O_2] \mu M')
ylabel('C burial')