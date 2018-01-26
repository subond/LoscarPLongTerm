%-------------------------------------
% function f = dafunPE(dic,alk,TC,S,P)
%
% returns co2,pco2,co3,ph,kh,o2,kspc,kspa
% 
% updates: 
%
% 03/23/06 Effect of Ca/Mg on K's included
% 12/26/05 P-Error in dafunPE corrected
%
%-------------------------------------
function [co2,pco2,co3,ph,kh,o2,kspc,kspa] = ...
    dafunPECA(dic,alk,TC,S,P,Ca,Mg)


global phflag k1k2flag Cam Mgm;

bor = 1.*(416.*(S./35.))*1.e-6;

tk = 273.15;           % [K] (for conversion [deg C] <-> [K])
T = TC + tk;           % TC [C]; T[K]
Cl = S / 1.80655;      % Cl = chlorinity; S = salinity (per mille)
cl3 = Cl.^(1/3);   
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl .* Cl;   % ionic strength
iom0 = 19.924*S/(1000.-1.005*S);
S_T = 0.14/96.062/1.80655*S;   % [mol/kg soln] total sulfate
                               %  Dickson and Goyet (1994) Ch.5 p.11


%-------------- solubility of O2 ---------------------
%  Weiss (1970) DSR, 17, p. 721.
%

A = [-177.7888 255.5907 146.4813 -22.2040];
B = [-0.037362 0.016504 -0.0020564];


lno2  = A(1)+A(2)*100./T+A(3)*log(T/100)+A(4)*(T/100) ...
	+S*(B(1)+B(2)*T/100+B(3)*(T/100).^2);

o2  = exp(lno2)/22.4;


% --------------------- kwater -----------------------------------
%
%       Millero (1995)(in Dickson and Goyet (1994, Chapter 5, p.18))
%       $K_w$ in mol/kg-soln.
%       pH-scale: pH$_{Hansson}$ ('total` scale).
                                                     

tmp1 = -13847.26./T + 148.96502 - 23.6521 .* log(T);
tmp2 = + (118.67./T - 5.977 + 1.0495.*log(T)).*sqrt(S) - 0.01615.*S;

lnkw =  tmp1 + tmp2;

if phflag == 0;
        kw  = exp(lnkw);
end;
if phflag == 1;
        lnkw = lnkw-log(total2free);
        kw  = exp(lnkw);
end;


%---------------------- kh (K Henry) ----------------------------
%
%               CO2(g) <-> CO2(aq.)
%               kh      = [CO2]/ p CO2
%
%   Weiss (1974)   [mol/kg/atm]
%
%                             
%
tmp = 9345.17 ./ T - 60.2409 + 23.3585 * log(T/100.);
nkhwe74 = tmp + S.*(0.023517-0.00023656*T+0.0047036e-4*T.*T);

tmpkh1= 9050.69 ./ T - 58.0931 + 22.2940 * log(T/100.);
nkhwe74l = tmpkh1 + S .* (0.027766-0.00025888*T+0.0050578e-4 * T .* T);

%kh= exp(nkhwe74l);
kh= exp(nkhwe74);


% --------------------- k1 ---------------------------------------
%   first acidity constant:
%   [H^+] [HCO_3^-] / [H_2CO_3] = K_1
%
%   Mehrbach et al (1973) refit by Lueker et al. (2000).
%
%   pH-scale: 'total'. mol/kg-soln

pk1mehr = 3633.86./T - 61.2172 + 9.6777.*log(T) - 0.011555.*S + 0.0001152.*S.*S;

if phflag == 0;
        k1mehr  = 10^(-pk1mehr);
end;
if phflag == 1;
	lnk1mehr = log(10^(-pk1mehr))-log(total2free);
        k1mehr   = exp(lnk1mehr);
end;


% --------------------- k2 ----------------------------------------
%
%   second acidity constant:
%   [H^+] [CO_3^--] / [HCO_3^-] = K_2
%
%   Mehrbach et al. (1973) refit by Lueker et al. (2000).
%
%   pH-scale: 'total'. mol/kg-soln

pk2mehr = 471.78./T + 25.9290 - 3.16967.*log(T) - 0.01781.*S + 0.0001122.*S.*S;

if phflag == 0;
        k2mehr  = 10^(-pk2mehr);
end;
if phflag == 1;
	lnk2mehr = log(10^(-pk2mehr))-log(total2free);
        k2mehr   = exp(lnk2mehr);
end;


%----------- Mehrbach.

k1 = k1mehr;
k2 = k2mehr;


% --------------------- kb  --------------------------------------------
%  Kbor = [H+][B(OH)4-]/[B(OH)3] = kp7 / km7
%
%   (Dickson, 1990 in Dickson and Goyet, 1994, Chapter 5, p. 14)
%   pH-scale: 'total'. mol/kg-soln


tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S.^(3./2.)-0.0996*S.*S);
tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S).*log(T);

lnkb = tmp1 ./ T + tmp2 + tmp3 + 0.053105*sqrt(S).*T;

if phflag == 0;
        kb  = exp(lnkb);
end;
if phflag == 1;
        lnkb = lnkb-log(total2free);
        kb  = exp(lnkb);
end;

% --------------------- Kspc (calcite) ----------------------------
%
% apparent solubility product of calcite
%
%  Kspc = [Ca2+]T [CO32-]T
%
%  where $[]_T$ refers to the equilibrium total 
% (free + complexed) ion concentration.
%
%  Mucci 1983 mol/kg-soln

tmp1 = -171.9065-0.077993.*T+2839.319./T+71.595.*log10(T);
tmp2 = +(-0.77712+0.0028426.*T+178.34./T).*sqrt(S);
tmp3 = -0.07711.*S+0.0041249.*S.^1.5;
log10kspc = tmp1 + tmp2 + tmp3;

kspc = 10.^(log10kspc);

% --------------------- Kspa (aragonite) ----------------------------
%
% apparent solubility product of aragonite
%
%  Kspa = [Ca2+]T [CO32-]T
%
%  where $[]_T$ refers to the equilibrium total 
% (free + complexed) ion concentration.
%
%  Mucci 1983 mol/kg-soln

tmp1 = -171.945-0.077993.*T+2903.293./T+71.595.*log10(T);
tmp2 = +(-0.068393+0.0017276.*T+88.135./T).*sqrt(S);
tmp3 = -0.10018.*S+0.0059415.*S.^1.5;
log10kspa = tmp1 + tmp2 + tmp3;

kspa = 10.^(log10kspa);


%-------- Pressure effect on K's (Millero, 95) -------%
if(P > 0.0)

RGAS = 8.314510;        % J mol-1 deg-1 (perfect Gas)  
R = 83.131;             % mol bar deg-1 (imperfect Gas?)
                        % conversion cm3 -> m3          *1.e-6
                        %            bar -> Pa = N m-2  *1.e+5
                        %                => *1.e-1 or *1/10

% index: k1 1, k2 2, kb 3, kw 4, ks 5, kf 6, kspc 7, kspa 8
%        k1p 9, k2p 10, k3p 11

%----- note: there is an error in Table 9 of Millero, 1995.
%----- The coefficients -b0 and b1
%----- have to be multiplied by 1.e-3!

%----- there are some more errors! 
%----- the signs (+,-) of coefficients in Millero 95 do not
%----- agree with Millero 79


a0 = -[25.5   15.82  29.48  25.60  18.03    9.78  48.76   46. ...
	14.51 23.12 26.57];
a1 =  [0.1271 -0.0219 0.1622 0.2324 0.0466 -0.0090 0.5304  0.5304 ...
	0.1211 0.1758 0.2020];
a2 =  [0.0     0.0    2.608 -3.6246 0.316  -0.942  0.0     0.0 ...
	-0.321 -2.647 -3.042]*1.e-3;
b0 = -[3.08   -1.13   2.84   5.13   4.53    3.91  11.76   11.76 ...
	2.67 5.15 4.08]*1.e-3;
b1 =  [0.0877 -0.1475 0.0    0.0794 0.09    0.054  0.3692  0.3692 ...
	0.0427 0.09 0.0714]*1.e-3;
b2 =  [0.0     0.0    0.0    0.0    0.0     0.0    0.0     0.0 ...
	0.0 0.0 0.0];

for ipc=1:length(a0);
  deltav(ipc)  =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
  deltak(ipc)  = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
  lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
end;

k1 = k1*exp(lnkpok0(1));
k2 = k2*exp(lnkpok0(2));
kb = kb*exp(lnkpok0(3));
kw = kw*exp(lnkpok0(4));
kspc = kspc*exp(lnkpok0(7));
kspa = kspa*exp(lnkpok0(8));

end;
%---------------------- END Pressure effect

if(1)
%-------- Mg/Ca effect on K's (Tyrrell Zeebe 04) -----%

xm   = Mgm/Cam;
xt   = Mg /Ca;
% sensitivity parameters for Mg effect on K*
sk1mg = 155.e-3;
sk2mg = 442.e-3;
% sensitivity parameters for Ca effect on K*
sk1ca = 33.73e-3;
sk2ca = 38.85e-3;

% add Mg correction K* (Ben-Yaakov & Goldhaber, 1973)
delk1   = sk1mg*k1*(Mg/Mgm-1.);
delk2   = sk2mg*k2*(Mg/Mgm-1.);
% add Ca correction K* (Ben-Yaakov & Goldhaber, 1973)
delk1ca = sk1ca*k1*(Ca/Cam-1.);
delk2ca = sk2ca*k2*(Ca/Cam-1.);
%

k1  = k1+delk1+delk1ca;
k2  = k2+delk2+delk2ca;

% add Mg correction Ksp* (Mucci & Morse, 1984)
% about 35% decrease for Mg/Ca from 5 to 1.

if(1) % new 03/23/06
alpha = 0.0833;    % slope of Ksp* vs. Mg/Ca (MM84)
kspc  = kspc*(1 - alpha*(xm-xt));
else  % old
alpha = 0.3655e-7; % slope of Ksp* vs. Mg/Ca (MM84)
kspc  = kspc    - alpha*(xm-xt);
end;


%-----------------------------------------------------%
end;


p5  = -1.;        
p4  = -alk-kb-k1;
p3  = dic*k1-alk*(kb+k1)+kb*bor+kw-kb*k1-k1*k2;
tmp = dic*(kb*k1+2.*k1*k2)-alk*(kb*k1+k1*k2)+kb*bor*k1;
p2  = tmp+(kw*kb+kw*k1-kb*k1*k2);
tmp = 2.*dic*kb*k1*k2-alk*kb*k1*k2+kb*bor*k1*k2;
p1  = tmp+(+kw*kb*k1+kw*k1*k2);
p0  = kw*kb*k1*k2;
p   = [p5 p4 p3 p2 p1 p0];
r   = roots(p);
h   = max(real(r));
co2 = dic/(1.+k1/h+k1*k2/h/h);
pco2 = co2*1e6/kh;
co3 = dic/(1+h/k2+h*h/k1/k2);
ph = -log10(h);

return;

