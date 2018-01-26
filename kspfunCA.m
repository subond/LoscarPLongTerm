%-------------------------------------
% function f = kspfun(TC,S,P,Ca,Mg)
%
% returns kspc,kspa
% 
% updates: 
%
% 10/20/07 new
%
%-------------------------------------
function [kspc,kspa] = kspfunCA(TC,S,P,Ca,Mg)


global phflag k1k2flag Cam Mgm;

tk = 273.15;           % [K] (for conversion [deg C] <-> [K])
T = TC + tk;           % TC [C]; T[K]


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
kspc = kspc*exp(lnkpok0(7));
kspa = kspa*exp(lnkpok0(8));

end;
%---------------------- END Pressure effect

if(1)
%-------- Mg/Ca effect on K's (Tyrrell Zeebe 04) -----%

xm   = Mgm/Cam;
xt   = Mg /Ca ;
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

return;



