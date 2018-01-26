%-------------------------------------
% function yp = LoscarDif9(t,y)
% provides difequations 
%  
% file: LoscarDif9.m
%
% LOSCAR Model: Long-term Ocean-atmosphere-Sediment 
%               CArbon cycle Reservoir Model
%
% ocean box model (+1 for atmosphere)
%
%
% updates: 
%
% 07/06/08 A few comments added
%
% 02/06/08 ffflag results are slightly different 
%          for MM-O2 from those saved Oct 2007.
%          Set KMMOX = 0.0.
%
% 10/27/07 Michaelis-Menton for Oxygen included.
% 
% 10/20/07 Indices of initial f0 (sediments, load)
%          corrected in Y0 when oxygen is included
%          (caused f>1 during erosion)
% 
% 10/20/07 TDflag: derivs with dYflag = 1 now called with
%          TCvt because TCv can't be updated (odeXX not called). 
%          Shouldn't be an issue for TDflag = 0 as
%          TCv IS updated.
%
% 05/13/07 constraints on CBl and d13CBl input
%          ccdA(1)-min(ccdA)
%          kk=150; min(d13c(kk,1))-d13c(1,1)
%
% 01/29/07 Tethys mv=3.5, rrain=
% 01/28/07 Oxygen included. Lots of changes:
%          hs, TH, TT, rrain 
% 01/22/07 [CO3=] grad: co3tv(jt,9)/co3tv(jt,7)
%
% 12/00/06 sedrate/phi (WRONG!) in old sed model
%          new is OK!
% 07/13/06 Temp for co3sat corrected
% 07/05/06 rcak for Ca/Cam corrected (Tethys)
% 07/01/06 Version B runs. Nearly same results.
% 06/26/06 New sediment model (Version B)
% 05/27/06 New erosion 
% 05/03/06 Pac CCD shallowed (x in THmfun)
% 03/23/06 Effect of Ca/Mg on K's included
% 03/12/06 Adjust dissolution to Ca. nc 2.4
% 03/11/06 shelf/deep rain           C13 missing. done
% 03/08/06 mix/biopump changed
% 03/04/06 Millero sat, Ca = 20, 1000 ppmv
% 03/01/06 H-Lat mix & 13Cin changed for basin-d13C
% 01/18/06 switch TH SO-NP-SO
% 12/29/05 Tethys ocean+sed complete
% 12/26/05 P-Error in dafunPE corrected
% 12/24/05 Tethys
% 12/10/05 Tuning good with water column diss
% 12/08/05 Tuning
% 11/26/05 10-box: rhos, phi = phi(fc)
% 11/22/05 03-box: DEQ in df/dt
% 11/05/05 porosity, phi = phi(fc)
% 10/27/05 error fixed (rsed could become < 0)
% 10/25/05 C13 10-box ocean + sediment
% 10/23/05 C13 10-box ocean
% 10/21/05 C13  3-box sediment
% 10/18/05 resumed
% 05/26/05 10-box (c,a,p) + sed complete
% 05/11/05 new file
%
%-------------------------------------
function yp = LoscarDif9(t,y)

global myflag Fem20 kt tst Dtst yst dYst stflag ...
       Cam Ca Mgm Mg y2s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%			10 BOX + sediment + Tethys
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if myflag == 01;

%========================================
% global variables:
% known to Loscar9 and LoscarDif9
%========================================

global Aoc rho Nb V A TH0 TH TS mv mhd kasv kkv TCv Sv Pv gp tA tI ...
    phflag k1k2flag EPH fEPL rrain ep ept ec ect REDPC REDNC REDO2C eI ...
    Ns asvA asvI asvP zv m2kg rhos frrf phic hs FiN ...
    phivtA phivtI phivtP phi0 phi1 gam ...
    phiiA phiiI phiiP phiiT dsv ...
    Fpr frain rsedv rburvtA rburvtI rburvtP dissvtA ...
    dissvtI dissvtP FprtA FprtI FprtP it co3satv ...
    Kd KS cst dsflag nc dYflag nu fsed ...
    Fpr13tA Fpr13tI Fpr13tP Fint Fin13t FSit FSi13t ...
    fc0A fc0I fc0P f13c0A f13c0I f13c0P ...
    dissv13tA dissv13tI dissv13tP ... 
    co3s0 as zs0 klid nli Rst epsp Rin FiN13 ...
    BlFlag CBl RBl kb FVC FVC13 Rvc pCSi nSi nCC Fkg Fkg13 ...
    ftys nOC TT asvT kliT fc0T f13c0T rsedvtT dissvtT dissv13tT ...
    FprtT Fpr13tT phivtT fcon swcon ...
    TCv0 TCvt ntL ntH fsh fdpv fshT nshT ...
    ffflag tem em RlsCtv ...
    fdox KMMOX vask DTS DTS2 ts3 DTS3 DTS4 ...
    k1st tfinal TDflag omegCSvt omegASvt ...
    THt mv0 mhd0 oxA CAvflag FSichck Finchck kspCHCK ...
    CHECK1 CHECK2 CHECK3 CHKFin CHKFSi CHCKbioL CHCKbioI CHCKbioD...
    CHKcarB1 CHKcarB2 CHKcarB3 CHKcar1 CHKcar2 CHKcar3 CHKcar4 CAv;    
    


ERR = 1; % report fc < 0 error on/off    
    
% shorthand for box volume 
V1 = V(1); V2 = V(2); V3 = V(3); V4 = V(4); V5  = V(05);
V6 = V(6); V7 = V(7); V8 = V(8); V9 = V(9); V10 = V(10);

VsvA = asvA*A(1)*hs;
VsvI = asvI*A(2)*hs;
VsvP = asvP*A(3)*hs;

if(ftys)
V11  = V(11); V12 = V(12); V13 = V(13);
VsvT = asvT*A(11)*hs;
end;

if(fsed)
n1 = nli(1); n2 = nli(2);
end;

%=========== CO2 emissions 
if(ffflag)
 Fem    = PEmisfun(t)*1e15/1.; % 16.6667
 Fem    = Fem/12; % (g C/y)/12-> mol/y
 d13Cem = -25.0;  % d13C emissions -25.
 Rem    = Rst*(d13Cem/1e3+1);
 Fem13  = Rem*Fem;
end;

%=========== Blast: Circ, Temp, Shelf Prod.
tcon = 0.0e5;
t2bl = 72.e3; % 72 80
if(swcon)
if(t > tcon & t < tcon+t2bl)
    fcon = 2;             % 2 = North Pac 
      TH = 0.25*TH0;      % 0.3(1,1) 0.40/.25(2,3) n0.32
      TS = 0.30*TH0;      % additional SO 0.6 0.4 0.5 n0.3
%      TH = flin(0.,3.e3,0.0,.25*TH0,t);
%      TS = flin(0.,3.e3,TH0,.30*TH0,t);
%%     TH = flin(2e4,t2bl,0.5*TH0,0.3*TH0,t);
%%     TS = flin(2e4,t2bl,0.2*TH0,0.4*TH0,t);
%     TCv = TCv0 + 4.0;    % +3.0 +4.0
%    fsh  = 7.30;          % 6.5(1,1) 7.50/6.30(2,3) n7.3
%   nshT  = 0.30;          % 0.5(1,1) 0.30/0.55(2,3).25
else
    fcon = 3;              % 3 = all SO
      TH = 1.0*TH0;         
      TS = 0.0*TH0;        % additional SO
%     TCv = flin(t2bl,4*t2bl,TCv0+4.0,TCv0,t);
%    fsh  = 4.50;           % 4.5(1,1) 4.50/4.0(2,3)
%   nshT  = 0.40;           % 0.6(1,1) 0.35/0.6(2,3).3
end;
end;



%=========== Blast: Temp
if(BlFlag == 2)
if(t > tcon & t < tcon+t2bl)
     TCv = TCv0 + 4.0;     % +3.0 +4.0
%if(t >= tcon & t < tcon+t2bl)
%     TCv = flin(0.,3.e3,TCv0,TCv0+4.0,t);
    fsh  = 7.30;          % 6.5(1,1) 7.50/6.30(2,3) n7.3
   nshT  = 0.30;          % 0.5(1,1) 0.30/0.55(2,3).25
else
     TCv = flin(t2bl,4*t2bl,TCv0+4.0,TCv0,t);
    fsh  = 4.50;           % 4.5(1,1) 4.50/4.0(2,3)
   nshT  = 0.40;           % 0.6(1,1) 0.35/0.6(2,3).3
end;
end;
    

% Tracer -> DIC ALK PO4 O2 DIC-13 Catm Catm-13 sediments ...
%
% Boxes:
% 
% Low, Interm, Deep, High
% Atlantic, Pacific, Indic
%
% AL   1      
% IL   2    
% PL   3     
% AI   4    
% II   5     
% PI   6     
% AD   7   
% ID   8   
% PD   9   
% H   10      

%      Nb 2Nb 3Nb 4Nb +1  +2 +2+Ns +2+2Ns +2+3Ns
%      10  20  30  40 41  42    58     74     90
%Y  =  [c   a   p  cc  C  CC  mcvA   mcvI   mcvP]';

%      Nb 2Nb 3Nb 4Nb +1  +2 +2+Ns  +2+2Ns   +2+3Ns
%      10  20  30  40 41  42    55      68       81
%Y0 = [c0  a0  p0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
%                            f13c0A  f13c0I   f13c0P]';
%                           +2+4Ns  +2+5Ns   +2+6Ns 
%                               94     107      120

if(fdox)
c    = y(     1:  Nb     ); % TCO2    (mol/m3)
a    = y(  Nb+1:2*Nb     ); % TA      (mol/m3)
p    = y(2*Nb+1:3*Nb     ); % PO4     (mol/m3)
CA   = y(3*Nb+1:4*Nb   );
dox  = y(4*Nb+1:5*Nb     ); % O2      (mol/m3)
cc   = y(5*Nb+1:6*Nb     ); % T13CO2  (mol/m3)
C    = y(6*Nb+1          ); % Catm    (mol/m2)
CC   = y(6*Nb+2          ); % 13Catm  (mol/m2)
if(fsed)
fcvA  = y(6*Nb+3     :6*Nb+2+1*Ns); % %calc 
fcvI  = y(6*Nb+3+1*Ns:6*Nb+2+2*Ns); % %calc
fcvP  = y(6*Nb+3+2*Ns:6*Nb+2+3*Ns); % %calc
%====== Carbon-13
f13cvA= y(6*Nb+3+3*Ns:6*Nb+2+4*Ns); % %calc 
f13cvI= y(6*Nb+3+4*Ns:6*Nb+2+5*Ns); % %calc 
f13cvP= y(6*Nb+3+5*Ns:6*Nb+2+6*Ns); % %calc
if(ftys)
fcvT  = y(6*Nb+3+6*Ns:6*Nb+2+7*Ns); % %calc
f13cvT= y(6*Nb+3+7*Ns:6*Nb+2+8*Ns); % %calc
end;
end;
else % fdox
c    = y(     1:  Nb     ); % TCO2    (mol/m3)
a    = y(  Nb+1:2*Nb     ); % TA      (mol/m3)
p    = y(2*Nb+1:3*Nb     ); % PO4     (mol/m3)
CA   = y(3*Nb+1:4*Nb     );
cc   = y(4*Nb+1:5*Nb     ); % T13CO2  (mol/m3)
C    = y(5*Nb+1          ); % Catm    (mol/m2)
CC   = y(5*Nb+2          ); % 13Catm  (mol/m2)
if(fsed)
fcvA  = y(5*Nb+3     :5*Nb+2+1*Ns); % %calc 
fcvI  = y(5*Nb+3+1*Ns:5*Nb+2+2*Ns); % %calc
fcvP  = y(5*Nb+3+2*Ns:5*Nb+2+3*Ns); % %calc
%====== Carbon-13
f13cvA= y(5*Nb+3+3*Ns:5*Nb+2+4*Ns); % %calc 
f13cvI= y(5*Nb+3+4*Ns:5*Nb+2+5*Ns); % %calc 
f13cvP= y(5*Nb+3+5*Ns:5*Nb+2+6*Ns); % %calc
if(ftys)
fcvT  = y(5*Nb+3+6*Ns:5*Nb+2+7*Ns); % %calc
f13cvT= y(5*Nb+3+7*Ns:5*Nb+2+8*Ns); % %calc
end;
end;
end; % fdox


if(fsed)
fcVV  = [fcvA fcvI fcvP f13cvA f13cvI f13cvP];
if(ERR & ~isempty(find(fcVV < 0.)) )
 disp('+++ WARNING: fc is negative +++'); 
if(0)
 myeps = 1.e-6;
 kzA   = find(fcvA < 0.);
 kzI   = find(fcvI < 0.);
 kzP   = find(fcvP < 0.); 
 fcvA(kzA)  = myeps;
 fcvI(kzI)  = myeps;
 fcvP(kzP)  = myeps; 
 k13zA   = find(f13cvA < 0.);
 k13zI   = find(f13cvI < 0.);
 k13zP   = find(f13cvP < 0.); 
 f13cvA(kzA)  = myeps;
 f13cvI(kzI)  = myeps;
 f13cvP(kzP)  = myeps; 
end;
end;  
end;

pco2a  =   C/(2.2e15/12/Aoc);
pcco2a =  CC/(2.2e15/12/Aoc);

% Deep-sea temperature sensitivity (Archer, 2005)

if(TDflag == 1 & dYflag == 0 & ftys == 0)
 DTeq = 3.0*log(pco2a/pCSi)/log(2); % T sens 1.5,3,4.5
 
 if(kt > k1st) % first call with new kt
     f1st = 1;
 else 
     f1st = 0;    
 end
 if(f1st == 1) % update TCv
      kvi   = [4 5 6];
      kvd   = [7 8 9];
     tlag   = [10. 200. 1000.];
    %tlag   = [1     1      1]*1000.;
  TCv(kkv)  = TCv(kkv)+...
      ((TCv0(kkv)+DTeq)-TCv(kkv))*Dtst/tlag(1); % 10   surf
  TCv(kvi)  = TCv(kvi)+...
      ((TCv0(kvi)+DTeq)-TCv(kvi))*Dtst/tlag(2); % 200  intmd
  TCv(kvd)  = TCv(kvd)+...
      ((TCv0(kvd)+DTeq)-TCv(kvd))*Dtst/tlag(3); % 1000 deep
  f1st = 0;
  k1st = kt;
 end
  TCvt(kt+1,:) = TCv; % store TCv in TCvt
end

if(0)
% Temp-CO2 feedback
TCv     = TCv0    *(pco2a/pCSi).^ntL;
TCv(10) = TCv0(10)*(pco2a/pCSi).^ntH;
if(ftys)
TCv(13) = TCv0(13)*(pco2a/pCSi).^ntH;
end;
end;

%====== TDflag: use TX to get fluxes !!!
if(TDflag == 1 & dYflag == 1)    
TX    = TCvt(it,:);
else
TX    = TCv;
end

% Conveyor slowdown
if(0)    
ftm = 1.5*(TX(1)-TCv0(1))/TCv0(1); % 2.0 2.2 6.2
%ftm = (t-1700)/700;
%ftm = 0.9;
TH  =  TH0*(1.-ftm);
mv  =  mv0*(1.-ftm);
mhd = mhd0*(1.-ftm/1.); % /1. /2.
end

% CO2 system & O2 of boxes (surface: 1=LA, 2=LI, 3=LP, 10=H)
% requires mol/kg
if(CAvflag == 2)
for k=1:Nb
[co2(k),pco2(k),co3(k),ph(k),kh(k),o2(k)] = ...
   dafunPECA(c(k)/rho,a(k)/rho,TX(k),Sv(k),Pv(k),CA(k)*1e-3,Mg);
end;
else
  for k=1:Nb
[co2(k),pco2(k),co3(k),ph(k),kh(k),o2(k)] = ...
   dafunPECA(c(k)/rho,a(k)/rho,TX(k),Sv(k),Pv(k),Ca,Mg);
  end;
end;


if(1)
%====== Carbon-13
Rb  = cc./c;

Tk  = TX + 273.15;
edb = -9866./Tk + 24.12;
edg = - 373./Tk +  0.19;
adb = edb/1e3+1;
adg = edg/1e3+1;
au  = 0.9995;

pcco2(kkv) = adb(kkv).*Rb(kkv)'.*pco2(kkv); % (adg^-1 dropped)
end;

%============== Biological Pump ========%
%
% Low Lat Export Corg
EPLv    = fEPL*mv(1:3)'.*p(4:6)/REDPC; % (m3/y*mol/m3 = mol/y)
if(ftys)
EPLv(4) = fEPL*mv(004)'.*p(012)/REDPC; % (m3/y*mol/m3 = mol/y)
end;


fcb = 1;
% biological calc response to low CO3L
if(0)
cb = 224.e-6;

%if(co3(1) < cb)
%    fcb = (co3(1)/cb)^1;
%else
%    fcb = 1.;
%end;  

ab  =  0.0935; 
bb  = -8.4805;
fcb =     (ab*co3(1)*1.e6+bb);
fcb = fcb/(ab    *cb*1.e6+bb);
%if(fcb < 1.e-5 | t > 2500)
if(fcb < 1.e-4)
fcb = 1.e-4;
end;

end;

 
EALv   = 2*EPLv/rrain*fcb;    % ALK export
PPLv   = EPLv*REDPC;      % PO4
ENLv   = EPLv*REDNC;      % NO3
EOLv   = EPLv*REDO2C;     % O2
ECALv  = EALv/2;          % Ca THIS ADDED!!!!! Ca transport
%total carbon: Corg+CaCO3
ECLv   = EPLv+EALv/2;

% High Lat Export
EAH   = 2*EPH/rrain;    % ALK export
EAH   = 0.000000000000; % !!!!!!!!!!!!!!!!!!!!!!!
PPH   = EPH*REDPC;      % PO4
ENH   = EPH*REDNC;      % NO3
EOH   = EPH*REDO2C;     % O2
ECAH  = EAH/2;          % Ca THIS ADDED!!!!! Ca transport
%total carbon: Corg+CaCO3
ECH   = EPH+EAH/2;

% total Corg export (g C/y)
ep = (sum(EPLv)+EPH)*12;
ec = (sum(EALv)+EAH)*12/2;


% fraction EPL, remineralized in I boxes
oI = 1-eI;

kL = [1:1:3];    
if(ftys)
kL = [kL 11]; 
end;
%======  Carbon-13
alp    = epsp/1e3+1;
% Corg
EPLvCC = alp*Rb(kL).*EPLv;
EPHCC  = alp*Rb(10) *EPH;
EALvCC =     Rb(kL).*EALv;
EAHCC  =     Rb(10) *EAH;
% Ctotal: Corg+CaCO3
ECLvCC = EPLvCC+EALvCC/2;
ECHCC  = EPHCC +EAHCC /2;

% calcite export
EAA   = EALv(1)+EAH*gp(7); % +EAH (no H seds, see below)
EAI   = EALv(2)+EAH*gp(8);
EAP   = EALv(3)+EAH*gp(9);


% Fpr is rain to sediments, not export
% Fpr = export - water column diss


FprA  = (1-nu)*EAA/2/A(1); % mol C/m2/y
FprI  = (1-nu)*EAI/2/A(2); % mol C/m2/y
FprP  = (1-nu)*EAP/2/A(3); % mol C/m2/y
%====== Carbon-13
EAACC   = EALvCC(1)+EAHCC*gp(7); % +EAH (no H seds, above)
EAICC   = EALvCC(2)+EAHCC*gp(8);
EAPCC   = EALvCC(3)+EAHCC*gp(9);
Fpr13A  = (1-nu)*EAACC/2/A(1); % mol C/m2/y
Fpr13I  = (1-nu)*EAICC/2/A(2); % mol C/m2/y
Fpr13P  = (1-nu)*EAPCC/2/A(3); % mol C/m2/y


if(fsed)
%====== shelf/deep rain
shlf  = 1;
if(shlf)
% clay rain to sediments AIP (Atl, Ind, Pac)   
%frem  = [1.00 1.00 1.00];       % 
frem  = [1.00 1.00 0.50];        % P: 0.1 0.05 !#!#!#!#!#!
frrfA = frrf*ones(1,Ns)*frem(1); %  kg  /m2/y
frrfI = frrf*ones(1,Ns)*frem(2); %  kg  /m2/y
frrfP = frrf*ones(1,Ns)*frem(3); %  kg  /m2/y

% carbonate rain to sediments AIP    
FprvA = FprA*ones(1,Ns); % mol C/m2/y
FprvI = FprI*ones(1,Ns); % mol C/m2/y
FprvP = FprP*ones(1,Ns); % mol C/m2/y

Fpr13vA = Fpr13A*ones(1,Ns); 
Fpr13vI = Fpr13I*ones(1,Ns);
Fpr13vP = Fpr13P*ones(1,Ns); 

% total C rain AIP
FAA    = FprA*A(1);
FII    = FprI*A(2);
FPP    = FprP*A(3);
jj     = 2;
fdpv(1)= ( FAA ...
   -fsh*sum(FprvA(1:jj   ).*asvA(1:jj   )*A(1)) ) ...
  ./    sum(FprvA(jj+1:Ns).*asvA(jj+1:Ns)*A(1))  ;
fdpv(2)= ( FII ...
   -fsh*sum(FprvI(1:jj   ).*asvI(1:jj   )*A(2)) ) ...
  ./    sum(FprvI(jj+1:Ns).*asvI(jj+1:Ns)*A(2))  ;
fdpv(3)= ( FPP ...
   -fsh*sum(FprvP(1:jj   ).*asvP(1:jj   )*A(3)) ) ...
  ./    sum(FprvP(jj+1:Ns).*asvP(jj+1:Ns)*A(3))  ;

fshvA = [fsh*ones(1,jj) fdpv(1)*ones(1,Ns-jj)];
fshvI = [fsh*ones(1,jj) fdpv(2)*ones(1,Ns-jj)];
fshvP = [fsh*ones(1,jj) fdpv(3)*ones(1,Ns-jj)];

frrfA = fshvA.*frrfA;
frrfI = fshvI.*frrfI;
frrfP = fshvP.*frrfP;

FprvA = fshvA.*FprvA;
FprvI = fshvI.*FprvI;
FprvP = fshvP.*FprvP;


if(0)
%FprvA(4:13) = 0.; 
%frrfA(4:13) = 0.;
FprvA = 0.*FprvA; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FprvI = 0.*FprvI; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FprvP = 0.*FprvP; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
frrfA = 0.*frrfA;
frrfI = 0.*frrfI;
frrfP = 0.*frrfP;
end;

Fpr13vA = fshvA.*Fpr13vA;
Fpr13vI = fshvI.*Fpr13vI;
Fpr13vP = fshvP.*Fpr13vP;

% test: FprvA; sum(FprvA.*asvA*A(1)); FAA;
end;

if(ftys)
EAT    = EALv  (4);
EATCC  = EALvCC(4);
FprT   = (1-nu)*EAT  /2/A(11); % mol C/m2/y
Fpr13T = (1-nu)*EATCC/2/A(11); % mol C/m2/y
FprvT  = FprT*ones(1,Ns);      % mol C/m2/y
if(shlf)
fremT  = 6.0;                  % 4. !#!#!#!#!#! clay
frrfT  = frrf*ones(1,Ns)*fremT;%  kg  /m2/y
FprvT  = FprT*ones(1,Ns);      % mol C/m2/y
Fpr13vT= Fpr13T*ones(1,Ns);    
fshT   = (fsh)^nshT;           % fshT < fsh !!! 0.6
FTT    = FprT*A(11);
fdpv(4)= ( FTT ...
  -fshT*sum(FprvT(1:jj   ).*asvT(1:jj   )*A(11)) ) ...
  ./    sum(FprvT(jj+1:Ns).*asvT(jj+1:Ns)*A(11))  ;

fshvT  = [fsh*ones(1,jj) fdpv(4)*ones(1,Ns-jj)];
frrfT  = fshvT.*frrfT;
FprvT  = fshvT.*FprvT;
Fpr13vT= fshvT.*Fpr13vT;
end;
end;

if( ~isempty(find(fdpv < 0.))  )
disp('+++ ERROR: fdpv is negative +++');    
fdpv    
disp('+++ Shelf rain too large    +++');    
return;
end;  

end; % fsed
% CA(k);
% CAv(kt)=CA(k)*1e-3;
if(fsed) % fsed: sediments
%%%%=================== Sediment Stuff =============%%%%

% update saturation state

if(1)
Tdv   = TX(klid);
Sdv   = Sv(klid);
zsatv = dsv;
%==== saturation at sediment levels
if(CAvflag == 2)
for k=1:Ns
[kspcSed(k),x] = ...
     kspfunCA(Tdv(k),Sdv(k),zsatv(k)/10.,CA(k)*1e-3,Mg);
 %kspcSed(k);
 %kspCHCK(kt)=kspcSed(k);
% co3satv = kspcSed(k)/(CA(k)*1e-3);
end;
co3satv = kspcSed./(CA(k)*1e-3);
else
    for k=1:Ns
[kspcSed(k),x] = ...
     kspfunCA(Tdv(k),Sdv(k),zsatv(k)/10.,Ca,Mg);
 kspcSed(k);
 kspCHCK(kt)=kspcSed(k);
 co3satv = kspcSed/Ca;
    end;
end;
 
%==== saturation of surface ocean boxes, needed???

if(0)
if(CAvflag == 2)
    for k=kkv
[kspcS(k),kspaS(k)] = ...
     kspfunCA(TX(k),Sv(k),Pv(k),CA(k)*1e-3,Mg);
end;
omegCSv = co3(kkv)*(CA(kkv)*1e-3)/kspcS(kkv);
omegASv = co3(kkv)*(CA(kkv)*1e-3)/kspaS(kkv);
    
end    
else
    for k=kkv
[kspcS(k),kspaS(k)] = ...
     kspfunCA(TX(k),Sv(k),Pv(k),Ca,Mg);
    end;
omegCSv = co3(kkv)*Ca./kspcS(kkv);
omegASv = co3(kkv)*Ca./kspaS(kkv);
    
end 
end;

% set omega and find supersat indices
omvA = co3(klid  )./co3satv; % Atl
omvI = co3(klid+1)./co3satv; % Ind
omvP = co3(klid+2)./co3satv; % Pac
kdA  = find(omvA >= 1.);
kdI  = find(omvI >= 1.);
kdP  = find(omvP >= 1.);
   

%====== Carbon-13
RsA  = f13cvA./fcvA;
RsI  = f13cvI./fcvI;
RsP  = f13cvP./fcvP;

% Ca, Mg/Ca correction for Sigman dissolution
alpha = 0.0833;
xm    = Mgm/Cam;
if(CAvflag == 2)
    for k=1:Ns
xtv(k)    = Mg /(CA(k)*1e-3) ;
rcak(k)  = ((CA(k)*1e-3)/Cam)*(1/(1-alpha*(xm-xtv(k))));
    end
else
    xt    = Mg /Ca ;
    rcak  = (Ca./Cam)*(1/(1-alpha*(xm-xt)));
end
    
% porosities as function of fc    
FF   = (phi1-phi0)/(1-phi1);
phiA = (phi0+FF*fcvA)./(1+FF*fcvA); 
phiI = (phi0+FF*fcvI)./(1+FF*fcvI); 
phiP = (phi0+FF*fcvP)./(1+FF*fcvP); 

% sed rate, m/y (kg/m2/y / kg*m3 = m/y)
rscvA = FprvA*m2kg/rhos/(1-phi1); 
rscvI = FprvI*m2kg/rhos/(1-phi1); 
rscvP = FprvP*m2kg/rhos/(1-phi1); 
rsrvA = frrfA     /rhos/(1-phi0); 
rsrvI = frrfI     /rhos/(1-phi0); 
rsrvP = frrfP     /rhos/(1-phi0); 
rsvA  = rscvA+rsrvA;
rsvI  = rscvI+rsrvI;
rsvP  = rscvP+rsrvP;

%====== Carbon-13
r13scvA = Fpr13vA*m2kg/rhos/(1-phi1); 
r13scvI = Fpr13vI*m2kg/rhos/(1-phi1); 
r13scvP = Fpr13vP*m2kg/rhos/(1-phi1); 

% dissolution
if(CAvflag == 1)
dKA   = KS*((co3satv-co3(klid  ))*rcak).^nc; % mol/m2/y
dKI   = KS*((co3satv-co3(klid+1))*rcak).^nc; % mol/m2/y
dKP   = KS*((co3satv-co3(klid+2))*rcak).^nc; % mol/m2/y
elseif(CAvflag == 2)
dKA   = KS*((co3satv-co3(klid  )).*rcak(klid)).^nc; % mol/m2/y
dKI   = KS*((co3satv-co3(klid+1)).*rcak(klid+1)).^nc; % mol/m2/y
dKP   = KS*((co3satv-co3(klid+2)).*rcak(klid+2)).^nc; % mol/m2/y
end
% fc^0.5
dissA = fcvA.^0.5.*dKA'; % mol/m2/y
dissI = fcvI.^0.5.*dKI'; % mol/m2/y
dissP = fcvP.^0.5.*dKP'; % mol/m2/y


% numerics: square drop in fc, as fc -> 0
if(0)
Bf = 1.e3;    
Bfs= 0.01; % 0.01
nf = 2.0; 
jA = find(fcvA < Bfs); 
jI = find(fcvI < Bfs);
jP = find(fcvP < Bfs);
if(ftys)
jT = find(fcvT < Bfs);
end;    
else %      linear drop in fc, as fc -> 0
Bf = sqrt(1./0.2); % 1./0.2
nf = 1.0; 
jA = find(fcvA < 1./Bf^2);
jI = find(fcvI < 1./Bf^2);
jP = find(fcvP < 1./Bf^2);
if(ftys)
jT = find(fcvT < 1./Bf^2);
end;    
end;

dissA(jA) = fcvA(jA).^nf.*dKA(jA)'*Bf; % mol/m2/y
dissI(jI) = fcvI(jI).^nf.*dKI(jI)'*Bf; % mol/m2/y
dissP(jP) = fcvP(jP).^nf.*dKP(jP)'*Bf; % mol/m2/y

% diss = 0, for omega > 1
dissA(kdA) = 0.;
dissI(kdI) = 0.;
dissP(kdP) = 0.;

% diss rate, m/y [(mol/m2/y)*kg/mol / kg*m3 = m/y]
% pure calcite/(1-phi1) = Delta h1
rdvA  = dissA'*m2kg/rhos/(1-phi1);     % m/y
rdvI  = dissI'*m2kg/rhos/(1-phi1);     % m/y
rdvP  = dissP'*m2kg/rhos/(1-phi1);     % m/y

%====== Carbon-13
r13dvA = RsA.*rdvA';
r13dvI = RsI.*rdvI';
r13dvP = RsP.*rdvP';

% burial rate, m/y
wvA   = rsvA-rdvA; 
wvI   = rsvI-rdvI; 
wvP   = rsvP-rdvP; 

% find erosion indices
lA = find(wvA < 0.);
lI = find(wvI < 0.);
lP = find(wvP < 0.);

% calcite burial (w>0)
wcvA  = fcvA.*wvA'.*(1-phiA)/(1-phi1);
wcvI  = fcvI.*wvI'.*(1-phiI)/(1-phi1);
wcvP  = fcvP.*wvP'.*(1-phiP)/(1-phi1);
fbA   = 1*ones(size(rdvA));
fbI   = 1*ones(size(rdvI));
fbP   = 1*ones(size(rdvP));

% calcite erosion (w<0)
wcvA(lA) = -(1-fc0A(lA)).*wvA(lA)'.*(1-phiiA(lA))/(1-phi0);
wcvI(lI) = -(1-fc0I(lI)).*wvI(lI)'.*(1-phiiI(lI))/(1-phi0);
wcvP(lP) = -(1-fc0P(lP)).*wvP(lP)'.*(1-phiiP(lP))/(1-phi0);
wcvA(lA) =   wcvA(lA) + rsrvA(lA)';
wcvI(lI) =   wcvI(lI) + rsrvI(lI)';
wcvP(lP) =   wcvP(lP) + rsrvP(lP)';
fbA (lA) = 0.;
fbI (lI) = 0.;
fbP (lP) = 0.;


%====== Carbon-13
w13cvA = RsA.*wcvA;
w13cvI = RsI.*wcvI;
w13cvP = RsP.*wcvP;

% dphi/dfc
dphiA = FF*(1-phi0)./(1+FF*fcvA).^2;
dphiI = FF*(1-phi0)./(1+FF*fcvI).^2;
dphiP = FF*(1-phi0)./(1+FF*fcvP).^2;
% G's
GA    = hs*(1-phiA-fcvA.*dphiA)/(1-phi1);
GI    = hs*(1-phiI-fcvI.*dphiI)/(1-phi1);
GP    = hs*(1-phiP-fcvP.*dphiP)/(1-phi1);

% dissolution (w>0) in mol/y, see above    

% dissolution (w<0) in mol/y
dissA(lA) = (1*rsvA(lA)-wvA(lA))*(1-phi1)*rhos/m2kg; % mol/m2/y
dissI(lI) = (1*rsvI(lI)-wvI(lI))*(1-phi1)*rhos/m2kg; % mol/m2/y
dissP(lP) = (1*rsvP(lP)-wvP(lP))*(1-phi1)*rhos/m2kg; % mol/m2/y


%-------------- Tethys --------------------------%
if(ftys)
omvT = co3(kliT  )./co3satv; % Tys
kdT  = find(omvT >= 1.);
RsT  = f13cvT./fcvT;
phiT = (phi0+FF*fcvT)./(1+FF*fcvT); 
rscvT = FprvT*m2kg/rhos/(1-phi1); 
rsrvT = frrfT     /rhos/(1-phi0); 
rsvT  = rscvT+rsrvT;
r13scvT = Fpr13vT*m2kg/rhos/(1-phi1); 
if(CAvflag == 1)
dKT   = KS*((co3satv-co3(kliT))*rcak).^nc; % mol/m2/y
elseif(CAvflag == 2)
 dKT   = KS*((co3satv-co3(kliT)).*rcak(kliT)).^nc; % mol/m2/y  
end
dissT = fcvT.^0.5.*dKT';                   % mol/m2/y
dissT(jT) = fcvT(jT).^nf.*dKT(jT)'*Bf;     % mol/m2/y
dissT(kdT) = 0.;
rdvT  = dissT'*m2kg/rhos/(1-phi1);     % m/y
r13dvT = RsT.*rdvT';
wvT   = rsvT-rdvT; 
lT = find(wvT < 0.);
wcvT  = fcvT.*wvT'.*(1-phiT)/(1-phi1);
fbT   = 1*ones(size(rdvT));
wcvT(lT) = -(1-fc0T(lT)).*wvT(lT)'.*(1-phiiT(lT))/(1-phi0);
fbT (lT) = 0.;
w13cvT = RsT.*wcvT;
dphiT = FF*(1-phi0)./(1+FF*fcvT).^2;
GT    = hs*(1-phiT-fcvT.*dphiT)/(1-phi1);
dissT(lT) = (1*rsvT(lT)-wvT(lT))*(1-phi1)*rhos/m2kg; % mol/m2/y

dissCv(:,:) = [dissA'.*asvA*A(01);  ...
               dissI'.*asvI*A(02);  ...
               dissP'.*asvP*A(03);  ...
               dissT'.*asvT*A(11)]; % -> mol C/y
           
dissC13v(:,:) = ...
        [RsA'.*dissA'.*asvA*A(01);  ...
         RsI'.*dissI'.*asvI*A(02);  ...
         RsP'.*dissP'.*asvP*A(03);  ...
         RsT'.*dissT'.*asvT*A(11)];
else  
% dissolution in mol/y
dissCv(:,:) = [dissA'.*asvA*A(01);  ...
               dissI'.*asvI*A(02);  ...
               dissP'.*asvP*A(03)]; % -> mol C/y


%====== Carbon-13
dissC13v(:,:) = [RsA'.*dissA'.*asvA*A(01);  ...
                 RsI'.*dissI'.*asvI*A(02);  ...
                 RsP'.*dissP'.*asvP*A(03)]; % -> mol C/y   
end;
% burCv(:,:) = [wcvA'.*asvA*A(01);  ...            %
%                wcvI'.*asvI*A(02);  ...           % !!!!!!THIS ADDED
%                wcvP'.*asvP*A(03);  ...           %       FOR BURIAL
%                wcvT'.*asvT*A(11)]; % -> mol C/y


% mass balance vars (dYflag)
if(1)
rbvA = wvA;
rbvI = wvI;
rbvP = wvP;
dissvA = dissA'.*asvA*A(01)*m2kg; % kg/y
dissvI = dissI'.*asvI*A(02)*m2kg; % kg/y
dissvP = dissP'.*asvP*A(03)*m2kg; % kg/y
dissv13A = RsA'.*dissvA;
dissv13I = RsI'.*dissvI;
dissv13P = RsP'.*dissvP;
if(ftys)
rbvT = wvT;
dissvT = dissT'.*asvT*A(11)*m2kg; % kg/y
dissv13T = RsT'.*dissvT;
end;    
end;


fcvAIP = [fcvA fcvI fcvP];
if(ERR &  ~isempty(find(fcvAIP < 0.))  )
disp('+++ ERROR: fc is negative +++');    
%fcvAIP    
%return;
end;  

if(ERR &  ~isempty(find(fcvAIP > 1.))  )
disp('+++ ERROR: fc > 1        +++');    
end;  


%%%%================ END Sediment Stuff =============%%%%
end; % fsed: sediments


% TH branches
gA = (1-tA);    %        export into     Deep Ind 
gI = (1-tA-tI); %        export into     Deep Pac


% volcanic degassing
Fvc   = FVC;

% CaSiO3 weathering
FSi   = FVC*(pco2a/pCSi)^nSi; % CaSiO3
FSichck(kt) = FSi;
% CaCO3  weathering
Fin   = FiN*(pco2a/pCSi)^nCC; % CaC O3
Finchck (kt) = Fin;
% All C13 Fluxes
FSi13 = FSi*Rin;
Fin13 = Fin*Rin;
Fvc13 = Fvc*Rvc;


if(~fsed)
Fin   = 0*FiN;
Fin13 = 0*FiN13;
FSi   = 0.;
FSi13 = 0.;
Fvc   = 0.;
Fvc13 = 0.;
Fkg   = 0.;
Fkg13 = 0.;
end;

if(dYflag == 1)
if(TDflag == 0)
    TCvt(it,:)    = TCv;    
end    
     THt(it  )    = TH;    
     ept(it  )    = ep;    
     ect(it  )    = ec;    
if(fsed) % sediments
    rburvtA(it,:) = rbvA;
    rburvtI(it,:) = rbvI;
    rburvtP(it,:) = rbvP;
    if(ftys)
    rsedvtT(it,:) = rbvT;
    end;        
    dissvtA(it,:) = dissvA;
    dissvtI(it,:) = dissvI;
    dissvtP(it,:) = dissvP;
  dissv13tA(it,:) = dissv13A;
  dissv13tI(it,:) = dissv13I;
  dissv13tP(it,:) = dissv13P;
    phivtA (it,:) = phiA';       
    phivtI (it,:) = phiI';       
    phivtP (it,:) = phiP';       
    FprtA  (it  ) = FprA;       
    FprtI  (it  ) = FprI;       
    FprtP  (it  ) = FprP;       
    Fpr13tA(it  ) = Fpr13A;       
    Fpr13tI(it  ) = Fpr13I;       
    Fpr13tP(it  ) = Fpr13P;       
    if(ftys)
    dissvtT(it,:) = dissvT;
  dissv13tT(it,:) = dissv13T;
    phivtT (it,:) = phiT';       
    FprtT  (it  ) = FprT;       
    Fpr13tT(it  ) = Fpr13T;       
    end;        
    Fint   (it  ) = Fin;       
    Fin13t (it  ) = Fin13;       
    FSit   (it  ) = FSi;       
    FSi13t (it  ) = FSi13; 
end; % sediments
%%%    it = it+1; % increased below (blflag)
end; % dYflag

% set all derivs to zero
%yp(1:3*Nb+1) = 0;
%yp(1:3*Nb+1+Ns) = 0;
%yp(1:3*Nb+1+3*Ns) = 0;
%yp(1:4*Nb+2+3*Ns) = 0;
if(fsed)
if(ftys)    
yp(1:5*Nb+2+8*Ns) = 0;
else
yp(1:5*Nb+2+6*Ns) = 0;
end;    
EXLv   = EPLv;
EXLvCC = EPLvCC;
EXH    = EPH;
EXHCC  = EPHCC;
x      = 0;
else % fsed
yp(1:5*Nb+2     ) = 0;
EXLv   = ECLv;
EXLvCC = ECLvCC;
EXH    = ECH;
EXHCC  = ECHCC;
x      = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Here is the right-hand side of the DEQ
%
%     Units ocean tracer: [mol/m3/y]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=================== TCO2 ======================%
% TH & mixing
cp = THmfun(c,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
% air-sea
for k=kkv
cp(k ) = cp(k ) + ...                           %
           kasv(k)*(pco2a-pco2(k))/V(k);        %
end;  
% bio pump Corg
for k=1:3
cp(k  ) = cp(k )  -    ECLv(k)/V(k  );    % L 
%tmpCAR=cp;
%save                    cp.DAT  tmpCAR -ASCII -DOUBLE -TABS;
%CHKcarB1(kt) = ECLv(k)/V(k  );
cp(k+3) = cp(k+3) + eI*EXLv(k)/V(k+3);    % I #!  EC or EP
%CHKcarB2(kt) = eI*EXLv(k)/V(k+3);
cp(k+6) = cp(k+6) + oI*EXLv(k)/V(k+6) ... % D #!(tot or Corg)
                  + nu*EALv(k)/V(k+6)/2;  % D ClmnDiss
%CHKcarB3(kt) = oI*EXLv(k)/V(k+6) ... 
%                  + nu*EALv(k)/V(k+6)/2;
end;    
cp(10)  = cp(10)  -    ECH/V10;           % H
for k=7:9                                 % DA,DI,DP
cp(k )  = cp(k )  +    EXH/V(k)*gp(k)...  % #!  EC or EP
                  + nu*EAH/V(k)*gp(k)/2;  %   
end;
if(ftys)
k = 11;
cp(k  ) = cp(k  ) -    ECLv(4)/V(k  );
cp(k+1) = cp(k+1) + eI*EXLv(4)/V(k+1);
cp(k+2) = cp(k+2) + oI*EXLv(4)/V(k+2) ...
                  + nu*EALv(4)/V(k+2)/2;  % D ClmnDiss
end;
% riverine & sediment fluxes
if(fsed)                                           %   #!
for k=1:3
cp(k  ) = cp(k  ) + 2*Fin*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                  + 2*FSi*Aoc/V(k)/nOC         ... % #! Si
                  - 1*Fkg*Aoc/V(k)/nOC;            % #! krgn
% CHKcar1(kt) = 2*Fin*Aoc/V(k)/nOC         ... % 
%                   + 2*FSi*Aoc/V(k)/nOC         ... 
%                   - 1*Fkg*Aoc/V(k)/nOC;           
cp(k  ) = cp(k  ) + sum(dissCv(k,1   :n1))/V(k  ); % diss L
%CHKcar2(kt) =sum(dissCv(k,1   :n1))/V(k  );
cp(k+3) = cp(k+3) + sum(dissCv(k,n1+1:n2))/V(k+3); % diss I
%CHKcar3(kt) = sum(dissCv(k,n1+1:n2))/V(k+3);
cp(k+6) = cp(k+6) + sum(dissCv(k,n2+1:Ns))/V(k+6); % diss D 
%CHKcar4(kt) = sum(dissCv(k,n2+1:Ns))/V(k+6);
end;
if(ftys)
k=11;    
cp(k  ) = cp(k  ) + 2*Fin*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                  + 2*FSi*Aoc/V(k)/nOC         ... % #! Si
                  - 1*Fkg*Aoc/V(k)/nOC;            % #! krgn
cp(k  ) = cp(k  ) + sum(dissCv(4,1   :n1))/V(k  ); % diss L
cp(k+1) = cp(k+1) + sum(dissCv(4,n1+1:n2))/V(k+1); % diss I
cp(k+2) = cp(k+2) + sum(dissCv(4,n2+1:Ns))/V(k+2); % diss D
end;    
end;



%=================== TA   ======================%
% TH & mixing
ap = THmfun(a,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
% bio pump CaCO3, Aorg
for k=1:3
ap(k  ) = ap(k ) -      EALv(k)/V(k  )+ENLv(k)/V(k  ) ;
ap(k+3) = ap(k+3)+eI*(x*EALv(k)/V(k+3)-ENLv(k)/V(k+3));   % 1#!
ap(k+6) = ap(k+6)+oI*(x*EALv(k)/V(k+6)-ENLv(k)/V(k+6))... % 1#!
                 +   nu*EALv(k)/V(k+6);                   % D ClmnDiss
end;    
ap(10)  = ap(10)  - EAH/V10        + ENH/V10;
for k=7:9                                              % DA,DI,DP
ap(k )  = ap(k )  +(x*EAH/V(k)     - ENH/V(k))*gp(k)...% 1#!  
                  +nu*EAH/V(k)                *gp(k);  %   
end;
if(ftys)
k = 11;
ap(k  ) = ap(k  )-      EALv(4)/V(k  )+ENLv(4)/V(k  ) ;
ap(k+1) = ap(k+1)+eI*(x*EALv(4)/V(k+1)-ENLv(4)/V(k+1));
ap(k+2) = ap(k+2)+oI*(x*EALv(4)/V(k+2)-ENLv(4)/V(k+2))... %
                 +   nu*EALv(4)/V(k+2);
end;
% riverine & sediment fluxes
if(fsed)                                           % #!
for k=1:3
ap(k  ) = ap(k  )+2*Fin*Aoc/V(k)/nOC          ...  % Fin  L
                 +2*FSi*Aoc/V(k)/nOC;              % #! Si
ap(k  ) = ap(k  )+2*sum(dissCv(k,1   :n1))/V(k  ); % diss L
ap(k+3) = ap(k+3)+2*sum(dissCv(k,n1+1:n2))/V(k+3); % diss I
ap(k+6) = ap(k+6)+2*sum(dissCv(k,n2+1:Ns))/V(k+6); % diss D
end;
if(ftys)
k=11;
ap(k  ) = ap(k  )+2*Fin*Aoc/V(k)/nOC          ...  % Fin  L
                 +2*FSi*Aoc/V(k)/nOC;              % #! Si
ap(k  ) = ap(k  )+2*sum(dissCv(4,1   :n1))/V(k  ); % diss L
ap(k+1) = ap(k+1)+2*sum(dissCv(4,n1+1:n2))/V(k+1); % diss I
ap(k+2) = ap(k+2)+2*sum(dissCv(4,n2+1:Ns))/V(k+2); % diss D
end;
end;
%=================== PO4   ======================%
% TH & mixing
pp = THmfun(p,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
% bio pump Porg
for k=1:3
pp(k  ) = pp(k )  -    PPLv(k)/V(k  );
pp(k+3) = pp(k+3) + eI*PPLv(k)/V(k+3);
pp(k+6) = pp(k+6) + oI*PPLv(k)/V(k+6);
end;    
pp(10)  = pp(10)  - PPH/V10;
for k=7:9
pp(k )  = pp(k )  + PPH/V(k)*gp(k);             % DA,DI,DP
end;
if(ftys)
k = 11;
pp(k  ) = pp(k  ) -    PPLv(4)/V(k  );
pp(k+1) = pp(k+1) + eI*PPLv(4)/V(k+1);
pp(k+2) = pp(k+2) + oI*PPLv(4)/V(k+2);
end;

% %===============================================
% %=================== CAlcium   ======================%
% %================================================
% % TH & mixing
% CAp = THmfun(CA,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
% % %bio pump CaCO3, Aorg
%  for k=1:3
% CAp(k  ) = CAp(k ) -      ECALv(k)/V(k  )+ENLv(k)/V(k  ) ;
% %CHCKbioL(kt) = ECALv(k)/V(k  );
% CAp(k+3) = CAp(k+3)+(eI*ECALv(k)/V(k+3)-ENLv(k)/V(k+3));
% %CHCKbioI(kt) = eI*ECALv(k)/V(k+3);
% CAp(k+6) = CAp(k+6)+(oI*ECALv(k)/V(k+6)-ENLv(k)/V(k+6))...
%                 +   nu*ECALv(k)/V(k+6)/2;                 % D ClmnDiss
% %CHCKbioD(kt) = oI*ECALv(k)/V(k+6) +nu*ECALv(k)/V(k+6)/2;
% end;    
% %CAp(10)  = CAp(10)  - ECAH/V10        + ENH/V10;
% CAp(10)  =0;% CAp(10)  - ECAH/V10;
% for k=7:9                                              % DA,DI,DP
% CAp(k )  = CAp(k )  +(ECAH/V(k)- ENH/V(k))*gp(k)...    
%                   +nu*ECAH/V(k)                *gp(k)/2;  %   
% end;
% if(ftys)
% k = 11;
% CAp(k  ) = CAp(k  )-      ECALv(4)/V(k  )+ENLv(4)/V(k  );
% CAp(k+1) = CAp(k+1)+(eI*ECALv(4)/V(k+1)-ENLv(4)/V(k+1));
% CAp(k+2) = CAp(k+2)+(oI*ECALv(4)/V(k+2)-ENLv(4)/V(k+2))...
%                  +   nu*ECALv(4)/V(k+2);
% end;
% % riverine & sediment fluxes
% if(fsed)                                           % #!
% for k=1:3
% CAp(k  ) = CAp(k  ) + 2*Fin*Aoc/V(k)/nOC          ...  % Fin  L
%                     + 2*FSi*Aoc/V(k)/nOC;              % #! Si
% %             %- 1*Fkg*Aoc/V(k)/nOC;            % #! krgn
% CHKFin(kt) = Fin;
% CHKFSi(kt) = FSi;
% %         CAp(k) = 0;
% CAp(k  ) = CAp(k  ) + sum(dissCv(k,1   :n1))/V(k  ); % diss L
% CHECK1(kt)=sum(dissCv(k,1   :n1))/V(k  );
% CAp(k+3) = CAp(k+3) + sum(dissCv(k,n1+1:n2))/V(k+3); % diss I
% CAp(k+6) = CAp(k+6) + sum(dissCv(k,n2+1:Ns))/V(k+6); % diss D
% CHECK3(kt)=sum(dissCv(k,n2+1:Ns))/V(k+6);
% %CHECK(kt)=dissCv(k,1   :n1);
% end;
% if(ftys)
% k=11;
% CAp(k  ) = CAp(k  )+2*Fin*Aoc/V(k)/nOC          ...  % Fin  L
%                  +2*FSi*Aoc/V(k)/nOC;              % #! Si
%              
% %             % - 1*Fkg*Aoc/V(k)/nOC;            % #! krgn
% %        CAp(k) = 0;
% CAp(k  ) = CAp(k  )+sum(dissCv(4,1   :n1))/V(k  ); % diss L
% CAp(k+1) = CAp(k+1)+sum(dissCv(4,n1+1:n2))/V(k+1); % diss I
% CAp(k+2) = CAp(k+2)+sum(dissCv(4,n2+1:Ns))/V(k+2); % diss D
% 
% end;
% end;
%=================== CALCIUM  ======================%
% TH & mixing
CAp = THmfun(CA,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
% bio pump CaCO3, Aorg
for k=1:3
CAp(k  ) = CAp(k ) -      ECALv(k)/V(k  );%+ENLv(k)/V(k  ) ;
CAp(k+3) = CAp(k+3)+eI*(x*ECALv(k)/V(k+3));%-ENLv(k)/V(k+3));   % 1#!
CAp(k+6) = CAp(k+6)+oI*(x*ECALv(k)/V(k+6))...%-ENLv(k)/V(k+6))... % 1#!
                 +   nu*ECALv(k)/V(k+6);                   % D ClmnDiss
end;    
CAp(10)  = CAp(10)  - EAH/V10;%        + ENH/V10;
for k=7:9                                              % DA,DI,DP
CAp(k )  = CAp(k )  +(x*EAH/V(k))...%     - ENH/V(k))*gp(k)...% 1#!  
                  +nu*EAH/V(k)                *gp(k);  %   
end;
if(ftys)
k = 11;
CAp(k  ) = CAp(k  )-      ECALv(4)/V(k  );%+ENLv(4)/V(k  ) ;
CAp(k+1) = CAp(k+1)+eI*(x*ECALv(4)/V(k+1));%-ENLv(4)/V(k+1));
CAp(k+2) = CAp(k+2)+oI*(x*ECALv(4)/V(k+2))...%-ENLv(4)/V(k+2))... %
                 +   nu*ECALv(4)/V(k+2);
end;
% riverine & sediment fluxes
if(fsed)                                           % #!
for k=1:3
CAp(k  ) = CAp(k  )+1*Fin*Aoc/V(k)/nOC          ...  % Fin  L
                 +1*FSi*Aoc/V(k)/nOC;                 % #! Si
           
CAp(k  ) = CAp(k  )+1*sum(dissCv(k,1   :n1))/V(k  ); % diss L
CAp(k+3) = CAp(k+3)+1*sum(dissCv(k,n1+1:n2))/V(k+3); % diss I
CAp(k+6) = CAp(k+6)+1*sum(dissCv(k,n2+1:Ns))/V(k+6); % diss D
end;
if(ftys)
k=11;
CAp(k  ) = CAp(k  )+1*Fin*Aoc/V(k)/nOC          ...  % Fin  L
                 +1*FSi*Aoc/V(k)/nOC;       % #! Si
             
CAp(k  ) = CAp(k  )+1*sum(dissCv(4,1   :n1))/V(k  ); % diss L
CAp(k+1) = CAp(k+1)+1*sum(dissCv(4,n1+1:n2))/V(k+1); % diss I
CAp(k+2) = CAp(k+2)+1*sum(dissCv(4,n2+1:Ns))/V(k+2); % diss D
end;
end;



%=================== O2     ======================%
if(fdox)
% TH & mixing
doxp = THmfun(dox,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
% air-sea (O2(k) = saturation)
for k=kkv
doxp(k) = doxp(k) + vask(k)*(o2(k)-dox(k))/V(k);
end;
% bio pump Corg
for k=1:3
KMMOX = 0.001; % 0.001
if(ffflag)
KMMOX = 0.0;
end
mmoxI     = dox(k+3)/(dox(k+3)+KMMOX)*sign(sign(dox(k+3))+1);
mmoxD     = dox(k+6)/(dox(k+6)+KMMOX)*sign(sign(dox(k+6))+1);
doxp(k  ) = doxp(k  ) +    ECLv(k)/V(k  )*REDO2C;       % L 
doxp(k+3) = doxp(k+3) - eI*EXLv(k)/V(k+3)*REDO2C*mmoxI; % I #!  EC or EP
doxp(k+6) = doxp(k+6) - oI*EXLv(k)/V(k+6)*REDO2C*mmoxD; % D #!(tot or Corg)
end;    
for k=7:9 % DA,DI,DP
mmoxD     = dox(k  )/(dox(k  )+KMMOX)*sign(sign(dox(k  ))+1); 
doxp(k )  = doxp(k  ) -    EXH/V(k)*gp(k)*REDO2C*mmoxD; % #!  EC or EP
end;
doxp(10)  = doxp(10 ) +    ECH/V10       *REDO2C;       % H
if(ftys)
k = 11;
mmoxI     = dox(k+1)/(dox(k+1)+KMMOX)*sign(sign(dox(k+1))+1);
mmoxD     = dox(k+2)/(dox(k+2)+KMMOX)*sign(sign(dox(k+2))+1);
doxp(k  ) = doxp(k  ) +    ECLv(4)/V(k  )*REDO2C;
doxp(k+1) = doxp(k+1) - eI*EXLv(4)/V(k+1)*REDO2C*mmoxI;
doxp(k+2) = doxp(k+2) - oI*EXLv(4)/V(k+2)*REDO2C*mmoxD;
end;
end;
%=================== T13CO2 ======================%
% TH & mixing
ccp = THmfun(cc,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
% air-sea
for k=kkv
ccp(k ) = ccp(k ) + ...                         %
   kasv(k)*au*(adg(k)*pcco2a-pcco2(k))/V(k);    %
end;  
% bio pump Corg
for k=1:3
ccp(k  ) = ccp(k )  -    ECLvCC(k)/V(k  );      % L 
ccp(k+3) = ccp(k+3) + eI*EXLvCC(k)/V(k+3);      % I #!  EC or EP
ccp(k+6) = ccp(k+6) + oI*EXLvCC(k)/V(k+6)...    % D #!(tot or Corg)
                    + nu*EALvCC(k)/V(k+6)/2;    % D ClmnDiss
end;    
ccp(10)  = ccp(10)  -    ECHCC/V10;             % H
for k=7:9                                       % DA,DI,DP
ccp(k )  = ccp(k )  +    EXHCC/V(k)*gp(k)...    % #!  EC or EP
                    + nu*EAHCC/V(k)*gp(k)/2;    %   
end;                                            %   
if(ftys)
k = 11;
ccp(k  ) = ccp(k  ) -    ECLvCC(4)/V(k  );
ccp(k+1) = ccp(k+1) + eI*EXLvCC(4)/V(k+1);
ccp(k+2) = ccp(k+2) + oI*EXLvCC(4)/V(k+2) ...
                    + nu*EALvCC(4)/V(k+2)/2;  % D ClmnDiss
end;
% riverine & sediment fluxes
if(fsed)                                        %   #!
for k=1:3
ccp(k  ) = ccp(k  ) + 2*Fin13*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                    + 2*FSi13*Aoc/V(k)/nOC         ... % #! Si
                    - 1*Fkg  *Aoc/V(k)/nOC*alp*Rb(k);  % #! krgn
ccp(k  ) = ccp(k  ) + sum(dissC13v(k,1   :n1))/V(k  ); % diss L
ccp(k+3) = ccp(k+3) + sum(dissC13v(k,n1+1:n2))/V(k+3); % diss I
ccp(k+6) = ccp(k+6) + sum(dissC13v(k,n2+1:Ns))/V(k+6); % diss D
end;
if(ftys)
k=11;    
ccp(k  ) = ccp(k  ) + 2*Fin13*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                    + 2*FSi13*Aoc/V(k)/nOC         ... % #! Si
                    - 1*Fkg  *Aoc/V(k)/nOC*alp*Rb(k);  % #! krgn
ccp(k  ) = ccp(k  ) + sum(dissC13v(4,1   :n1))/V(k  ); % diss L
ccp(k+1) = ccp(k+1) + sum(dissC13v(4,n1+1:n2))/V(k+1); % diss I
ccp(k+2) = ccp(k+2) + sum(dissC13v(4,n2+1:Ns))/V(k+2); % diss D
end;    
end;




%=================== C  atm ====================%
Cp    = sum(                                 ...%
      -kasv(kkv).*(pco2a-pco2(kkv))          ...% mol/m2/y
        )/Aoc                                ...%
      -1*Fin   + Fvc   - 2*FSi   + Fkg;         % wthr #!
%=================== C13 atm ===================%
CCp   = sum(                                 ...%
-kasv(kkv).*au.*(adg(kkv)*pcco2a-pcco2(kkv)) ...% mol/m2/y
        )/Aoc                                ...%  
      -1*Fin13 + Fvc13 - 2*FSi13 + Fkg13;       % wthr #!
% Anthropogenic CO2
if(ffflag)
Cp    = Cp + Fem  /Aoc; % input(t) (mol/y)-> mol/m2/y    
CCp   = CCp+ Fem13/Aoc; % input(t) (mol/y)-> mol/m2/y    
end;

%============== Sediment Boxes =================%
%
if(fsed) % sediments                                  % #!
for l=1:Ns
fcvAp(l)   = ( fbA(l)*rscvA(l) - fbA(l)*rdvA(l) - wcvA(l) )/GA(l);
fcvIp(l)   = ( fbI(l)*rscvI(l) - fbI(l)*rdvI(l) - wcvI(l) )/GI(l);
fcvPp(l)   = ( fbP(l)*rscvP(l) - fbP(l)*rdvP(l) - wcvP(l) )/GP(l);
f13cvAp(l) = ( fbA(l)*r13scvA(l) - fbA(l)*r13dvA(l) - w13cvA(l) )/GA(l);
f13cvIp(l) = ( fbI(l)*r13scvI(l) - fbI(l)*r13dvI(l) - w13cvI(l) )/GI(l);
f13cvPp(l) = ( fbP(l)*r13scvP(l) - fbP(l)*r13dvP(l) - w13cvP(l) )/GP(l);
if(ftys)
fcvTp(l)   = ( fbT(l)*rscvT(l) - fbT(l)*rdvT(l) - wcvT(l) )/GT(l);
f13cvTp(l) = ( fbT(l)*r13scvT(l) - fbT(l)*r13dvT(l) - w13cvT(l) )/GT(l);
end;    
end; % for
else % sediments
 for l=1:Ns
 fcvAp(l)  = 0;           
 fcvIp(l)  = 0;           
 fcvPp(l)  = 0;  
 mcAp(l)   = 0;           
 mcIp(l)   = 0;           
 mcPp(l)   = 0;  
 f13cvAp(l)= 0;           
 f13cvIp(l)= 0;           
 f13cvPp(l)= 0;           
 m13cAp(l) = 0;           
 m13cIp(l) = 0;           
 m13cPp(l) = 0;           
 if(ftys)
 fcvTp(l)  = 0;  
 mcTp(l)   = 0;  
 f13cvTp(l)= 0;           
 m13cTp(l) = 0;           
 end;     
end;   
end; % sediments 


%tmpCO(1) = cp(kb);
tmpCA    = Cp;

% Tethys
if(ftys)
kT = 13;
tmpCO(2) = cp(kT);
else;
kT = 1;
tmpCO(2) = 0.0;
end;

if(BlFlag == 2)
% inital shot (e.g. 2x < 1,000y, Roehl) 1111111111111111111
DTS = 06.e3;   % n06 [10 1] e3 years 04
fB  = oxA;     % Atl fraction here  0.460 .4 n.32
fx  = oxA;     % Atl fraction below 1 .4
fY  = 0.000;   % Tys fraction 0.000
fM  = 1-fB-fY; % Atm fraction
if(t >= tcon & t <= tcon + DTS)
  cp(kb) =  cp(kb)+   fB *    CBl/12/V(kb)/DTS;  % add X Gt   C ocn:/V2
 ccp(kb) = ccp(kb)+   fB *RBl*CBl/12/V(kb)/DTS;  % add X Gt 13C ocn:/V2
 
  Cp     =  Cp    +   fM *    CBl/12/Aoc  /DTS;  % add X Gt C atm:/Aoc
 CCp     = CCp    +   fM *RBl*CBl/12/Aoc  /DTS;  % add X Gt C atm:/Aoc
%
% oxygen
doxp(kb) =doxp(kb)-   fB*   2*CBl/12/V(kb)/DTS;  %
% Tethys
  cp(kT) =  cp(kT)+   fY*     CBl/12/V(kT)/DTS;  % add X Gt   C ocn:/V2
  
 ccp(kT) = ccp(kT)+   fY* RBl*CBl/12/V(kT)/DTS;  % add X Gt   C ocn:/V2
doxp(kT) =doxp(kT)-   fY*   2*CBl/12/V(kT)/DTS;  %
CAp(kb)=0;%%% THIS ADDED
CAp(kT)=0; %%%THIS ADDED

end;


if(1)
% continuous release                  22222222222222222222
DTS2 = 14.e3; % 65 55 75e3 years 20 n14
fB  = 0.0*(  fx);  % Atl fraction 0 n.15 .2 .1
fM  = 0.0*(1-fx);  % Atm fraction 0 n.15 .2 .1
%if(t >= tcon + 1*DTS & t < tcon + 1*DTS + DTS2) % 5*
if(t >= tcon + DTS & t < tcon + DTS + DTS2) % 5*
  cp(kb) =  cp(kb)+fB*    CBl/12/V(kb)/DTS2;  % add X Gt   C ocn:/V2
 ccp(kb) = ccp(kb)+fB*RBl*CBl/12/V(kb)/DTS2;  % add X Gt 13C ocn:/V2
   Cp    =  Cp    +fM*    CBl/12/Aoc  /DTS2;  % add X Gt C atm:/Aoc
  CCp    = CCp    +fM*RBl*CBl/12/Aoc  /DTS2;  % add X Gt C atm:/Aoc
%
% oxygen
doxp(kb) =doxp(kb)-fB*  2*CBl/12/V(kb)/DTS2;  %
end;
end;


if(1)
% another bump                        333333333333333333333
ts3 = 15.e3;  % onset
DTS3= 05.e3;  % 06 [10 1] e3 years n05
fB  = .1*(  fx);  % .07 Atl fraction .6 .047 n.1
fM  = .1*(1-fx);  % .07 Atm fraction .6 .047 n.1
if(t >= tcon + ts3 & t <= tcon + ts3 + DTS3)
%if(t >= tcon + 16.e3 & t <= tcon + 26.e3)
  cp(kb) =  cp(kb)+   fB *    CBl/12/V(kb)/DTS3;  % add X Gt   C ocn:/V2
 ccp(kb) = ccp(kb)+   fB *RBl*CBl/12/V(kb)/DTS3;  % add X Gt 13C ocn:/V2
  Cp     =  Cp    +   fM *    CBl/12/Aoc  /DTS3;  % add X Gt C atm:/Aoc
 CCp     = CCp    +   fM *RBl*CBl/12/Aoc  /DTS3;  % add X Gt C atm:/Aoc
%
% oxygen
doxp(kb) =doxp(kb)-   fB*   2*CBl/12/V(kb)/DTS3;  %
end;
end;


if(1)
% continuous release                  44444444444444444444444
DTS4 = 42.e3; % 65 55 75e3 years 42 50
fB  = .4*(  fx);  % .75 Atl fraction .9 .7 .47 .6 n.5
fM  = .4*(1-fx);  % .75 Atm fraction .9 .7 .47 .6 n.5
%if(t >= tcon + 1*DTS & t < tcon + 1*DTS + DTS4) % 5*
%if(t >= tcon + DTS & t < tcon + 75.e3) % 
if(t >= tcon + ts3 + DTS3 & t < tcon + ts3 + DTS3 + DTS4) % 5*
  cp(kb) =  cp(kb)+fB*    CBl/12/V(kb)/DTS4;  % add X Gt   C ocn:/V2
 ccp(kb) = ccp(kb)+fB*RBl*CBl/12/V(kb)/DTS4;  % add X Gt 13C ocn:/V2
   Cp    =  Cp    +fM*    CBl/12/Aoc  /DTS4;  % add X Gt C atm:/Aoc
  CCp    = CCp    +fM*RBl*CBl/12/Aoc  /DTS4;  % add X Gt C atm:/Aoc
%
% oxygen
doxp(kb) =doxp(kb)-fB*  2*CBl/12/V(kb)/DTS4;  %
end;
end;
end; % BlFlag


if(dYflag == 1)
 if(kb)
 RlsCtv(it) = (cp(kb) - tmpCO(1))*12*V(kb) ...  % g C
            + (cp(kT) - tmpCO(2))*12*V(kT) ...  % g C
            + (Cp     - tmpCA)*12*Aoc;          % g C
 end;
it = it+1;
end;



if(0)
 cp(11:13) =          0;         % T
 ap(11:13) =          0;         % T
 pp(11:13) =          0;         % T
 CAp(11:13)=          0;
ccp(11:13) =          0;         % T
end;

% all in one  

if(fdox)  
yp(     1:  Nb     ) = cp;
yp(  Nb+1:2*Nb     ) = ap;
yp(2*Nb+1:3*Nb     ) = pp;
yp(3*Nb+1:4*Nb     ) = CAp ;
yp(4*Nb+1:5*Nb     ) = doxp;
yp(5*Nb+1:6*Nb     ) = ccp;
yp(6*Nb+1          ) = Cp;
yp(6*Nb+2          ) = CCp;
if(fsed)
yp(6*Nb+3     :6*Nb+2+1*Ns) = fcvAp; 
yp(6*Nb+3+1*Ns:6*Nb+2+2*Ns) = fcvIp; 
yp(6*Nb+3+2*Ns:6*Nb+2+3*Ns) = fcvPp; 
yp(6*Nb+3+3*Ns:6*Nb+2+4*Ns) = f13cvAp; 
yp(6*Nb+3+4*Ns:6*Nb+2+5*Ns) = f13cvIp; 
yp(6*Nb+3+5*Ns:6*Nb+2+6*Ns) = f13cvPp;
if(ftys)
yp(6*Nb+3+6*Ns:6*Nb+2+7*Ns) = fcvTp; 
yp(6*Nb+3+7*Ns:6*Nb+2+8*Ns) = f13cvTp;    
end;    
end;
else % fdox
yp(     1:  Nb     ) = cp;
yp(  Nb+1:2*Nb     ) = ap;
yp(2*Nb+1:3*Nb     ) = pp;
yp(3*Nb+1:4*Nb     ) = CAp ;
yp(4*Nb+1:5*Nb     ) = ccp;
yp(5*Nb+1          ) = Cp;
yp(5*Nb+2          ) = CCp;
if(fsed)
yp(5*Nb+3     :5*Nb+2+1*Ns) = fcvAp; 
yp(5*Nb+3+1*Ns:5*Nb+2+2*Ns) = fcvIp; 
yp(5*Nb+3+2*Ns:5*Nb+2+3*Ns) = fcvPp; 
yp(5*Nb+3+3*Ns:5*Nb+2+4*Ns) = f13cvAp; 
yp(5*Nb+3+4*Ns:5*Nb+2+5*Ns) = f13cvIp; 
yp(5*Nb+3+5*Ns:5*Nb+2+6*Ns) = f13cvPp;
if(ftys)
yp(5*Nb+3+6*Ns:5*Nb+2+7*Ns) = fcvTp; 
yp(5*Nb+3+7*Ns:5*Nb+2+8*Ns) = f13cvTp;    
end;    
end;
end; % fdox


yp = yp';


end; % myflag == 01



%============================%
% function flin
%============================%
function YI = flin(ts,te,Ys,Ye,t)

dt = te-ts;
dY = Ye-Ys;

if    (t <= ts)
  YI = Ys;
elseif(t >= te)
  YI = Ye;
else
  YI = Ys + dY.*(t-ts)/dt;
end;  

return;
