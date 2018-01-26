%--------------------------------------
%
% file: Loscar.m
%
% LOSCAR Model: Long-term Ocean-atmosphere-Sediment 
%               CArbon cycle Reservoir Model
%
% ocean box model (+1 for atmosphere)
%
% 
%    
%
% updates: 
%
% 01/20/09 (Malte) List of Variable description (begin); a few changes
% regarding the flags in the 'save' section
%
% 01/15/09 (Malte) LineWidth control at the beginning of plot section;
% comments
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
% 
%-------------------------------------------
%    VARIABLE AND PARAMETER DESCRIPTIONS   %      
%-------------------------------------------
%
% t0 and tfinal are defined somewhere around line 750
%
%
% Mass Balance:
% DXs = Change of [Mass_Carbon Alkalinity CCD_depth] between initial cond. and last timestep
% Mci = initial Mass Carbon
% Mc = Mass Carbon at last timestep
%
%
% output variables: (in the diff.equations the SI-Units (mol, kg, m3...)
% are used!
%
% c = TCO2 (Total CO2) [mmol/kg]
% cc = T13CO2  [mol/m3]
% ccdA, ccdI, ccdP = CCD (Carbonate compensation depth [m]
% a = TA (Total Alkalinity)
% p = Phosphate [mikromol/kg]
% dox = [% Oxic Respiration]
% co3tv = carbonate ion conz. [mikromol/kg]
% pco2v = partial pressure CO2 [mikro atm]
% pHtv = pH
% d13c = c13/c12 ratio
% fc = % CaCO3 amount in sediments (dry weight)
% mc = CaCO3 in sediment box [kg/m3]
%
% FiN = CaCO3 in-flux [mol/(m2*y)]
% 
% C =  Catm    [mol/m2]
% CC = 13Catm  [mol/m2]
% 
%----------------------------------------


solflag = 1; % 0: skip solver/load
if (solflag == 1)
    clear all all;
    solflag = 1;
end; 
logax = 0;               % plot: log axes on/off
axx   = [-0.5e5 2.0e5];  % x-axes limits (time)
%axx   = [000 1000];

global myflag kasflag Fem20 kt tst Dtst yst dYst stflag ...
       Cam Ca Mgm Mg y2s;

kt  = 1;       % counter time step

Cam = 10.3e-3; % 10.3 (mol/kg) Calcium   modern
Mgm = 55.2050e-3; % 53.0 (mol/kg) Magnesium modern
Ca  = Cam;
Mg  = Mgm;





y2s = 3600.*24.*365.; % year to seconds

myflag = 01; 
             % 01: 10-box model + N sediment + Tethys
             
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%			10 BOX + sediment + Tethys
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if myflag  == 01;
if solflag == 01;


%========================================
% global variables:
% known to Loscar and LoscarDif
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
    THt mv0 mhd0 oxA;  


% +++++ Edit this line to add correct path to solver ! +++++ %

addpath('~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/myode');

plotflag = 0; % plot results

loadf   = 1; % load initial steady-state; check if pRef and pCsi are set right according to the initial steady state loaded
BlFlag  = 0; % Blast; 0: no carbon blast, 1: cont. release over time, 2: cont release with pulses, 4: test  

fcon    = 1; % 1: NADW, 2: NPDW, 3: SO (PETM)
swcon   = 0; % 1: SO -> NP -> SO switch

fpaleo  = 0; % use fpaleo when using bath =>5 to overwrite the ftys flag at some points in the code!
bath    = 2; % 1,2,3 bathymetry (use 2 for modern and petm runs);(4 = 15 Ma Bathymetry, 5 = 34 Ma Bathymetry, 6 = 56 Ma Bathymetry, 7 = 67 Ma Bathymetry)
dsflag  = 3; % 1,2,3 dissolution parameter
   
ftys    = 0; % Tethys
fsed    = 1; % include sediments
parflag = 1; % write parameter to file

ffflag  = 0; % Anthropogenic CO2 (fossil fuel)
fdox    = 1; % include dissolved oxygen

TDflag  = 0; % Temp sens to doubling CO2

ccdrun  = 0;   % parameter run: CBl, oxA
oxA     = 0; % fraction released in deep Atl 0 0.4

fmalte  = 1;  % pRef and pCSi = xxxx, otherwise pRef and pCsi = 280 (modern ocean)
fsaveendstate = 1;  % save Y_end of the run; only working if parflag=1
contourflag = 1;     % see Section; saves variable X

% other changes (no flags):
Tempchange = 6.94;           % change temp. of each box with this value (normal = 0)


% CHECKLIST time series runs:
% 1.) turn Blastflag off
% 2.) adjust bathymety (bath + fpaleo flag) and circulation mode (fcon
% flag); adjust fA3, fT and VD !
% 3.) adjust basin temp. offset (Tempchange)
% 4.) adjust initial pCO2 (pRef and PCSi and the ratio of both (weathering rate), fmalte)
% 5.) adjust Mg and Ca
% 6.) adjust CCD (fsh and nshT (shelf to open ocean production ratio))
% 7.) save steady state
% 8.) turn Blastflag on / load steady state
% 9.) RUN


disp('    ');
disp('@==================== RUN Loscar ====================@');
disp('    ');
load_Blast_fcon_bathym_diss_TDflag_fmalte_fsaveendstate_contourflag = ...
    sprintf('  %d    %d    %d     %d      %d     %d     %d      %d      %d',...
    loadf,BlFlag,fcon,bath,dsflag,TDflag,fmalte,fsaveendstate,contourflag);
load_Blast_fcon_bathym_diss_TDflag_fmalte_fsaveendstate_contourflag


% Tracer -> DIC ALK PO4 O2 DIC-13 Catm Catm-13 sediments ...
%
% Boxes:
% 
% Low, Interm, Deep, High
% Atlantic, Pacific, Indic, (Tethys)
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
%(TL  11
% TI  12
% TD  13)

Voc  = 1.2918235e18;    % (m3) volume ocean
Aoc  = 3.49e14;         % (m2) area ocean
Hav  = Voc/Aoc;         % (m)  average depth

rho  = 1.025e3;         % kg/m3 (1.025: Toggweiler)
rhos = 2.50e3;          % kg/m3 sed. density 2.50
m2kg = 100/1e3;         % mol C -> kg CaCO3

REDPC  =   1/130; % 130     Redfield  P:C
REDNC  =  15/130; % 15/130  Redfield  N:C
REDO2C = 165/130; % 165/130 Redfield O2:C 169

% flags CO2 system
phflag   = 0;     
k1k2flag = 1;	

% Number of oceans
if(ftys)
nOC = 4;
else
nOC = 3;
end;    

on3  = ones(1,3);




%----------------------------------------%
%----------------------------------------%
%-- discrete 5 Ma bathymetry timesteps --%
%----------- basin sizes ----------------%
%----------------------------------------%

% O  Ma  = bath 11 (modern)
% 5  Ma  = bath 12 (AIP interp. between 0 and 15)
% 10 Ma  = bath 13 (AIP interp. between 0 and 15)
% 15 Ma  = bath 14 (Herold)
% 20 Ma  = bath 15 (AIP interp. between 15 and 56)
% ----- (23 Ma: Tethys closure) ----- %
% 25 Ma  = bath 16 (AIP interp. between 15 and 56, Tethys between 23 and 56)
% 30 Ma  = bath 17 (AIP interp. between 15 and 56, Tethys between 23 and 56)
% ----- (33 Ma: circulation switch NADW --> SO/NADW)
% 35 Ma  = bath 18 (AIP interp. between 15 and 56, Tethys between 23 and 56)
% 40 Ma  = bath 19 (AIP interp. between 15 and 56, Tethys between 23 and 56)
% 45 Ma  = bath 20 (AIP interp. between 15 and 56, Tethys between 23 and 56)
% 50 Ma  = bath 21 (AIP interp. between 15 and 56, Tethys between 23 and 56)
% 56 Ma  = bath 22 (Bice)
% 60 Ma  = bath 23 (AIP interp. between 56 and 67, Tethys between 56 and 67)
% 67 Ma  = bath 24 (Sewall)


if(11 <= bath) & (bath <= 15) % 0 Ma - 20 Ma (no Tethys)
%--------------------- Ocean Boxes
%        A   I   P 
%fA3  = [.26 .18 .46    ];    % Area fraction Modern 0 Ma;     bath = 11
% fA3  = [.25 .18 .47  ];    % Area fraction 5  Ma;           bath = 12
% fA3  = [.24 .18 .48  ];    % Area fraction 10  Ma;          bath = 13
% fA3  = [.23 .19 .48  ];    % Area fraction 15 Ma (Herold);  bath = 14
 fA3  = [.22 .19 .49  ];    % Area fraction 20  Ma;          bath = 15

fA   = [fA3 fA3 fA3 .10];  % 
A    = fA*Aoc;
HLI  = [100. 900.];         % (m) height L I boxes
HD   = Hav-sum(HLI);
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 ...
        HLID(3)*on3 250.];  % (m) height of boxes
V    = fA.*H*Aoc;           % (m) Volume of boxes
% volume below H box = A(10)*(Hav-H(10))
% add to deep boxes
V(7:9) = V(7:9) + A(10)*(Hav-H(10))/3;
H = V./A;    
end

if(16 <= bath) & (bath <= 24)
    
% fA3  = [.21 .19 .49  ];    % Area fraction 25  Ma;           bath = 16
% fA3  = [.20 .19 .49  ];    % Area fraction 30  Ma;           bath = 17
% fA3  = [.19 .18 .50  ];    % Area fraction 35  Ma;           bath = 18
% fA3  = [.18 .17 .50  ];    % Area fraction 40  Ma;           bath = 19
% fA3  = [.17 .16 .51  ];    % Area fraction 45  Ma;           bath = 20
% fA3  = [.16 .16 .51  ];    % Area fraction 50  Ma;           bath = 21
% fA3  = [.15 .14 .52];        % Area fraction 56  Ma; (Bice)    bath = 22
% fA3  = [.15 .14 .53  ];    % Area fraction 60  Ma;           bath = 23
 fA3  = [.15 .13 .55  ];    % Area fraction 67  Ma; (Sewall)  bath = 24


fH   = 0.10;                 % Area fraction H box


% fT = 0.01;                 % Area fraction 25  Ma;           bath = 16
% fT = 0.02;                 % Area fraction 30  Ma;           bath = 17
% fT = 0.03;                 % Area fraction 35  Ma;           bath = 18
% fT = 0.05;                 % Area fraction 40  Ma;           bath = 19
% fT = 0.06;                 % Area fraction 45  Ma;           bath = 20
% fT = 0.07;                 % Area fraction 50  Ma;           bath = 21
% fT   = 0.09;                 % Area fraction 56  Ma;           bath = 22
% fT = 0.08;                 % Area fraction 50  Ma;           bath = 23
 fT = 0.07;                 % Area fraction 67  Ma;           bath = 24



fA   = [fA3 fA3 fA3 fH fT*on3]; %
A    = fA*Aoc;
HLI  = [100. 900.];             % (m) height L I boxes
DTM  = sum(HLI);                % (m) depth thermocline
HH   = 250.;                    % (m) depth H box
%HDT = [1000.];                 % (m) height Deep Tethys
HDT  = [ 200.];                 % (m) height Deep Tethys
Vres = Voc-(DTM*(1-fH)+HH*fH+HDT*fT)*Aoc;


% VD   = Vres*[23.0 18.0 59.]/100;% (m3) Vol Deep AIP; 25 Ma;    bath = 16
% VD   = Vres*[22.0 18.0 60.]/100;% (m3) Vol Deep AIP; 30 Ma;    bath = 17
% VD   = Vres*[21.0 18.0 61.]/100;% (m3) Vol Deep AIP; 35 Ma;    bath = 18
% VD   = Vres*[20.0 17.0 63.]/100;% (m3) Vol Deep AIP; 40 Ma;    bath = 19
% VD   = Vres*[19.0 16.0 65.]/100;% (m3) Vol Deep AIP; 45 Ma;    bath = 20
% VD   = Vres*[17.0 17.0 66.]/100;% (m3) Vol Deep AIP; 50 Ma;    bath = 21
% VD   = Vres*[16.0 16.0 68.]/100;% (m3) Vol Deep AIP; 56 Ma;      bath = 22
% VD   = Vres*[15.0 15.0 70.]/100;% (m3) Vol Deep AIP; 60 Ma;    bath = 23
 VD   = Vres*[13.0 14.0 73.]/100;% (m3) Vol Deep AIP; 67 Ma;    bath = 24




HD   = VD./(fA3*Aoc);           % (m)  H   Deep AIP
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 HD ...
        HH HLI HDT];            % (m) height of boxes
V    = A.*H;  
end


%----------------------------------------%
%----------------------------------------%
%------ END discrete timesteps END ------%
%----------------------------------------%
%----------------------------------------%






% Middle Miocene Bathymetry (15 Ma)
if(bath == 4)
 %--------------------- Ocean Boxes
%        A   I   P 
fA3  = [.23 .18 .49    ];  % Area fraction
fA   = [fA3 fA3 fA3 .10];  % 0.05
A    = fA*Aoc;
HLI  = [100. 900.];         % (m) height L I boxes
HD   = Hav-sum(HLI);
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 ...
        HLID(3)*on3 250.];  % (m) height of boxes
V    = fA.*H*Aoc;           % (m) Volume of boxes
% volume below H box = A(10)*(Hav-H(10))
% add to deep boxes
V(7:9) = V(7:9) + A(10)*(Hav-H(10))/3;
H = V./A;
 
end


if(bath == 2)
%--------------------- Ocean Boxes
%        A   I   P 
if(ftys)
fA3  = [.15 .14 .52]; %03/31/06 % Area fraction AIP
%fA3 = [.17 .18 .46];           % Area fraction AIP
fH   = 0.10;                    % Area fraction H box
fT   = 0.09;                    % Area fraction Tethys
fA   = [fA3 fA3 fA3 fH fT*on3]; %
A    = fA*Aoc;
HLI  = [100. 900.];             % (m) height L I boxes
DTM  = sum(HLI);                % (m) depth thermocline
HH   = 250.;                    % (m) depth H box
%HDT = [1000.];                 % (m) height Deep Tethys
HDT  = [ 200.];                 % (m) height Deep Tethys
Vres = Voc-(DTM*(1-fH)+HH*fH+HDT*fT)*Aoc;
%VD  = Vres*fA3/(sum(fA3));     % (m3) Vol Deep AIP
VD   = Vres*[16.0 16.0 68.]/100;% (m3) Vol Deep AIP
%      04/02/06


HD   = VD./(fA3*Aoc);           % (m)  H   Deep AIP
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 HD ...
        HH HLI HDT];            % (m) height of boxes
V    = A.*H;    
else
%--------------------- Ocean Boxes
%        A   I   P 
fA3  = [.26 .18 .46    ];  % Area fraction
fA   = [fA3 fA3 fA3 .10];  % 0.05
A    = fA*Aoc;
HLI  = [100. 900.];         % (m) height L I boxes
HD   = Hav-sum(HLI);
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 ...
        HLID(3)*on3 250.];  % (m) height of boxes
V    = fA.*H*Aoc;           % (m) Volume of boxes
% volume below H box = A(10)*(Hav-H(10))
% add to deep boxes
V(7:9) = V(7:9) + A(10)*(Hav-H(10))/3;
H = V./A;
end;
end

% 34 Ma Bathymetry
if(bath == 5)
%--------------------- Ocean Boxes
%        A   I   P 
if(ftys)
fA3  = [.19 .17 .51]; %03/31/06 % Area fraction AIP

fH   = 0.10;                    % Area fraction H box
fT   = 0.03;                    % Area fraction Tethys
fA   = [fA3 fA3 fA3 fH fT*on3]; %
A    = fA*Aoc;
HLI  = [100. 900.];             % (m) height L I boxes
DTM  = sum(HLI);                % (m) depth thermocline
HH   = 250.;                    % (m) depth H box
%HDT = [1000.];                 % (m) height Deep Tethys
HDT  = [ 200.];                 % (m) height Deep Tethys
Vres = Voc-(DTM*(1-fH)+HH*fH+HDT*fT)*Aoc;
%VD  = Vres*fA3/(sum(fA3));     % (m3) Vol Deep AIP
VD   = Vres*[23. 19. 58.]/100;% (m3) Vol Deep AIP
%      04/02/06


HD   = VD./(fA3*Aoc);           % (m)  H   Deep AIP
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 HD ...
        HH HLI HDT];            % (m) height of boxes
V    = A.*H;    
end;
end

% 56 Ma Bathymetry
if(bath == 6)
%--------------------- Ocean Boxes
%        A   I   P 
if(ftys)
fA3  = [.15 .14 .52]; %03/31/06 % Area fraction AIP

fH   = 0.10;                    % Area fraction H box
fT   = 0.09;                    % Area fraction Tethys
fA   = [fA3 fA3 fA3 fH fT*on3]; %
A    = fA*Aoc;
HLI  = [100. 900.];             % (m) height L I boxes
DTM  = sum(HLI);                % (m) depth thermocline
HH   = 250.;                    % (m) depth H box
%HDT = [1000.];                 % (m) height Deep Tethys
HDT  = [ 200.];                 % (m) height Deep Tethys
Vres = Voc-(DTM*(1-fH)+HH*fH+HDT*fT)*Aoc;
%VD  = Vres*fA3/(sum(fA3));     % (m3) Vol Deep AIP
VD   = Vres*[16.0 16.0 68.]/100;% (m3) Vol Deep AIP
%      04/02/06


HD   = VD./(fA3*Aoc);           % (m)  H   Deep AIP
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 HD ...
        HH HLI HDT];            % (m) height of boxes
V    = A.*H;    
end;
end

% 67 Ma Bathymetry
if(bath == 7)
%--------------------- Ocean Boxes
%        A   I   P 
if(ftys)
fA3  = [.15 .13 .55]; %03/31/06 % Area fraction AIP

fH   = 0.10;                    % Area fraction H box
fT   = 0.07;                    % Area fraction Tethys
fA   = [fA3 fA3 fA3 fH fT*on3]; %
A    = fA*Aoc;
HLI  = [100. 900.];             % (m) height L I boxes
DTM  = sum(HLI);                % (m) depth thermocline
HH   = 250.;                    % (m) depth H box
%HDT = [1000.];                 % (m) height Deep Tethys
HDT  = [ 1200.];                 % (m) height Deep Tethys
Vres = Voc-(DTM*(1-fH)+HH*fH+HDT*fT)*Aoc;
%VD  = Vres*fA3/(sum(fA3));     % (m3) Vol Deep AIP
VD   = Vres*[13.0 14.0 73.0]/100;% (m3) Vol Deep AIP
%      04/02/06


HD   = VD./(fA3*Aoc);           % (m)  H   Deep AIP
HLID = [HLI HD];
H    = [HLID(1)*on3 HLID(2)*on3 HD ...
        HH HLI HDT];            % (m) height of boxes
V    = A.*H;    
end;
end




%   Number of ocean boxes
Nb  = length(V);           
onV = ones(1,length(V));

% volume of basins
for i=1:3
VO(i) = sum(V([i i+3 i+6]));
end;
VO(4) = V(10);
if (ftys)
VO(5) = sum(V(11:13));
end


% Temperature and feedback
%ntL   = 0.4; % 0.3 Low Lat sensitivity
%ntH   = 0.5; % 0.4 HighLat sensitivity

TC3  = [20. 10. 2.];      % (degC) temp. of boxes
TCv0 = [TC3(1)*on3 TC3(2)*on3 TC3(3)*on3 2.0]+Tempchange;
if(ftys) 
TC3  = [25. 16. 12.];     % (degC) temp. of boxes
TCv0 = [TC3(1)*on3 TC3(2)*on3 TC3(3)*on3 12.0+0];
TCT  = [18. 14. 12.+0];   % 18/25 16/14 12
TCv0 = [TCv0 TCT]+0;
if(fpaleo)
TC3  = [20. 10. 2.];     % (degC) temp. of boxes
TCv0 = [TC3(1)*on3 TC3(2)*on3 TC3(3)*on3 2.0+0];
TCT  = [20. 10. 2.];   % 18/25 16/14 12
TCv0 = [TCv0 TCT]+Tempchange;

% PETM Values:
% TC3  = [25. 16. 12.];     % (degC) temp. of boxes
% TCv0 = [TC3(1)*on3 TC3(2)*on3 TC3(3)*on3 12.0+0];
% TCT  = [18. 14. 12.+0];   % 18/25 16/14 12
% TCv0 = [TCv0 TCT]+0;

end
end;


TCv   = TCv0;
Soc   = 34.72;             % Sal whole ocean 
Sv    = onV*Soc;           % Salinity vector

% Pressure vector. Note: H(k+6) all different
for k=1:3
Hv3(k,:)= [H(k)/2 H(k)+H(k+3)/2 H(k)+H(k+3)+H(k+6)/2]; 
end;
Pv      = [Hv3(:,1)' Hv3(:,2)' Hv3(:,3)' H(10)/2]/10; 
if(ftys) 
k = 11;
HTv     = [H(k)/2 H(k)+H(k+1)/2 H(k)+H(k+1)+H(k+2)/2];
Pv      = [Pv HTv/10];
end;

% Overturning
TH0  = 20.e6*3600*24*365;  % (m3/y) 25 20 Sv conveyor transport 
if(ftys)
TH0  = 25.e6*3600*24*365;  % (m3/y) 25 20 Sv conveyor transport 
if(fpaleo)
TH0  = 20.e6*3600*24*365;  % (m3/y) 25 20 Sv conveyor transport 

% PETM:
%TH0  = 25.e6*3600*24*365;  % (m3/y) 25 20 Sv conveyor transport
end;
end




TT   = 02.e6*3600*24*365;  % (m3/y) 03 02 Sv conveyor transport 
TH   = TH0;
TS   = 0.0;
% TH branches
tA = 0.20;      % 0.20 upwelled into intermdt Atl 0.27 0.15 
tI = 0.20;      % 0.20 upwelled into intermdt Ind 0.29 0.30 


% mixing A   I   P   TLI TII      %     3.5 3.5 8.5
mv0   = [5.5 4.5 6.5 2.5 2.]*1e6; %  Sv 5.5 4.5 6.5 2.5 2
% high-deep
mhd0  = [03. 02. 8.0 1.0]*1e6;    %  Sv 3 2 8 | 4 4 6
if(ftys)
mv0   = [3.5 3.5 7.0 3.2 2.]*1e6; %  Sv 5.5 4.5 8.5 3.0 2
mhd0  = [04. 04. 6.0 0.7]*1e6;    %  Sv 3 2 8 | 4 4 6
if(fpaleo)
mv0   = [5.5 4.5 6.5 2.5 2.]*1e6; %  Sv 5.5 4.5 6.5 2.5 2
mhd0  = [03. 02. 8.0 1.0]*1e6;    %  Sv 3 2 8 | 4 4 6

% PETM:
%mv0   = [3.5 3.5 7.0 3.2 2.]*1e6; %  Sv 5.5 4.5 8.5 3.0 2
%mhd0  = [04. 04. 6.0 0.7]*1e6;    %  Sv 3 2 8 | 4 4 6
end
end;
mv0   = 3.8*mv0 *365*24*3600;     % (m3/y) 3.8 4.0
mhd0  = 1.3*mhd0*365*24*3600;     % (m3/y) 1.3
mv    = mv0;
mhd   = mhd0;

% air-sea CO2/O2
kasv      = NaN*onV;
vask      = NaN*onV;
if(ftys)
kkv       = [1 2 3 10 11];
else 
kkv       = [1 2 3 10];
end;
xkh       = 1.*0.06;    %  Wally's CO2 exchange coeff.
kasv(kkv) = xkh*A(kkv); % (mol/uatm/y) air sea exch coeff Llat                               
pv        = 3.*365;     % (m/day) -> (m/y) piston velocity
vask(kkv) = pv*A(kkv);  %  m3/y

%============== Biological Pump ============%
%
EPH   = 1.8*A(10); % (mol/y) 1.8 1.6 H Export, mol C/m2/y*A = mol/y
rrain = 6.1;       % 6.1 6.2 6.7 export rain ratio (Corg/CaCO3)
                   % 5.9(?) 6.1(2,3)

if(ftys)
rrain = 6.7;       % 6.7 7 4.2 export rain ratio (Corg/CaCO3)
                    % 8.0(1,1) 6.3(1,3) 6.6/6.2/6.0?(2,3) 
if(fpaleo)
 rrain = 6.1; 
 
 % PETM:
% rrain = 6.7;
end
end;               
nu    = 0.31;      % 0.31 water column dissolution
fEPL  = 0.80;      % 0.80 0.9 LL utilization
eI    = 0.78;      % 0.78 0.8 fraction EPL, remineralized in I boxes
if(~fsed)
nu    = 0.;
end;

% fraction EPH, remineralized in deep A,I,P boxes
gp      = 0.*ones(1,Nb);
gp(7:9) = [.3 .3 .4];   % .3 .3 .4
%gp(7:9) = [1   1  1]/3; % .7 .3 0


% malte
% modern ocean pRef and pCSi = 280

%============== silicate weathering: volc degass
pRef  = 280.;        % uatm, weathering  ref  280
pCSi  = 280.;        % uatm, std-stt atm pCO2 280
if(fmalte)
pRef  = 849.;        % uatm, weathering  ref  1000
pCSi  = 849.;        % uatm, std-stt atm pCO2 1000
end

if(ftys)
pRef  = 0500.*1.;    % uatm, weathering ref   500  574 350  750 /2
pCSi  = 1000.*1.;    % uatm, std-stt atm pCO2 560 1000 700 1000 /2
    if(fpaleo)
     pRef  = 280.;        % uatm, weathering  ref  280
     pCSi  = 280.;        % uatm, std-stt atm pCO2 280
        if(fmalte)        % CHANGE here!!!!
        pRef  = 1949./2;        % uatm, weathering  ref  1000
        pCSi  = 1949.;        % uatm, std-stt atm pCO2 1000  
       % pRef  = 848.6/2;        % uatm, weathering  ref  1000
       % pCSi  = 848.6;        % uatm, std-stt atm pCO2 1000   
        end
    end
end; 
FVC   = 1*5.e12/Aoc;  % mol C, degassing /m2/y @280 uatm
nSi   = 0.2;          % 0.2 0.3
FVC   = FVC*(pCSi/pRef)^nSi; % initial 

%============== kerogen oxidation
Fkg   = 1*09.e12/Aoc;        % mol C    /m2/y 09
if(ftys)
Fkg   =   05.e12/Aoc;        % mol C    /m2/y 05

    if(fpaleo)
    Fkg   = 1*09.e12/Aoc;        % mol C    /m2/y 09
    
    % PETM:
    %Fkg   =   05.e12/Aoc;        % mol C    /m2/y 05
    end
end

%============== CaCO3 in-flux ===============%
%
FiN   = 1.0*12.e12/Aoc;      % mol C    /m2/y riverine flux 1.3
nC0   = 0.40;                % 0.4 
nCC   = 0.40;                % 0.4 0.3 1.0 0.5
FiN   = FiN*(pCSi/pRef)^nC0; % River Influx
Fpr   = 3.*FiN;              % mol C    /m2/y production 3.6
% rain of 'remainder'
%frrf  = 1.5*1.15*0.180;     %  g/cm2/ky remainder 1.5*1.15
frrf  = 0.35;                %  g/cm2/ky remainder .311
frrf  = frrf*1e4/1e3/1e3;    % -> kg/ m2/ y
frain = Fpr*m2kg/(Fpr*m2kg+frrf);

%======= Carbon-13
Rst    = 0.011;   % 13C: R standard (value irrelevant)
epsp   = -27.;    % -27 fractionation Corg
if(ftys)
epsp   = -33.;    % -33 fractionation Corg
    if(fpaleo)
    epsp   = -27.;    % -27 fractionation Corg
    end
end;
d13Cin = +3.0;    % d13C of riverine flux 3.0 2.0 2.6
Rin    = Rst*(d13Cin/1e3+1);
FiN13  = Rin*FiN; % mol C    /m2/y riverine flux
% silicate weathering: volc degass
d13Cvc = -3.0;    % d13C -5 +0.3 -0.7 +2.0 
if(ftys)
d13Cvc = -5.0;    % d13C -5 +0.3 -0.7 +2.0 
    if(fpaleo)
    d13Cvc = -3.0;    % d13C -5 +0.3 -0.7 +2.0 
    end
end;
Rvc    = Rst*(d13Cvc/1e3+1);
FVC13  = Rvc*FVC; % mol C    /m2/y
% kerogen oxidation
d13Ckg = -22.3;   % d13C -22.3 -28.3
Rkg    = Rst*(d13Ckg/1e3+1);
Fkg13  = Rkg*Fkg; % mol C    /m2/y


if(fsed) %========= sediments ============== fsed
% dissolution parameter
%
% ~fc*(1-om)
if    (dsflag == 1)
Kd    = 1.*365/100   % 1/d -> 1/y
nc    = 4.5;         % 4.50 calc dissolution order
elseif(dsflag == 2)
% ~fc^0.5*(1-om)
nc    = 2.35;        % 2.35 calc dissolution order
Kd    = 0.12;
elseif(dsflag == 3)
% ~fc^0.5*(cs-c)     % Kd defnd in DEQ
nc    = 2.40;        % 2.40 2.35 calc dissolution order
cst   = 100.e-6;     %
KS    = 20.36e10;    % mol/m2/y
end;

%====== shelf/deep rain
% fsh and nshT sets the fraction of shelf production vs. open ocean
% production (modern ocean = 1)
fsh    = 1.00;            %  STANDARD
nshT   = 1.00;            %  STANDARD

% increase shelf rain:
fsh = 5.02;
nshT = 0.3;



% at PETM the fraction of shelf production to open ocean prod. was higher
if(ftys)
fsh    = 4.5;            % 4.5 ORIGINAL VALUE
nshT   = 0.4;            % 0.4 ORIGINAL VALUE


    if(fpaleo)
    % increase shelf rain:    % Change here:
    fsh    = 2.72;           
    nshT   = 0.3;           
    end
end;

%----------------------------------------%
%----------------------------------------%
%-- discrete 5 Ma bathymetry timesteps --%
%---------- depth profiles --------------%
%----------------------------------------%

% PROCEDURE (example interp. between 0 and 56 Ma at 15 Ma):
% the following is done for each basin (AIP (and T))!
% A0   = asvA at 0  Ma (%)
% A56  = asvA at 56 Ma (%)
% A0m  = A0*atlantic_fraction_at_0*Aoc = A0*0.26*Aoc (m)
% A56m = A56*0.15*Aoc (m)
%
% for i = 1:13                             % length(dsv) = 13
% A15m(i) = interp1([0 56],[A0m(i) A56m(i)],15,'linear')
% A15(i)  = A34m(i)/(sum(A15m))             % 0.23 = fraction Atlantic at 15 Ma
% end

% INDICES
% O  Ma  = bath 11 (modern)
% 5  Ma  = bath 12 (AIP interp. between 0 and 56)
% 10 Ma  = bath 13 (AIP interp. between 0 and 56)
% 15 Ma  = bath 14 (AIP interp. between 0 and 56)
% 20 Ma  = bath 15 (AIP interp. between 0 and 56)
% ----- (23 Ma: Tethys closure) ----- %
% 25 Ma  = bath 16 (AIP interp. between 0 and 56, Tethys same as in 56)
% 30 Ma  = bath 17 (AIP interp. between 0 and 56, Tethys same as in 56)
% ----- (33 Ma: circulation switch NADW --> SO/NADW)
% 35 Ma  = bath 18 (AIP interp. between 0 and 56, Tethys same as in 56)
% 40 Ma  = bath 19 (AIP interp. between 0 and 56, Tethys same as in 56)
% 45 Ma  = bath 20 (AIP interp. between 0 and 56, Tethys same as in 56)
% 50 Ma  = bath 21 (AIP interp. between 0 and 56, Tethys same as in 56)
% 56 Ma  = bath 22 (Bice)
% 60 Ma  = bath 23 (AIP interp. between 56 and 67, Tethys between 56 and 67)
% 67 Ma  = bath 24 (Sewall)

if(bath == 11) % 0 Ma
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
asvA  = [7.0297  5.1729  1.9106  2.3882  4.2988  4.2988  9.6711  ...
         9.6711 16.2389 16.2389 11.1712 11.1712  0.7388]/100;
asvI  = [3.5710  2.6844  1.5907  1.9884  5.0146  5.0146 12.6299 ...
        12.6299 18.3221 18.3221  8.4957  8.4957  1.2407]/100; 
asvP  = [1.6358  2.5901  1.4484  1.8105  3.4372  3.4372 10.9275 ...
        10.9275 17.5411 17.5411 13.4784 13.4784  1.7468]/100;    
end

if(bath == 12) % 5 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m) 
asvA = [0.067144 0.054324 0.022749 0.026014 0.043251 0.042558 0.098114 0.097285 0.158244 0.159805 0.111670 0.111843 0.006992];
asvI = [0.033357 0.028863 0.018687 0.024819 0.049872 0.051224 0.125656 0.126281 0.180605 0.180044 0.084663 0.084397 0.011527];
asvP = [0.014892 0.026143 0.016091 0.018831 0.032361 0.035946 0.108152 0.112352 0.168084 0.177866 0.135020 0.138531 0.015725];
end

if(bath == 13) % 10 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m) 
asvA = [0.063734 0.057132 0.026689 0.028321 0.043536 0.042094 0.099633 0.097906 0.153762 0.157011 0.111626 0.111985 0.006564];
asvI = [0.030907 0.030966 0.021581 0.029959 0.049587 0.052346 0.124986 0.126263 0.177880 0.176736 0.084356 0.083814 0.010612];
asvP = [0.013461 0.026380 0.017662 0.019541 0.030396 0.037485 0.107055 0.115359 0.160924 0.180265 0.135251 0.142193 0.014022];
end

if(bath == 14) % 15 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)    
asvA = [0.060033 0.060179 0.030966 0.030825 0.043846 0.041591 0.101281 0.098580 0.148898 0.153978 0.111578 0.112139 0.006100];
asvI = [0.028353 0.033157 0.024599 0.035316 0.049290 0.053516 0.124288 0.126244 0.175041 0.173289 0.084037 0.083207 0.009658];
asvP = [0.012061 0.026612 0.019198 0.020236 0.028475 0.038989 0.105982 0.118298 0.153925 0.182610 0.135477 0.145773 0.012357];
end

if(bath == 15) % 20 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m) 
asvA = [0.056003 0.063497 0.035623 0.033552 0.044183 0.041042 0.103075 0.099314 0.143600 0.150676 0.111526 0.112307 0.005594];
asvI = [0.025690 0.035443 0.027746 0.040903 0.048980 0.054736 0.123560 0.126224 0.172079 0.169692 0.083704 0.082573 0.008663];
asvP = [0.010693 0.026839 0.020700 0.020914 0.026597 0.040459 0.104934 0.121172 0.147082 0.184903 0.135698 0.149273 0.010729];
end

if(bath == 16) % 25 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)   
asvA = [0.051597 0.067124 0.040714 0.036532 0.044552 0.040442 0.105037 0.100116 0.137810 0.147066 0.111469 0.112491 0.005042];
asvI = [0.022909 0.037830 0.031031 0.046737 0.048657 0.056010 0.122799 0.126203 0.168987 0.165938 0.083357 0.081911 0.007624];
asvP = [0.009354 0.027060 0.022169 0.021578 0.024760 0.041897 0.103909 0.123983 0.140390 0.187146 0.135914 0.152696 0.009137];
asvT = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;  
end

if(bath == 17) % 30 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)    
asvA = [0.046761 0.071106 0.046303 0.039804 0.044956 0.039784 0.107191 0.100997 0.131453 0.143103 0.111406 0.112692 0.004435];
asvI = [0.020003 0.040324 0.034465 0.052834 0.048318 0.057342 0.122005 0.126181 0.165756 0.162014 0.082994 0.081220 0.006538];
asvP = [0.008045 0.027277 0.023605 0.022228 0.022964 0.043304 0.102905 0.126732 0.133843 0.189340 0.136126 0.156045 0.007580];
asvT = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;  
end

if(bath == 18) % 35 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)   
asvA = [0.041429 0.075497 0.052465 0.043412 0.045403 0.039059 0.109565 0.101969 0.124444 0.138734 0.111337 0.112915 0.003766];
asvI = [0.016963 0.042933 0.038057 0.059211 0.047965 0.058735 0.121174 0.126158 0.162375 0.157910 0.082614 0.080497 0.005403];
asvP = [0.006764 0.027489 0.025011 0.022863 0.021206 0.044681 0.101924 0.129423 0.127437 0.191486 0.136332 0.159321 0.006056];
asvT = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;  
end

if(bath == 19) % 40 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)  
asvA = [0.035519 0.080362 0.059295 0.047410 0.045897 0.038254 0.112197 0.103045 0.116676 0.133891 0.111260 0.113161 0.003025];
asvI = [0.013779 0.045665 0.041818 0.065889 0.047595 0.060193 0.120304 0.126134 0.158836 0.153611 0.082216 0.079740 0.004213];
asvP = [0.005511 0.027697 0.026387 0.023485 0.019485 0.046028 0.100963 0.132056 0.121168 0.193587 0.136535 0.162528 0.004565];
asvT = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;  
end

if(bath == 20) % 45 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
asvA = [0.028933 0.085785 0.066906 0.051866 0.046448 0.037358 0.115130 0.104245 0.108019 0.128494 0.111175 0.113436 0.002198];
asvI = [0.010442 0.048529 0.045762 0.072890 0.047206 0.061722 0.119391 0.126109 0.155125 0.149106 0.081799 0.078946 0.002966];
asvP = [0.004283 0.027900 0.027733 0.024093 0.017801 0.047347 0.100023 0.134633 0.115031 0.195643 0.136733 0.165667 0.003105];
asvT = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;  
end

if(bath == 21) % 50 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)   
asvA = [0.021546 0.091867 0.075442 0.056863 0.047066 0.036353 0.118419 0.105590 0.098311 0.122442 0.111079 0.113744 0.001272];
asvI = [0.006939 0.051535 0.049900 0.080238 0.046799 0.063327 0.118434 0.126083 0.151230 0.144377 0.081361 0.078113 0.001658];
asvP = [0.003082 0.028099 0.029052 0.024689 0.016152 0.048638 0.099102 0.137157 0.109022 0.197657 0.136927 0.168741 0.001676];
asvT = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;  
end

if(bath == 22) % 56 Ma
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
asvA  = [1.1407 10.0216  8.7160  6.3724  4.7915  3.4973 12.2935 ...
        10.7438  8.4983 11.4134 11.0948 11.4167 1e-6]/100;
asvI  = [0.2501  5.5345  5.5145  8.9550  4.6283  6.5361 11.7221 ...
        12.6050 14.6295 13.8384  8.0807  7.7057 1e-6]/100; 
asvP  = [0.1673  2.8333  3.0599  2.5389  1.4218  5.0153  9.8023 ...
        14.0117 10.1975 20.0019 13.7155 17.2346 1e-6]/100; 
asvT  = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;   
end

if(bath == 23) % 60 Ma
dsv  = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)  
asvA = [0.007258 0.133449 0.064087 0.053504 0.051560 0.043129 0.110058 0.117849 0.107462 0.119188 0.114726 0.077722 6.363636e-09];
asvI = [0.001633 0.075778 0.043798 0.072831 0.046132 0.062388 0.103513 0.103176 0.119629 0.112961 0.202319 0.055836 6.533337e-09];
asvP = [0.001042 0.037944 0.021983 0.019107 0.013230 0.035436 0.072793 0.108565 0.087615 0.160192 0.313064 0.123068 0.005955];
asvT = [0.048831 0.457820 0.166704 0.057408 0.021598 0.036698 0.030085 0.016683 0.018780 0.029225 0.113026 0.003138 6.923079e-09];
end

if(bath == 24) % 67 Ma
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
asvA  = [3.e-15 19.1609 2.3711 3.5620 5.7941 5.7403 8.7525 13.6071 14.6801 12.8035 12.1339 1.3945 0]/100;
asvI  = [8.e-16 11.4288 2.2414 4.1322 4.5850 5.6785 7.7681 6.0069 6.9374 6.5049 43.1324 1.5844 0]/100; 
asvP  = [7.e-16 5.3848 0.7729 0.8714 1.1596 1.1087 3.1049 5.6362 6.3856 9.4298 60.4115 4.1537 1.5810]/100; 
asvT  = [2.e-15 44.0848 3.7635 2.5700 1.5607 2.9391 3.6885 3.5025 3.8955 2.5679 31.4274 0 0]/100;    
end

%----------------------------------------%
%----------------------------------------%
%------ END discrete timesteps END ------%
%----------------------------------------%
%----------------------------------------%










%----------------- sediment boxes (bathymetry) ------------%
% use bath == 2 for modern or PETM
% use bath == 4 for 15 Ma
% use bath == 5 for 34 Ma
% use bath == 6 for 56 Ma
% use bath == 7 for 67 Ma

if(bath == 4)
% % 15 Ma Bathymetry (Herold) (DO NOT USE!)
% dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% % area fraction A I P
% asvA  = [1.e-14 13.6308 1.3753 2.5557 1.7192 1.6522 2.4777 5.5467 8.1215 13.4460 19.3860 30.0889 0]/100;
% asvI  = [8.e-15 8.4714 0.9512 1.9060 1.6199 1.0306 2.2203 6.0971 9.4091 13.8967 26.9934 27.4043 0]/100; 
% asvP  = [5.e-15 6.5732 1.0424 2.5070 1.2601 1.3069 3.5827 7.1425 11.2253 14.0526 16.4338 34.8735 0]/100; 
 
% 15 Ma Bathymetry (Modern (0Ma)-Bice(56Ma) Interpolation)
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% area fraction A I P
asvA  = [0.0600 0.0601 0.0309 0.0308 0.0438 0.0415 0.1012 0.0985 0.1488 0.1539 0.1115 0.1121 0.0060];
asvI  = [0.0283 0.0331 0.0246 0.0353 0.0492 0.0535 0.1242 0.1262 0.1750 0.1732 0.0840 0.0832 0.0096]; 
asvP  = [0.0120 0.0266 0.0192 0.0202 0.0284 0.0389 0.1059 0.1183 0.1539 0.1826 0.1354 0.1457 0.0123];
end

% old
if    (bath == 1)
dsv   = [ .1 .6 1.5 2.5 3.5 4.5 5.5 6.5]*1000;
asvA  = [ 7.0297  5.1729  4.2988  8.5975 19.3421 32.4777 22.3425  0.7388]/100.;
asvI  = [ 3.5710  2.6844  3.5792 10.0293 25.2598 36.6442 16.9915  1.2407]/100.; 
asvP  = [ 1.6358  2.5901  3.2590  6.8744 21.8550 35.0822 26.9567  1.7468]/100.;
if(ftys)
asvT  = [16.3934 16.3934 16.3934 20.4918 20.4918  8.1967  1.6393  0.0001]/100.;
end;
end



% modern and petm
if(bath == 2)
if(ftys)  % 03/31/06, 2x2\deg Bice  
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% area fraction A I P
asvA  = [1.1407 10.0216  8.7160  6.3724  4.7915  3.4973 12.2935 ...
        10.7438  8.4983 11.4134 11.0948 11.4167 1e-6]/100;
asvI  = [0.2501  5.5345  5.5145  8.9550  4.6283  6.5361 11.7221 ...
        12.6050 14.6295 13.8384  8.0807  7.7057 1e-6]/100; 
asvP  = [0.1673  2.8333  3.0599  2.5389  1.4218  5.0153  9.8023 ...
        14.0117 10.1975 20.0019 13.7155 17.2346 1e-6]/100; 
asvT  = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;
else %#! old 03/22/06
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% area fraction A I P
asvA  = [7.0297  5.1729  1.9106  2.3882  4.2988  4.2988  9.6711  ...
         9.6711 16.2389 16.2389 11.1712 11.1712  0.7388]/100;
asvI  = [3.5710  2.6844  1.5907  1.9884  5.0146  5.0146 12.6299 ...
        12.6299 18.3221 18.3221  8.4957  8.4957  1.2407]/100; 
asvP  = [1.6358  2.5901  1.4484  1.8105  3.4372  3.4372 10.9275 ...
        10.9275 17.5411 17.5411 13.4784 13.4784  1.7468]/100; 
%if(ftys)
%asvT  = [16.3934 16.3934  7.2860  9.1075 10.2459 10.2459 10.2459 ...
%         10.2459  4.0984  4.0984  0.8197  0.8197  0.0001]/100;   
%end;
end; % ftys old/new
end;


% OLD ONE: interpolated between Herold (15Ma) and Bice (56Ma) (Do not
% use!!!)
% % 34 Ma Bathymetry 
% if(bath == 5)
% dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% % area fraction A I P
% asvA  = [0.0041 0.1231 0.0404 0.0394 0.0283 0.0232 0.0604 0.0743 0.0825 0.1270 0.1637 0.2329 0];
% asvI  = [0.0009 0.0733 0.0272 0.0464 0.0278 0.0316 0.0591 0.0862 0.1143 0.1387 0.1964 0.1975 3.e-9]; 
% asvP  = [0.0008 0.0477 0.0201 0.0252 0.0133 0.0308 0.0657 0.1044 0.1073 0.1691 0.1512 0.2639 0]; 
% asvT  = [0.0705 0.4653 0.2240 0.0715 0.0242 0.0399 0.0270 0.0085 0.0098 0.0308 0.0235 0.0045 9.e-9;]; 
% end


% NEW ONE: interpolated between modern and Bice
% 34 Ma Bathymetry 
if(bath == 5)
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% area fraction A I P
asvA  = [0.0425 0.0746 0.0512 0.0426 0.0453 0.0392 0.1090 0.1017 0.1258 0.1396 0.1113 0.1128 0.0039];
asvI  = [0.0175 0.0424 0.0373 0.0579 0.0480 0.0584 0.1213 0.1261 0.1630 0.1587 0.0826 0.0806 0.0056]; 
asvP  = [0.0070 0.0274 0.0247 0.0227 0.0215 0.0444 0.1021 0.1288 0.1287 0.1910 136292 0.1586 0.0063]; 
asvT  = [0.0705 0.4653 0.2240 0.0715 0.0242 0.0399 0.0270 0.0085 0.0098 0.0308 0.0235 0.0045 9.e-9;]; 
end

% 56 Ma Bathymetry (same as in 2 (Bice))
if(bath == 6)
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% area fraction A I P
asvA  = [1.1407 10.0216  8.7160  6.3724  4.7915  3.4973 12.2935 ...
        10.7438  8.4983 11.4134 11.0948 11.4167 1e-6]/100;
asvI  = [0.2501  5.5345  5.5145  8.9550  4.6283  6.5361 11.7221 ...
        12.6050 14.6295 13.8384  8.0807  7.7057 1e-6]/100; 
asvP  = [0.1673  2.8333  3.0599  2.5389  1.4218  5.0153  9.8023 ...
        14.0117 10.1975 20.0019 13.7155 17.2346 1e-6]/100; 
asvT  = [7.0534 46.5363 22.4068  7.1501  2.4261  3.9946  2.7063 ...
         0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100; 
end

% 67 Ma Bathymetry (Sewall)
if(bath == 7)
dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
% area fraction A I P
asvA  = [3.e-15 19.1609 2.3711 3.5620 5.7941 5.7403 8.7525 13.6071 14.6801 12.8035 12.1339 1.3945 0]/100;
asvI  = [8.e-16 11.4288 2.2414 4.1322 4.5850 5.6785 7.7681 6.0069 6.9374 6.5049 43.1324 1.5844 0]/100; 
asvP  = [7.e-16 5.3848 0.7729 0.8714 1.1596 1.1087 3.1049 5.6362 6.3856 9.4298 60.4115 4.1537 1.5810]/100; 
asvT  = [2.e-15 44.0848 3.7635 2.5700 1.5607 2.9391 3.6885 3.5025 3.8955 2.5679 31.4274 0 0]/100;    
end



% old
if(bath == 3)
dsv   = [.1 .6 1 1.5 2 2.5 2.75 3 3.25 3.5 3.75 4 ...
        4.25 4.5 4.75 5. 5.25 5.5 6.0 6.5]*1000; % depth (m)
% area fraction A I P
asvA  = [07.0297  5.1729  1.9106  2.3882  4.2988  4.2988  4.8355  4.8355  4.8355  4.8355  8.1194  8.1194  8.1194  8.1194  5.5856  5.5856  5.5856  5.5856  0.3694  0.3694]/100;
asvI  = [03.5710  2.6844  1.5907  1.9884  5.0146  5.0146  6.3149  6.3149  6.3149  6.3149  9.1610  9.1610  9.1610  9.1610  4.2479  4.2479  4.2479  4.2479  0.6204  0.6204]/100; 
asvP  = [01.6358  2.5901  1.4484  1.8105  3.4372  3.4372  5.4638  5.4638  5.4638  5.4638  8.7705  8.7705  8.7705  8.7705  6.7392  6.7392  6.7392  6.7392  0.8734  0.8734]/100; 
if(ftys)
asvT  = [16.3934 16.3934  7.2860  9.1075 10.2459 10.2459  5.1229  5.1229  5.1229  5.1229  2.0492  2.0492  2.0492  2.0492  0.4098  0.4098  0.4098  0.4098  0.0001  0.0001]/100;
end;
end








% Number of sediment boxes
Ns    = length(dsv);   
onNs  = ones(1,Ns);

% kLID: assign sediment to ocean boxes
%       Low, Interm, or Deep
kl = find(dsv <= HLI(1));
ki = find(dsv >  HLI(1) & dsv <= (HLI(1)+HLI(2)));
kd = find(                dsv >  (HLI(1)+HLI(2)));
klid(kl) = 01; klid(ki) = 04; klid(kd) = 07;
if(ftys)
kliT(kl) = 11; kliT(ki) = 12; kliT(kd) = 13;
end
nli      = [length(kl) (length(kl)+length(ki))];

 zv   = [000.:10:6000.]; % z, continuous (1 or 10 m)
lzv   = length(zv);


% calculate calcite saturation at depth of sediment boxes

satflg = 2; % 1/2, 1:Wally 2:Millero

zsatv   = dsv;
if     (satflg == 1)
as      = 0.189/1.e3;
zs0     = 3.82e3;     % 3.82 m
co3s0   = 88.7e-6;    % 88.7 mol/kg
co3satv = co3s0*exp(as*(zsatv-zs0));
elseif(satflg == 2)
% 2. Millero
Cam     = 10.3e-3;    % 10.3 (mol/kg) modern
Mgm     = 53.0e-3;    % 53.0 (mol/kg)
Ca      = Cam;    
Mg      = Mgm; 


% CHANGE HERE:
Ca = 19.34e-3; % 10.3 (mol/kg) Calcium 
Mg = 31.35e-3; % 53.0 (mol/kg) Magnesium 





if(ftys)
Ca      = 20.0e-3;    % 10.3 20.0
Mg      = 30.0e-3;    % 53.0 30.0 
    if(fpaleo)    
    % CHANGE HERE:  
   % Ca = Cam; % 10.3 (mol/kg) Calcium 
   % Mg = Mgm; % 53.0 (mol/kg) Magnesium 
    Ca = 19.34e-3; % 10.3 (mol/kg) Calcium 
    Mg = 31.35e-3; % 53.0 (mol/kg) Magnesium 
    end
end;
% calculate ratio
Mgca_ratio = Mg/Ca



% Temp for co3sat corrected 07/13/06
Tdv = TCv(klid);
Sdv =  Sv(klid);
for k=1:Ns
[kspc(k),x] = ...
     kspfun(Tdv(k),Sdv(k),zsatv(k)/10.,Ca,Mg);
end;
co3satv  = kspc/Ca;

end; % satflag



%-------------- Porosity --------------------%
%phic    = 0.78;      % porosity 0.75 0.78
if(phic) % do NOT use exist! phic does exist (global!)
phiiA   = ones(1,Ns)*phic;
phiiI   = ones(1,Ns)*phic;
phiiP   = ones(1,Ns)*phic;
if(ftys)
phiiT   = ones(1,Ns)*phic;
end;
else % phic
phi0    = 0.85;      % porosity max  0.85 0.88
gam     = 0.23;      %               0.23 0.28
phi1    = phi0-gam;  % porosity min
end;

hs   = 0.08;         % (m ) bioturbated layer 0.08 0.1
VsvA = asvA*A(1)*hs; % (m3) Volume sediment boxes A
VsvI = asvI*A(2)*hs; % (m3) Volume sediment boxes I
VsvP = asvP*A(3)*hs; % (m3) Volume sediment boxes P
VsA  = sum(VsvA);
VsI  = sum(VsvI);
VsP  = sum(VsvP);
if(ftys)
VsvT = asvT*A(11)*hs;% (m3) Volume sediment boxes T
VsT  = sum(VsvT);
end;

% set initial calcite fraction
fc0A = 0.46*ones(1,Ns);
fc0I = 0.46*ones(1,Ns);
fc0P = 0.46*ones(1,Ns);
if(ftys)
fc0T = 0.46*ones(1,Ns);
end;
% calc initial phi
if(isempty(phic)) % phi = phi(fc)
FF    = (phi1-phi0)/(1-phi1);
phiiA = (phi0+FF*fc0A)./(1+FF*fc0A); 
phiiI = (phi0+FF*fc0I)./(1+FF*fc0I); 
phiiP = (phi0+FF*fc0P)./(1+FF*fc0P); 
if(ftys)
phiiT = (phi0+FF*fc0T)./(1+FF*fc0T); 
end;
end;
% calc initial calcite mass
mc0vA  = (fc0A.*rhos.*(1-phiiA));
mc0vI  = (fc0I.*rhos.*(1-phiiI));
mc0vP  = (fc0P.*rhos.*(1-phiiP));
if(ftys)
mc0vT  = (fc0T.*rhos.*(1-phiiT));
end;    


%====== Carbon-13
m13c0vA = Rin*mc0vA;     % -> kg CaCO3/m3 Rin
m13c0vI = Rin*mc0vI;     % -> kg CaCO3/m3 Rin
m13c0vP = Rin*mc0vP;     % -> kg CaCO3/m3 Rin
f13c0A  = fc0A.*m13c0vA./mc0vA;
f13c0I  = fc0I.*m13c0vI./mc0vI;
f13c0P  = fc0P.*m13c0vP./mc0vP;
if(ftys)
m13c0vT = Rin*mc0vT;     % -> kg CaCO3/m3 Rin
f13c0T  = fc0T.*m13c0vT./mc0vT;
end;

end; %============================== fsed END

%===============================================%
%
% Define initial conditions
%
%===============================================%


c0   = [2.30*onV];        % mmol/kg (DIC at t=0)
%c0   = [1:1:10]*.32;
a0   = [2.40*onV];        % mmol/kg (ALK at t=0)
%a0   = [1:1:10]*.33;
p0   = [2.50*onV]*1e-3;   % mmol/kg (PO4 at t=0)
p0   = p0*0.87;
%p0   = [1:1:10]*.25*1e-3;

if(fdox)
dox0 = [0.20*onV];         % mol/m3 (O2 t=0)
%surface 
for k=kkv
[x,x,x,x,x,dox0(k)] = ...
   dafunPE(1,1,TCv(k),Sv(k),Pv(k),1,1);
end;
end;

%====== Carbon-13
d13c0 = [2.35*on3 0.5*on3 0.67*on3 1.63];
if(ftys)
d13c0 = [[d13c0] 2.35 0.5 0.67];
end;    
R0    = (d13c0/1e3+1.)*Rst;
cc0   = R0.*c0;

c0    = c0  *1e-3.*rho; % mol/m3
a0    = a0  *1e-3.*rho; % mol/m3
p0    = p0  *1e-3.*rho; % mol/m3
cc0   = cc0 *1e-3.*rho; % mol/m3

% Modern Ocean pCO2 initial cond.
% 280 ppm
C0 = 280*2.2e15/12/Aoc; % 280 
						% (mol/m2) atm. CO2 inventory / m2
                        % 1 ppmv = 2.2 Pg C
%====== Carbon-13
d13C0 = -6.45;
CC0   = C0.*(d13C0/1e3+1.)*Rst;
                        
%Y0 = [c0 a0 p0 C0 mc0vA]';
%Y0 = [c0 a0 p0 C0 mc0vA mc0vI mc0vP]';
%      Nb 2Nb 3Nb 4Nb +1  +2 +2+Ns +2+2Ns +2+3Ns
%      10  20  30  40 41  42    58     74     90
%Y0 = [c0  a0  p0 cc0 C0 CC0 mc0vA  mc0vI  mc0vP]';

if(fsed) %<<<<<<<<<<<<<<<<<<<<<<< fsed
if(fdox) % fdox
%     Nb 2Nb 3Nb 4Nb  5Nb +1  +2 +2+Ns  +2+2Ns   +2+3Ns
%     10  20  30   40  50 51  52    65      78       91
Y0 = [c0  a0  p0 dox0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                           f13c0A  f13c0I   f13c0P]';
%                          +2+4Ns  +2+5Ns   +2+6Ns 
%                             104    117       130
else % fdox 
%     Nb 2Nb 3Nb 4Nb +1  +2 +2+Ns  +2+2Ns   +2+3Ns
%     10  20  30  40 41  42    55      68       81
Y0 = [c0  a0  p0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                           f13c0A  f13c0I   f13c0P]';
%                          +2+4Ns  +2+5Ns   +2+6Ns 
%                              94     107      120
end; % fdox
if(ftys) % ftys
if(fdox) % fdox
Y0 = [c0  a0  p0 dox0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                                f13c0A  f13c0I   f13c0P ...
                                  fc0T  f13c0T]';
else % fdox
Y0 = [c0  a0  p0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                           f13c0A  f13c0I   f13c0P ...
                             fc0T  f13c0T]';
end; % fdox
end; % ftys   
else
if(fdox) % fdox
Y0 = [c0 a0 p0 dox0 cc0 C0 CC0]';
else % fdox
Y0 = [c0 a0 p0 cc0 C0 CC0]';
end; % fdox
end; %<<<<<<<<<<<<<<<<<<<<<<< fsed

% Number of ocean tracer (not atm, not sediment)
if(fdox)
 NO = 5; % c a p ox cc
else
 NO = 4; % c a p    cc
end;


% set integration time                      

t0     =   0e5; % (y) time (start)
if    (BlFlag == 1)
tfinal =   20000; % (y) time (end) 2e5
%tfinal = 1e7;
elseif(BlFlag == 2)
tfinal =   2e5; % (y) time (end) 2e5 7e5
else 
tfinal =   1e8; % (y) time (end) 1e7 --> used for steady state calibr.run
end;


%==== Anthropogenic CO2 ======================%
if(ffflag)
t0     =   1700.; % 1700.
tfinal =   20000.; % 2100 3000 10000

%load 'dat\co2PEm.tex';
dirstr = '5000_0500.dat';
co2PEm = load(['dat\Emss\EmssScen\' dirstr]);

tem   = co2PEm(:,1);
em    = co2PEm(:,2);
% deep sea temp
Dtst = 0.0;
k1st =  -1;
TCvt(1,:) = TCv0;
end;

Dt = tfinal-t0;

%===============================================%
%
% Alternatively, load initial conditions 
%
%===============================================%


if (loadf == 1)
%     load YPE10Sed.DAT;
%     load dat\Modern\NoSed0\YPE10Sed.DAT;
  if    (bath == 1)
%     load dat\Modern\B1D1BL0\YPE10Sed;
%     load dat\B1D1BL0\YPE10Sed;
%     load dat\B1D3BL0\YPE10Sed;
  elseif(bath == 2)
% load ~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/dat/Modern/B2D3BL0/YPE10Sed.DAT;  % STANDARD!
%load ~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/output/results/2ndsubmission/no_bathy_constant_mgca/highpco2/steadystates/Y0_67Ma.DAT
 load ~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/output/results/2ndsubmission/no_bathy_var_mgca/geocarb/steadystates/Y0_67Ma.DAT
   %  load ~/Desktop/malte_ubuntu/myLoscar/steady_states/june2_all_change_but_bathy/Y0_15Ma_mgca_var_lowandhigh_pco2.DAT
    % load D:\malte\myLoscar\output\runs_feb27_all_change_but_bathy_and_mgca\steadystates\highpco2\Y0_67.dat
%load D:\malte\myLoscar\output\runs_feb_25_all_change_but_bathymetrie\steadystates\lowpco2\Y0_10.dat
      %load  D:\malte\myLoscar\steady_states\higher_ppm\2000ppm_ocean\Y0_pRef2000.DAT     % steady state ocean with pRef = xxxx ppm
                                                                                % don't forget to
                                                                                % change also pRef and PCSi 
                                                                                
      % load D:\malte\myLoscar\steady_states\lower_ccd\Y0_fsh5dot8.dat              % change also fsh and nshT !!!
      %load D:\malte\myLoscar\steady_states\MgCa\Y0_mgca_1_6.dat                                    % CaMg
                                                                                
%     load dat\Modern\CO2XLS\YPE10Sed.DAT;     
%     load dat\B2D1BL0\YPE10Sed;
%%      load dat\B2D3BL0\YPE10Sed.DAT;
%     load dat\B2D3BL0\Lpco2\YPE10Sed.DAT;
%     load dat\B2D3BL0\Hpco2\YPE10Sed.DAT;
%load ~/Desktop/malte_ubuntu/myLoscar/dat/PETM/YPE10Sed.DAT
  elseif(bath == 3)
%     load dat\B3D3BL0\YPE10Sed;

  elseif(bath == 4)
     % load ~/Desktop/malte_ubuntu/myLoscar/steady_states/other_bathymetry/26may_only_bath_changed_nadw/Y0_15Ma_nadw.DAT;
     % load ~/Desktop/malte_ubuntu/myLoscar/steady_states/all_change_incl_bathy_28may/Y0_15Ma_allchange.DAT
      load ~/Desktop/malte_ubuntu/myLoscar/output/run5june_all_change_incl_bathy/steady_states/Y0_15Ma.DAT   
   elseif(bath >= 5) & (bath <= 10)
      %load ~/Desktop/malte_ubuntu/myLoscar/dat/PETM/YPE10Sed.DAT    % PETM steady state
      %load ~/Desktop/malte_ubuntu/myLoscar/steady_states/other_bathymetry/13may_only_bath_changed_nadw/Y0_67Ma.DAT
     % load ~/Desktop/malte_ubuntu/myLoscar/steady_states/other_bathymetry/18may_only_bath_changend_nadw/Y0_34Ma.DAT
    %  load ~/Desktop/malte_ubuntu/myLoscar/steady_states/other_bathymetry/18may_only_bath_changed_so/Y0_67Ma_so.DAT
      %load ~/Desktop/malte_ubuntu/myLoscar/steady_states/other_bathymetry/26may_only_bath_changed_nadw/Y0_67Ma_nadw.DAT
      %load ~/Desktop/malte_ubuntu/myLoscar/steady_states/other_bathymetry/26may_only_bath_changed_so/Y0_34Ma_so.DAT
     % load ~/Desktop/malte_ubuntu/myLoscar/steady_states/other_bathymetry/26may_only_bath_changed_nadw/Y0_34Ma_nadw.DAT
      %load ~/Desktop/malte_ubuntu/myLoscar/output/run5june_all_change_incl_bathy/steady_states/Y0_56Ma_low_so.DAT
      % load ~/Desktop/Y067.DAT
      load ~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/output/results/2ndsubmission/bathy/geocarb/steadystates/Y0_67Ma.DAT  

  elseif(11 <= bath) & (bath <= 24)
    load ~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/output/results/2ndsubmission/bathy/highpco2/steadystates/Y0_67Ma.DAT  
   %load ~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/steady_states/other_bathymetry/26may_only_bath_changed_nadw/Y0_34Ma_nadw.DAT
end
    
   % Y0 = YPE10Sed;
  Y0 = Y0_67Ma;
  %Y0 = Y0_34Ma_nadw;
  %Y0 = Y0_162_plus0;
   %Y0 = Y0_280_plus5;
   % Y0 = Y0_67Ma;
  %Y0 = Y0_34Ma_nadw;
   %Y0 = Y0_15Ma;
   %Y0 = Y0_67;
   %Y0 = YPE10Sed;
   % Y0 = Y0_fsh5dot8
   %Y0 = Y0_mgca_1_6;
  % Y0 = Y0_15Ma_mgca_var_lowandhigh_pco2;
   %Y0 = Y0_pRef2000;
   
    % Cbl = Amount of carbon in the Blast
   CBl = 3000.e15; % Blast       2200  2000  3000 n2500
d13CBl = -50.;     % Blast d13C -55   -60   -34   n-50
   RBl = Rst*(d13CBl/1.e3+1.);
 C13Bl = RBl*CBl;    
    kb = 07;       % 07
   kbb = kb+3*Nb;
   
   % !!!!!!!!!!!!!
   % not used:
   % here the carbonblast gets added to the initial cond.
if(BlFlag == 4)  
  Y0(5*Nb+1)  = Y0(5*Nb+1) +  CBl/12/Aoc;   % add X Gt C atm:/Aoc (Aoc = AreaOcean)
  Y0(5*Nb+2)  = Y0(5*Nb+2) +C13Bl/12/Aoc;   % add X Gt C atm:/Aoc
%   Y0(kb)  = Y0(kb) +  CBl/12/V(kb); % add X Gt C ocn:/V2
%   Y0(kbb) = Y0(kbb)+C13Bl/12/V(kb); % add X Gt C ocn:/V2
end;


% adjust PO4
%Y0(2*Nb+1:3*Nb) = Y0(2*Nb+1:3*Nb)*0.87;
%Y0(2*Nb+1:3*Nb) = Y0(2*Nb+1:3*Nb)*1.2;
%Y0(2*Nb+1:3*Nb) = Y0(2*Nb+1:3*Nb)*1.2;
% set initial Oxygen
%Y0(3*Nb+1:4*Nb) = dox0;
% adjust TCO2
%Y0(   1:  Nb) = Y0(   1:  Nb)/1.03;
%Y0(Nb+1:2*Nb) = Y0(Nb+1:2*Nb)/1.03;
end; % loadf


% initial inventory
Mci  = sum(Y0(     1:   Nb).*V')/Voc ... % TC inventory/V
          +Y0(NO*Nb+1     ) *Aoc/Voc;    % atmosphere					
Mai  = sum(Y0(   Nb+1:2*Nb).*V')/Voc;    % TA inventory/V
Mpi  = sum(Y0( 2*Nb+1:3*Nb).*V')/Voc;    % P  inventory/V
Mcci = sum(Y0( 3*Nb+1:4*Nb).*V')/Voc ... % T13C inventory/V
          +Y0(NO*Nb+2     ) *Aoc/Voc;    % atmosphere					

      
% initial total C inventory ocean+atm (Pg C)
Mci*Voc*12/1e15;      
      
if(fsed) % fsed   
if(fdox) % fdox
fc0A   =  Y0(5*Nb+3     :5*Nb+2+1*Ns);    
fc0I   =  Y0(5*Nb+3+1*Ns:5*Nb+2+2*Ns);    
fc0P   =  Y0(5*Nb+3+2*Ns:5*Nb+2+3*Ns);    
f13c0A =  Y0(5*Nb+3+3*Ns:5*Nb+2+4*Ns);    
f13c0I =  Y0(5*Nb+3+4*Ns:5*Nb+2+5*Ns);    
f13c0P =  Y0(5*Nb+3+5*Ns:5*Nb+2+6*Ns); 
if(ftys)
fc0T   =  Y0(5*Nb+3+6*Ns:5*Nb+2+7*Ns);    
f13c0T =  Y0(5*Nb+3+7*Ns:5*Nb+2+8*Ns); 
end   
else % fdox
fc0A   =  Y0(4*Nb+3     :4*Nb+2+1*Ns)    
fc0I   =  Y0(4*Nb+3+1*Ns:4*Nb+2+2*Ns)    
fc0P   =  Y0(4*Nb+3+2*Ns:4*Nb+2+3*Ns);    
f13c0A =  Y0(4*Nb+3+3*Ns:4*Nb+2+4*Ns);    
f13c0I =  Y0(4*Nb+3+4*Ns:4*Nb+2+5*Ns);    
f13c0P =  Y0(4*Nb+3+5*Ns:4*Nb+2+6*Ns); 
if(ftys)
fc0T   =  Y0(4*Nb+3+6*Ns:4*Nb+2+7*Ns);    
f13c0T =  Y0(4*Nb+3+7*Ns:4*Nb+2+8*Ns); 
end    
end % fdox   
% calc initial phi
if(isempty(phic)) % phi = phi(fc)
FF    = (phi1-phi0)/(1-phi1);
phiiA = (phi0+FF*fc0A)./(1+FF*fc0A); 
phiiI = (phi0+FF*fc0I)./(1+FF*fc0I); 
phiiP = (phi0+FF*fc0P)./(1+FF*fc0P); 
if(ftys) % ftys
phiiT = (phi0+FF*fc0T)./(1+FF*fc0T); 
end; % ftys
end; % phi = phi(fc)
end; % fsed


% get derivatives at t0
dYflag = 0;
dY0dt  = LoscarDif(t0,Y0);

% output format
format compact;



%++++++++++++++++++++++ SOLVER ++++++++++++++++++++++++++%
%
% solve differential equations.
% use matlab routine ode (Runge-Kutta)
% with function '....dif' containing the difeq.
%
% set options, 'RelTol' (1e-3 by default) 
%              'AbsTol' (all components 1e-6 by default).
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%options = odeset('RelTol',1e-3,'AbsTol',1e-4);
options = odeset('RelTol',1e-3,'AbsTol',1e-3);
%options = odeset('RelTol',1e-5,'AbsTol',1e-5);
if(ffflag) % 0.05*Dt pH-Contour: -3,-4
options = odeset('RelTol',1e-3,'AbsTol',1e-4,'Maxstep',0.05*Dt);
%options = odeset('RelTol',1e-4,'AbsTol',1e-4,'Maxstep',0.05*Dt);
end;


tic % start clock
[tv,Y] = myode15s('LoscarDif',[t0 tfinal],Y0',options); % 23t 15s
toc % stop  clock
lt = length(tv);


% save solution Y(end)
savf = 0;
if (savf == 1)
	YPE10Sed = Y(lt,:)';
    
    save                    YPE10Sed.DAT  YPE10Sed -ASCII -DOUBLE -TABS;
        
   %save dat\Modern\B2D3BL0\YPE10Sed.DAT  YPE10Sed -ASCII -DOUBLE -TABS;
   %save dat\Modern\NoSed0\YPE10Sed.DAT   YPE10Sed -ASCII -DOUBLE -TABS;
   %save dat\B2D3BL0\YPE10Sed.DAT YPE10Sed -ASCII -DOUBLE -TABS; % <------
   %save dat\B2D3BL0\Lpco2\YPE10Sed.DAT  YPE10Sed -ASCII -DOUBLE -TABS;
   %save dat\B2D3BL0\Hpco2\YPE10Sed.DAT  YPE10Sed -ASCII -DOUBLE -TABS;
   %save dat\Modern\CO2XLS\YPE10Sed.DAT YPE10Sed -ASCII -DOUBLE -TABS;
end;

axx = [t0    tfinal];
if(ftys)
axx = [-0.5e5 tfinal];
    if(fpaleo)
    axx = [t0    tfinal];
    end
end;
end; %%%%%%% solflag


% call dif again at all t to get dYdt(t)
% plus other vars & fluxes
% ############### check if temp changes !!!!!!!!
dYflag = 1;
if(dYflag == 1)
it   = 1;
dYdt = ones(size(Y'));
for i=1:lt
dYdt(:,i) = LoscarDif(tv(i),Y(i,:)');
end;
end;

%===============================================%
% Rename variables (solution was stored in Y)
%===============================================%


if(fsed)   
% values sed-boxes
for l=1:Ns
fcA  (:,l) = Y( :,l+NO*Nb+2     );   %    
fcI  (:,l) = Y( :,l+NO*Nb+2+  Ns);   %    
fcP  (:,l) = Y( :,l+NO*Nb+2+2*Ns);   %    
mcalA(:,l) = Y( :,l+NO*Nb+2     ).*(rhos*(1-phivtA(:, l)));
mcalI(:,l) = Y( :,l+NO*Nb+2+  Ns).*(rhos*(1-phivtI(:, l)));
mcalP(:,l) = Y( :,l+NO*Nb+2+2*Ns).*(rhos*(1-phivtP(:, l)));
McalvA( l) = Y(lt,l+NO*Nb+2     ).*(rhos*(1-phivtA(lt,l)));
McalvI( l) = Y(lt,l+NO*Nb+2+  Ns).*(rhos*(1-phivtI(lt,l)));
McalvP( l) = Y(lt,l+NO*Nb+2+2*Ns).*(rhos*(1-phivtP(lt,l)));
%====== Carbon 13
f13cA  (:,l) = Y( :,l+NO*Nb+2+3*Ns);  %    
f13cI  (:,l) = Y( :,l+NO*Nb+2+4*Ns);  %    
f13cP  (:,l) = Y( :,l+NO*Nb+2+5*Ns);  %    
m13calA(:,l) = mcalA( :,l).*f13cA( :,l)./fcA( :,l);
m13calI(:,l) = mcalI( :,l).*f13cI( :,l)./fcI( :,l);
m13calP(:,l) = mcalP( :,l).*f13cP( :,l)./fcP( :,l);
M13calvA( l) = mcalA(lt,l).*f13cA(lt,l)./fcA(lt,l);
M13calvI( l) = mcalI(lt,l).*f13cI(lt,l)./fcI(lt,l);
M13calvP( l) = mcalP(lt,l).*f13cP(lt,l)./fcP(lt,l);
if(ftys)
fcT    (:,l) = Y( :,l+NO*Nb+2+6*Ns);  %    
f13cT  (:,l) = Y( :,l+NO*Nb+2+7*Ns);  %    
mcalT  (:,l) = Y( :,l+NO*Nb+2+6*Ns).*(rhos*(1-phivtT(:, l)));
McalvT (  l) = Y(lt,l+NO*Nb+2+6*Ns).*(rhos*(1-phivtT(lt,l)));
m13calT(:,l) = mcalT( :,l).*f13cT( :,l)./fcT( :,l);
M13calvT( l) = mcalT(lt,l).*f13cT(lt,l)./fcT(lt,l);
end;    
end;

% final Mcal and fc
McalA = sum(McalvA.*VsvA)/VsA;
McalI = sum(McalvI.*VsvI)/VsI;
McalP = sum(McalvP.*VsvP)/VsP;
fcfA  = fcA(lt,:);
fcfI  = fcI(lt,:);
fcfP  = fcP(lt,:);
%====== Carbon 13
M13calA = sum(M13calvA.*VsvA)/VsA;
M13calI = sum(M13calvI.*VsvI)/VsI;
M13calP = sum(M13calvP.*VsvP)/VsP;
f13cfA  = f13cA(lt,:);
f13cfI  = f13cI(lt,:);
f13cfP  = f13cP(lt,:);

% test: Average must equal Rin
RsAf = (f13cfA./fcfA/Rst-1)*1e3;
RsIf = (f13cfI./fcfI/Rst-1)*1e3;
RsPf = (f13cfP./fcfP/Rst-1)*1e3;

if(ftys)
McalT = sum(McalvT.*VsvT)/VsT;
fcfT  = fcT(lt,:);
M13calT = sum(M13calvT.*VsvT)/VsT;
f13cfT  = f13cT(lt,:);
RsTf = (f13cfT./fcfT/Rst-1)*1e3;
end;

end; %----------- end sediments


% inventories/average concentrations ocean
MCO_init = sum(Y(1, 1:Nb).*V)    *12/1e15;      %initial CO2 ocean [Pg]
MCA_init = Y(1,NO*Nb+1  )*Aoc*12/1e15;          %initial CO2 atmosph [Pg]
total_C_init = MCO_init + MCA_init;             %initial total carbon [Pg]
Ma_init  = sum(Y(1,   Nb+1:2*Nb).*V)/Voc;       %initial Alkalinity/V


Mc  = sum(Y(lt,     1:   Nb).*V)/Voc ... % TC inventory/V
         +Y(lt,NO*Nb+1     )*Aoc/Voc;    % atmosphere
Ma  = sum(Y(lt,   Nb+1:2*Nb).*V)/Voc;    % TA inventory/V
Mp  = sum(Y(lt, 2*Nb+1:3*Nb).*V)/Voc;    % PO4 inventory/V
Mcc = sum(Y(lt, 3*Nb+1:4*Nb).*V)/Voc ... % TC inventory/V
         +Y(lt,NO*Nb+2     )*Aoc/Voc;    % atmosphere

DX0 = [Mc-Mci Ma-Mai Mp-Mpi Mcc-Mcci];
DX0

% total C inventory ocean     (Pg C)
MCO=sum(Y(lt, 1:Nb).*V)    *12/1e15;
% total C inventory atm       (Pg C)
MCA=    Y(lt,NO*Nb+1  )*Aoc*12/1e15;
% total C inventory ocean+atm (Pg C)
Mc*Voc*12/1e15;


for i=1:lt-1
tvb(i)     = tv(i);    
dtv(i)     = tv(i+1)-tv(i);
end

if(fsed)   
%----------------- MASS BALANCE CHECK -------------%
%
% mass balance c and a. At final t
% Dc  = total influx - rain + dissolution
%
% mass balance mcal.    At final t
% Dmc = total rain - burial - dissolution

% calculate burial and dissolution at all t.
% (kg/y). (times dt) -> kg

for i=1:lt-1
for k=1:Ns    
burialtA(i,k) = (  mcalA(i+1,k)+  mcalA(i,k)) ...
               *(rburvtA(i+1,k)+rburvtA(i,k)) ...
                       *asvA(k)*A(1)*dtv(i)/4; % kg
burialtI(i,k) = (  mcalI(i+1,k)+  mcalI(i,k)) ...
               *(rburvtI(i+1,k)+rburvtI(i,k)) ...
                       *asvI(k)*A(2)*dtv(i)/4; % kg
burialtP(i,k) = (  mcalP(i+1,k)+  mcalP(i,k)) ...
               *(rburvtP(i+1,k)+rburvtP(i,k)) ...
                       *asvP(k)*A(3)*dtv(i)/4; % kg
DisstA(i,k)   = ...
       (dissvtA(i+1,k)+dissvtA(i,k))*dtv(i)/2; % kg
DisstI(i,k)   = ...
       (dissvtI(i+1,k)+dissvtI(i,k))*dtv(i)/2; % kg
DisstP(i,k)   = ...
       (dissvtP(i+1,k)+dissvtP(i,k))*dtv(i)/2; % kg
Diss13tA(i,k)   = ...
       (dissv13tA(i+1,k)+dissv13tA(i,k))*dtv(i)/2; % kg
Diss13tI(i,k)   = ...
       (dissv13tI(i+1,k)+dissv13tI(i,k))*dtv(i)/2; % kg
Diss13tP(i,k)   = ...
       (dissv13tP(i+1,k)+dissv13tP(i,k))*dtv(i)/2; % kg
if(ftys)
burialtT(i,k) = (  mcalT(i+1,k)+  mcalT(i,k)) ...
               *(rsedvtT(i+1,k)+rsedvtT(i,k)) ...
                       *asvT(k)*A(11)*dtv(i)/4; % kg
DisstT(i,k)   = ...
       (dissvtT(i+1,k)+dissvtT(i,k))*dtv(i)/2; % kg
Diss13tT(i,k)   = ...
       (dissv13tT(i+1,k)+dissv13tT(i,k))*dtv(i)/2; % kg
end;    
end;
FprdtA  (i) = (FprtA  (i+1)+FprtA(  i))*dtv(i)/2;
FprdtI  (i) = (FprtI  (i+1)+FprtI  (i))*dtv(i)/2;
FprdtP  (i) = (FprtP  (i+1)+FprtP  (i))*dtv(i)/2;
Fpr13dtA(i) = (Fpr13tA(i+1)+Fpr13tA(i))*dtv(i)/2;
Fpr13dtI(i) = (Fpr13tI(i+1)+Fpr13tI(i))*dtv(i)/2;
Fpr13dtP(i) = (Fpr13tP(i+1)+Fpr13tP(i))*dtv(i)/2;
Findt   (i) = (Fint   (i+1)+Fint   (i))*dtv(i)/2; % mol C/m2
Fin13dt (i) = (Fin13t (i+1)+Fin13t (i))*dtv(i)/2; % mol C/m2
FSidt   (i) = (FSit   (i+1)+FSit   (i))*dtv(i)/2; % mol C/m2
FSi13dt (i) = (FSi13t (i+1)+FSi13t (i))*dtv(i)/2; % mol C/m2
if(ftys)
FprdtT  (i) = (FprtT  (i+1)+FprtT  (i))*dtv(i)/2;
Fpr13dtT(i) = (Fpr13tT(i+1)+Fpr13tT(i))*dtv(i)/2;
end;    
end;

% Atlantic deep dissolution
DisstAd = DisstA(:,4:Ns);

% "sediment core depth above base", base at t=0

for k=1:Ns    
rdtA(:,k) = rburvtA(1:lt,k).*[0 dtv]';
rdtI(:,k) = rburvtI(1:lt,k).*[0 dtv]';
rdtP(:,k) = rburvtP(1:lt,k).*[0 dtv]';
if (ftys)
rdtT(:,k) = rsedvtT(1:lt,k).*[0 dtv]';
end;
end;

for i=1:lt
for k=1:Ns    
dabA(i,k) = sum(rdtA(1:i,k));
dabI(i,k) = sum(rdtI(1:i,k));
dabP(i,k) = sum(rdtP(1:i,k));
if (ftys)
dabT(i,k) = sum(rdtT(1:i,k));
end;
end;
end;

% integral (sum) over time/boxes
BurialA  = sum(sum(burialtA))/VsA; % kg /m3 sed
BurialI  = sum(sum(burialtI))/VsI; % kg /m3 sed
BurialP  = sum(sum(burialtP))/VsP; % kg /m3 sed
BurialCA = BurialA*VsA/Voc/m2kg;   % mol/m3 ocean
BurialCI = BurialI*VsI/Voc/m2kg;   % mol/m3 ocean
BurialCP = BurialP*VsP/Voc/m2kg;   % mol/m3 ocean
BBurialC = sum(BurialCA+BurialCI+BurialCP);
% ovet time (per sed box)
DissAk   = sum(DisstA,1)/m2kg*12/1e15; % mol -> PgC
DissIk   = sum(DisstI,1)/m2kg*12/1e15; % mol -> PgC
DissPk   = sum(DisstP,1)/m2kg*12/1e15; % mol -> PgC
% time and sed boxes
DissA    = sum(sum(DisstA))  /VsA; % kg /m3 sed
DissI    = sum(sum(DisstI))  /VsI; % kg /m3 sed
DissP    = sum(sum(DisstP))  /VsP; % kg /m3 sed
DissCA   = DissA*VsA/Voc/m2kg;     % mol/m3 ocean
DissCI   = DissI*VsI/Voc/m2kg;     % mol/m3 ocean
DissCP   = DissP*VsP/Voc/m2kg;     % mol/m3 ocean
DDissC   = sum(DissCA+DissCI+DissCP);
%== Carbon-13
Diss13A    = sum(sum(Diss13tA))  /VsA; % kg /m3 sed
Diss13I    = sum(sum(Diss13tI))  /VsI; % kg /m3 sed
Diss13P    = sum(sum(Diss13tP))  /VsP; % kg /m3 sed
Diss13CA   = Diss13A*VsA/Voc/m2kg;     % mol/m3 ocean
Diss13CI   = Diss13I*VsI/Voc/m2kg;     % mol/m3 ocean
Diss13CP   = Diss13P*VsP/Voc/m2kg;     % mol/m3 ocean
DDiss13C   = sum(Diss13CA+Diss13CI+Diss13CP);
%- deep Atl
DissCAd    = sum(sum(DisstAd))/Voc/m2kg; % 

if(ftys)
BurialT  = sum(sum(burialtT))/VsT; % kg /m3 sed
BurialCT = BurialT*VsT/Voc/m2kg;   % mol/m3 ocean
BBurialC = sum(BurialCA+BurialCI+BurialCP+BurialCT);
DissT    = sum(sum(DisstT))  /VsT; % kg /m3 sed
DissCT   = DissT*VsT/Voc/m2kg;     % mol/m3 ocean
DDissC   = sum(DissCA+DissCI+DissCP+DissCT);
Diss13T    = sum(sum(Diss13tT))  /VsT; % kg /m3 sed
Diss13CT   = Diss13T*VsT/Voc/m2kg;     % mol/m3 ocean
DDiss13C   = sum(Diss13CA+Diss13CI+Diss13CP+Diss13CT);
end;    


% right-hand sides
FinT   =  sum(Findt)/Dt;                    % mol/m2/y
FSiT   =  sum(FSidt)/Dt;                    % mol/m2/y
FprT   =  sum(  sum(FprdtA*A(1)) ...        % mol/m2/y
               +sum(FprdtI*A(2)) ...
               +sum(FprdtP*A(3)) )/Aoc/Dt;
% right-hand sides
Fin13T   =  sum(Fin13dt)/Dt;                % mol/m2/y
FSi13T   =  sum(FSi13dt)/Dt;                % mol/m2/y
Fpr13T   =  sum(sum(Fpr13dtA*A(1)) ...      % mol/m2/y
               +sum(Fpr13dtI*A(2)) ...
               +sum(Fpr13dtP*A(3)) )/Aoc/Dt;
           
if(ftys)           
FprT   = FprT   + sum( sum(FprdtT  *A(11)) )/Aoc/Dt;
Fpr13T = Fpr13T + sum( sum(Fpr13dtT*A(11)) )/Aoc/Dt;
end;    
           
rhc    = (FinT+FSiT-FprT)*Aoc*Dt/Voc+DDissC; % mol/m3 ocean
%rhMc   =  FprT*m2kg*Aoc*Dt/Vs-Burial-Diss;  % kg /m3 sed

rhc13  = (Fin13T+FSi13T-Fpr13T)*Aoc*Dt/Voc+DDiss13C; % mol/m3 ocean

% DXs = Change of [Mass_Carbon Alkalinity CCD_depth] between initial cond. and last timestep
% Mci = initial Mass 
% Mc = Mass at last timestep
DXs = [(Mc-Mci)-rhc (Ma-Mai)-2*rhc Mcc-Mcci-rhc13];


end; % fsed ------- MASS BALANCE CHECK END -------%



%==== Anthropogenic CO2 uptake ===============%
if(ffflag)
FCO2t   = dYdt(NO*Nb+1,:)*Aoc*12/1e15; % Gt C/y, dC/dt atmosphere
% CO2 emissions, see PEmisfun
for i=1:lt
Femt(i) = PEmisfun(tv(i));
end;
% uptake is diff between em. and dC/dt atm
Fupt = Femt-FCO2t;
% cumulative em/uptake
for i=1:lt-1
dMFemt(i) = (Femt(i+1)+Femt(i))*dtv(i)/2;
dMFupt(i) = (Fupt(i+1)+Fupt(i))*dtv(i)/2;
end;
MFemt = sum(dMFemt);
MFupt = sum(dMFupt);

% [MFemt MFupt]
% [max(pco2t) pco2t(end)]

end;%=========================================%


% ocean tracer and pCO2
if(fdox)
for k=1:Nb
c  (:,k)   = Y(:,     k)/(1e-3*rho); % mol/m3 -> mmol/kg
a  (:,k)   = Y(:,  Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
p  (:,k)   = Y(:,2*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
dox(:,k)   = Y(:,3*Nb+k)           ; % mol/m3
cc (:,k)   = Y(:,4*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
%====== Carbon-13
d13c(:,k) = ((cc(:,k)./c(:,k))/Rst-1.)*1e3; % Ocean
end;
C        = Y(:,5*Nb+1); 
CC       = Y(:,5*Nb+2); 
pco2t    =  C/(2.2e15/12/Aoc);
pcco2t   = CC/(2.2e15/12/Aoc);
pco2a    =  pco2t(lt);
pcco2a   = pcco2t(lt);
pco2a
doxf = dox(end,:);
else  % fdox
for k=1:Nb
c  (:,k)   = Y(:,     k)/(1e-3*rho); % mol/m3 -> mmol/kg
a  (:,k)   = Y(:,  Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
p  (:,k)   = Y(:,2*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
cc (:,k)   = Y(:,3*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
%====== Carbon-13
d13c(:,k) = ((cc(:,k)./c(:,k))/Rst-1.)*1e3; % Ocean
end;
C        = Y(:,4*Nb+1); 
CC       = Y(:,4*Nb+2); 
pco2t    =  C/(2.2e15/12/Aoc);
pcco2t   = CC/(2.2e15/12/Aoc);
pco2a    =  pco2t(lt);
pcco2a   = pcco2t(lt);
pco2a
end; % fdox



%====== Carbon-13
d13cv = d13c(lt,:);                  % Ocean
d13C  = ((CC./C       )/Rst-1.)*1e3; % Atmosphere
d13Ca = ((CC(lt)/C(lt))/Rst-1.)*1e3; % Atmosphere

d13CV = [d13cv d13Ca];
d13CV

% 13C surf-deep diff & excursion
k = 1;
Dd13cv  = d13c(1,1:3)-d13c(1,7:9);
d13cex = max(d13c(:,k))-min(d13c(:,k));

% Final carbonate system (needs mol/kg)

for k=1:Nb
dic = c(lt,k)*1e-3;
alk = a(lt,k)*1e-3; 
TC = TCvt(lt,k);
S  = Sv(k);
P  = Pv(k);          
[co2,pco2,co3,ph,kh,o2] = dafunPE(dic,alk,TC,S,P,Ca,Mg);
khv(k)   = kh;
co2v(k)  = co2;
co3v(k)  = co3;
phv(k)   = ph;
o2v(k)   = o2;
pco2v(k) = pco2;
end;

% carbonate ion, ph over time (needs mol/kg)
for i=1:lt
for k=1:Nb
dic = c(i,k)*1e-3;
alk = a(i,k)*1e-3; 
TC = TCvt(i,k);
S  = Sv(k);
P  = Pv(k);          
[co2,pco2,co3,ph,kh,o2] = dafunPE(dic,alk,TC,S,P,Ca,Mg);
co3tv(i,k)  = co3;
 phtv(i,k)  = ph;
end;
end;

%==== calculate omega(t)
%==== saturation of surface ocean boxes

for i=1:lt
l = 1;
for k=kkv
[kspcS,kspaS] = ...
     kspfun(TCvt(i,k),Sv(k),Pv(k),Ca,Mg);
omegCSvt(i,l) = co3tv(i,k)*Ca./kspcS;
omegASvt(i,l) = co3tv(i,k)*Ca./kspaS;
l = l + 1;
end
end

% final deep CO3 and sediment values
co3Af = co3tv(lt,07)*1e6;
co3If = co3tv(lt,08)*1e6;
co3Pf = co3tv(lt,09)*1e6;
co3ALLf = [co3Af co3If co3Pf]';
if(ftys)
co3Tf = co3tv(lt,13)*1e6;
co3ALLf = [co3Af co3If co3Pf co3Tf]';
end;
if(fsed) % sediments
for l=1:Ns
mcalAi(l) = mcalA(01,l)*VsvA(l)/(1e15*100/12/1e3); % Gt C
mcalIi(l) = mcalI(01,l)*VsvI(l)/(1e15*100/12/1e3); % Gt C
mcalPi(l) = mcalP(01,l)*VsvP(l)/(1e15*100/12/1e3); % Gt C
mcalfA(l) = mcalA(lt,l)*VsvA(l)/(1e15*100/12/1e3); % Gt C
mcalfI(l) = mcalI(lt,l)*VsvI(l)/(1e15*100/12/1e3); % Gt C
mcalfP(l) = mcalP(lt,l)*VsvP(l)/(1e15*100/12/1e3); % Gt C
if(ftys)
mcalTi(l) = mcalT(01,l)*VsvT(l)/(1e15*100/12/1e3); % Gt C
mcalfT(l) = mcalT(lt,l)*VsvT(l)/(1e15*100/12/1e3); % Gt C
end;
end;
  fcALLf = [fcfA; fcfI; fcfP]';
McalALLi = [sum(mcalAi) sum(mcalIi) sum(mcalPi)]';
McalALLf = [sum(mcalfA) sum(mcalfI) sum(mcalfP)]';
if(ftys)
  fcALLf = [fcfA; fcfI; fcfP; fcfT]';
McalALLi = [sum(mcalAi) sum(mcalIi) sum(mcalPi) sum(mcalTi)]';
McalALLf = [sum(mcalfA) sum(mcalfI) sum(mcalfP) sum(mcalfT)]';
end;    
McALL_fi = [McalALLf McalALLi];
MMcali   = sum(McalALLi);
MMcalf   = sum(McalALLf);
MM_fi    = [MMcalf MMcali];

co3ALLf
 fcALLf
McALL_fi
MM_fi



% saturation horizon deep boxes
if    (satflg == 1)
co3szv = co3s0*exp(as*(zv-zs0));
elseif(satflg == 2)
for k=1:lzv
[x,x,x,x,x,x,kspc(k),x] = dafunPE(1e-3,1e-3,2,Soc,zv(k)/10.,Ca,Mg);
end;
co3szv = kspc/Ca;
end;

    
for i=1:lt
[tmp, jshA] = min(abs(co3tv(i,07)-co3szv));
[tmp, jshI] = min(abs(co3tv(i,08)-co3szv));
[tmp, jshP] = min(abs(co3tv(i,09)-co3szv));
shtA(i) = zv(jshA);
shtI(i) = zv(jshI);
shtP(i) = zv(jshP);
if(ftys)
[tmp, jshT] = min(abs(co3tv(i,13)-co3szv));
shtT(i) = zv(jshT);
end;
end;

shA   = shtA(lt);
shI   = shtI(lt);
shP   = shtP(lt);

sh_ALL = [shA shI shP];

if(ftys)
shT   = shtT(lt);
sh_ALL = [sh_ALL shT];
end;

sh_ALL


%------------ CCD --------------------%
fccd  = 0.10; % 0.10 0.05
dd    = 0.01;
nd    = 0;
fccdv = [fccd-nd*dd:dd:fccd+nd*dd];
ld    = length(fccdv);

clear jccdA jccdI jccdP;
if(ftys) clear jccdT; end;
for i=1:lt
yA(i,:)         = interp1(dsv,fcA(i,:),zv,'cubic');
yI(i,:)         = interp1(dsv,fcI(i,:),zv,'cubic');
yP(i,:)         = interp1(dsv,fcP(i,:),zv,'cubic');
if(ftys)
yT(i,:)         = interp1(dsv,fcT(i,:),zv,'cubic');
end;    
for k=1:ld 
[tmp, jccdA(i,k)] = min(abs(yA(i,:)-fccdv(k)));
[tmp, jccdI(i,k)] = min(abs(yI(i,:)-fccdv(k)));
[tmp, jccdP(i,k)] = min(abs(yP(i,:)-fccdv(k)));
if(ftys)
[tmp, jccdT(i,k)] = min(abs(yT(i,:)-fccdv(k)));
end;    
end;
end;

if(nd == 0)
ccdA = zv(jccdA);
ccdI = zv(jccdI);
ccdP = zv(jccdP);
if(ftys)
ccdT = zv(jccdT);    
end;
else
ccdA = sum(zv(jccdA),2)/ld;
ccdI = sum(zv(jccdI),2)/ld;
ccdP = sum(zv(jccdP),2)/ld;
if(ftys)
ccdT = sum(zv(jccdT),2)/ld;
end;    
end;


% initial sediment inventory vs. dissolution
MciGtC    = MMcali;                % Gt C
DissGtC   = DDissC *12*Voc/1e15;   % Gt C
DissGtCA  = DissCA *12*Voc/1e15;   % Gt C
DissGtCAd = DissCAd*12*Voc/1e15;   % Gt C
MciGtC;
DissGtC
DissGtCA;
rers    = DissGtC   /MciGtC;
rersA   = DissGtCA /sum(mcalAi);
rersAd  = DissGtCAd/sum(mcalAi(kd));
rersGAd = [rers rersAd];
rers;
rersA;
% estimate
rest1 = 1+fc0A(1)/(1-fc0A(1));
%rest2 = 1+fc0A(1)/(1-fc0A(1))*(1-phi0)/(1-phi1);
% average of all sed boxes
rest2 = 1+fc0A./(1-fc0A)*(1-phi0)/(1-phi1);
rest2 = sum(asvA.*rest2');
r_estA12 = [rest1 rest2];
r_estA12;

%>>>>>>>>>>>>>>>>>>>>>> NEW 07/01/06


RdissA = DissAk./mcalAi;
RdissI = DissIk./mcalIi;
RdissP = DissPk./mcalPi;
RestA  = 1+fc0A./(1-fc0A)*(1-phi0)/(1-phi1);
RestI  = 1+fc0I./(1-fc0I)*(1-phi0)/(1-phi1);
RestP  = 1+fc0P./(1-fc0P)*(1-phi0)/(1-phi1);

Rdiss_Rest_A = [RdissA; RestA']';
Rdiss_Rest_A

end; % sediments


if(fsed)
if(1)
%====== shelf/deep rain
% global rain
% asv(1:2):  <600m
% asv(3:Ns): >600m
rainsh(1) = sum(fsh    *FprtA(lt)*asvA(1: 2)*A(1));
rainds(1) = sum(fdpv(1)*FprtA(lt)*asvA(3:Ns)*A(1));

rainsh(2) = sum(fsh    *FprtI(lt)*asvI(1: 2)*A(2));
rainds(2) = sum(fdpv(2)*FprtI(lt)*asvI(3:Ns)*A(2));

rainsh(3) = sum(fsh    *FprtP(lt)*asvP(1: 2)*A(3));
rainds(3) = sum(fdpv(3)*FprtP(lt)*asvP(3:Ns)*A(3));

if(ftys)
rainsh(4) = sum(fshT   *FprtT(lt)*asvT(1: 2)*A(11)); 
rainds(4) = sum(fdpv(4)*FprtT(lt)*asvT(3:Ns)*A(11));
end;


%
RainSH = sum(rainsh); % mol C y
RainDS = sum(rainds); % mol C y
RainT  = RainSH+RainDS;

Rdssh = RainDS/RainSH;

% areas
Ash(1) = sum(asvA(1: 2)*A(1));
Ads(1) = sum(asvA(3:Ns)*A(1));
Ash(2) = sum(asvI(1: 2)*A(2));
Ads(2) = sum(asvI(3:Ns)*A(2));
Ash(3) = sum(asvP(1: 2)*A(3));
Ads(3) = sum(asvP(3:Ns)*A(3));
if(ftys)
Ash(4) = sum(asvT(1: 2)*A(11));
Ads(4) = sum(asvT(3:Ns)*A(11));
end;

% total
TAsh = sum(Ash);
TAds = sum(Ads);
TRA  = (TAsh+TAds);
%(1-fH)*Aoc

% rain per square meter in basins

rshsqm = rainsh./Ash; % mol/m2/y
rdssqm = rainds./Ads; % mol/m2/y

% average
rSHsqm = sum(rainsh)/TAsh;
rDSsqm = sum(rainds)/TAds;

rsqm = rSHsqm/rDSsqm;

Rdssh_rsqm = [Rdssh rsqm];
Rdssh_rsqm

% [fsh*[1 1 1] fshT]./fdpv
% rshsqm./rdssqm

% remainder flux per sqm
% 
frrfsh = frrf*(RainSH/TAsh)/(RainT/TRA);
frrfds = frrf*(RainDS/TAds)/(RainT/TRA);

frrf_frrfsh_frrfds = [frrf frrfsh frrfds]*1e2;
frrf_frrfsh_frrfds

end;
end;


if(ccdrun)
  ccdv0   = [  ccdA(1)   ccdI(1)   ccdP(1)   ccdT(1)]'; 
  ccdvmin = [min(ccdA) min(ccdI) min(ccdP) min(ccdT)]'; 
  pCO2Max = max(pco2t);
  [tmp,k1]=min(abs(tv-10.e3));
  [tmp,k2]=min(abs(tv-72.e3));
  pCO2Av  = mean(pco2t(k1:k2));
  oxminA  = [min(dox(:,4)) min(dox(:,7))]';
  ccdv0_ccdvmin = [ccdv0 ccdvmin];
  ccdv0_ccdvmin
  pCO2Max
  oxminA
  pCO2Av
  
  PgC     = num2str(CBl/1.e15);
  Aox     = num2str(oxA*100);
  PgC_Aox = [PgC '_' Aox];
  PgC_Aox
  
  % write to file
  X = [CBl/1.e15; oxA; ccdv0; ccdvmin; pCO2Max; oxminA; pCO2Av];
  outstr = ['dat\B2D3BL2SW\ccdrun\out' PgC_Aox '.DAT'];
  %outstr = ['dat\B2D3BL2SW\ccdrun\NoSwcon\out' PgC_Aox '.DAT'];
  save(outstr,'X','-ascii','-double','-tabs');
return
end    


if(ffflag)
  % pCO2, pH
  pCO2Max = max(pco2t);
  pHAInit = phtv(1,1);
  pHAMin  = min(phtv(:,1));
  DeltapH = phtv(1,1)-min(phtv(:,1));
  pCO2Max_pHAInit_pHAMin_DeltapH = ...
      [pCO2Max/1.e3 pHAInit pHAMin DeltapH];
  pCO2Max_pHAInit_pHAMin_DeltapH
  % saturation, calcite
  omCInitL =     omegCSvt(1,1);
  omCInitH =     omegCSvt(1,4);
  omCMinL  = min(omegCSvt(:,1));
  omCMinH  = min(omegCSvt(:,4));
  % saturation, aragonite
  omAInitL =     omegASvt(1,1);
  omAInitH =     omegASvt(1,4);
  omAMinL  = min(omegASvt(:,1));
  omAMinH  = min(omegASvt(:,4));
  % CO3=
  co3InitL =     co3tv(1,01);
  co3InitH =     co3tv(1,10);
  co3MinL  = min(co3tv(:,01));
  co3MinH  = min(co3tv(:,10));
  if(exist('dirstr','var'))
  PgC     = str2num(dirstr(1:4));
  RelTime = str2num(dirstr(6:9));
  PgC_RelTime = [PgC RelTime];
  PgC_RelTime
  
  % write to file
  X = [PgC RelTime pCO2Max pHAInit pHAMin DeltapH]';
  X2= [PgC RelTime omCInitL omCInitH omCMinL omCMinH ...
                   omAInitL omAInitH omAMinL omAMinH ...
                   co3InitL co3InitH co3MinL co3MinH]';
  %outstr = ['dat\Emss\EmssScen\out' dirstr ];
  %outstr = ['dat\Emss\EmssScen\sat\out' dirstr ];
  %save(outstr,'X' ,'-ascii','-double','-tabs');
  %save(outstr,'X2','-ascii','-double','-tabs');
  end
  %return
end

  if(contourflag)
  % pCO2, pH
  pCO2Init = pco2t(1);
  pCO2Max = max(pco2t);
  DeltapCO2 =max(pco2t)-pco2t(1);
  
  % pH Atlantik low lat. surface box
  pH_A_Init = phtv(1,1);
  pH_A_Min  = min(phtv(:,1));    
  DeltaApH = phtv(1,1)-min(phtv(:,1));
  
  % save initial and minimum calcite saturation state (Atlantic surface box)
  Omega_init = omegCSvt(1,1);
  Omega_min  = min(omegCSvt(:,1));
  
  %pCO2Init_pCO2Max_DeltapCO2_pH_A_Init_pH_A_Min_DeltaApH = ...
   %   [pCO2Init pCO2Max DeltapCO2 pH_A_Init pH_A_Min DeltaApH];
  
 % pCO2Init_pCO2Max_DeltapCO2_pH_A_Init_pH_A_Min_DeltaApH
  

  
  % write to file 
  
  %pco2initial runs:
  %X = [pCO2Init pCO2Max DeltapCO2 pH_A_Init pH_A_Min DeltaApH]';
  
  % ccd runs:
  %X = [pCO2Init pCO2Max pH_A_Init pH_A_Min DeltaApH ccdA(lt) ccdI(lt) ccdP(lt)];
  
  % Mg/Ca runs (saves also the bath flag index):
  X = [pCO2Init pCO2Max DeltapCO2 pH_A_Init pH_A_Min DeltaApH Mg Ca Mgca_ratio Omega_init Omega_min bath]';
  
  end
  
 

% save the last timestep results for Y
Y_end = Y(length(tv),:);
Y_end = Y_end';             % transform to column vektor


DIC_of_all_boxes_init = sum(Y(1,     1:   Nb).*V)/Voc    % [mmol/kg]
TA_of_all_boxes_init  = sum(Y(1,   Nb+1:2*Nb).*V)/Voc   % [mmol/kg]

DIC_of_all_boxes_end = sum(Y(lt,     1:   Nb).*V)/Voc    % [mmol/kg]
TA_of_all_boxes_end  = sum(Y(lt,   Nb+1:2*Nb).*V)/Voc   % [mmol/kg]

%=====================================================%
%=================== plot results ====================%
%=====================================================%
%if(~plotflag)
%return
%end    
if(plotflag)

LineWidth = 1.3;
fs = 14;
%cs = 'gggkkkrrrb';
%sstr = '- ---.- ---.- ---. -';    
%lstr = ' HLAIADALIIIDILPIPDP';
%lstr0 = 'LALILPIAIIIPDADIDP H';

cs = 'gggkkkrrrbgkr';
sstr = '- ---.- ---.- ---. -: : : ';    
lstr0 = 'LALILPIAIIIPDADIDP HLTITDT';
for k=1:Nb
lstr(k,:) = sprintf('%s',lstr0(2*k-1:2*k));
end;


if(logax == 1 & axx(1) < 0)
  axx(1) = 1.e-1;
end;  

% tv11 is the time-vector tv with axx(1) added infront of it
tv11   = [axx(1)    tv(2)]; % tv(2)
c1     = [c(1,:)' c(1,:)'];
a1     = [a(1,:)' a(1,:)'];
p1     = [p(1,:)' p(1,:)'];
co3tv1 = [co3tv(1,:)' co3tv(1,:)'];
 phtv1 = [ phtv(1,:)'  phtv(1,:)'];

%---------- TCO2 ---------------%
figure(1)
clf; 
box  on;
hold on;
for k=1:Nb
plot(tv  ,c  (:,k),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
for k=1:Nb
plot(tv11,c1 (k,:),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('TCO_2 (mmol kg^{-1})');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);
refresh;

%return;

%-------------- TA --------------------%
figure(2) 
clf; 
box  on;
hold on;
for k=1:Nb
plot(tv  ,a  (:,k),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
for k=1:Nb
plot(tv11,a1 (k,:),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('TA (mmol kg^{-1})');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);


%-------------- PO4 --------------------%
figure(3)
clf; 
box  on;
hold on;
for k=1:Nb
plot(tv  ,p (:,k)*1e3,sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
for k=1:Nb
plot(tv11,p1(k,:)*1e3,sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('PO_4 (\mumol kg^{-1})');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);


%-------------- O2 --------------------%
if(fdox)
dox1   = [dox(1,:)' dox(1,:)'];

figure(35)
clf; 
box  on;
hold on;
for k=1:Nb
plot(tv  ,dox (:,k),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
for k=1:Nb
plot(tv11,dox1(k,:),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('O_2 (mol m^{-3})');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);

figure(36)
clf;
box  on;
hold on;
for k=1:Nb
plot(tv  ,dox (:,k)./(dox (:,k)+KMMOX)*100,...
    sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
for k=1:Nb
plot(tv11,dox1(k,:)./(dox1(k,:)+KMMOX)*100,...
    sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('%Oxic Respiration');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);


end; % fdox


%-------------- CO3 --------------------%
figure(4)
clf; 
box  on;
hold on;
for k=1:Nb
plot(tv  ,co3tv (:,k)*1e6,sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
for k=1:Nb
plot(tv11,co3tv1(k,:)*1e6,sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('[CO_3^{2-}] (\mumol kg^{-1})');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);

if(1)
%-------------- ph --------------------%
figure(50)
clf; 
box  on;
hold on;
for k=1:Nb
plot(tv  ,phtv (:,k),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
for k=1:Nb
plot(tv11,phtv1(k,:),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('pH');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);
end;

%-------------- pCO2 --------------------%
systr = 'rdrsrx            bdro';    
pco2t1  = [pco2t(1)' pco2t(1)'];
figure(5)
plot(tv  ,pco2t ,'b-','LineWidth',LineWidth);
hold on;
for k=kkv
plot(tv(lt),pco2v(k),systr(2*k-1:2*k));
end;
plot(tv11,pco2t1,'b-','LineWidth',LineWidth);
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('pCO_2 (\muatm)');

pstr0 = 'a LALILP HLT';
for k=1:6
pstr(k,:) = sprintf('%s',pstr0(2*k-1:2*k));
end;
Hp=legend(pstr,1);
set(Hp,'FontSize',10);


%---------- d13C ----------------%
d13c1  = [d13c(1,:)' d13c(1,:)'];
d13C1  = [d13C(1)'   d13C(1)'  ];

MS  = 'MarkerSize';
ms  = 4;
MFC = 'MarkerFaceColor';
Dam = 7.;

%-------------- load PETM data
if(ftys)
%== Zachos data
ZchF3a = load('~/Desktop/malte_ubuntu/myLoscar/dat/PETM/ZchF3aNEW.DAT');
age3  = ZchF3a(:,3)*1e3; 
c13c3 = ZchF3a(:,04);
%== Roehl data
% new age model, G^3, 2007
NEW690     = load('~/Desktop/malte_ubuntu/myLoscar/dat/PETM/NEW690.txt');
d690n      = NEW690(:,2);
a690n      = NEW690(:,3);
% c13 data
Roehl690   = load('~/Desktop/malte_ubuntu/myLoscar/dat/PETM/Roehl690.txt');
d690       = Roehl690(:,1);
% interpolate 
age690blk  = interp1(d690n,a690n,d690)*1e3;
c13c690blk = Roehl690(:,5);
end % ftys


figure(6)
clf; 
box  on;
hold on;
for k=1:Nb
plot(tv  ,d13c (:,k),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
plot(tv  ,d13C +Dam,'m-','LineWidth',LineWidth);
for k=1:Nb
plot(tv11,d13c1(k,:),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',LineWidth);
end;
plot(tv11,d13C1+Dam,'m-','LineWidth',LineWidth);
if(ftys)
plot(age3,     c13c3     ,'b--s',MFC,'b',MS,ms,'LineWidth',LineWidth);
plot(age690blk,c13c690blk,'r- d',MFC,'r',MS,ms,'LineWidth',LineWidth);
end
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
set(gca,'YDir','reverse');
xlabel('Time (y)');
ylabel('\delta^{13}C (???)');
Hl=legend(lstr,1);
set(Hl,'FontSize',10);


if(fsed)    
%---------- fc vs. depth----------%
figure(7)
clf; 
plot(fcfA*1e2,dsv,'b-d','LineWidth',LineWidth);
hold on;
plot(fcfI*1e2,dsv,'m-s','LineWidth',LineWidth);
plot(fcfP*1e2,dsv,'k-o','LineWidth',LineWidth);
if(ftys)
plot(fcfT*1e2,dsv,'g-p','LineWidth',LineWidth);
end;
hold off;
set(gca,'YDir','reverse');
set(gca,'FontSize',fs);
set(gca,'YDir','reverse');
xlabel('CaCO_3 (wt %)');
ylabel('Depth (m)');
%axis([0 100 0 5.5]);
legend('Atlantic','Indic','Pacific',2);
if(ftys)
legend('Atlantic','Indic','Pacific','Tethys',2);
end;


%---------- CO3, SatHor -----------%
if(1)
shtA1  = [shtA(1)' shtA(1)'];    
shtI1  = [shtI(1)' shtI(1)'];    
shtP1  = [shtP(1)' shtP(1)'];    
ccdA1  = [ccdA(1)' ccdA(1)'];    
ccdI1  = [ccdI(1)' ccdI(1)'];    
ccdP1  = [ccdP(1)' ccdP(1)'];    
    
pstr8 = '';    
yl = [0 120];
yr = [6.2 0];

figure(8)
clf; 
f2  = (yl(2)-yl(1))/(yr(2)-yr(1));
f1  =  yl(2)-f2*yr(2);
b1  = [axx yl(1) yl(2)];
axis off;
a0 = get(gca,'Position');
%a = [a0(1) a0(2) a0(3) .7*a0(4)];
h = axes('position',a0);
axis(b1);
if(logax)
set(gca,'XScale','log');
end;
hold on;
for k=7:9
pl(k-6)= ...
plot(tv  ,co3tv (:,k)*1e6   ,sstr(2*k-1:2*k),'LineWidth',LineWidth);
plot(tv11,co3tv1(k,:)*1e6   ,sstr(2*k-1:2*k),'LineWidth',LineWidth);
end;
pl(2)=  ...
plot(tv  ,f1+shtA *f2/1000,'r- ','LineWidth',LineWidth);
plot(tv  ,f1+shtI *f2/1000,'r--','LineWidth',LineWidth);
plot(tv  ,f1+shtP *f2/1000,'r-.','LineWidth',LineWidth);
plot(tv11,f1+shtA1*f2/1000,'r- ','LineWidth',LineWidth);
plot(tv11,f1+shtI1*f2/1000,'r--','LineWidth',LineWidth);
plot(tv11,f1+shtP1*f2/1000,'r-.','LineWidth',LineWidth);
pl(3)=  ...
plot(tv  ,f1+f2*ccdA /1e3,'k- ','LineWidth',LineWidth);
plot(tv  ,f1+f2*ccdI /1e3,'k--','LineWidth',LineWidth);
plot(tv  ,f1+f2*ccdP /1e3,'k-.','LineWidth',LineWidth);
plot(tv11,f1+f2*ccdA1/1e3,'k- ','LineWidth',LineWidth);
plot(tv11,f1+f2*ccdI1/1e3,'k--','LineWidth',LineWidth);
plot(tv11,f1+f2*ccdP1/1e3,'k-.','LineWidth',LineWidth);
if(ftys)
k=13;   
shtT1  = [shtT(1)' shtT(1)'];    
ccdT1  = [ccdT(1)' ccdT(1)'];    

plot(tv  ,co3tv (:,k)*1e6   ,sstr(2*k-1:2*k),'LineWidth',LineWidth);
plot(tv11,co3tv1(k,:)*1e6   ,sstr(2*k-1:2*k),'LineWidth',LineWidth);
plot(tv  ,f1+shtT *f2/1000,'r:','LineWidth',LineWidth);
plot(tv11,f1+shtT1*f2/1000,'r:','LineWidth',LineWidth);
plot(tv  ,f1+f2*ccdT /1e3,'k:','LineWidth',LineWidth);
plot(tv11,f1+f2*ccdT1/1e3,'k:','LineWidth',LineWidth);
end;

hold off;
set(gca,'FontSize',fs);
xlabel('Time (y)');
ylabel('[CO_3^{=}]');
%-------------------- right and top axes
% note: b1(3) <-> b1(4) exch. if yr(2) < yr(1)
b2 = [b1(1) b1(2) (b1(4)-f1)/f2 (b1(3)-f1)/f2];
ap = get(gca,'Position');
Hx  = axes('position',ap);
set(Hx,'Color','none');
set(Hx,'YAxisLocation','right');
set(Hx,'XAxisLocation','top');
set(gca,'XTickLabel',[]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fs);
Hx=ylabel('Sat. Horizon/CCD (km)');
set(Hx,'VerticalAlig','bottom','Rotation',-90);
axis(b2);
Hl=legend([pl(1:3)],'[CO_3^{2-}]','SH','CCD',2);
set(Hl,'FontSize',8);
end;

fcA1     = [fcA(1,:)' fcA(1,:)'];
fcI1     = [fcI(1,:)' fcI(1,:)'];
fcP1     = [fcP(1,:)' fcP(1,:)'];

%----------------------- fc -------------------%
% axis
axfc = [tv(1) tv(lt) -10. 100.];
%-------------- fc A --------------------%
cs = 'bgrmckybgrmckybgrmcky';
figure(9)
clf; 
box  on;
hold on;
for l=1:Ns
plot(tv  ,fcA (:,l)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
for l=1:Ns
plot(tv11,fcA1(l,:)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
%plot([tv(1) tv(lt)],[fccd fccd],'k');
hold off;
set(gca,'FontSize',fs);
axis(axfc);
set(gca,'YDIr','reverse');
set(gca,'XLim',axx);
title('Atlantic');
xlabel('Time (y)');
ylabel('CaCO_3 (wt %)');

clear dstr;
for l=1:Ns
dstr(l,:) = sprintf('%4.0f m',dsv(l));
end;
Hl=legend(dstr,4);
set(Hl,'FontSize',10);


%-------------- fc I --------------------%
figure(10)
clf; 
box  on;
hold on;
for l=1:Ns
plot(tv  ,fcI (:,l)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
for l=1:Ns
plot(tv11,fcI1(l,:)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
%plot([tv(1) tv(lt)],[fccd fccd],'k');
hold off;
set(gca,'FontSize',fs);
axis(axfc);
set(gca,'YDIr','reverse');
set(gca,'XLim',axx);
title('Indic');
xlabel('Time (y)');
ylabel('CaCO_3 (wt %)');

Hl=legend(dstr,4);
set(Hl,'FontSize',10);


%-------------- fc P --------------------%
figure(11)
clf; 
box  on;
hold on;
for l=1:Ns
plot(tv  ,fcP (:,l)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
for l=1:Ns
plot(tv11,fcP1(l,:)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
%plot([tv(1) tv(lt)],[fccd fccd],'k');
hold off;
set(gca,'FontSize',fs);
axis(axfc);
set(gca,'YDIr','reverse');
set(gca,'XLim',axx);
title('Pacific');
xlabel('Time (y)');
ylabel('CaCO_3 (wt %)');

Hl=legend(dstr,4);
set(Hl,'FontSize',10);


if(ftys)
fcT1     = [fcT(1,:)' fcT(1,:)'];
%-------------- fc T --------------------%
figure(12)
clf; 
box  on;
hold on;
for l=1:Ns
plot(tv  ,fcT (:,l)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
for l=1:Ns
plot(tv11,fcT1(l,:)*1e2,'--','Color',cs(l),'LineWidth',LineWidth);
end;
plot([tv(1) tv(lt)],[fccd fccd],'k');
hold off;
set(gca,'FontSize',fs);
axis(axfc);
set(gca,'YDIr','reverse');
set(gca,'XLim',axx);
title('Tethys');
xlabel('Time (y)');
ylabel('CaCO_3 (wt %)');


Hl=legend(dstr,4);
set(Hl,'FontSize',10);
end;

if(0)
d13fcA = (f13cA./fcA/Rst-1)*1e3;
d13fcI = (f13cI./fcI/Rst-1)*1e3;
d13fcP = (f13cP./fcP/Rst-1)*1e3;
d13fcT = (f13cT./fcT/Rst-1)*1e3;
%-------------- f13c A --------------------%
figure(13)
clf; 
plot(tv  ,d13fcA(:,1),'-','Color',cs(1),'LineWidth',LineWidth);
%plot(tv  ,f13cA(:,1),'-','Color',cs(1));
%plot(tv  ,m13calA(:,1),'-','Color',cs(1));
hold on;
for l=2:Ns
plot(tv  ,d13fcA(:,l),'--','Color',cs(l),'LineWidth',LineWidth);
plot(tv  ,d13fcI(:,l),'--','Color','r','LineWidth',LineWidth);
plot(tv  ,d13fcP(:,l),'--','Color','m','LineWidth',LineWidth);
plot(tv  ,d13fcT(:,l),'--','Color','g','LineWidth',LineWidth);
%plot(tv  ,f13cA(:,l),'--','Color',cs(l));
%plot(tv  ,m13calA(:,l),'--','Color',cs(l));
end;
hold off;
set(gca,'FontSize',fs);
%set(gca,'YDIr','reverse');
set(gca,'XLim',axx);
set(gca,'XScale','log');
title('Atlantic');
xlabel('Time (y)');
ylabel('^{13}f_c');

Hl=legend(dstr,4);
set(Hl,'FontSize',10);
end;

end; % sediments

%--------------Emission plot---------------------%
% edited: malte


if (BlFlag == 1)
    Emission_per_year = zeros(t0+1,tfinal-1);
    timespan = t0+1:tfinal-1;
    
    for i = t0+1:DTS
        Emission_per_year(i) = (CBl/1.e15)/DTS;  % conversion in Pg C
    end
    
    for i = DTS+1:tfinal-1
        Emission_per_year(i) = 0;
    end
    
    figure(55)
    %semilogx(timespan,Emission_per_year,'Linewidth',LineWidth)
    plot(timespan,Emission_per_year,'Linewidth',LineWidth)
    
    title('Carbon Emissions per year')
    xlabel('t [years]')
    ylabel('Emission Carbon [Pg C]')
    axis([t0 tfinal 0 2*((CBl/1.e15)/DTS)])
end

%--------------------------------------------%

%----------- Compare model vs. data --------------------%
%
% SALINITY Atl  Ind  Pac 
Sal =   [35.88 35.05 34.82 ...
         35.15 35.00 34.55 ...
         34.92 34.73 34.65 34.20]; 

% TOTAL CO2   Atl  Ind  Pac Atl  Ind  Pac
%Ctun =      [2.00 2.00 2.00 NaN NaN NaN 2.20 2.30 2.33 NaN];
Ctun =     [[2.068 2.006 2.006]-0.050 ...
             2.170 2.195 2.218        ...
             2.191 2.292 2.327 2.115-0.050];
Ctun = Ctun.*Soc./Sal;         
% TOTAL ALK
%Atun =      [2.35 2.30 2.30 NaN NaN NaN 2.33 2.42 2.42 NaN];
Atun =      [2.369 2.310 2.295 ...
             2.322 2.328 2.314 ...
             2.334 2.393 2.405 2.294];
Atun = Atun.*Soc./Sal;         
% PO4
%Ptun =      [0.50 0.50 0.50 NaN NaN NaN 2.30 2.30 2.30 NaN];
Ptun =      [0.347 0.351 0.481 ...
             1.431 1.779 2.159 ...
             1.521 2.423 2.657 1.421]/1.025; % (umol/l -> umol/kg)
Ptun = Ptun.*Soc./Sal;         
% O2
O2tun =     [0.2278 0.1971 0.2217 ...
             0.1821 0.1293 0.1184 ...
             0.2434 0.1596 0.1368 0.3097]; % (mol/m3)
% d13C
d13Ctun =   [2.5 2.5 2.5 ...
             0.9 0.5 0.5 ...
             0.8 0.4 0.0 1.6];
% CARBONATE ION
%cartun = [100 81 67]; ZW03
cartun = [101 80 70]; % calc.

cs = 'gggkkkrrrbgkr';
sy = 'ds^ds^ds^vppp';  
fs = 10;
lw = 1.5;
MS = 8;
figure(14)
clf;
%---- TCO2 TA
%subplot(221)
subplot(231)
plot(c(lt,1),a(lt,1),sy(1),'Color',cs(1),'LineW',lw,'MarkerSi',MS);
hold on;
plot(Ctun(1),Atun(1),sy(1),'Color',cs(1),'MarkerSi',MS, ...
    'MarkerFaceColor',cs(1),'MarkerEdgeColor','k');
for k=2:Nb
pl(k) = ...
plot(c(lt,k),a(lt,k),sy(k),'Color',cs(k),'LineW',lw,'MarkerSi',MS);
end;
for k=2:10
Pl(k) = ...
plot(Ctun(k),Atun(k),sy(k),'Color',cs(k),'MarkerSi',MS, ...
    'MarkerFaceColor',cs(k),'MarkerEdgeColor','k');
end;
hold off;
%axis([1.8 2.41 2.2 2.5]);
set(gca,'FontSize',fs);
xlabel('TCO_2 (mmol/kg)');
ylabel('TA (mmol/kg)');

legend([pl(4) Pl(4)],'Model','Obsrv',2);


if(1)
%---- PO4 d13C
%subplot(222)
subplot(232)
plot(p(lt,1)*1e3,d13c(lt,1),sy(1),'Color',cs(1),'LineW',lw,'MarkerSi',MS);
hold on;
plot(Ptun(1)    ,d13Ctun(1),sy(1),'Color',cs(1),'MarkerSi',MS, ...
    'MarkerFaceColor',cs(1),'MarkerEdgeColor','k');
for k=2:Nb
plot(p(lt,k)*1e3,d13c(lt,k),sy(k),'Color',cs(k),'LineW',lw,'MarkerSi',MS);
end;
for k=2:10
plot(Ptun(k)    ,d13Ctun(k),sy(k),'Color',cs(k),'MarkerSi',MS, ...
    'MarkerFaceColor',cs(k),'MarkerEdgeColor','k');
end;
hold off;
%axis([1.9 2.4 2.3 2.45]);
set(gca,'YDir','r');
set(gca,'FontSize',fs);
xlabel('PO4 (\mumol/kg)');
ylabel('\delta^{13}C (???)');

if(fdox)
%---- PO4 O2
%subplot(222)
subplot(233)
plot(p(lt,1)*1e3,dox(lt,1),sy(1),'Color',cs(1),'LineW',lw,'MarkerSi',MS);
hold on;
plot(Ptun(1)    ,O2tun(1),sy(1),'Color',cs(1),'MarkerSi',MS, ...
    'MarkerFaceColor',cs(1),'MarkerEdgeColor','k');
for k=2:Nb
plot(p(lt,k)*1e3,dox(lt,k),sy(k),'Color',cs(k),'LineW',lw,'MarkerSi',MS);
end;
for k=2:10
plot(Ptun(k)    ,O2tun(k),sy(k),'Color',cs(k),'MarkerSi',MS, ...
    'MarkerFaceColor',cs(k),'MarkerEdgeColor','k');
end;
hold off;
%axis([1.9 2.4 2.3 2.45]);
set(gca,'FontSize',fs);
xlabel('PO4 (\mumol/kg)');
ylabel('O_2 (mol/m3)');
end;


%if(0)
%subplot(223)
subplot(234)
plot(1,co3ALLf(1),sy(7),'Color',cs(7),'LineW',lw,'MarkerSi',MS);
hold on;
plot(1,cartun (1),sy(7),'Color',cs(7),'MarkerSi',MS, ...
    'MarkerFaceColor',cs(7),'MarkerEdgeColor','k');
for k=2:3
plot(k,co3ALLf(k),sy(k+6),'Color',cs(k+6),'LineW',lw,'MarkerSi',MS);
plot(k,cartun (k),sy(k+6),'Color',cs(k+6),'MarkerSi',MS, ...
'MarkerFaceColor',cs(k+6),'MarkerEdgeColor','k');
end;
if(ftys)
k=4;
plot(k,co3Tf     ,sy(k+9),'Color',cs(k+9),'LineW',lw,'MarkerSi',MS);
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XTickLabel',[]);
ylm = get(gca,'YLim');
text([1 2 3],[1 1 1]*ylm(1)-5,['Atl';'Ind';'Pac'],...
    'FontSize',fs,'Hor','c');
ylabel('Deep [CO_3^{2-}] (\mumol/kg)');

end;

%subplot(224)
subplot(236)
plot(0,0,'dk',0,0,'sk',0,0,'^k',0,0,'vk');
Hl=legend('Atl','Ind','Pac','HL',1);
if(ftys)
plot(0,0,'dk',0,0,'sk',0,0,'^k',0,0,'vk',0,0,'pk');
Hl=legend('Atl','Ind','Pac','HL','Tet',1);
end;
set(Hl,'FontSize',fs);
set(gca,'XTick',[],'YTick',[]);
axis([1 10 0 15]);
text(2,13,'Surf','Color','g','FontSize',fs)
text(2,11,'Intm','Color','k','FontSize',fs)
text(2,09,'Deep','Color','r','FontSize',fs)
text(2,07,'High','Color','b','FontSize',fs)

% inventory/d13Catm
MCstr   = sprintf('%5.0f',MCO);
d13Cstr = sprintf('%3.2f',d13Ca);
text(01,-1.5,['MC = ' MCstr ' (35760) Pg C']);
text(01,-3.0,['d13Ca = ' d13Cstr ' (-6.30) %^{_o}']);


if(logax == 1)
if(fsed)    
if(ftys)    
kf = [1 2 3 4 5 6 9 10 11 12];  
else 
kf = [1 2 3 4 5 6 9 10 11];  
end;    
else % fsed
kf = [1 2 3 4 5 6];  
end;
for k=kf
figure(k)
set(gca,'XScale','log');
end;
end;

end % end plotflag

% write parameter to file
if(parflag == 1)
fid = fopen('parms.tex','w');
fprintf(fid,  '%2.2f %2.2f %2.2f %2.2f fA\n',[fA3 fA(10)]);
fprintf(fid,  '%6.4e    TH\n',TH/1e6/(3600*24*365));
fprintf(fid,  '%2.2f           tA\n',tA);
fprintf(fid,  '%2.2f           tI\n',tI);
for i=1:length(mv)
fprintf(fid,  '%4.2f ',mv(i)/1e6/(3600*24*365)/3.8);
end;
fprintf(fid,  'mv\n');
for i=1:length(mhd)
fprintf(fid,  '%4.2f ',mhd(i)/1e6/(3600*24*365)/1.3);
end;
fprintf(fid,  'mhd\n');
for i=1:Nb
fprintf(fid,  '%2.2f ',TCv0(i));
end;
fprintf(fid,  'TCv0\n');
fprintf(fid,  '%6.4e   EPH\n',EPH/A(10));
fprintf(fid,  '%2.2f %2.2f %2.2f gp\n',gp(7:9));
fprintf(fid,  '%6.4e  fEPL\n',fEPL);
fprintf(fid,  '%6.4e rrain\n',rrain);
fprintf(fid,  '%6.4e    nu\n',nu);
fprintf(fid,  '%6.4e    eI\n',eI);
fprintf(fid,  '%2.2f          Ns\n',Ns);
fprintf(fid,  '%6.4e   FiN\n',FiN*Aoc);
fprintf(fid,  '%6.4e  frrf\n',frrf*1e2);
fprintf(fid,  '%6.4e    hs\n',hs);
fprintf(fid,  '%6.4e  phi0\n',phi0);
fprintf(fid,  '%6.4e   gam\n',gam);
fprintf(fid,  '%6.4e    Kd\n',Kd);
fprintf(fid,  '%6.4e    nc\n',nc);
fprintf(fid,  '%6.4e   nSi\n',nSi);
fprintf(fid,  '%6.4e   nCC\n',nCC);
fclose(fid);
end;

%-------------- print figures

% print -depsc -f1  L:\Matlab\box\PETM\eps\TCOT.eps
% print -depsc -f2  L:\Matlab\box\PETM\eps\TA.eps
% print -depsc -f3  L:\Matlab\box\PETM\eps\PO4.eps
% print -depsc -f4  L:\Matlab\box\PETM\eps\CO3.eps
% print -depsc -f6  L:\Matlab\box\PETM\eps\d13C.eps
% print -depsc -f7  L:\Matlab\box\PETM\eps\fcz.eps
% print -depsc -f14 L:\Matlab\box\PETM\eps\tune.eps

% print -dpng  -f14 L:\Matlab\box\PETM\eps\tune.png


%-------------- save output
% create a directory for the dirstr, the date is first saved
% in the working directory and then moved to the assigned directory

%uncomment if(0) to not save the results, also uncomment end statement!
%if(0)
%if(BlFlag == 2)

dirstr = ' ~/Documents/Studium/PHD_Hawaii/thesis/myLoscar/output/test'
if(parflag == 1)


%dirstr = ' dat\Modern\NoSed';
%dirstr = ' dat\Modern\NoWth';
%dirstr = ' dat\Modern\B2D3';
    
%dirstr = ' dat\B1D1BL0';
%dirstr = ' dat\B1D1BL2SW';

%dirstr = ' dat\B2D3BL2SW';       % <---------
%dirstr = ' dat\B2D3BL2SW\new';   % <---------   
%dirstr = ' dat\B2D3BL2SW\Hpco2';   
%dirstr = ' dat\B2D3BLXPG';
%dirstr = ' dat\B2D3BL2ky10';


%dirstr = ' dat\B2D3BL2ky10\Atm2000';
%dirstr = ' dat\B2D3BL2ky10\Atm6000';

%dirstr = ' dat\B3D3BL0';
%dirstr = ' dat\B3D3BL2SW';

%--- fossil fuel, bio calc
%dirstr = ' dat\FFYnk';
%dirstr = ' dat\FFbcYnk';

%--- fossil fuel NatYEarth

%dirstr = ' dat\Modern\NatYEarth\SiYes';
%dirstr = ' dat\Modern\NatYEarth\SiNo';
%dirstr = ' dat\Modern\NatYEarth\TS45';
%dirstr = ' dat\Modern\NatYEarth\TS15';

%--- Emss pH
%dirstr = ' dat\Emss\5000_0500';
%dirstr = ' dat\Emss\3500_0500';
%dirstr = ' dat\Emss\3500_0500\T15';
%dirstr = ' dat\Emss\3500_0500\T45';
%dirstr = ' dat\Emss\1500_0280';
%dirstr = ' dat\Emss\1000_0500';


%dirstr = ' dat\Emss\SlowCirc\3500_0500';
%dirstr = ' dat\Emss\SlowCirc\3500_0500SC';


dsvp = dsv';
save vars.tex  Nb Ns d13CBl -ASCII -DOUBLE -TABS;
save dsv.DAT   dsvp  -ASCII -DOUBLE -TABS;
save tv.DAT    tv    -ASCII -DOUBLE -TABS;
save c.DAT     c     -ASCII -DOUBLE -TABS;
save a.DAT     a     -ASCII -DOUBLE -TABS;
save kkv.DAT   kkv   -ASCII -DOUBLE -TABS;
save pco2t.DAT pco2t -ASCII -DOUBLE -TABS;
save pco2v.DAT pco2v -ASCII -DOUBLE -TABS;
save co3tv.DAT co3tv -ASCII -DOUBLE -TABS;
save phtv.DAT  phtv  -ASCII -DOUBLE -TABS;
save d13c.DAT  d13c  -ASCII -DOUBLE -TABS;
save d13CA.DAT d13C  -ASCII -DOUBLE -TABS;

save TCvt.DAT  TCvt  -ASCII -DOUBLE -TABS;

if(0)
save THt.DAT   THt   -ASCII -DOUBLE -TABS;
end

save OmegCS.DAT omegCSvt  -ASCII -DOUBLE -TABS;
save OmegAS.DAT omegASvt  -ASCII -DOUBLE -TABS;

save fDTS.DAT DTS DTS2 ts3 DTS3 DTS4 -ASCII -DOUBLE -TABS;
save RlsCtv.DAT RlsCtv -ASCII -DOUBLE -TABS;
save dox.DAT    dox    -ASCII -DOUBLE -TABS;

save shtA.DAT  shtA  -ASCII -DOUBLE -TABS;
save shtI.DAT  shtI  -ASCII -DOUBLE -TABS;
save shtP.DAT  shtP  -ASCII -DOUBLE -TABS;
if(ftys)
save shtT.DAT  shtT  -ASCII -DOUBLE -TABS;
end;

if(fsed)
save fcA.DAT   fcA   -ASCII -DOUBLE -TABS;
save fcI.DAT   fcI   -ASCII -DOUBLE -TABS;
save fcP.DAT   fcP   -ASCII -DOUBLE -TABS;

save ccdA.DAT  ccdA  -ASCII -DOUBLE -TABS;
save ccdI.DAT  ccdI  -ASCII -DOUBLE -TABS;
save ccdP.DAT  ccdP  -ASCII -DOUBLE -TABS;
if(ftys)
save fcT.DAT   fcT   -ASCII -DOUBLE -TABS;
save ccdT.DAT  ccdT  -ASCII -DOUBLE -TABS;
end;
end

if(plotflag)
if(BlFlag) % malte
save emission_time.DAT timespan -ASCII -DOUBLE -TABS;
save emission_carbon.DAT Emission_per_year -ASCII -DOUBLE -TABS;
end
end

if(fsaveendstate) % save end state of all variables Y(:,1:130) 
save Y0_pRef.DAT    Y_end   -ASCII  -DOUBLE -TABS;    
end





% save mass balance
% DXs = Change of [Mass_Carbon Alkalinity CCD_depth] between initial cond. and last timestep
save Massbalance.DAT DXs -ASCII -DOUBLE -TABS;
end

if(contourflag)
save X.DAT      X   -ASCII -DOUBLE -TABS;    
end

if any([parflag contourflag])
cmdstr(1,:) = 'mv parms.tex';
cmdstr(2,:) = 'mv  vars.tex';
cmdstr(3,:) = 'mv     *.DAT';

for k=1:3
  dostr(k,:)     = [cmdstr(k,:) dirstr]; % moving is done here
 [status,result] = dos(dostr(k,:))
end;    
 
%         dostr2 = ['copy' dirstr '\YPE10Sed.DAT'];
%[status,result] = dos(dostr2)
end % moving



end; % myflag


return;

