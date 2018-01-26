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
%    0
%
% updates:
% 05/19/10 44Ca tracer added but only for CAflag=0 &
%          fsed && ftys && fdox =1
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
%--------------------------------------

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
Mgm = 53.0e-3; % 53.0 (mol/kg) Magnesium modern
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
        % known to Loscar10 and LoscarDif10
        %========================================
        
        global Aoc rho Nb V A TH0 TH TS mv mhd kasv kkv TCv Sv Pv gp tA tI ...
            phflag k1k2flag EPH fEPL rrain ep ept ec ect REDPC REDNC REDO2C eI ...
            Ns asvA asvI asvP zv m2kg m2kgca rhos frrf frrfca phic hs FiN ...
            phivtA phivtI phivtP phi0 phi1 gam ...
            phiiA phiiI phiiP phiiAca phiiIca phiiPca phiiT phiiTca dsv ...
            Fpr frain rsedv rburvtA rburvtI rburvtP dissvtA ...
            dissvtI dissvtP FprtA FprtI FprtP it co3satv ...
            Kd KS cst dsflag nc dYflag nu fsed ...
            Fpr13tA Fpr13tI Fpr13tP Fpr44tA Fpr44tI Fpr44tP Fint Fin13t Fin44t Fpr44tT FSit FSi13t FSi44t...
            fc0A fc0I fc0P fca0A fca0I fca0P f13c0A f13c0I f13c0P ...
            dissv13tA dissv13tI dissv13tP ...
            dissv44tA dissv44tI dissv44tP...
            co3s0 as zs0 klid nli Rst epsp epspca Rin FiN13 ...
            BlFlag CBl RBl kb FVC FVC13 Rvc pCSi nSi nCC Fkg Fkg13 ...
            ftys nOC TT asvT kliT fc0T fca0T f13c0T rsedvtT dissvtT dissv13tT dissv44tT ...
            FprtT Fpr13tT Fpr44tT phivtT fcon swcon ...
            TCv0 TCvt ntL ntH fsh fshI fshP fdpv fshT nshT ...
            ffflag tem em RlsCtv ...
            fdox KMMOX vask DTS DTS2 ts3 DTS3 DTS4 ...
            k1st tfinal TDflag omegCSvt omegASvt ...
            THt mv0 mhd0 oxA CAvflag FSichck Finchck kspCHCK...
            CHECK1 CHECK2 CHECK3 CHKFin CHKFSi CHCKbioL CHCKbioI CHCKbioD ...
            CHKcarB1 CHKcarB2 CHKcarB3 CHKcar1 CHKcar2 CHKcar3 CHKcar4 CAv g Cinpc Cadv slj avCALC AJDE Voc...
            RinCA RcasAcheck RsAcheck f13cvAchck f44cavAchck fcvAchck CHKcarC1 CHKcarC2 CHKcarC3...
            CHKcarD1 CHKcarD2 CHKcarD3 CHKcarD4 Rbchck Rbcachck...
            rburvtAca rburvtIca rburvtPca rsedvtTca phivtAca phivtIca phivtPca...
            FprtAca FprtIca FprtPca dissvtTca phivtTca FprtTca EPLv EXLv ECLv...
            fwcv fwsv...
            fwsv1 fwgv fmgv pco2gca tgc kgc fbgv fmcv FSi0 fbcv pco20 Fkgs fbch fbbv...
            frkc fekc fdkc flakc FERT ACT RUN gamma fGG FiN0 Rkg LTflag acv ybbv ybv ycv epspv wcvtA wcvtI wcvtP wcvtT...
            dissvtAca dissvtIca dissvtPca dissvtTca epspcaV inorgF epspF biopF;
        
        
        
        
        
        % Loscar-GEOCARB long-term stuff
        LTflag = 1; % Flag to switch Long-Term, LOSCAR-GEOCARB run
        
        kgc=0;
        pco20=300;  %default 300 ppm
        FERT=0.40;
        ACT=0.05;
                 [ybbv,ybv,ycv,aav,bbv,ccv,fGG,fgkc,frkc,fekc,fdkc,flakc,fbch,fbbv,fwcv,fwsv1,fwsv2,fbcv,fbgv,rco2,pco2gca,fwgv,fmcv,fmgv,ggc,cgc,dg,dc]=myGeoCarbMODULE1(1);
        
        % Initial year in Ma
        tgc=59;
        if(LTflag)
            % The length of the run in years
            runtime=0.5e6;
        else
            runtime=4.e5;
        end
        
        %%%%Payne data set for d44Ca and d13C end-Permian
        pd =  csvread('PayneData.csv');
        stel = pd(:,1); % stratigraphic elevation in meters
        d44cap = pd(:,2);
        d13cp = pd(:,4);
        
        depth = [0 160]; % used for calculating time
        for i=1:length(stel)
            time(i) = interp1(depth, [0 1] , stel(i));
        end
        
        
        % +++++ Edit this line to add correct path to solver ! +++++ %
        % *****LINES CHANGED TO GET PETM set-up with Ca= 10.3; 292 540 659 and save
        %         addpath('C:\Users\Nemanja\Dropbox\geocarb\Loscar d44Ca\myode');
        addpath('myode');
        CBl = 12000.e15; % Blast end=Perm      13200, 43200, 10000
        d13CBl = -5.;     % Blast d13C -55   -60   -34   n-50
        % set DTS - duration of the initial carbon input pulse
        % this for BlFlag = 2;
        DTS = 0.5e5;    % 6e3 years for PETM,1e5 and 0.5e5 end-Permian
        
        
        nSi0  = 0.2;          % 0.2
        nC0   = 0.4;                % 0.4
        nSi   = 0.1;          % 0.2 0.3
        nCC   = 0.2;                % 0.4 0.3 1.0 0.5
        plotflag = 1; % plot results
        CAvflag  =2; % 2:calcium is a tracer and changes,1:calcium is a tracer but stays constant
        % 0: Ca is not a tracer
        savf     = 0; % save end state
        loadf    = 1; % load initial steady-state
        
        BlFlag   = 0; % Blast 1: shot0, 2: cont release
        inorgF   = 0; % inorganic CaCO3 precip. should be 1 when runninig Permian
        % for any other time period should be 0
        
        epspF = 0;   %Ca fractionation depends on CO32-;
        biopF = 0;      %increases intial bio pump
        
        fcon     = 3; % 1: NADW, 2: NPDW, 3: SO
        swcon    = 0; % 1: SO -> NP -> SO switch
        
        bath     = 2; % 1,2,3 bathymetry; 4-made up bathym for end perm
        dsflag   = 3; % 1,2,3 dissolution parameter
        
        ftys     = 1; % Tethys
        fsed     = 1; % include sediments
        parflag  = 0; % write parameter to file
        
        ffflag   = 0; % Anthropogenic CO2 (fossil fuel)
        fdox     = 1; % include dissolved oxygen
        
        TDflag   = 0; % Temp sens to doubling CO2
        
        ccdrun   = 0;   % parameter run: CBl, oxA
        oxA      = 0.4; % fraction released in deep Atl
        
        
        disp('    ');
        disp('@==================== RUN Loscar ====================@');
        disp('    ');
        load_Blast_fcon_bathym_diss_TDflag = ...
            sprintf('  %d    %d    %d     %d      %d     %d',...
            loadf,BlFlag,fcon,bath,dsflag,TDflag);
        load_Blast_fcon_bathym_diss_TDflag
        
        
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
        
        
        Voc  = 1.2918235e18;    % (m3) volume ocean
        Aoc  = 3.49e14;         % (m2) area ocean
        Hav  = Voc/Aoc;         % (m)  average depth
        
        rho  = 1.025e3;         % kg/m3 (1.025: Toggweiler)
        rhos = 2.50e3;          % kg/m3 sed. density 2.50
        m2kg = 100/1e3;         % mol C -> kg CaCO3
        m2kgca = 100/1e3;         % mol Ca -> kg CaCO3
        
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
        TCv0 = [TC3(1)*on3 TC3(2)*on3 TC3(3)*on3 2.0];
        CA3  = [10.3 10.3 10.3]*1e-3; % mol/kg Ca of boxes
        CAv0 = [CA3(1)*on3 CA3(2)*on3 CA3(3)*on3 10.3];
        if(ftys)
            TC3  = [25. 16. 12.];     % (degC) temp. of boxes
            TCv0 = [TC3(1)*on3 TC3(2)*on3 TC3(3)*on3 12.0+0];
            TCT  = [18. 14. 12.+0];   % 18/25 16/14 12
            TCv0 = [TCv0 TCT]+0;
            %CA3  = [20 20 20]*1e-3; % mol/kg Ca of boxes
            if(inorgF)
                CA3  = [10 10 10]*1e-3; % mol/kg Ca of boxes
            else
                CA3  = [20 20 20]*1e-3; % mol/kg Ca of boxes
            end
            CAv0 = [CA3(1)*on3 CA3(2)*on3 CA3(3)*on3 CA3];
        end;
        TCv   = TCv0;
        Soc   = 34.72;             % Sal whole ocean
        Sv    = onV*Soc;           % Salinity vector
        CAv= CAv0;
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
        end;
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
        
        if(biopF)
            EPH   = 2*1.8*A(10); % (mol/y) 1.8 1.6 H Export, mol C/m2/y*A = mol/y
            fEPL  = 0.95;      % 0.80 0.9 LL utilization
        else
            EPH   = 1.8*A(10); % (mol/y) 1.8 1.6 H Export, mol C/m2/y*A = mol/y
            fEPL  = 0.80;      % 0.80 0.9 LL utilization
        end
        rrain = 6.1;       % 6.1 6.2 6.7 export rain ratio (Corg/CaCO3)
        % 5.9(?) 6.1(2,3)
        
        if(ftys)
            if(inorgF)
                rrain = 6.7*1e20;       % 6.7 7 4.2 export rain ratio (Corg/CaCO3)
            else
                rrain = 6.7;
            end
        end;               % 8.0(1,1) 6.3(1,3) 6.6/6.2/6.0?(2,3)
        nu    = 0.31;      % 0.31 water column dissolution
        
        eI    = 0.78;      % 0.78 0.8 fraction EPL, remineralized in I boxes
        if(~fsed)
            nu    = 0.;
        end;
        
        % fraction EPH, remineralized in deep A,I,P boxes
        gp      = 0.*ones(1,Nb);
        gp(7:9) = [.3 .3 .4];   % .3 .3 .4
        %gp(7:9) = [1   1  1]/3; % .7 .3 0
        
        
        %============== silicate weathering: volc degass
        pRef  = 280.;        % uatm, weathering  ref  280
        pCSi  = 280.;        % uatm, std-stt atm pCO2 280
        if(ftys)
            if(~inorgF)
                pRef  = 0500.*1.;    % uatm, weathering ref   500  574 350  750 /2
                pCSi = 1000.*1;     % uatm, std-stt atm pCO2 560 1000 700 1000 /2 887.8625
            else
                
                if(biopF)
                    pRef  = 0157.25*1.;    % uatm, weathering ref   370 500  574 350  750 /2
                    pCSi = 850.*1;     % uatm, std-stt atm pCO2 560 1000 700 1000 /2 887.8625
                else    
                    pRef  = 0157.25*1.;    % uatm, weathering ref   500  574 350  750 /2
                    pCSi = 850.*1;     % uatm, std-stt atm pCO2 2000
                end
                
            end
            if(LTflag)
                pCSi  = pco2gca(tgc);
            end;
        end;
        if(inorgF)
            FVC   = 1*3.e12/Aoc;  % 5.e12 mol C, degassing /m2/y @280 uatm
        else
            FVC   = 1*5.e12/Aoc;  % 5.e12 mol C, degassing /m2/y @280 uatm
        end
        
        %nSi   = 1.5;          % 0.2 0.3
        FVC   = FVC*(pCSi/pRef)^nSi0; % initial 9.2071
        if(LTflag)
            FVC   = (fmcv(tgc)*1.e12)/Aoc;
            FSi0   = ((fbcv(tgc)-fwcv(tgc))*1.e12)/Aoc;
        end;
        
        %============== kerogen oxidation
        Fkg   = 1*09.e12/Aoc;        % mol C    /m2/y 09
        Fkgs  = 1*09.e12/Aoc;
        if(ftys)
            Fkg   =   05.e12/Aoc;        % mol C    /m2/y 05 3.5425
            Fkgs  = 1*05.e12/Aoc;
            if(LTflag)
                Fkg   =   (fbgv(tgc)*1.e12)/Aoc;        % mol C    /m2/y 05
                Fkgs  =((fwgv(tgc)+fmgv(tgc))*1.e12)/Aoc;
            end
        end
        
        %============== CaCO3 in-flux ===============%
        %
        if(inorgF)
            FiN   = 1.0*5.e12/Aoc;      % 12e12 mol C    /m2/y riverine flux 1.3
        else
            FiN   = 1.0*12.e12/Aoc;      % 12e12 mol C    /m2/y riverine flux 1.3
        end
        
        %         nCC   = 0.40;                % 0.4 0.3 1.0 0.5
        FiN   = FiN*(pCSi/pRef)^nC0; % 12.7732
        Fpr   = 3.*FiN;              % mol C    /m2/y production 3.6
        if(LTflag)
            FiN   = (fwcv(tgc)*1.e12)/Aoc; % 12.7732
            FiN0  = (fwcv(tgc)*1.e12)/Aoc; % 12.7732
        end
        
        
        % rain of 'remainder'
        %frrf  = 1.5*1.15*0.180;     %  g/cm2/ky remainder 1.5*1.15
        frrf  = 0.35;                %  g/cm2/ky remainder .311
        frrf  = frrf*1e4/1e3/1e3;    % -> kg/ m2/ y
        frain = Fpr*m2kg/(Fpr*m2kg+frrf);
        
        frrfca  = 0.35;                %  g/cm2/ky remainder .311
        frrfca  = frrfca*1e4/1e3/1e3;    % -> kg/ m2/ y
        
        %======= Carbon-13
        Rst    = 0.011;   % 13C: R standard (value irrelevant)
        epsp   = -27.;    % -27 fractionation Corg
        if(ftys && ~LTflag)
            epsp   = -33.;    % -33 fractionation Corg
        end;
        if(ftys && LTflag)
            epsp   = -acv(tgc);   % fractionation Corg
        end;
        d13Cin = +2.0;    % d13C of riverine flux 3.0 2.0 2.6
        Rin    = Rst*(d13Cin/1e3+1);
        FiN13  = Rin*FiN; % mol C    /m2/y riverine flux
        % silicate weathering: volc degass
        d13Cvc = -3.0;    % d13C -5 +0.3 -0.7 +2.0
        if(ftys)
            d13Cvc = -5.0;    % d13C -5 +0.3 -0.7 +2.0
        end;
        Rvc    = Rst*(d13Cvc/1e3+1);
        FVC13  = Rvc*FVC; % mol C    /m2/y
        % kerogen oxidation
        d13Ckg = -22.3;   % d13C -22.3 -28.3
        Rkg    = Rst*(d13Ckg/1e3+1);
        %         if(LTflag)
        Fkg13  = Rkg*Fkgs; % mol C    /m2/y
        %         else
        %             Fkg13  = Rkg*Fkg;
        %         end
        
        %======= Calcium-44
        Rstca    = 0.0208;   % 13C: R standard (value irrelevant)
        epspca   = -1.4;    % -1.3 fractionation between seawater and carbonate minerals
        if(ftys)
            if(epspF)
                if(biopF)
                    epspca = -0.9249;
                else
                    epspca   = -0.9983; %1.0111, 1.0305
                end
            else
                epspca   = -1.4;    % -1.3 fractionation between seawater and carbonate minerals
            end
        end;
        d44CAin = -0.6;    % d44Ca of riverine flux 3.0 2.0 2.6
        RinCA    = Rstca*(d44CAin/1e3+1);
        FiN44ca  = RinCA*FiN; % mol C    /m2/y riverine flux
        % silicate weathering: volc degass
        %d44CAvc = -3.0;    % d13C -5 +0.3 -0.7 +2.0
        if(ftys)
            %d44CAvc = -5.0;    % d13C -5 +0.3 -0.7 +2.0
        end;
        % RvcCA    = Rstca*(d44CAvc/1e3+1);
        % FVC13ca  = RvcCA*FVC; % mol C    /m2/y
        % kerogen oxidation
        % d44CAkg = -22.3;   % d13C -22.3 -28.3
        % Rkgca    = Rstca*(d44CAkg/1e3+1);
        % Fkg13ca  = Rkgca*Fkg; % mol C    /m2/y
        
        
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
            fsh    = 1.00;            %  increase shelf rain
            fshI   = 1;
            fshP   = 1;
            nshT   = 1.00;            %
            if(ftys)
                if(inorgF)
                    fsh    = 1;             % 8.9
                    fshI   = 1;            % 17.2
                    fshP   = 1;            % 33.3
                    nshT   = 0.27;            % 0.6(1,1) 0.35/0.6(2,3).3
                else
                    fsh    = 4.5;             % 8.9
                    fshI   = 4.5;            % 17.2
                    fshP   = 4.5;            % 33.3
                    nshT   = 0.4;            % 0.6(1,1) 0.35/0.6(2,3).3
                end
            end;
            
            %----------------- sediment boxes (bathymetry) ------------%
            %
            if    (bath == 1)
                dsv   = [ .1 .6 1.5 2.5 3.5 4.5 5.5 6.5]*1000;
                asvA  = [ 7.0297  5.1729  4.2988  8.5975 19.3421 32.4777 22.3425  0.7388]/100.;
                asvI  = [ 3.5710  2.6844  3.5792 10.0293 25.2598 36.6442 16.9915  1.2407]/100.;
                asvP  = [ 1.6358  2.5901  3.2590  6.8744 21.8550 35.0822 26.9567  1.7468]/100.;
                if(ftys)
                    asvT  = [16.3934 16.3934 16.3934 20.4918 20.4918  8.1967  1.6393  0.0001]/100.;
                end;
            elseif(bath == 2)
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
            elseif(bath == 3)
                dsv   = [.1 .6 1 1.5 2 2.5 2.75 3 3.25 3.5 3.75 4 ...
                    4.25 4.5 4.75 5. 5.25 5.5 6.0 6.5]*1000; % depth (m)
                % area fraction A I P
                asvA  = [07.0297  5.1729  1.9106  2.3882  4.2988  4.2988  4.8355  4.8355  4.8355  4.8355  8.1194  8.1194  8.1194  8.1194  5.5856  5.5856  5.5856  5.5856  0.3694  0.3694]/100;
                asvI  = [03.5710  2.6844  1.5907  1.9884  5.0146  5.0146  6.3149  6.3149  6.3149  6.3149  9.1610  9.1610  9.1610  9.1610  4.2479  4.2479  4.2479  4.2479  0.6204  0.6204]/100;
                asvP  = [01.6358  2.5901  1.4484  1.8105  3.4372  3.4372  5.4638  5.4638  5.4638  5.4638  8.7705  8.7705  8.7705  8.7705  6.7392  6.7392  6.7392  6.7392  0.8734  0.8734]/100;
                if(ftys)
                    asvT  = [16.3934 16.3934  7.2860  9.1075 10.2459 10.2459  5.1229  5.1229  5.1229  5.1229  2.0492  2.0492  2.0492  2.0492  0.4098  0.4098  0.4098  0.4098  0.0001  0.0001]/100;
                end;
                
            elseif(bath == 4)
                % 03/31/06, 2x2\deg Bice
                dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000; % depth (m)
                % area fraction A I P
                asvA  = [1.1407*2 10.0216*2  8.7160  6.3724  4.7915  3.4973 12.2935 ...
                    10.7438  8.4983 11.4134 11.0948-1.1407 11.4167-10.0216 1e-6]/100;
                asvI  = [0.2501*2  5.5345*2  5.5145  8.9550  4.6283  6.5361 11.7221 ...
                    12.6050 14.6295 13.8384  8.0807-0.2501  7.7057-5.5345 1e-6]/100;
                asvP  = [0.1673*2  2.8333*2  3.0599  2.5389  1.4218  5.0153  9.8023 ...
                    14.0117 10.1975 20.0019 13.7155-0.1673 17.2346-2.8333 1e-6]/100;
                asvT  = [7.0534*3 46.5363 22.4068-2*7.0534  7.1501  2.4261  3.9946  2.7063 ...
                    0.8532  0.9814  3.0802  2.3583  0.4533 1e-6]/100;
                
            end;
            
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
                
                if(ftys)
                    if(inorgF)
                        Ca      = 10.00e-3;    % Modern 10.3 PETM 20.0
                    else
                        Ca      = 20.00e-3;    % Modern 10.3 PETM 20.0
                    end
                    Mg      = 30.0e-3;    % 53.0 30.0
                end;
                
                if(CAvflag == 0)
                    Tdv = TCv(klid);
                    Sdv =  Sv(klid);
                    for k=1:Ns
                        [kspc(k),x] = ...
                            kspfun(Tdv(k),Sdv(k),zsatv(k)/10.,Ca,Mg);
                    end;
                    co3satv  = kspc/Ca;
                elseif(CAvflag == 1)
                    Tdv = TCv(klid);
                    Sdv =  Sv(klid);
                    for k=1:Ns
                        [kspc(k),x] = ...
                            kspfun(Tdv(k),Sdv(k),zsatv(k)/10.,Ca,Mg);
                    end;
                    co3satv  = kspc/Ca;
                else
                    % Temp for co3sat corrected 07/13/06
                    Tdv = TCv(klid);
                    Sdv =  Sv(klid);
                    Cadv = CAv(klid);
                    CALCv=Ca;
                    clear kspc;
                    for k=1:Ns
                        [kspc(k),x] = ...
                            kspfun(Tdv(k),Sdv(k),zsatv(k)/10.,CALCv,Mg);
                    end;
                    co3satv  = kspc./(CALCv*1e-3);
                    
                end
            end;% satflag
            kspc
            
            
            %-------------- Porosity --------------------%
            %phic    = 0.78;      % porosity 0.75 0.78
            if(phic) % do NOT use exist! phic does exist (global!)
                phiiA   = ones(1,Ns)*phic;
                phiiI   = ones(1,Ns)*phic;
                phiiP   = ones(1,Ns)*phic;
                phiiAca   = ones(1,Ns)*phic;
                phiiIca   = ones(1,Ns)*phic;
                phiiPca   = ones(1,Ns)*phic;
                if(ftys)
                    phiiT   = ones(1,Ns)*phic;
                    phiiTca   = ones(1,Ns)*phic;
                    
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
            fca0A = 0.46*ones(1,Ns);
            fca0I = 0.46*ones(1,Ns);
            fca0P = 0.46*ones(1,Ns);
            if(ftys)
                fc0T = 0.46*ones(1,Ns);
                fca0T = 0.46*ones(1,Ns);
            end;
            % calc initial phi
            if(isempty(phic)) % phi = phi(fc)
                FF    = (phi1-phi0)/(1-phi1);
                phiiA = (phi0+FF*fc0A)./(1+FF*fc0A);
                phiiI = (phi0+FF*fc0I)./(1+FF*fc0I);
                phiiP = (phi0+FF*fc0P)./(1+FF*fc0P);
                
                FFca    = (phi1-phi0)/(1-phi1);
                phiiAca = (phi0+FFca*fca0A)./(1+FFca*fca0A);
                phiiIca = (phi0+FFca*fca0I)./(1+FFca*fca0I);
                phiiPca = (phi0+FFca*fca0P)./(1+FFca*fca0P);
                if(ftys)
                    phiiT = (phi0+FF*fc0T)./(1+FF*fc0T);
                    phiiTca = (phi0+FFca*fca0T)./(1+FFca*fca0T);
                end;
            end;
            % calc initial calcite mass
            mc0vA  = (fc0A.*rhos.*(1-phiiA));
            mc0vI  = (fc0I.*rhos.*(1-phiiI));
            mc0vP  = (fc0P.*rhos.*(1-phiiP));
            mca0vA  = (fca0A.*rhos.*(1-phiiAca));
            mca0vI  = (fca0I.*rhos.*(1-phiiIca));
            mca0vP  = (fca0P.*rhos.*(1-phiiPca));
            if(ftys)
                mc0vT  = (fc0T.*rhos.*(1-phiiT));
                mca0vT  = (fca0T.*rhos.*(1-phiiTca));
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
            
            %====== Calcium-44
            m44ca0vA = RinCA*mca0vA;     % -> kg CaCO3/m3 Rin
            m44ca0vI = RinCA*mca0vI;     % -> kg CaCO3/m3 Rin
            m44ca0vP = RinCA*mca0vP;     % -> kg CaCO3/m3 Rin
            f44ca0A  = fca0A.*m44ca0vA./mca0vA;
            f44ca0I  = fca0I.*m44ca0vI./mca0vI;
            f44ca0P  = fca0P.*m44ca0vP./mca0vP;
            if(ftys)
                m44ca0vT = RinCA*mca0vT;     % -> kg CaCO3/m3 Rin
                f44ca0T  = fca0T.*m44ca0vT./mca0vT;
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
        if(biopF)
            p0   = [4.20*onV]*1e-3;   % mmol/kg (PO4 at t=0)
        else
            p0   = [2.50*onV]*1e-3;   % mmol/kg (PO4 at t=0)
        end
        p0   = p0*0.87;
        if(ftys == 1)
            if(inorgF)
                CA0  = [10.00*onV];      % PETM 20
            else
                CA0  = [20.00*onV];      % PETM 20
            end
        else
            CA0  = [10.30*onV];
        end;
        %p0   = [1:1:10]*.25*1e-3;
        
        if(fdox)
            dox0 = [0.20*onV];         % mol/m3 (O2 t=0)
            %surface
            if(CAvflag>0)
                for k=kkv
                    [x,x,x,x,x,dox0(k)] = ...
                        dafunPECA(1,1,TCv(k),Sv(k),Pv(k),1,1);
                end;
            else
                for k=kkv
                    [x,x,x,x,x,dox0(k)] = ...
                        dafunPECA(1,1,TCv(k),Sv(k),Pv(k),1,1);
                end;
            end
        end;
        
        %====== Carbon-13
        d13c0 = [2.35*on3 0.5*on3 0.67*on3 1.63];
        if(ftys)
            d13c0 = [[d13c0] 2.35 0.5 0.67];
        end;
        R0    = (d13c0/1e3+1.)*Rst;
        cc0   = R0.*c0;
        
        %====== Calcium-44
        d44ca0 = [2.35*on3 0.5*on3 0.67*on3 1.63];
        if(ftys)
            d44ca0 = [[d44ca0] 2.35 0.5 0.67];
        end;
        Rca0    = (d44ca0/1e3+1.)*Rstca;
        cca0   = Rca0.*CA0;
        
        c0    = c0  *1e-3.*rho; % mol/m3
        a0    = a0  *1e-3.*rho; % mol/m3
        p0    = p0  *1e-3.*rho; % mol/m3
        CA0   = CA0 *1e-3.*rho;
        cc0   = cc0 *1e-3.*rho; % mol/m3
        cca0   = cca0 *1e-3.*rho; % mol/m3
        if(LTflag)
            C0 = pco2gca(tgc)*2.2e15/12/Aoc; % 280
            % (mol/m2) atm. CO2 inventory / m2
            % 1 ppmv = 2.2 Pg C
        else
            C0 = 280*2.2e15/12/Aoc; % 280
            % (mol/m2) atm. CO2 inventory / m2
            % 1 ppmv = 2.2 Pg C
        end
        %====== Carbon-13
        d13C0 = -6.45;
        CC0   = C0.*(d13C0/1e3+1.)*Rst;
        
        %====== Calcium-44
        d44CA0 = -6.45; %%%%%%%%%%%%% PROLLY DONT NEED THIS CUZ THIS IS ATMOSPH
        CAC0   = C0.*(d44CA0/1e3+1.)*Rst;
        
        %Y0 = [c0 a0 p0 C0 mc0vA]';
        %Y0 = [c0 a0 p0 C0 mc0vA mc0vI mc0vP]';
        %      Nb 2Nb 3Nb 4Nb +1  +2 +2+Ns +2+2Ns +2+3Ns
        %      10  20  30  40 41  42    58     74     90
        %Y0 = [c0  a0  p0 cc0 C0 CC0 mc0vA  mc0vI  mc0vP]';
        
        if(fsed) %<<<<<<<<<<<<<<<<<<<<<<< fsed
            if(fdox) % fdox
                %     Nb 2Nb 3Nb 4Nb  5Nb +1  +2 +2+Ns  +2+2Ns   +2+3Ns
                %     10  20  30   40  50 51  52    65      78       91
                if(CAvflag == 0)
                    Y0 = [c0  a0  p0 dox0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                        f13c0A  f13c0I   f13c0P]';
                else
                    Y0 = [c0  a0  p0 CA0 dox0 cca0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                        f13c0A  f13c0I   f13c0P fca0A fca0I fca0P...
                        f44ca0A f44ca0I f44ca0P]';
                    %                          +2+4Ns  +2+5Ns   +2+6Ns
                    %                             104    117       130
                end;
            else % fdox
                %     Nb 2Nb 3Nb 4Nb +1  +2 +2+Ns  +2+2Ns   +2+3Ns
                %     10  20  30  40 41  42    55      68       81
                if(CAvflag == 0)
                    Y0 = [c0  a0  p0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                        f13c0A  f13c0I   f13c0P ...
                        f44ca0A f44ca0I f44ca0P]';
                else
                    Y0 = [c0  a0  p0 CA0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                        f13c0A  f13c0I   f13c0P ...
                        f44ca0A f44ca0I f44ca0P]';
                    %                          +2+4Ns  +2+5Ns   +2+6Ns
                    %                              94     107      120
                end;
            end; % fdox
            if(ftys) % ftys
                if(fdox) % fdox
                    if(CAvflag == 0)
                        Y0 = [c0  a0  p0 dox0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                            f13c0A  f13c0I   f13c0P ...
                            fc0T  f13c0T]';
                        
                    else
                        Y0 = [c0  a0  p0 CA0 dox0 cca0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                            f13c0A  f13c0I   f13c0P fca0A fca0I fca0P...
                            f44ca0A f44ca0I f44ca0P ...
                            fc0T  f13c0T ...
                            fca0T f44ca0T]';
                    end;
                else % fdox
                    if(CAvflag == 0)
                        Y0 = [c0  a0  p0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                            f13c0A  f13c0I   f13c0P f44ca0A f44ca0I f44ca0P ...
                            fc0T  f13c0T ...
                            f44ca0T]';
                    else
                        Y0 = [c0  a0  p0 CA0 cc0 C0 CC0  fc0A    fc0I     fc0P ...
                            f13c0A  f13c0I   f13c0P f44ca0A f44ca0I f44ca0P ...
                            fc0T  f13c0T ...
                            f44ca0T]';
                    end;
                end; % fdox
            end; % ftys
        else
            if(fdox) % fdox
                if(CAvflag == 0)
                    Y0 = [c0 a0 p0 dox0 cc0 C0 CC0]';
                else
                    Y0 = [c0 a0 p0 CA0 dox0 cca0 cc0 C0 CC0]';
                end;
            else % fdox
                if(CAvflag == 0)
                    Y0 = [c0 a0 p0 cc0 C0 CC0]';
                else
                    Y0 = [c0 a0 p0 CA0 cc0 C0 CC0]';
                end;
            end; % fdox
        end; %<<<<<<<<<<<<<<<<<<<<<<< fsed
        
        % Number of ocean tracer (not atm, not sediment)
        if(CAvflag == 0)
            if(fdox)
                NO = 5; % c a p ox cc
            else
                NO = 4; % c a p    cc
            end;
        end;
        if(CAvflag == 1 || CAvflag == 2)
            if(fdox)
                NO = 7; % c a p ox cc
            else
                NO = 6; % c a p    cc
            end;
        end;
        
        
        % set integration time
        
        t0     =   0e5; % (y) time (start)
        if    (BlFlag == 1)
            tfinal =   2e5; % (y) time (end) 2e5
        elseif(BlFlag == 2)
            tfinal =   runtime; % (y) time (end) 2e5 7e5
        else
            tfinal =   runtime; % (y) time (end) 1e7
        end;
        
        
        %==== Anthropogenic CO2 ======================%
        if(ffflag)
            t0     =   1700.; % 1700.
            tfinal =   3000.; % 2100 3000 10000
            
            %load 'dat\co2PEm.tex';
            dirstr = '1000_0500.dat';
            co2PEm = load(['dat\Emss\EmssScen\' dirstr]);
            
            tem   = co2PEm(:,1);
            em    = co2PEm(:,2);
            % deep sea temp
            Dtst = 0.0;
            k1st =  -1;
            TCvt(1,:) = TCv0;
        end;
        Dtst = 0.0;
        k1st =  -1;
        TCvt(1,:) = TCv0;
        Dt = tfinal-t0;
        
        %===============================================%
        %
        % Alternatively, load initial conditions
        %
        %===============================================%
        
        if(ftys == 1)
            if (loadf == 1)
                %   load YPE10Sed.DAT;
                %     load dat\Modern\NoSed0\YPE10Sed.DAT;
                if    (bath == 1)
                    %     load dat\Modern\B1D1BL0\YPE10Sed;
                    %     load dat\B1D1BL0\YPE10Sed;
                    %     load dat\B1D3BL0\YPE10Sed;
                elseif(bath == 2 && CAvflag == 2 && ~LTflag && inorgF)
                    if(epspF)
                        if(biopF)
                            load dat/PERM44Ca/YPE10SedCAv1epspBioPump.DAT
                        else
                            load dat/PERM44Ca/YPE10SedCAv1epsp.DAT
                        end
                        %                         load dat/PERM44Ca/YPE10NoSED.DAT
                    else
                        load dat/PERM44Ca/YPE10SedCAv1.DAT %for PETM set-up but Ca=20.00
                    end
                    %%load dat\PETMCa\YPE10SedCAv10.DAT  %for PETM set-up but Ca=10.30
                elseif(bath == 2 && CAvflag == 2 && LTflag && ~inorgF)
                    load dat/PETM44Ca/YPE10SedCAv1LT.DAT %for PETM set-up but Ca=20.00
                    %%load dat\PETMCa\YPE10SedCAv10.DAT  %for PETM set-up but Ca=10.30
                elseif(bath == 2 && CAvflag == 1 && ~LTflag && inorgF)
                    load dat\PERM44Ca\YPE10SedCA1.DAT
                elseif(bath == 2 && CAvflag == 0 && ~LTflag && ~inorgF)
                    load dat\PETM44Ca\YPE10Sed.DAT
                elseif(bath == 2 && CAvflag == 0 && LTflag && ~inorgF)
                    load dat\PETM44Ca\YPE10SedLT.DAT
                elseif(bath == 2 && CAvflag == 2 && ~LTflag && ~inorgF)
                    load dat/PETM44Ca/YPE10SedCa.DAT
                    %%     load dat\Modern\B2D3BL0\YPE10Sed.DAT;
                    %     load dat\Modern\CO2XLS\YPE10Sed.DAT;
                    %     load dat\B2D1BL0\YPE10Sed;
                    %%     load dat\B2D3BL0\YPE10Sed.DAT;
                    %     load dat\B2D3BL0\Lpco2\YPE10Sed.DAT;
                    %     load dat\B2D3BL0\Hpco2\YPE10Sed.DAT;
                    %                 elseif(bath == 3)
                    %     load dat\B3D3BL0\YPE10Sed;
                elseif(bath == 4 && CAvflag == 2 && ~LTflag)
                    load dat/PERM44Ca/YPE10SedCAbsBath.DAT
                    
                    
                end;
                if(CAvflag == 2 && ~LTflag && inorgF)
                    if(bath == 4)
                        Y0=YPE10SedCAbsBath;
                    else
                        if(epspF)
                            if(biopF)
                                Y0 = YPE10SedCAv1epspBioPump;
                            else
                                Y0 = YPE10SedCAv1epsp;
                                
                            end
                            %                             Y0 = YPE10NoSED;
                        else
                            Y0 = YPE10SedCAv1; %for PETM set-up but Ca=20.00
                        end
                        %Y0 = YPE10SedCAv10; %for PETM set-up but Ca=10.30
                    end
                elseif(CAvflag == 2 && ~LTflag && ~inorgF)
                    Y0 = YPE10SedCa;
                elseif(CAvflag == 2 && LTflag)
                    Y0 = YPE10SedCAv1LT; %for PETM set-up but Ca=20.00
                    %Y0 = YPE10SedCAv10; %for PETM set-up but Ca=10.30
                elseif(CAvflag == 1 && ~LTflag)
                    Y0 = YPE10SedCA1;
                elseif(CAvflag == 0 && ~LTflag)
                    Y0 = YPE10Sed;
                elseif(CAvflag == 0 && LTflag)
                    Y0 = YPE10SedLT;
                end;
                % CBl = 3000.e15; % Blast       2200  2000  3000 n2500
                %                 d13CBl = -15.;     % Blast d13C -55   -60   -34   n-50
                RBl = Rst*(d13CBl/1.e3+1.);
                C13Bl = RBl*CBl;
                kb = 07;       % 07
                kbb = kb+3*Nb;
                if(BlFlag == 1)
                    if(CAvflag == 0)
                        Y0(6*Nb+1)  = Y0(6*Nb+1) +  CBl/12/Aoc;   % add X Gt C atm:/Aoc
                        Y0(6*Nb+2)  = Y0(6*Nb+2) +C13Bl/12/Aoc;   % add X Gt C atm:/Aoc
                    else
                        Y0(7*Nb+1)  = Y0(7*Nb+1) +  CBl/12/Aoc;   % add X Gt C atm:/Aoc
                        Y0(7*Nb+2)  = Y0(7*Nb+2) +C13Bl/12/Aoc;   % add X Gt C atm:/Aoc
                        %   Y0(kb)  = Y0(kb) +  CBl/12/V(kb); % add X Gt C ocn:/V2
                        %   Y0(kbb) = Y0(kbb)+C13Bl/12/V(kb); % add X Gt C ocn:/V2
                    end;
                end;
            end;
        else
            
            if (loadf == 1)
                %   load YPE10Sed.DAT;
                %     load dat\Modern\NoSed0\YPE10Sed.DAT;
                if    (bath == 1)
                    %     load dat\Modern\B1D1BL0\YPE10Sed;
                    %     load dat\B1D1BL0\YPE10Sed;
                    %     load dat\B1D3BL0\YPE10Sed;
                elseif(bath == 2 && CAvflag == 2)
                    load dat/Modern/YPE10SedMcav.DAT
                elseif(bath == 2 && CAvflag == 1)
                    load dat\Modern\YPE10SedMca.DAT
                elseif(bath == 2 && CAvflag == 0)
                    load dat/Modern/YPE10SedM.DAT
                    %%     load dat\Modern\B2D3BL0\YPE10Sed.DAT;
                    %     load dat\Modern\CO2XLS\YPE10Sed.DAT;
                    %     load dat\B2D1BL0\YPE10Sed;
                    %%     load dat\B2D3BL0\YPE10Sed.DAT;
                    %     load dat\B2D3BL0\Lpco2\YPE10Sed.DAT;
                    %     load dat\B2D3BL0\Hpco2\YPE10Sed.DAT;
                elseif(bath == 3)
                    %     load dat\B3D3BL0\YPE10Sed;
                    
                    
                end;
                if(CAvflag == 2)
                    Y0 = YPE10SedMcav;
                elseif(CAvflag == 1)
                    Y0 = YPE10SedMca;
                elseif(CAvflag == 0)
                    Y0 = YPE10SedM;
                end;
                %CBl = 3000.e15; % Blast       2200  2000  3000 n2500
                % d13CBl = -50.;     % Blast d13C -55   -60   -34   n-50
                %    RBl = Rst*(d13CBl/1.e3+1.);
                %  C13Bl = RBl*CBl;
                %     kb = 07;       % 07
                %    kbb = kb+3*Nb;
                % if(BlFlag == 1)
                %   if(CAvflag == 0)
                %     Y0(5*Nb+1)  = Y0(5*Nb+1) +  CBl/12/Aoc;   % add X Gt C atm:/Aoc
                %   Y0(5*Nb+2)  = Y0(5*Nb+2) +C13Bl/12/Aoc;   % add X Gt C atm:/Aoc
                %   else
                %     Y0(6*Nb+1)  = Y0(6*Nb+1) +  CBl/12/Aoc;   % add X Gt C atm:/Aoc
                %   Y0(6*Nb+2)  = Y0(6*Nb+2) +C13Bl/12/Aoc;   % add X Gt C atm:/Aoc
                % %   Y0(kb)  = Y0(kb) +  CBl/12/V(kb); % add X Gt C ocn:/V2
                % %   Y0(kbb) = Y0(kbb)+C13Bl/12/V(kb); % add X Gt C ocn:/V2
                %   end;
                % end;
            end;
        end; %ftys
        % adjust PO4
        %Y0(2*Nb+1:3*Nb) = Y0(2*Nb+1:3*Nb)*0.87;
        if(savf == 1 && ftys == 1)
            Y0(2*Nb+1:3*Nb) = Y0(2*Nb+1:3*Nb)*1.0133; %THIS ONLY WHEN CALC ST. ST
            
        elseif(savf == 1 && ftys == 0)
            Y0(2*Nb+1:3*Nb) = Y0(2*Nb+1:3*Nb);
        end;
        %Y0(5*Nb+1:6*Nb ) = cca0;
        %Y0(3*Nb+1:4*Nb) = Y0(3*Nb+1:4*Nb)*1.1;
        %Y0(2*Nb+1:3*Nb) = Y0(2*Nb+1:3*Nb)*1.2;
        % set initial Oxygen
        %Y0(3*Nb+1:4*Nb) = dox0;
        % adjust TCO2
        %Y0(   1:  Nb) = Y0(   1:  Nb)/1.03;
        %Y0(Nb+1:2*Nb) = Y0(Nb+1:2*Nb)/1.03;
        % loadf
        
        
        % initial inventory
        if(CAvflag == 0)
            Mci  = sum(Y0(     1:   Nb).*V')/Voc ... % TC inventory/V
                +Y0(NO*Nb+1     ) *Aoc/Voc;    % atmosphere
            Mai  = sum(Y0(   Nb+1:2*Nb).*V')/Voc;    % TA inventory/V
            Mpi  = sum(Y0( 2*Nb+1:3*Nb).*V')/Voc;    % P  inventory/V
            Mcci = sum(Y0( 4*Nb+1:5*Nb).*V')/Voc ... % T13C inventory/V
                +Y0(NO*Nb+2     ) *Aoc/Voc;    % atmosphere
        else
            Mci  = sum(Y0(     1:   Nb).*V')/Voc ... % TC inventory/V
                +Y0(NO*Nb+1     ) *Aoc/Voc;    % atmosphere
            Mai  = sum(Y0(   Nb+1:2*Nb).*V')/Voc;    % TA inventory/V
            Mpi  = sum(Y0( 2*Nb+1:3*Nb).*V')/Voc;    % P  inventory/V
            MCAi  = sum(Y0( 3*Nb+1:4*Nb).*V')/Voc;    % Ca inventory
            Mcci = sum(Y0( 6*Nb+1:7*Nb).*V')/Voc ... % T13C inventory/V
                +Y0(NO*Nb+2     ) *Aoc/Voc;    % atmosphere
            MCA44i  = sum(Y0( 5*Nb+1:6*Nb).*V')/Voc;    % Ca44 inventory
            
        end;
        % initial total C inventory ocean+atm (Pg C)
        Mci*Voc*12/1e15;
        if(CAvflag > 0)
            % initial total Ca inventory ocean (Pg Ca)
            MCAi*Voc*40/1e15;
        end;
        if(CAvflag == 0)
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
        else
            if(fsed) % fsed
                if(fdox) % fdox
                    fc0A   =  Y0(7*Nb+3     :7*Nb+2+1*Ns);
                    fc0I   =  Y0(7*Nb+3+1*Ns:7*Nb+2+2*Ns);
                    fc0P   =  Y0(7*Nb+3+2*Ns:7*Nb+2+3*Ns);
                    f13c0A =  Y0(7*Nb+3+3*Ns:7*Nb+2+4*Ns);
                    f13c0I =  Y0(7*Nb+3+4*Ns:7*Nb+2+5*Ns);
                    f13c0P =  Y0(7*Nb+3+5*Ns:7*Nb+2+6*Ns);
                    fca0A  =  Y0(7*Nb+3+6*Ns:7*Nb+2+7*Ns);
                    fca0I  =  Y0(7*Nb+3+7*Ns:7*Nb+2+8*Ns);
                    fca0P  =  Y0(7*Nb+3+8*Ns:7*Nb+2+9*Ns);
                    f44ca0A = Y0(7*Nb+3+9*Ns:7*Nb+2+10*Ns);
                    f44ca0I = Y0(7*Nb+3+10*Ns:7*Nb+2+11*Ns);
                    f44ca0P = Y0(7*Nb+3+11*Ns:7*Nb+2+12*Ns);
                    if(ftys)
                        fc0T   =  Y0(7*Nb+3+12*Ns:7*Nb+2+13*Ns);
                        f13c0T =  Y0(7*Nb+3+13*Ns:7*Nb+2+14*Ns);
                        fca0T  =  Y0(7*Nb+3+14*Ns:7*Nb+2+15*Ns);
                        f44ca0T=  Y0(7*Nb+3+15*Ns:7*Nb+2+16*Ns);
                    end
                else % fdox
                    fc0A   =  Y0(6*Nb+3     :6*Nb+2+1*Ns)
                    fc0I   =  Y0(6*Nb+3+1*Ns:6*Nb+2+2*Ns)
                    fc0P   =  Y0(6*Nb+3+2*Ns:6*Nb+2+3*Ns);
                    f13c0A =  Y0(6*Nb+3+3*Ns:6*Nb+2+4*Ns);
                    f13c0I =  Y0(6*Nb+3+4*Ns:6*Nb+2+5*Ns);
                    f13c0P =  Y0(6*Nb+3+5*Ns:6*Nb+2+6*Ns);
                    f44ca0A = Y0(6*Nb+3+6*Ns:6*Nb+2+7*Ns);
                    f44ca0I = Y0(6*Nb+3+7*Ns:6*Nb+2+8*Ns);
                    f44ca0P = Y0(6*Nb+3+8*Ns:6*Nb+2+9*Ns);
                    if(ftys)
                        fc0T   =  Y0(6*Nb+3+9*Ns:6*Nb+2+10*Ns);
                        f13c0T =  Y0(6*Nb+3+10*Ns:6*Nb+2+11*Ns);
                        f44ca0T = Y0(6*Nb+3+11*Ns:6*Nb+2+12*Ns);
                    end
                end % fdox
                % calc initial phi
                if(isempty(phic)) % phi = phi(fc)
                    FF    = (phi1-phi0)/(1-phi1);
                    phiiA = (phi0+FF*fc0A)./(1+FF*fc0A);
                    phiiI = (phi0+FF*fc0I)./(1+FF*fc0I);
                    phiiP = (phi0+FF*fc0P)./(1+FF*fc0P);
                    
                    FFca    = (phi1-phi0)/(1-phi1);
                    phiiAca = (phi0+FFca*fca0A)./(1+FFca*fca0A);
                    phiiIca = (phi0+FFca*fca0I)./(1+FFca*fca0I);
                    phiiPca = (phi0+FFca*fca0P)./(1+FFca*fca0P);
                    if(ftys) % ftys
                        phiiT = (phi0+FF*fc0T)./(1+FF*fc0T);
                        phiiTca = (phi0+FFca*fca0T)./(1+FFca*fca0T);
                    end; % ftys
                end; % phi = phi(fc)
            end; % fsed
        end; %CAvflag
        
        % get derivatives at t0
        dYflag = 0;
        dY0dt  = LoscarDifdCAPerm1(t0,Y0);
        
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
        [tv,Y] = myode15s('LoscarDifdCAPerm1',[t0 tfinal],Y0',options); % 23t 15s
        toc % stop  clock
        lt = length(tv);
        
        
        % save solution Y(end)
        if(ftys == 1)
            if(bath == 4)
                if (savf == 1 && CAvflag == 2 && ~LTflag)
                    YPE10SedCAbsBath = Y(lt,:)';
                    
                    save                    dat/PERM44Ca/YPE10SedCAbsBath.DAT  YPE10SedCAbsBath -ASCII -DOUBLE -TABS;
                end
                
            else
                if (savf == 1 && CAvflag == 2 && ~LTflag && inorgF)
                    if(epspF)
                        if(biopF)
                            YPE10SedCAv1epspBioPump = Y(lt,:)';
                            save                    dat/PERM44Ca/YPE10SedCAv1epspBioPump.DAT  YPE10SedCAv1epspBioPump -ASCII -DOUBLE -TABS;
                        else
                            YPE10SedCAv1epsp = Y(lt,:)';
                            save                    dat/PERM44Ca/YPE10SedCAv1epsp.DAT  YPE10SedCAv1epsp -ASCII -DOUBLE -TABS;
                        end
                    else
                        YPE10SedCAv1 = Y(lt,:)';
                        
                        save                    dat/PERM44Ca/YPE10SedCAv1.DAT  YPE10SedCAv1 -ASCII -DOUBLE -TABS;
                    end
                end
                if (savf == 1 && CAvflag == 2 && LTflag && ~inorgF)
                    YPE10SedCAv1LT = Y(lt,:)';
                    
                    
                    save                    dat/PETM44Ca/YPE10SedCAv1LT.DAT  YPE10SedCAv1LT -ASCII -DOUBLE -TABS;
                end
                if (savf == 1 && CAvflag == 1 && ~LTflag && inorgF)
                    YPE10SedCA1 = Y(lt,:)';
                    
                    save                    dat/PERM44Ca/YPE10SedCA1.DAT  YPE10SedCA1 -ASCII -DOUBLE -TABS;
                end
                if(savf == 1 && CAvflag == 0 && ~LTflag && inorgF)
                    YPE10Sed = Y(lt,:)';
                    
                    save                    dat/PERM44Ca/YPE10Sed.DAT YPE10Sed -ASCII -DOUBLE -TABS;
                end
                if(savf == 1 && CAvflag == 0 && LTflag)
                    YPE10SedLT = Y(lt,:)';
                    
                    save                    dat/PERM44Ca/YPE10SedLT.DAT YPE10SedLT -ASCII -DOUBLE -TABS;
                end
                if (savf == 1 && CAvflag == 2 && ~LTflag && ~inorgF)
                    YPE10SedCa = Y(lt,:)';
                    
                    save                    dat/PETM44Ca/YPE10SedCa.DAT  YPE10SedCa -ASCII -DOUBLE -TABS;
                end
            end
        else
            if(savf == 1 && CAvflag == 0)
                YPE10SedM = Y(lt,:)';
                
                save                    dat/Modern/YPE10SedM.DAT  YPE10SedM -ASCII -DOUBLE -TABS;
            end
            if(savf == 1 && CAvflag == 1)
                YPE10SedMca = Y(lt,:)';
                
                save                    dat\Modern\YPE10SedMca.DAT  YPE10SedMca -ASCII -DOUBLE -TABS;
            end
            if(savf == 1 && CAvflag == 2)
                YPE10SedMcav = Y(lt,:)';
                
                save                    dat/Modern/YPE10SedMcav.DAT  YPE10SedMcav -ASCII -DOUBLE -TABS;
            end
            %save dat\Modern\B2D3BL0\YPE10Sed.DAT  YPE10Sed -ASCII -DOUBLE -TABS;
            %save dat\Modern\NoSed0\YPE10Sed.DAT   YPE10Sed -ASCII -DOUBLE -TABS;
            %save dat\B2D3BL0\YPE10Sed.DAT YPE10Sed -ASCII -DOUBLE -TABS; % <------
            %save dat\B2D3BL0\Lpco2\YPE10Sed.DAT  YPE10Sed -ASCII -DOUBLE -TABS;
            %save dat\B2D3BL0\Hpco2\YPE10Sed.DAT  YPE10Sed -ASCII -DOUBLE -TABS;
            %save dat\Modern\CO2XLS\YPE10Sed.DAT YPE10Sed -ASCII -DOUBLE -TABS;
        end;
        
        axx = [t0    tfinal];
        %         if(ftys)
        axx = [-0.5e5 tfinal];
        %         end;
    end; %%%%%%% solflag
    
    
    % call dif again at all t to get dYdt(t)
    % plus other vars & fluxes
    % ############### check if temp changes !!!!!!!!
    dYflag = 1;
    if(dYflag == 1)
        it   = 1;
        dYdt = ones(size(Y'));
        for i=1:lt
            dYdt(:,i) = LoscarDifdCAPerm1(tv(i),Y(i,:)');
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
            
            if(CAvflag>0)
                fcaA  (:,l) = Y( :,l+NO*Nb+2+6*Ns);   %
                fcaI  (:,l) = Y( :,l+NO*Nb+2+7*Ns);   %
                fcaP  (:,l) = Y( :,l+NO*Nb+2+8*Ns);   %
                
                mCAlA(:,l) = Y(:,l+NO*Nb+2+6*Ns).*(rhos*(1-phivtA(:, l)));
                mCAlI(:,l) = Y(:,l+NO*Nb+2+7*Ns).*(rhos*(1-phivtI(:, l)));
                mCAlP(:,l) = Y(:,l+NO*Nb+2+8*Ns).*(rhos*(1-phivtP(:, l)));
                MCAlvA( l) = Y(lt,l+NO*Nb+2+6*Ns).*(rhos*(1-phivtA(lt,l)));
                MCAlvI( l) = Y(lt,l+NO*Nb+2+7*Ns).*(rhos*(1-phivtI(lt,l)));
                MCAlvP( l) = Y(lt,l+NO*Nb+2+8*Ns).*(rhos*(1-phivtP(lt,l)));
                %====== Calcium 44
                f44caA (:,l) = Y( :,l+NO*Nb+2+9*Ns);
                f44caI (:,l) = Y( :,l+NO*Nb+2+10*Ns);
                f44caP (:,l) = Y( :,l+NO*Nb+2+11*Ns);
                
                m44calA(:,l) = mCAlA( :,l).*f44caA( :,l)./fcaA( :,l);
                m44calI(:,l) = mCAlI( :,l).*f44caI( :,l)./fcaI( :,l);
                m44calP(:,l) = mCAlP( :,l).*f44caP( :,l)./fcaP( :,l);
                
                M44calvA( l) = mCAlA(lt,l).*f44caA(lt,l)./fcaA(lt,l);
                M44calvI( l) = mCAlI(lt,l).*f44caI(lt,l)./fcaI(lt,l);
                M44calvP( l) = mCAlP(lt,l).*f44caP(lt,l)./fcaP(lt,l);
            end;
            if(ftys)
                if(CAvflag>0)
                    % Carbon
                    fcT    (:,l) = Y( :,l+NO*Nb+2+12*Ns);
                    f13cT  (:,l) = Y( :,l+NO*Nb+2+13*Ns);
                    mcalT  (:,l) = Y( :,l+NO*Nb+2+12*Ns).*(rhos*(1-phivtT(:, l)));
                    McalvT (  l) = Y(lt,l+NO*Nb+2+12*Ns).*(rhos*(1-phivtT(lt,l)));
                    m13calT(:,l) = mcalT( :,l).*f13cT( :,l)./fcT( :,l);
                    M13calvT( l) = mcalT(lt,l).*f13cT(lt,l)./fcT(lt,l);
                    
                    % Calcium
                    fcaT    (:,l) = Y( :,l+NO*Nb+2+14*Ns);
                    mCAlT  (:,l) = Y( :,l+NO*Nb+2+14*Ns).*(rhos*(1-phivtT(:, l)));
                    MCAlvT (  l) = Y(lt,l+NO*Nb+2+14*Ns).*(rhos*(1-phivtT(lt,l)));
                    f44caT  (:,l) = Y( :,l+NO*Nb+2+15*Ns);
                    m44calT(:,l) = mCAlT( :,l).*f44caT( :,l)./fcaT( :,l);
                    M44calvT( l) = mCAlT(lt,l).*f44caT(lt,l)./fcaT(lt,l);
                    
                else
                    
                    fcT    (:,l) = Y( :,l+NO*Nb+2+6*Ns);  % 9*Ns
                    f13cT  (:,l) = Y( :,l+NO*Nb+2+7*Ns);  % 10*Ns
                    mcalT  (:,l) = Y( :,l+NO*Nb+2+6*Ns).*(rhos*(1-phivtT(:, l)));
                    McalvT (  l) = Y(lt,l+NO*Nb+2+6*Ns).*(rhos*(1-phivtT(lt,l)));
                    m13calT(:,l) = mcalT( :,l).*f13cT( :,l)./fcT( :,l);
                    M13calvT( l) = mcalT(lt,l).*f13cT(lt,l)./fcT(lt,l);
                    
                end;
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
        
        
        if(CAvflag>0)
            % final MCAl and fca
            MCAlA = sum(MCAlvA.*VsvA)/VsA;
            MCAlI = sum(MCAlvI.*VsvI)/VsI;
            MCAlP = sum(MCAlvP.*VsvP)/VsP;
            fcafA  = fcaA(lt,:);
            fcafI  = fcaI(lt,:);
            fcafP  = fcaP(lt,:);
            
            %====== Calcium 44
            M44calA = sum(M44calvA.*VsvA)/VsA;
            M44calI = sum(M44calvI.*VsvI)/VsI;
            M44calP = sum(M44calvP.*VsvP)/VsP;
            f44cafA  = f44caA(lt,:);
            f44cafI  = f44caI(lt,:);
            f44cafP  = f44caP(lt,:);
        end;
        
        % test: Average must equal Rin
        RsAf = (f13cfA./fcfA/Rst-1)*1e3;
        RsIf = (f13cfI./fcfI/Rst-1)*1e3;
        RsPf = (f13cfP./fcfP/Rst-1)*1e3;
        
        % test: Average must equal RinCA
        if(CAvflag>0)
            RsAfca = (f44cafA./fcafA/Rstca-1)*1e3;
            RsIfca = (f44cafI./fcafI/Rstca-1)*1e3;
            RsPfca = (f44cafP./fcafP/Rstca-1)*1e3;
        end;
        
        if(ftys)
            McalT = sum(McalvT.*VsvT)/VsT;
            fcfT  = fcT(lt,:);
            M13calT = sum(M13calvT.*VsvT)/VsT;
            f13cfT  = f13cT(lt,:);
            RsTf = (f13cfT./fcfT/Rst-1)*1e3;
            if(CAvflag>0)
                MCAlT = sum(MCAlvT.*VsvT)/VsT;
                fcafT  = fcaT(lt,:);
                M44calT = sum(M44calvT.*VsvT)/VsT;
                f44cafT  = f44caT(lt,:);
                RsTfca = (f44cafT./fcafT/Rstca-1)*1e3;
            end;
        end;
        
    end; %----------- end sediments
    
    
    % inventories/average concentrations ocean
    if(CAvflag == 0)
        Mc  = sum(Y(lt,     1:   Nb).*V)/Voc ... % TC inventory/V
            +Y(lt,NO*Nb+1     )*Aoc/Voc;    % atmosphere
        Ma  = sum(Y(lt,   Nb+1:2*Nb).*V)/Voc;    % TA inventory/V
        Mp  = sum(Y(lt, 2*Nb+1:3*Nb).*V)/Voc;    % TA inventory/V
        Mcc = sum(Y(lt, 4*Nb+1:5*Nb).*V)/Voc ... % TC inventory/V
            +Y(lt,NO*Nb+2     )*Aoc/Voc;    % atmosphere
        
        DX0 = [Mc-Mci Ma-Mai Mp-Mpi Mcc-Mcci];
        DX0
    else
        Mc  = sum(Y(lt,     1:   Nb).*V)/Voc ... % TC inventory/V
            +Y(lt,NO*Nb+1     )*Aoc/Voc;    % atmosphere
        Ma  = sum(Y(lt,   Nb+1:2*Nb).*V)/Voc;    % TA inventory/V
        Mp  = sum(Y(lt, 2*Nb+1:3*Nb).*V)/Voc;    % TA inventory/V
        MCA  = sum(Y(lt, 3*Nb+1:4*Nb).*V)/Voc;   % Ca inventory
        Mcc = sum(Y(lt, 6*Nb+1:7*Nb).*V)/Voc ... % TC inventory/V
            +Y(lt,NO*Nb+2     )*Aoc/Voc;    % atmosphere
        DX0 = [Mc-Mci Ma-Mai Mp-Mpi Mcc-Mcci MCA-MCAi];
        DX0
    end;
    
    
    % total C inventory ocean     (Pg C)
    MCO=sum(Y(lt, 1:Nb).*V)    *12/1e15;
    % total C inventory atm       (Pg C)
    MCatm=    Y(lt,NO*Nb+1  )*Aoc*12/1e15;
    % total C inventory ocean+atm (Pg C)
    Mctot=Mc*Voc*12/1e15
    
    if(CAvflag>0)
        %total Ca inventory ocean     (Pg Ca)
        MCAL=sum(Y(lt, 3*Nb+1:4*Nb).*V)    *40/1e15;
        % total Ca inventory ocean (Pg Ca)
        MCatot=MCA*Voc*40/1e15
    end;
    
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
                
                if(CAvflag>0)
                    DisstACA(i,k)   = ...
                        (dissvtAca(i+1,k)+dissvtAca(i,k))*dtv(i)/2; % kg
                    DisstICA(i,k)   = ...
                        (dissvtIca(i+1,k)+dissvtIca(i,k))*dtv(i)/2; % kg
                    DisstPCA(i,k)   = ...
                        (dissvtPca(i+1,k)+dissvtPca(i,k))*dtv(i)/2; % kg
                    Diss44tA(i,k)   = ...
                        (dissv44tA(i+1,k)+dissv44tA(i,k))*dtv(i)/2; % kg
                    Diss44tI(i,k)   = ...
                        (dissv44tI(i+1,k)+dissv44tI(i,k))*dtv(i)/2; % kg
                    Diss44tP(i,k)   = ...
                        (dissv44tP(i+1,k)+dissv44tP(i,k))*dtv(i)/2; % kg
                end;
                
                if(ftys)
                    burialtT(i,k) = (  mcalT(i+1,k)+  mcalT(i,k)) ...
                        *(rsedvtT(i+1,k)+rsedvtT(i,k)) ...
                        *asvT(k)*A(11)*dtv(i)/4; % kg
                    DisstT(i,k)   = ...
                        (dissvtT(i+1,k)+dissvtT(i,k))*dtv(i)/2; % kg
                    Diss13tT(i,k)   = ...
                        (dissv13tT(i+1,k)+dissv13tT(i,k))*dtv(i)/2; % kg
                    if(CAvflag>0)
                        DisstTCA(i,k)   = ...
                            (dissvtTca(i+1,k)+dissvtTca(i,k))*dtv(i)/2; % kg
                        Diss44tT(i,k)   = ...
                            (dissv44tT(i+1,k)+dissv44tT(i,k))*dtv(i)/2; % kg
                    end;
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
            if(CAvflag>0)
                FprdtACA  (i) = (FprtAca  (i+1)+FprtAca(  i))*dtv(i)/2;
                FprdtICA  (i) = (FprtIca  (i+1)+FprtIca  (i))*dtv(i)/2;
                FprdtPCA  (i) = (FprtPca  (i+1)+FprtPca  (i))*dtv(i)/2;
                Fpr44dtA(i) = (Fpr44tA(i+1)+Fpr44tA(i))*dtv(i)/2;
                Fpr44dtI(i) = (Fpr44tI(i+1)+Fpr44tI(i))*dtv(i)/2;
                Fpr44dtP(i) = (Fpr44tP(i+1)+Fpr44tP(i))*dtv(i)/2;
                Fin44dt (i) = (Fin44t (i+1)+Fin44t (i))*dtv(i)/2; % mol C/m2
                FSi44dt (i) = (FSi44t (i+1)+FSi44t (i))*dtv(i)/2; % mol C/m2
            end;
            
            
            if(ftys)
                FprdtT  (i) = (FprtT  (i+1)+FprtT  (i))*dtv(i)/2;
                Fpr13dtT(i) = (Fpr13tT(i+1)+Fpr13tT(i))*dtv(i)/2;
                if(CAvflag>0)
                    FprdtTCA  (i) = (FprtTca  (i+1)+FprtTca  (i))*dtv(i)/2;
                    Fpr44dtT(i) = (Fpr44tT(i+1)+Fpr44tT(i))*dtv(i)/2;
                end;
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
        if(CAvflag>0)
            DissACA    = sum(sum(DisstACA))  /VsA; % kg /m3 sed
            DissICA    = sum(sum(DisstICA))  /VsI; % kg /m3 sed
            DissPCA    = sum(sum(DisstPCA))  /VsP; % kg /m3 sed
            DissCAA   = DissACA*VsA/Voc/m2kg;     % mol/m3 ocean
            DissCAI   = DissICA*VsI/Voc/m2kg;     % mol/m3 ocean
            DissCAP   = DissPCA*VsP/Voc/m2kg;     % mol/m3 ocean
            DDissCA   = sum(DissCAA+DissCAI+DissCAP);
        end
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
        %== Calcium-44
        if(CAvflag>0)
            Diss44A    = sum(sum(Diss44tA))  /VsA; % kg /m3 sed
            Diss44I    = sum(sum(Diss44tI))  /VsI; % kg /m3 sed
            Diss44P    = sum(sum(Diss44tP))  /VsP; % kg /m3 sed
            Diss44CA   = Diss44A*VsA/Voc/m2kg;     % mol/m3 ocean
            Diss44CI   = Diss44I*VsI/Voc/m2kg;     % mol/m3 ocean
            Diss44CP   = Diss44P*VsP/Voc/m2kg;     % mol/m3 ocean
            DDiss44C   = sum(Diss44CA+Diss44CI+Diss44CP);
        end;
        
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
            
            if(CAvflag>0)
                DissTCA    = sum(sum(DisstTCA))  /VsT; % kg /m3 sed
                DissCAT   = DissTCA*VsT/Voc/m2kg;     % mol/m3 ocean
                DDissCA   = sum(DissCAA+DissCAI+DissCAP+DissCAT);
                Diss44T    = sum(sum(Diss44tT))  /VsT; % kg /m3 sed
                Diss44CT   = Diss44T*VsT/Voc/m2kg;     % mol/m3 ocean
                DDiss44C   = sum(Diss44CA+Diss44CI+Diss44CP+Diss44CT);
            end;
            
        end;
        
        
        % right-hand sides
        FinT   =  sum(Findt)/Dt;                    % mol/m2/y
        FSiT   =  sum(FSidt)/Dt;                    % mol/m2/y
        FprT   =  sum(  sum(FprdtA*A(1)) ...        % mol/m2/y
            +sum(FprdtI*A(2)) ...
            +sum(FprdtP*A(3)) )/Aoc/Dt;
        if(CAvflag>0)
            FprTCA   =  sum(  sum(FprdtACA*A(1)) ...        % mol/m2/y
                +sum(FprdtICA*A(2)) ...
                +sum(FprdtPCA*A(3)) )/Aoc/Dt;
        end
        % right-hand sides
        Fin13T   =  sum(Fin13dt)/Dt;                % mol/m2/y
        FSi13T   =  sum(FSi13dt)/Dt;                % mol/m2/y
        Fpr13T   =  sum(sum(Fpr13dtA*A(1)) ...      % mol/m2/y
            +sum(Fpr13dtI*A(2)) ...
            +sum(Fpr13dtP*A(3)) )/Aoc/Dt;
        % right-hand sides
        if(CAvflag>0)
            Fin44T   =  sum(Fin44dt)/Dt;                % mol/m2/y
            FSi44T   =  sum(FSi44dt)/Dt;                % mol/m2/y
            Fpr44T   =  sum(sum(Fpr44dtA*A(1)) ...      % mol/m2/y
                +sum(Fpr44dtI*A(2)) ...
                +sum(Fpr44dtP*A(3)) )/Aoc/Dt;
        end;
        
        if(ftys)
            FprT   = FprT   + sum( sum(FprdtT  *A(11)) )/Aoc/Dt;
            Fpr13T = Fpr13T + sum( sum(Fpr13dtT*A(11)) )/Aoc/Dt;
            if(CAvflag>0)
                FprTCA   = FprTCA   + sum( sum(FprdtTCA  *A(11)) )/Aoc/Dt;
                Fpr44T = Fpr44T + sum( sum(Fpr44dtT*A(11)) )/Aoc/Dt;
            end;
        end;
        
        rhc    = (FinT+FSiT-FprT)*Aoc*Dt/Voc+DDissC; % mol/m3 ocean
        if(CAvflag>0)
            rhCA    = (FinT+FSiT-FprTCA)*Aoc*Dt/Voc+DDissCA; % mol/m3 ocean
        end
        %rhMc   =  FprT*m2kg*Aoc*Dt/Vs-Burial-Diss;  % kg /m3 sed
        
        rhc13  = (Fin13T+FSi13T-Fpr13T)*Aoc*Dt/Voc+DDiss13C; % mol/m3 ocean
        
        if(CAvflag>0)
            rhCA44  = (Fin44T+FSi44T-Fpr44T)*Aoc*Dt/Voc+DDiss44C; % mol/m3 ocean
        end;
        
        DXs = [(Mc-Mci)-rhc (Ma-Mai)-2*rhc Mcc-Mcci-rhc13];
        DXs
        
        %====== steady state: Fin + FSi = burial ?
        
        % total influx. unit: [F*Aoc] = mol/y
        Fintot = (Fint(lt)+FSit(lt))*Aoc/1.e12;
        
        % total burial. unit  [wcvt] = m/y
        % (CaCO3+H2O)*(1-phi1) = pure CaCO3
        % *area = m3/y. *rhos = kg/y. *(100/1e3) = mol/y
        
        ke = [2 Ns]; % 2: shelf, Ns: all
        for i=1:2
            kk = [1:1:ke(i)];
            tmp =   sum(wcvtA(lt,kk).*asvA(kk))*A(1) ...
                + sum(wcvtI(lt,kk).*asvI(kk))*A(2) ...
                + sum(wcvtP(lt,kk).*asvP(kk))*A(3);
            if(ftys)
                tmp = tmp ...
                    + sum(wcvtT(lt,kk).*asvT(kk))*A (11);
            end
            ttmp(i) = tmp*(1-phi1)*rhos/m2kg/1.e12;
        end
        Fbursh  = ttmp(1);
        Fburtot = ttmp(2);
        
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
    if(CAvflag == 0)
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
    else
        if(fdox)
            for k=1:Nb
                c  (:,k)   = Y(:,     k)/(1e-3*rho); % mol/m3 -> mmol/kg
                a  (:,k)   = Y(:,  Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
                p  (:,k)   = Y(:,2*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
                CA (:,k)   = Y(:,3*Nb+k)/(1e-3*rho);
                dox(:,k)   = Y(:,4*Nb+k)           ; % mol/m3
                cac (:,k)   = Y(:,5*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
                cc (:,k)   = Y(:,6*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
                %====== Carbon-13
                d13c(:,k) = ((cc(:,k)./c(:,k))/Rst-1.)*1e3; % Ocean
                d44ca(:,k) = ((cac(:,k)./CA(:,k))/Rstca-1.)*1e3; % Ocean
                DIC(:,k)=c(:,k);
                DIC13(:,k)=cc(:,k);
                Cac(:,k)=CA(:,k);
                Ca44(:,k)=cac(:,k);
            end;
            C        = Y(:,7*Nb+1);
            CC       = Y(:,7*Nb+2);
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
                CA (:,k)   = Y(:,3*Nb+k)/(1e-3*rho);
                cac (:,k)   = Y(:,4*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
                cc (:,k)   = Y(:,5*Nb+k)/(1e-3*rho); % mol/m3 -> mmol/kg
                %====== Carbon-13
                d13c(:,k) = ((cc(:,k)./c(:,k))/Rst-1.)*1e3; % Ocean
                d44ca(:,k) = ((cac(:,k)./CA(:,k))/Rstca-1.)*1e3; % Ocean
            end;
            C        = Y(:,6*Nb+1);
            CC       = Y(:,6*Nb+2);
            pco2t    =  C/(2.2e15/12/Aoc);
            pcco2t   = CC/(2.2e15/12/Aoc);
            pco2a    =  pco2t(lt);
            pcco2a   = pcco2t(lt);
            pco2a
        end; % fdox
    end; %CAvflag
    
    
    %====== Carbon-13
    d13cv = d13c(lt,:);                  % Ocean
    d13C  = ((CC./C       )/Rst-1.)*1e3; % Atmosphere
    d13Ca = ((CC(lt)/C(lt))/Rstca-1.)*1e3; % Atmosphere
    
    %====== Calcium-44
    if(CAvflag>0)
        d44cav = d44ca(lt,:)                  % Ocean
    end;
    d13CV = [d13cv d13Ca];
    d13CV
    
    % 13C surf-deep diff & excursion
    k = 1;
    Dd13cv  = d13c(1,1:3)-d13c(1,7:9);
    d13cex = max(d13c(:,k))-min(d13c(:,k));
    
    
    if(CAvflag == 0)
        
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
    else
        for k=1:Nb
            dic = c(lt,k)*1e-3;
            alk = a(lt,k)*1e-3;
            Cal = CA(lt,k)*1e-3;%%%THIS ADDED !!!!!!!!!!!!!!
            TC = TCvt(lt,k);
            S  = Sv(k);
            P  = Pv(k);
            [co2,pco2,co3,ph,kh,o2] = dafunPECA(dic,alk,TC,S,P,Cal,Mg);
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
                Cal  = CA(i,k)*1e-3;%%%THIS ADDED !!!!!!!!!!!!!!
                TC = TCvt(i,k);
                S  = Sv(k);
                P  = Pv(k);
                [co2,pco2,co3,ph,kh,o2] = dafunPECA(dic,alk,TC,S,P,Cal,Mg);
                co3tv(i,k)  = co3;
                phtv(i,k)  = ph;
            end;
        end;
        
        %==== calculate omega(t)
        %         %==== saturation of surface ocean boxes
        %         for i=1:lt
        %             l = 1;
        %             for k=kkv
        %                 [kspcS(k),kspaS(k)] = ...
        %                     kspfun(TCvt(i,k),Sv(k),Pv(k),CA(i,k)*1e-3,Mg);
        %                 omegCSvt(i,l) = co3tv(i,k)*CA(i,k)/1000/kspcS(k);
        %                 omegASvt(i,l) = co3tv(i,k)*CA(i,k)/1000/kspaS(k);
        %                 l = l + 1;
        %             end
        %         end
        %==== saturation of all ocean boxes
        for i=1:lt
            l = 1;
            for k=1:Nb
                [kspcS(k),kspaS(k)] = ...
                    kspfun(TCvt(i,k),Sv(k),Pv(k),CA(i,k)*1e-3,Mg);
                omegCSvt(i,l) = co3tv(i,k)*CA(i,k)/1000/kspcS(k);
                omegASvt(i,l) = co3tv(i,k)*CA(i,k)/1000/kspaS(k);
                l = l + 1;
            end
        end
    end;%CAvflag
    
    
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
        if(CAvflag == 0)
            if    (satflg == 1)
                co3szv = co3s0*exp(as*(zv-zs0));
            elseif(satflg == 2)
                for k=1:lzv
                    [x,x,x,x,x,x,kspc(k),x] = dafunPE(1e-3,1e-3,2,Soc,zv(k)/10.,Ca,Mg);
                end;
                co3szv = kspc/Ca;
            end;
        elseif(CAvflag == 1)
            if    (satflg == 1)
                co3szv = co3s0*exp(as*(zv-zs0));
            elseif(satflg == 2)
                for k=1:lzv
                    [x,x,x,x,x,x,kspc(k),x] = dafunPE(1e-3,1e-3,2,Soc,zv(k)/10.,Ca,Mg);
                end;
                co3szv = kspc/Ca;
            end;
        elseif(CAvflag == 2)
            if    (satflg == 1)
                co3szv = co3s0*exp(as*(zv-zs0));
            elseif(satflg == 2)
                %      mnCA=mean(CA); %average value of calcium over the timesteps for every box
                %      CAr=mnCA.*V/Voc; %relative calcium concentration for each box
                %      avCALC=sum(CAr);
                %avCALC=sum(CALC.*V/Voc)
                for k=1:lzv
                    [x,x,x,x,x,x,kspc(k),x] = dafunPECA(1e-3,1e-3,2,Soc,zv(k)/10.,avCALC*1e-3,Mg);
                    %co3szv = kspc(k)/Cal(k);
                end;
                co3szv = kspc./(avCALC*1e-3);
            end;
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
            rainsh(1) = sum(fsh    *FprtA(lt)*asvA(1: 2)*A(1));
            rainds(1) = sum(fdpv(1)*FprtA(lt)*asvA(3:Ns)*A(1));
            
            rainsh(2) = sum(fshI    *FprtI(lt)*asvI(1: 2)*A(2));
            rainds(2) = sum(fdpv(2)*FprtI(lt)*asvI(3:Ns)*A(2));
            
            rainsh(3) = sum(fshP    *FprtP(lt)*asvP(1: 2)*A(3));
            rainds(3) = sum(fdpv(3)*FprtP(lt)*asvP(3:Ns)*A(3));
            
            if(ftys)
                rainsh(4) = sum(fshT   *FprtT(lt)*asvT(1: 2)*A(11));
                rainds(4) = sum(fdpv(4)*FprtT(lt)*asvT(3:Ns)*A(11));
            end;
            
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
    
    FSiinrease = (max(FSichck)/FSichck(1)-1)*100 %FSi increase in percent
    Finincrease= (max(Finchck)/Finchck(1)-1)*100 %Fin increase in percent
    
    
    if(CAvflag > 0)
        Ca_chng=((max(CA)./(CA(1,:)))-1)*100
        Ca_chngDP=((max(CA(:,9))./(CA(1,9)))-1)*100
    end;
    %Alk_chng=((max(a)./(a(1,:)))-1)*100
    % dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
    % save                    dCalcium.DAT  dCalcium -ASCII -DOUBLE -TABS;
    % load dCalcium.dat
    %  bDP(:,:)=dCalcium
    % % ovco=bDP
    % % dCADP=ovco(slj,1)'
    %        bDPstr = [ 'myDataFile' num2str(slj) '.dat' ];
    %         %save(dCalcium);
    % save(bDPstr,'bDP', '-ASCII')
    
    %clear all all
    %return
    %=====================================================%
    %=================== plot results ====================%
    %=====================================================%
    if(~plotflag)
        return
    end
    
    
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
    if(CAvflag == 0)
        tv11   = [axx(1)    tv(2)]; % tv(2)
        c1     = [c(1,:)' c(1,:)'];
        a1     = [a(1,:)' a(1,:)'];
        p1     = [p(1,:)' p(1,:)'];
        co3tv1 = [co3tv(1,:)' co3tv(1,:)'];
        phtv1 = [ phtv(1,:)'  phtv(1,:)'];
        FSichck1 = [FSichck(1,:)' FSichck(1,:)'];
        Finchck1 = [Finchck(1,:)' Finchck(1,:)'];
    else
        tv11   = [axx(1)    tv(2)]; % tv(2)
        c1     = [c(1,:)' c(1,:)'];
        a1     = [a(1,:)' a(1,:)'];
        p1     = [p(1,:)' p(1,:)'];
        CA1    = [CA(1,:)' CA(1,:)'];
        co3tv1 = [co3tv(1,:)' co3tv(1,:)'];
        phtv1 = [ phtv(1,:)'  phtv(1,:)'];
        FSichck1 = [FSichck(1,:)' FSichck(1,:)'];
        Finchck1 = [FiN*Aoc FiN*Aoc];
        
        omegCSvt1 = [omegCSvt(1,:)' omegCSvt(1,:)'];
    end
    %---------- TCO2 ---------------%
    figure(1)
    clf;
    box  on;
    hold on;
    for k=1:Nb
        plot(tv  ,c  (:,k),sstr(2*k-1:2*k),'Color',cs(k));
    end;
    for k=1:Nb
        plot(tv11,c1 (k,:),sstr(2*k-1:2*k),'Color',cs(k));
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
        plot(tv  ,a  (:,k),sstr(2*k-1:2*k),'Color',cs(k));
    end;
    for k=1:Nb
        plot(tv11,a1 (k,:),sstr(2*k-1:2*k),'Color',cs(k));
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
        plot(tv  ,p (:,k)*1e3,sstr(2*k-1:2*k),'Color',cs(k));
    end;
    for k=1:Nb
        plot(tv11,p1(k,:)*1e3,sstr(2*k-1:2*k),'Color',cs(k));
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
            plot(tv  ,dox (:,k),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        for k=1:Nb
            plot(tv11,dox1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
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
                sstr(2*k-1:2*k),'Color',cs(k));
        end;
        for k=1:Nb
            plot(tv11,dox1(k,:)./(dox1(k,:)+KMMOX)*100,...
                sstr(2*k-1:2*k),'Color',cs(k));
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
        plot(tv  ,co3tv (:,k)*1e6,sstr(2*k-1:2*k),'Color',cs(k));
    end;
    for k=1:Nb
        plot(tv11,co3tv1(k,:)*1e6,sstr(2*k-1:2*k),'Color',cs(k));
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
            plot(tv  ,phtv (:,k),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        for k=1:Nb
            plot(tv11,phtv1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
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
    plot(tv  ,pco2t ,'b-');
    hold on;
    for k=kkv
        plot(tv(lt),pco2v(k),systr(2*k-1:2*k));
    end;
    plot(tv11,pco2t1,'b-');
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
    
    if(CAvflag>0)
        d44ca1  = [d44ca(1,:)' d44ca(1,:)'];
        d44CA1  = [d13C(1)'   d13C(1)'  ];
    end;
    
    MS  = 'MarkerSize';
    ms  = 4;
    MFC = 'MarkerFaceColor';
    Dam = 7.;
    
    %-------------- load PETM data
    if(ftys)
        %== Zachos data
        ZchF3a = load('dat/Zachos/ZchF3aNEW.DAT');
        age3  = ZchF3a(:,3)*1e3;
        c13c3 = ZchF3a(:,04);
        %== Roehl data
        % new age model, G^3, 2007
        NEW690     = load('dat/Zachos/NEW690.txt');
        d690n      = NEW690(:,2);
        a690n      = NEW690(:,3);
        % c13 data
        Roehl690   = load('dat/Zachos/Roehl690.txt');
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
        plot(tv  ,d13c (:,k),sstr(2*k-1:2*k),'Color',cs(k));
    end;
    plot(tv  ,d13C +Dam,'m-');
    for k=1:Nb
        plot(tv11,d13c1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
    end;
    plot(tv11,d13C1+Dam,'m-');
    if(ftys)
        plot((time-0.26)*1e6,d13cp,'o');
        %         plot(age3,     c13c3     ,'b--s',MFC,'b',MS,ms);
        %         plot(age690blk,c13c690blk,'r- d',MFC,'r',MS,ms);
    end
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'YDir','reverse');
    xlabel('Time (y)');
    ylabel('\delta^{13}C (�)');
    Hl=legend(lstr,1);
    set(Hl,'FontSize',10);
    
    if(CAvflag>0)
        figure(60)
        clf;
        box  on;
        hold on;
        for k=1:Nb
            plot(tv  ,d44ca (:,k),sstr(2*k-1:2*k),'Color',cs(k));
            %plot(tv  ,d13c (:,k),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        %plot(tv  ,d13C +Dam,'m-');
        for k=1:Nb
            plot(tv11,d44ca1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
            %plot(tv11,d13c1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        %plot(tv11,d13C1+Dam,'m-');
        if(ftys)
            %plot(age3,     c13c3     ,'b--s',MFC,'b',MS,ms);
            %plot(age690blk,c13c690blk,'r- d',MFC,'r',MS,ms);
        end
        hold off;
        set(gca,'FontSize',fs);
        set(gca,'XLim',axx);
        %set(gca,'YDir','reverse');
        xlabel('Time (y)');
        ylabel('\delta^{44}Ca (�)');
        Hl=legend(lstr,1);
        set(Hl,'FontSize',10);
    end;
    
    if(fsed)
        %---------- fc vs. depth----------%
        figure(7)
        clf;
        plot(fcfA*1e2,dsv,'b-d');
        hold on;
        plot(fcfI*1e2,dsv,'m-s');
        plot(fcfP*1e2,dsv,'k-o');
        if(ftys)
            plot(fcfT*1e2,dsv,'g-p');
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
                    plot(tv  ,co3tv (:,k)*1e6   ,sstr(2*k-1:2*k));
                plot(tv11,co3tv1(k,:)*1e6   ,sstr(2*k-1:2*k));
            end;
            pl(2)=  ...
                plot(tv  ,f1+shtA *f2/1000,'r- ');
            plot(tv  ,f1+shtI *f2/1000,'r--');
            plot(tv  ,f1+shtP *f2/1000,'r-.');
            plot(tv11,f1+shtA1*f2/1000,'r- ');
            plot(tv11,f1+shtI1*f2/1000,'r--');
            plot(tv11,f1+shtP1*f2/1000,'r-.');
            pl(3)=  ...
                plot(tv  ,f1+f2*ccdA /1e3,'k- ');
            plot(tv  ,f1+f2*ccdI /1e3,'k--');
            plot(tv  ,f1+f2*ccdP /1e3,'k-.');
            plot(tv11,f1+f2*ccdA1/1e3,'k- ');
            plot(tv11,f1+f2*ccdI1/1e3,'k--');
            plot(tv11,f1+f2*ccdP1/1e3,'k-.');
            if(ftys)
                k=13;
                shtT1  = [shtT(1)' shtT(1)'];
                ccdT1  = [ccdT(1)' ccdT(1)'];
                
                plot(tv  ,co3tv (:,k)*1e6   ,sstr(2*k-1:2*k));
                plot(tv11,co3tv1(k,:)*1e6   ,sstr(2*k-1:2*k));
                plot(tv  ,f1+shtT *f2/1000,'r:');
                plot(tv11,f1+shtT1*f2/1000,'r:');
                plot(tv  ,f1+f2*ccdT /1e3,'k:');
                plot(tv11,f1+f2*ccdT1/1e3,'k:');
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
            plot(tv  ,fcA (:,l)*1e2,'--','Color',cs(l));
        end;
        for l=1:Ns
            plot(tv11,fcA1(l,:)*1e2,'--','Color',cs(l));
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
            plot(tv  ,fcI (:,l)*1e2,'--','Color',cs(l));
        end;
        for l=1:Ns
            plot(tv11,fcI1(l,:)*1e2,'--','Color',cs(l));
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
            plot(tv  ,fcP (:,l)*1e2,'--','Color',cs(l));
        end;
        for l=1:Ns
            plot(tv11,fcP1(l,:)*1e2,'--','Color',cs(l));
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
                plot(tv  ,fcT (:,l)*1e2,'--','Color',cs(l));
            end;
            for l=1:Ns
                plot(tv11,fcT1(l,:)*1e2,'--','Color',cs(l));
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
        
        figure(102)
        clf;
        subplot(221)
        
        box  on;
        hold on;
        for l=1:Ns
            plot(tv  ,fcA (:,l)*1e2,'--','Color',cs(l));
        end;
        for l=1:Ns
            plot(tv11,fcA1(l,:)*1e2,'--','Color',cs(l));
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
        subplot(222)
        
        box  on;
        hold on;
        for l=1:Ns
            plot(tv  ,fcI (:,l)*1e2,'--','Color',cs(l));
        end;
        for l=1:Ns
            plot(tv11,fcI1(l,:)*1e2,'--','Color',cs(l));
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
        subplot(223)
        
        box  on;
        hold on;
        for l=1:Ns
            plot(tv  ,fcP (:,l)*1e2,'--','Color',cs(l));
        end;
        for l=1:Ns
            plot(tv11,fcP1(l,:)*1e2,'--','Color',cs(l));
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
            subplot(224)
            
            box  on;
            hold on;
            for l=1:Ns
                plot(tv  ,fcT (:,l)*1e2,'--','Color',cs(l));
            end;
            for l=1:Ns
                plot(tv11,fcT1(l,:)*1e2,'--','Color',cs(l));
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
        
        if(1)
            d13fcA = (f13cA./fcA/Rst-1)*1e3;
            d13fcI = (f13cI./fcI/Rst-1)*1e3;
            d13fcP = (f13cP./fcP/Rst-1)*1e3;
            if(ftys)
                d13fcT = (f13cT./fcT/Rst-1)*1e3;
            end
            
            if(CAvflag~=0)
                d44fcA = (f44caA./fcaA/Rstca-1)*1e3;
                d44fcI = (f44caI./fcaI/Rstca-1)*1e3;
                d44fcP = (f44caP./fcaP/Rstca-1)*1e3;
                if(ftys)
                    d44fcT = (f44caT./fcaT/Rstca-1)*1e3;
                end
            end
            %-------------- f13c A --------------------%
            figure(13)
            clf;
            subplot(211)
            plot(tv  ,d13fcA(:,1),'-','Color',cs(1));
            %plot(tv  ,f13cA(:,1),'-','Color',cs(1));
            %plot(tv  ,m13calA(:,1),'-','Color',cs(1));
            hold on;
            for l=2:Ns
                plot(tv  ,d13fcA(:,l),'--','Color',cs(l));
                plot(tv  ,d13fcI(:,l),'--','Color','r');
                plot(tv  ,d13fcP(:,l),'--','Color','m');
                if(ftys)
                    plot(tv  ,d13fcT(:,l),'--','Color','g');
                end
                %plot(tv  ,f13cA(:,l),'--','Color',cs(l));
                %plot(tv  ,m13calA(:,l),'--','Color',cs(l));
            end;
            plot((time-0.26)*1e6,d13cp,'o');
            hold off;
            set(gca,'FontSize',fs);
            %set(gca,'YDIr','reverse');
            set(gca,'XLim',axx);
            %             set(gca,'XScale','log');
            title('Atlantic');
            xlabel('Time (y)');
            %             ylabel('^{13}f_c');
            ylabel('\delta^{13} C (sediments)');
            
            Hl=legend(dstr,4);
            set(Hl,'FontSize',10);
            if(CAvflag~=0)
                subplot(212)
                plot(tv  ,d44fcA(:,1),'-','Color',cs(1));
                
                hold on;
                for l=2:Ns
                    plot(tv  ,d44fcA(:,l),'--','Color',cs(l));
                    plot(tv  ,d44fcI(:,l),'--','Color','r');
                    plot(tv  ,d44fcP(:,l),'--','Color','m');
                    if(ftys)
                        plot(tv  ,d44fcT(:,l),'--','Color','g');
                    end
                end;
                hold off;
            end
            set(gca,'FontSize',fs);
            %set(gca,'YDIr','reverse');
            set(gca,'XLim',axx);
            %             set(gca,'XScale','log');
            title('Atlantic');
            xlabel('Time (y)');
            %             ylabel('^{44}f_c');
            ylabel('\delta^{44} Ca (sediments)');
            
            Hl=legend(dstr,4);
            set(Hl,'FontSize',10);
        end;
        
        %         cs = 'bgrmc';
        cs = 'gggkkkrrrbgkr';
        sstr = '- ---.- ---.- ---. -: : : ';
        if(inorgF)
            figure(40)
            clf;
            box  on;
            hold on;
            for lk=1:Nb
                plot(tv  ,omegCSvt (:,lk),sstr(2*lk-1:2*lk),'Color',cs(lk));
            end;
            hold off;
            set(gca,'FontSize',fs);
            %         axis([]);
            %         set(gca,'YDIr','reverse');
            set(gca,'XLim',axx);
            title('');
            xlabel('Time (y)');
            ylabel('\Omega (surface)');
            
            clear dstr;
            sk=1;
            for k=1:Nb
                
                lstr(sk,:) = sprintf('%s',lstr0(2*k-1:2*k));
                sk=sk+1;
            end;
            Hl=legend(lstr,4);
            set(Hl,'FontSize',10);
        end
        
    end; % sediments
    
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
        ylabel('\delta^{13}C (�)');
        
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
    
    figure(20)
    clf;
    box  on;
    hold on;
    if(CAvflag == 1)
        for k=1:Nb
            plot(tv  ,CA  (:,k),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        for k=1:Nb
            plot(tv11,CA1 (k,:),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        
        hold off;
        set(gca,'FontSize',fs);
        set(gca,'XLim',axx);
        xlabel('Time (y)');
        ylabel('Ca^2^+ (mmol kg^{-1})');
        Hl=legend(lstr,4);
        set(Hl,'FontSize',10);
        refresh;
    end;
    if(CAvflag == 2)
        for k=1:Nb
            plot(tv  ,CA  (:,k),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        for k=1:Nb
            plot(tv11,CA1 (k,:),sstr(2*k-1:2*k),'Color',cs(k));
        end;
        
        hold off;
        set(gca,'FontSize',fs);
        set(gca,'XLim',axx);
        xlabel('Time (y)');
        ylabel('Ca^2^+ (mmol kg^{-1})');
        Hl=legend(lstr,4);
        set(Hl,'FontSize',10);
        refresh;
    end;
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
    figure(21)
    plot(tv, FSichck,'m')
    xlabel('Time (y)');
    ylabel('Silicate flux increase');
    hold on
    plot(tv,Finchck,'k')
    hold off
    
    figure(100)
    clf;
    subplot(711)
    % bar(tv, Cinpc,'r')
    ylim([0 0.51])
    ylabel('Input (Pg C y^-^1) ')
    set(gca,'YTickLabel',[])
    set(gca,'XTickLabel',[])
    
    subplot(712)
    hold on
    hndl=plot(tv, FSichck,'m',tv,Finchck,'g');
    set(hndl,'LineWidth',1)
    legend('F_S_i','F_C')
    hndl1=plot(tv11, FSichck1(9,:),'m',tv11,Finchck1,'g');
    set(hndl1,'LineWidth',1)
    %xlabel('Time (y)');
    %ylabel('Silicate flux increase');
    %hold on
    %plot(tv,Finchck,'g','LineWidth',1)
    %plot(tv11,Finchck1(9,:),'g','LineWidth',1)
    hold off
    ylabel('F_C, F_S_i (mol C y^{-1})')
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[])
    
    subplot(713)
    box  on;
    hold on;
    hndl1=plot(tv  ,a  (:,9),'-b',tv  ,c  (:,9),'-r');
    set(hndl1,'LineWidth',1)
    legend('TA (DP)','TC (DP)')
    hndl=plot(tv11,a1 (9,:),'-b',tv11,c1 (9,:),'-r');
    set(hndl,'LineWidth',1)
    %legend('TA')
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[])
    % set(gca,'XTickLabel',[])
    % plot(tv,TCAR*1e3,'k',tv,TALK*1e3,'r');
    ylabel('TA, TC (mmol kg^{-1})')
    
    
    % clf;
    % box  on;
    % hold on;
    % plot(tv  ,a  (:,9),sstr(2*9-1:2*9),'Color',cs(9));
    % plot(tv11,a1 (9,:),sstr(2*9-1:2*9),'Color',cs(9));
    % hold off;
    % set(gca,'FontSize',fs);
    % set(gca,'XLim',axx);
    % xlabel('Time (y)');
    % ylabel('TA (mmol kg^{-1})');
    %
    %
    % clf;
    % box  on;
    % hold on;
    % plot(tv  ,c  (:,9),sstr(2*10-1:2*10),'Color',cs(10));
    % plot(tv11,c1 (9,:),sstr(2*10-1:2*10),'Color',cs(10));
    % hold off;
    % set(gca,'FontSize',fs);
    % set(gca,'XLim',axx);
    % ylabel('TA, TC (mmol/kg)')
    % set(gca,'XTickLabel',[])
    % Hl=legend(lstr,4);
    % set(Hl,'FontSize',10);
    % set(gca,'XTickLabel',[])
    % legend('TA','TC')
    % refresh;
    
    
    subplot(714);
    %clf;
    box  on;
    hold on;
    %for k=1:Nb
    plot(tv  ,co3tv (:,9)*1e6,sstr(2*9-1:2*9),'Color',cs(9));
    %end;
    %for k=1:Nb
    plot(tv11,co3tv1(9,:)*1e6,sstr(2*9-1:2*9),'Color',cs(9));
    %end;
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('[CO_3^{2-}] (\mumol kg^{-1})');
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    legend('CO_3^{2-}(DP)')
    set(gca,'XTickLabel',[])
    
    subplot(715);
    plot(tv  ,ccdP,'k-','LineWidth',2);
    hold on
    plot(tv11,ccdP1,'k-','LineWidth',2);
    hold off
    set(gca,'XDir','normal','YDir','reverse')
    ylabel('CCD (m)')
    set(gca,'XTickLabel',[])
    legend('Pacific')
    
    subplot(716);
    plot(tv  ,pco2t ,'r-','LineWidth',2);
    hold on;
    %for k=kkv
    %plot(tv(lt),pco2v(k),systr(2*k-1:2*k));
    %end;
    plot(tv11,pco2t1,'r-','LineWidth',2);
    hold off;
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('pCO_2 (\muatm)');
    
    %Ylim([-2*1e13 2*1e13])
    
    if(CAvflag~=0)
        subplot(717);plot(tv/(1e3),Ca*1e3)
        plot(tv  ,CA  (:,9),'-b','LineWidth',2);
        hold on
        plot(tv11,CA1 (9,:),'-b','LineWidth',2);
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[-50 0 50 100 150 200])
    xlabel('Time (ky)');
    ylabel('Ca^2^+ (mmol kg^{-1})');
    legend('Ca^2^+ (DP)')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    xlabel('Time (ky)')
    ylabel('[Ca^2^+] mmol/kg')
    refresh;
    
    figure(101)
    clf;
    
    subplot(711);
    plot(tv  ,pco2t ,'r-','LineWidth',2);
    hold on;
    %for k=kkv
    %plot(tv(lt),pco2v(k),systr(2*k-1:2*k));
    %end;
    plot(tv11,pco2t1,'r-','LineWidth',2);
    hold off;
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('pCO_2 (\muatm)');
    
    if(ftys)
        if(inorgF)
            title('end-Permian')
        else
            title('PETM')
        end
    else
        title('Modern')
    end
    
    subplot(712)
    box on
    hold on
    hndl=plot(tv, FSichck/.1e13,'m',tv,Finchck/.1e13,'g');
    set(hndl,'LineWidth',1)
    legend('F_S_i','F_C')
    hndl1=plot(tv11, FSichck1(9,:)/.1e13,'m',tv11,Finchck1/.1e13,'g');
    set(hndl1,'LineWidth',1)
    %xlabel('Time (y)');
    %ylabel('Silicate flux increase');
    %hold on
    %plot(tv,Finchck,'g','LineWidth',1)
    %plot(tv11,Finchck1(9,:),'g','LineWidth',1)
    hold off
    ylabel('F_C, F_S_i (10^{12} mol C y^{-1})')
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[])
    
    % MEAN SURFACE OCEAN alkalinity and total carbon over time
    if(ftys)
        aMean=(a(:,1).*V(1)+a(:,2).*V(2)+a(:,3)*V(3)+a(:,10)*V(10)+a(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
        cMean=(c(:,1).*V(1)+c(:,2).*V(2)+c(:,3)*V(3)+c(:,10)*V(10)+c(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
    else
        aMean=(a(:,1).*V(1)+a(:,2).*V(2)+a(:,3)*V(3)+a(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));
        cMean=(c(:,1).*V(1)+c(:,2).*V(2)+c(:,3)*V(3)+c(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));
    end
    
    aMean1 = [aMean(1)'  aMean(1)'];
    cMean1 = [cMean(1)'  cMean(1)'];
    
    subplot(713)
    box  on;
    hold on;
    %     hndl1=plot(tv  ,a  (:,3),'-b',tv  ,c  (:,3),'-r');
    hndl1=plot(tv  ,aMean,'-b',tv  ,cMean,'-r');
    set(hndl1,'LineWidth',1)
    legend('TA','TC ')
    %     hndl=plot(tv11,a1 (3,:),'-b',tv11,c1 (3,:),'-r');
    hndl=plot(tv11,aMean1,'-b',tv11,cMean1,'-r');
    set(hndl,'LineWidth',1)
    %legend('TA')
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[])
    % set(gca,'XTickLabel',[])
    % plot(tv,TCAR*1e3,'k',tv,TALK*1e3,'r');
    ylabel('TA, TC (mmol kg^{-1})')
    
    
    % MEAN SURFACE OCEAN carbonate ion conc. over time
    if(ftys)
        co3Mean=(co3tv(:,1).*V(1)+co3tv(:,2).*V(2)+co3tv(:,3)*V(3)+co3tv(:,10)*V(10)+co3tv(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
    else
        co3Mean=(co3tv(:,1).*V(1)+co3tv(:,2).*V(2)+co3tv(:,3)*V(3)+co3tv(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));
    end
    co3Mean1 = [co3Mean(1)'  co3Mean(1)'];
    
    subplot(714);
    %clf;
    box  on;
    hold on;
    %for k=1:Nb
    %     plot(tv  ,co3tv (:,3)*1e6,sstr(2*3-1:2*3),'Color',cs(3));
    plot(tv  ,co3Mean.*1e6,sstr(2*3-1:2*3),'Color',cs(3));
    %end;
    %for k=1:Nb
    %     plot(tv11,co3tv1(3,:)*1e6,sstr(2*3-1:2*3),'Color',cs(3));
    plot(tv11,co3Mean1*1e6,sstr(2*3-1:2*3),'Color',cs(3));
    %end;
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('[CO_3^{2-}] (\mumol kg^{-1})');
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    legend('mean surface ocean')
    set(gca,'XTickLabel',[])
    
    if(ftys)
        d13cbulkV = ((d13fcA(:,1)*BurialCA+d13fcP(:,1)*BurialCP+d13fcI(:,1)*BurialCI+d13fcT(:,1)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA)+...
            (d13fcA(:,2)*BurialCA+d13fcP(:,2)*BurialCP+d13fcI(:,2)*BurialCI+d13fcT(:,2)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA))/2;
        d44cabulkV = ((d44fcA(:,1)*BurialCA+d44fcP(:,1)*BurialCP+d44fcI(:,1)*BurialCI+d44fcT(:,1)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA)+...
            (d44fcA(:,2)*BurialCA+d44fcP(:,2)*BurialCP+d44fcI(:,2)*BurialCI+d44fcT(:,2)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA))/2;
    else
        d13cbulkV = ((d13fcA(:,1)*BurialCA+d13fcP(:,1)*BurialCP+d13fcI(:,1)*BurialCI)./(BurialCP+BurialCI+BurialCA)+...
            (d13fcA(:,2)*BurialCA+d13fcP(:,2)*BurialCP+d13fcI(:,2)*BurialCI)./(BurialCP+BurialCI+BurialCA))/2;
        d44cabulkV = ((d44fcA(:,1)*BurialCA+d44fcP(:,1)*BurialCP+d44fcI(:,1)*BurialCI)./(BurialCP+BurialCI+BurialCA)+...
            (d44fcA(:,2)*BurialCA+d44fcP(:,2)*BurialCP+d44fcI(:,2)*BurialCI)./(BurialCP+BurialCI+BurialCA))/2;
        
    end
    
    d13cbulkV1 = [d13cbulkV(1)'  d13cbulkV(1)'];
    d44cabulkV1 = [d44cabulkV(1)'  d44cabulkV(1)'];
    
    subplot(715);
    box on;
    hold on;
    myh1= plot(tv, d13cbulkV,'-b','LineWidth',2 );
    myh2= plot(tv11, d13cbulkV1,'-b','LineWidth',2 );
    myh3= plot((time-0.26)*1e6,d13cp,'o');
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('bulk \delta^{13}C');
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    legend([myh1 myh3], 'Mean shallow', 'data - Payne et al. 2010')
    set(gca,'XTickLabel',[])
    
    
    %     hold on
    %     plot(tv  ,ccdP,'k-','LineWidth',2);
    %     plot(tv  ,ccdA,'g-','LineWidth',2);
    %     plot(tv  ,ccdI,'r-','LineWidth',2);
    %
    %     if(ftys)
    %         plot(tv  ,ccdT,'b-','LineWidth',2);
    %     end
    %     plot(tv11,ccdP1,'k-','LineWidth',2);
    %     plot(tv11,ccdA1,'g-','LineWidth',2);
    %     plot(tv11,ccdI1,'r-','LineWidth',2);
    %     plot(tv11,ccdT1,'b-','LineWidth',2);
    %
    %     hold off
    %     set(gca,'XDir','normal','YDir','reverse')
    %     ylabel('CCD (m)')
    %     set(gca,'XTickLabel',[])
    %     legend('Pacific','Atlantic','Indic','Tethys')
    
    
    if(CAvflag~=0)
        subplot(716);
        box on;
        %         plot(tv  ,d44ca (:,9),sstr(2*9-1:2*9),'Color',cs(9));
        hold on
        myh1 = plot(tv, d44cabulkV,'-b','LineWidth',2);
        myh2 = plot(tv11,d44cabulkV1,'-b', 'LineWidth', 2);
        myh3 = plot((time-0.26)*1e6,d44cap,'o','Color','r');
        %         plot(tv11,d44ca1(9,:),sstr(2*9-1:2*9),'Color',cs(9));
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    ylabel('Ca^2^+ (mmol kg^{-1})');
    legend([myh1 myh3],'Mean shallow','data - Payne et al. 2010')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    set(gca,'XTickLabel',[])
    ylabel('bulk \delta^{44} Ca')
    
    
    
    
    if(CAvflag~=0)
        
        % MEAN SURFACE OCEAN calcium ion conc. over time
        if(ftys)
            caMean=(CA(:,1).*V(1)+CA(:,2).*V(2)+CA(:,3)*V(3)+CA(:,10)*V(10)+CA(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
        else
            caMean=(CA(:,1).*V(1)+CA(:,2).*V(2)+CA(:,3)*V(3)+CA(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));
        end
        caMean1 = [caMean(1)'  caMean(1)'];
        
        subplot(717);
        %         plot(tv  ,CA  (:,3),'-b','LineWidth',2);
        plot(tv  ,caMean,'-b','LineWidth',2);
        hold on
        %         plot(tv11,CA1 (3,:),'-b','LineWidth',2);
        plot(tv11,caMean1,'-b','LineWidth',2);
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %     set(gca,'XTickLabel',[-50 0 50 100 150 200 250 300 350 400])
    xlabel('Time (ky)');
    ylabel('Ca^2^+ (mmol kg^{-1})');
    legend('Ca^2^+')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    xlabel('Time (ky)')
    ylabel('[Ca^2^+] mmol/kg')
    
    
    refresh;
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
    
    if(0)
        %if(BlFlag == 2)
        dirstr = ' dat/FinalPERMsims/Sim1';
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
        save tv11.DAT    tv11    -ASCII -DOUBLE -TABS;
        save c.DAT     c     -ASCII -DOUBLE -TABS;
        save a.DAT     a     -ASCII -DOUBLE -TABS;
        save kkv.DAT   kkv   -ASCII -DOUBLE -TABS;
        save pco2t.DAT pco2t -ASCII -DOUBLE -TABS;
        save pco2v.DAT pco2v -ASCII -DOUBLE -TABS;
        save co3tv.DAT co3tv -ASCII -DOUBLE -TABS;
        save phtv.DAT  phtv  -ASCII -DOUBLE -TABS;
        save d13c.DAT  d13c  -ASCII -DOUBLE -TABS;
        save d13CA.DAT d13C  -ASCII -DOUBLE -TABS;
        
        
        save V.DAT V  -ASCII -DOUBLE -TABS;
        
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
        
        save FiN.DAT  Finchck  -ASCII -DOUBLE -TABS;
        save FSi.DAT  FSichck  -ASCII -DOUBLE -TABS;
        save FiN1.DAT  Finchck1  -ASCII -DOUBLE -TABS;
        save FSi1.DAT  FSichck1  -ASCII -DOUBLE -TABS;
        if(ftys)
            save shtT.DAT  shtT  -ASCII -DOUBLE -TABS;
        end;
        
        if(fsed)
            save fcA.DAT   fcA   -ASCII -DOUBLE -TABS;
            save fcI.DAT   fcI   -ASCII -DOUBLE -TABS;
            save fcP.DAT   fcP   -ASCII -DOUBLE -TABS;
            
            save d13fcA.DAT   d13fcA   -ASCII -DOUBLE -TABS;
            save d13fcI.DAT   d13fcI   -ASCII -DOUBLE -TABS;
            save d13fcP.DAT   d13fcP   -ASCII -DOUBLE -TABS;
            
            if(CAvflag>0)
                save d44ca.DAT  d44ca  -ASCII -DOUBLE -TABS;
                
                save d44fcA.DAT   d44fcA   -ASCII -DOUBLE -TABS;
                save d44fcI.DAT   d44fcI   -ASCII -DOUBLE -TABS;
                save d44fcP.DAT   d44fcP   -ASCII -DOUBLE -TABS;
                save CA.DAT CA  -ASCII -DOUBLE -TABS;
                
                save BurialCA.DAT BurialCA  -ASCII -DOUBLE -TABS;
                save BurialCI.DAT BurialCI  -ASCII -DOUBLE -TABS;
                save BurialCP.DAT BurialCP  -ASCII -DOUBLE -TABS;
                save BurialCT.DAT BurialCT  -ASCII -DOUBLE -TABS;
            end
            
            save ccdA.DAT  ccdA  -ASCII -DOUBLE -TABS;
            save ccdI.DAT  ccdI  -ASCII -DOUBLE -TABS;
            save ccdP.DAT  ccdP  -ASCII -DOUBLE -TABS;
            if(ftys)
                save fcT.DAT   fcT   -ASCII -DOUBLE -TABS;
                save d13fcT.DAT   d13fcT   -ASCII -DOUBLE -TABS;
                save ccdT.DAT  ccdT  -ASCII -DOUBLE -TABS;
                if(CAvflag>0)
                    save d44fcT.DAT   d44fcT   -ASCII -DOUBLE -TABS;
                end
            end;
        end;
        
        cmdstr(1,:) = 'move parms.tex';
        cmdstr(2,:) = 'move  vars.tex';
        cmdstr(3,:) = 'move     *.DAT';
        
        for k=1:3
            dostr(k,:)     = [cmdstr(k,:) dirstr];
            [status,result] = dos(dostr(k,:))
        end;
        
        dostr2 = ['copy' dirstr '/YPE10SedCAv1.DAT'];
        [status,result] = dos(dostr2)
        
    end; % save
    
end; % myflag



return;