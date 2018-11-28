%-------------------------------------
% function yp = LoscarDif10(t,y)
% provides difequations
%
% file: LoscarDif10.m
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
function yp = LoscarDifLT10(t,y)

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
    % known to Loscar10 and LoscarDif10
    %========================================

    global Aoc rho Nb V A TH0 TH TS mv mhd kasv kkv TCv Sv Pv gp tA tI ...
        phflag k1k2flag EPH fEPL rrain ep ept ec ect REDPC REDNC REDO2C eI ...
        Ns asvA asvI asvP zv m2kg m2kgca rhos frrf frrfca phic hs FiN ...
        phivtA phivtI phivtP phivtAca phivtIca phivtPca phi0 phi1 gam ...
        phiiA phiiI phiiP phiiAca phiiIca phiiPca phiiT phiiTca dsv ...
        Fpr frain rsedv rburvtA rburvtI rburvtP dissvtA ...
        dissvtI dissvtP FprtA FprtI FprtP FprtAca FprtIca FprtPca it co3satv ...
        Kd KS cst dsflag nc dYflag nu fsed ...
        Fpr13tA Fpr13tI Fpr13tP Fpr44tA Fpr44tI Fpr44tP ...
        Fint Fin13t Fin44t FSi44t FSit FSi13t ...
        fc0A fc0I fc0P fca0A fca0I fca0P f13c0A f13c0I f13c0P ...
        dissv13tA dissv13tI dissv13tP ...
        dissv44tA dissv44tI dissv44tP...
        co3s0 as zs0 klid nli Rst epsp epspca Rin FiN13 ...
        BlFlag CBl RBl kb FVC FVC13 Rvc pCSi nSi nCC Focb0 Fkg13 ...
        ftys nOC TT asvT kliT fc0T fca0T f13c0T rsedvtT rsedvtTca dissvtT dissvtTca dissv13tT dissv44tT ...
        FprtT FprtTca Fpr13tT Fpr44tT phivtT phivtTca fcon swcon ...
        TCv0 TCvt ntL ntH fsh fshI fshP fdpv fshT nshT ...
        ffflag tem em RlsCtv ...
        fdox KMMOX vask DTS DTS2 ts3 DTS3 DTS4 ...
        k1st tfinal TDflag omegCSvt omegASvt ...
        THt mv0 mhd0 oxA CAvflag FSichck Finchck kspCHCK ...
        CHECK1 CHECK2 CHECK3 CHKFin CHKFSi CHCKbioL CHCKbioI CHCKbioD...
        CHKcarB1 CHKcarB2 CHKcarB3 CHKcar1 CHKcar2 CHKcar3 CHKcar4 CAv Cadv Cinpc avCALC Voc...
        RinCA RcasAcheck RsAcheck f13cvAchck f44cavAchck fcvAchck CHKcarC1 CHKcarC2 CHKcarC3...
        CHKcarD1 CHKcarD2 CHKcarD3 CHKcarD4 Rbchck Rbcachck...
        rburvtAca rburvtIca rburvtPca dissvtAca dissvtIca dissvtPca EPLv EXLv ECLv...
        fwcv fwsv frkc fekc fdkc fGGi flakc...
        fwsv1 fwgv fmgv pco2gca tgc kgc fbgv fmcv FSi0 fbcv fbch fbbv...
         flakc ACT RUN gamma ws  pco20 pflag burCv FiN0 LTflag Rkg Focw0...
         ybbv ybv ycv runtime FERT epspv wcvtA wcvtI wcvtP wcvtT epspcaV inorgF epspF biopF Fpw0...
        fEPLv EPHv Rstca dincaV epsSensF Fopb0 Ffep0 Fcap0 oxicf0 anoxf0 Pfeed chck1st Focbv...
        EALvv ECALvv ECAHv EAHv dissCavv EPLvv oIv oIpv PPLvv Ffepv Fcapv counter capk PPHv Pscenario O0 Floegel totCexp0...
        rREG0 rREGv REDPCv REDPC0 meanSpco2v EPLv0 dbc alphagc Fpwv fT fT0 capk0...
        po4bf0 ocbf0 Q10 smoothcon beta eIv eIpv ocbf eIi Faom13 Fmeth13 Faom Fmeth;


    if(LTflag)
        % To run steady-state set the ones below to zero
        kgct=0; %t/1e6
        kgct1=0; %kgct
        DUM=0;      %1

    end


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
    if(LTflag && smoothcon)
        if(tgc>0 && tgc<23)
            TH = interp1([1 22],[TH0 0],tgc);
            TS = interp1([1 22],[0 TH0],tgc);
        elseif(tgc>=23)
            TH = 0.0;
            TS = TH0;
        else
%             TH = TH0;
%             TS = 0.0;
        end
    end



    %=========== Blast: Temp
%     if(LTflag)
        if(BlFlag == 2 & ftys == 1)
            if(t > tcon & t < tcon+t2bl)
%                 TCv = TCv0 + 4.0;     % +3.0 +4.0
                %if(t >= tcon & t < tcon+t2bl)
                %     TCv = flin(0.,3.e3,TCv0,TCv0+4.0,t);
                fsh  = 6.30;          % 6.5(1,1) 7.50/6.30(2,3) n7.3
                nshT  = 0.30;          % 0.5(1,1) 0.30/0.55(2,3).25
            else
%                 TCv = flin(t2bl,4*t2bl,TCv0+4.0,TCv0,t);
                fsh  = 4.50;           % 4.5(1,1) 4.50/4.0(2,3)
                nshT  = 0.40;           % 0.6(1,1) 0.35/0.6(2,3).3
            end;
        end;
%     end;

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
    if(CAvflag == 0)
        if(fdox)
            c    = y(     1:  Nb     ); % TCO2    (mol/m3)
            a    = y(  Nb+1:2*Nb     ); % TA      (mol/m3)
            p    = y(2*Nb+1:3*Nb     ); % PO4     (mol/m3)
            dox  = y(3*Nb+1:4*Nb     ); % O2      (mol/m3)
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
        else % fdox
            c    = y(     1:  Nb     ); % TCO2    (mol/m3)
            a    = y(  Nb+1:2*Nb     ); % TA      (mol/m3)
            p    = y(2*Nb+1:3*Nb     ); % PO4     (mol/m3)
            cc   = y(3*Nb+1:4*Nb     ); % T13CO2  (mol/m3)
            C    = y(4*Nb+1          ); % Catm    (mol/m2)
            CC   = y(4*Nb+2          ); % 13Catm  (mol/m2)
            if(fsed)
                fcvA  = y(4*Nb+3     :4*Nb+2+1*Ns); % %calc
                fcvI  = y(4*Nb+3+1*Ns:4*Nb+2+2*Ns); % %calc
                fcvP  = y(4*Nb+3+2*Ns:4*Nb+2+3*Ns); % %calc
                %====== Carbon-13
                f13cvA= y(4*Nb+3+3*Ns:4*Nb+2+4*Ns); % %calc
                f13cvI= y(4*Nb+3+4*Ns:4*Nb+2+5*Ns); % %calc
                f13cvP= y(4*Nb+3+5*Ns:4*Nb+2+6*Ns); % %calc
                if(ftys)
                    fcvT  = y(4*Nb+3+6*Ns:4*Nb+2+7*Ns); % %calc
                    f13cvT= y(4*Nb+3+7*Ns:4*Nb+2+8*Ns); % %calc
                end;
            end;
        end; % fdox
    else
        if(fdox)
            c    = y(     1:  Nb     ); % TCO2    (mol/m3)
            a    = y(  Nb+1:2*Nb     ); % TA      (mol/m3)
            p    = y(2*Nb+1:3*Nb     ); % PO4     (mol/m3)
            CA   = y(3*Nb+1:4*Nb   );
            dox  = y(4*Nb+1:5*Nb     ); % O2      (mol/m3)
            cac   = y(5*Nb+1:6*Nb     ); % 44Ca  (mol/m3)
            cc   = y(6*Nb+1:7*Nb     ); % T13CO2  (mol/m3)
            C    = y(7*Nb+1          ); % Catm    (mol/m2)
            CC   = y(7*Nb+2          ); % 13Catm  (mol/m2)
            if(fsed)
                fcvA  = y(7*Nb+3     :7*Nb+2+1*Ns); % %calc
                fcvI  = y(7*Nb+3+1*Ns:7*Nb+2+2*Ns); % %calc
                fcvP  = y(7*Nb+3+2*Ns:7*Nb+2+3*Ns); % %calc
                %====== Carbon-13
                f13cvA= y(7*Nb+3+3*Ns:7*Nb+2+4*Ns); % %calc
                f13cvI= y(7*Nb+3+4*Ns:7*Nb+2+5*Ns); % %calc
                f13cvP= y(7*Nb+3+5*Ns:7*Nb+2+6*Ns); % %calc
                fcavA=  y(7*Nb+3+6*Ns:7*Nb+2+7*Ns); % %calc
                fcavI=  y(7*Nb+3+7*Ns:7*Nb+2+8*Ns); % %calc
                fcavP=  y(7*Nb+3+8*Ns:7*Nb+2+9*Ns); % %calc
                f44cavA=y(7*Nb+3+9*Ns:7*Nb+2+10*Ns); % %calc
                f44cavI=y(7*Nb+3+10*Ns:7*Nb+2+11*Ns); % %calc
                f44cavP=y(7*Nb+3+11*Ns:7*Nb+2+12*Ns); % %calc
                if(ftys)
                    fcvT  = y(7*Nb+3+12*Ns:7*Nb+2+13*Ns); % %calc
                    f13cvT= y(7*Nb+3+13*Ns:7*Nb+2+14*Ns); % %calc
                    fcavT=  y(7*Nb+3+14*Ns:7*Nb+2+15*Ns); % %calc
                    f44cavT= y(7*Nb+3+15*Ns:7*Nb+2+16*Ns); % %calc
                end;
            end;
        else % fdox
            c    = y(     1:  Nb     ); % TCO2    (mol/m3)
            a    = y(  Nb+1:2*Nb     ); % TA      (mol/m3)
            p    = y(2*Nb+1:3*Nb     ); % PO4     (mol/m3)
            CA   = y(3*Nb+1:4*Nb     );
            cac   = y(4*Nb+1:5*Nb     ); % T13CO2  (mol/m3)
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
        end; % fdox
    end; %CAvflag

    if(fsed)
        if(ftys)
            fcVV  = [fcvA fcvI fcvP fcvT f13cvA f13cvI f13cvP f13cvT];% f44cavA f44cavI f44cavP ];
        else
            fcVV  = [fcvA fcvI fcvP f13cvA f13cvI f13cvP];% f44cavA f44cavI f44cavP ];
        end

        if(ERR & ~isempty(find(fcVV < 0.)) )
            disp('+++ WARNING: fc is negative +++');
            if(0)
                myeps = 1.e-6;
                kzA   = find(fcvA < 0.);
                kzI   = find(fcvI < 0.);
                kzP   = find(fcvP < 0.);
                kzT   = find(fcvT < 0.);
                fcvA(kzA)  = myeps;
                fcvI(kzI)  = myeps;
                fcvP(kzP)  = myeps;
                fcvT(kzT)  = myeps;
                k13zA   = find(f13cvA < 0.);
                k13zI   = find(f13cvI < 0.);
                k13zP   = find(f13cvP < 0.);
                k13zT   = find(f13cvT < 0.);
                f13cvA(kzA)  = myeps;
                f13cvI(kzI)  = myeps;
                f13cvP(kzP)  = myeps;
                f13cvT(kzT)  = myeps;
                k44zA   = find(f44cavA < 0.);
                k44zI   = find(f44cavI < 0.);
                k44zP   = find(f44cavP < 0.);
                k44zT   = find(f44cavT < 0.);
                f44cavA(kzA)  = myeps;
                f44cavI(kzI)  = myeps;
                f44cavP(kzP)  = myeps;
                f44cavT(kzT)  = myeps;
            end;
        end;
    end;

    pco2a  =   C/(2.2e15/12/Aoc);
    pcco2a =  CC/(2.2e15/12/Aoc);

    % Deep-sea temperature sensitivity (Archer, 2005)
    % Modern
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
    %LPEE
    if(TDflag == 1 & dYflag == 0 & ftys == 1)
        DTeq = 3.0*log(pco2a/pCSi)/log(2); % T sens 1.5,3,4.5

        if(kt > k1st) % first call with new kt
            f1st = 1;
        else
            f1st = 0;
        end
        if(f1st == 1) % update TCv
            if(ftys)
                kvi   = [4 5 6 12];
                kvd   = [7 8 9 13];
            else
                kvi   = [4 5 6];
                kvd   = [7 8 9];
            end
            %tlag   = [10. 200. 1000.];
            %             tlag   = [10     100      1000].*1e2;

            tlag   = [1     1      1].*1e3;
            TCv(kkv)  = TCv(kkv)+...
                ((TCv0(kkv)+DTeq)-TCv(kkv))*Dtst/tlag(1); % 10   surf
            TCv(kvi)  = TCv(kvi)+...
                ((TCv0(kvi)+DTeq)-TCv(kvi))*Dtst/tlag(2); % 200  intmd
            TCv(kvd)  = TCv(kvd)+...
                ((TCv0(kvd)+DTeq)-TCv(kvd))*Dtst/tlag(3); % 1000 deep
            f1st = 0;
            k1st = kt;
        end
        %         TCvt(1,:)=[25    25    25    16    16    16     12     12     12    12    18    14    12];
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
    % if(TDflag == 1 & dYflag == 1)
    % TX    = TCvt(it,:);
    % else
    TX    = TCv;
    % end

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
    elseif(CAvflag == 1)
        for k=1:Nb
            [co2(k),pco2(k),co3(k),ph(k),kh(k),o2(k)] = ...
                dafunPECA(c(k)/rho,a(k)/rho,TX(k),Sv(k),Pv(k),Ca,Mg);
        end;
    elseif(CAvflag == 0)
        for k=1:Nb
            [co2(k),pco2(k),co3(k),ph(k),kh(k),o2(k)] = ...
                dafunPE(c(k)/rho,a(k)/rho,TX(k),Sv(k),Pv(k),Ca,Mg);
        end;
    end;


    if(1)
        %====== Carbon-13
        Rb  = cc./c;
        if(CAvflag>0)
            Rbca=cac./CA;
        end;
%         Rbchck(kt+1,:)=Rb;
        if(CAvflag>0)
%             Rbcachck(kt+1,:)=Rbca;
        end;
        Tk  = TX + 273.15;
        edb = -9866./Tk + 24.12;
        edg = - 373./Tk +  0.19;
        adb = edb/1e3+1;
        adg = edg/1e3+1;
        au  = 0.9995;

        pcco2(kkv) = adb(kkv).*Rb(kkv)'.*pco2(kkv); % (adg^-1 dropped)
    end;
    % if running long term non steady-state
    if(LTflag && DUM)
        epsp=(interp1(1:length(alphagc),-alphagc,tgc-kgct));
%         epsp=(interp1(1:length(acv),-acv,tgc-kgct));
        %     epsp=flin(0,runtime,-acv(tgc),-acv(tgc),t);
    end
    epspv(kt)=epsp;


    % For PETM not Permian

    if(t==0)
        %        testA=0.0;
        %    testI=0.0;
        %    testP=0.0;
        %    testT=0.0;
        %     EPLv(1)=1e-40;
        %     EPLv(2)=1e-40;
        %     EPLv(3)=1e-40;
        %     EPLv(4)=1e-40;
        %     fEPL = 0;
        %     EPH = 0;
        %     rrain  = 4;
    end

    if(biopF)
        if(t<=DTS)
            fEPL = flin(0,DTS, 0.95, 0.0, t);
            EPH = flin (0,DTS, 2*1.8*A(10),0,t);
            %      Focb0 = flin (0,50e3, 1.0*5.e12/Aoc,0.0*5.e12/Aoc,t);
        end
        if(t>DTS)
            fEPL = 0 + 0.95*(t-DTS)/(4e5-DTS);
            EPH =  0 + 2*1.8*A(10)*(t-DTS)/(4e5-DTS);
            %      Focb0 = flin (0,50e3, 1.0*5.e12/Aoc,0.0*5.e12/Aoc,t);
        end

    end

    fEPLv(kt) = fEPL;
    EPHv(kt) = EPH;

    %============== Biological Pump ========%
    %
    % Low Lat Export Corg
    if(Floegel)
        meanSpco2=(pco2(1)*V1+pco2(2)*V2+pco2(3)*V3+pco2(10)*V10+pco2(11)*V11)/(V1+V2+V3+V10+V11);
        REDPC = 1/(0.06*meanSpco2+88);%1.0068e+003
%         REDPC = flin(0,3e5, 1/147.8592, 1/190,  t);
        meanSpco2v(kt) =meanSpco2;
    end
    
    REDPCv(kt) = REDPC;
    EPLv    = fEPL*mv(1:3)'.*p(4:6)/REDPC; % (m3/y*mol/m3 = mol/y)
    if(ftys)
        EPLv(4) = fEPL*mv(004)'.*p(012)/REDPC*fT/fT0;%*A(11)/2.7920e+13;%*flin(25,59,0.01,fT0,tgc)/fT0; % (m3/y*mol/m3 = mol/y)
    end;
    if(counter == 0)
       EPLv0 = EPLv; 
%        display('initialized')
    end
    if(1)
        lk = 1;
        if(CAvflag>0)
            for k=kkv
                [kspcS(k),kspaS(k)] = ...
                    kspfunCA(TX(k),Sv(k),Pv(k),CA(k)*1e-3,Mg);
                omegCSedv(lk) = co3(k)*(CA(k)*1e-3)/kspcS(k);
                %         omegASv = co3(k)*(CA(k)*1e-3)/kspaS(k);
                lk = lk + 1;
            end;
        else
            for k=kkv
                [kspcS(k),kspaS(k)] = ...
                    kspfunCA(TX(k),Sv(k),Pv(k),Ca,Mg);
                omegCSedv(lk) = co3(k)*(Ca)/kspcS(k);
                %         omegASv = co3(k)*(CA(k)*1e-3)/kspaS(k);
                lk = lk + 1;
            end;

        end
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

    if(ftys)
        if(inorgF)
            kcarb=0.031; % mol C /m2/y 0.031
            kcarb0=kcarb;
        else
            kcarb=0; %
            kcarb0=kcarb;
        end

    else
        kcarb=0; %
        kcarb0=kcarb;
    end
    etta = 2.0; %2.0, 1.7
    % omegCSedv(1)=2.1031;
    % omegCSedv(2)=2.1141;
    % omegCSedv(3)=2.1035;
    % omegCSedv(4)=1.1923;
    % omegCSedv(5)=1.6344;
    if(ftys)
        if(inorgF)
            if(omegCSedv(1)>=1)
                testA = kcarb*(omegCSedv(1)-1)^etta;
            else
                testA = kcarb0;
            end
            if(omegCSedv(2)>=1)
                testI = kcarb*(omegCSedv(2)-1)^etta;
            else
                testI = kcarb0;
            end
            if(omegCSedv(3)>=1)
                testP = kcarb*(omegCSedv(3)-1)^etta;
            else
                testP = kcarb0;
            end
            if(omegCSedv(5)>=1)
                testT = kcarb*(omegCSedv(5)-1)^etta;
            else
                testT = kcarb0;
            end

        else
            testA = 0;
            testI = 0;
            testP = 0;
            testT = 0;
        end
    else
        testA = 0;
        testI = 0;
        testP = 0;
        testT = 0;
    end

    if(ftys)

    end

% if(t>=60e3 && t<=100e3)
%    fcb=0.7; 
% else
%    fcb=1;  
% end


    EALv    = 2*EPLv/rrain*fcb+2*testA*Aoc;    % ALK export.

    EALv(1) = 2*EPLv(1)/rrain*fcb+2*testA*A(1);  %*Aoc/nOC
    EALv(2) = 2*EPLv(2)/rrain*fcb+2*testI*A(2);
    EALv(3) = 2*EPLv(3)/rrain*fcb+2*testP*A(3);
    if(ftys)
        EALv(4) = 2*EPLv(4)/rrain*fcb+2*testT*A(11);
    end
    
    if(Pscenario==13)
       EALv    = 2*EPLv0/rrain*fcb+2*testA*Aoc;    % ALK export.

    EALv(1) = 2*EPLv0(1)/rrain*fcb+2*testA*A(1);  %*Aoc/nOC
    EALv(2) = 2*EPLv0(2)/rrain*fcb+2*testI*A(2);
    EALv(3) = 2*EPLv0(3)/rrain*fcb+2*testP*A(3);
    if(ftys)
        EALv(4) = 2*EPLv0(4)/rrain*fcb+2*testT*A(11);
    end  
    end



    PPLv   = EPLv*REDPC;      % PO4
    ENLv   = EPLv*REDNC;      % NO3
    EOLv   = EPLv*REDO2C;     % O2
    ECALv  = EALv./2;         % Ca THIS ADDED!!!!! Ca transport

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




    kL = [1:1:3];
    if(ftys)
        kL = [kL 11];
    end;
    %======  Carbon-13
    alp    = epsp/1e3+1;
    if(LTflag)
        alpC   = -0.5/1e3+1;
    else
        alpC   = -0.0/1e3+1;
    end
    % Corg
    EPLvCC = alp*Rb(kL).*EPLv;
    EPHCC  = alp*Rb(10) *EPH;
    EALvCC =   alpC*  Rb(kL).*EALv;
    EAHCC  =   alpC*  Rb(10) *EAH;
    % Ctotal: Corg+CaCO3
    ECLvCC = EPLvCC+EALvCC/2;
    ECHCC  = EPHCC +EAHCC /2;




    % if(t<=4e4)
    %  epspca=flin(0,4e4,-1.0328,-1.13,t);
    % end
    % if(t>4e4 && t<=8e4)
    %  epspca=flin(4e4,8e4,-1.13,-1.0328,t);
    % else
    %     epspca = -1.0328;
    % end

    tStart=0e4; %5e4
    tEnd=5e4; %10e4
    % if(t>tStart && t<tEnd)
    %     dinca = flin(tStart,tEnd, -0.6, -0.2, t);
    % %     epspca = flin(tStart,tEnd, -1.4, -0.7, t);
    % elseif(t>=tEnd)
    %     dinca = -0.2;
    % %     epspca = -0.7;
    % else
    %     dinca = -0.6;
    % %     epspca = -1.4;
    % end

    % if(t<5e4)
    % dinca = -0.9;
    % else
    % dinca = -0.6;
    % end
    %  dincaV(kt) = dinca;
    % RinCA    = Rstca*(dinca/1e3+1);
    if(inorgF  )
        if(epsSensF==0)
            co3Mean=(co3(1).*V(1)+co3(2).*V(2)+co3(3)*V(3)+co3(10)*V(10)+co3(11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
            epspca=-1.31+3.69*co3Mean*1e3;
        elseif(epsSensF==1)
            if(t>tStart && t<=tEnd)
                epspca = flin(tStart,tEnd, -0.9937, -1.1611, t);
            elseif(t>tEnd)
                epspca = flin(tEnd,runtime, -1.1611, -0.9937, t);
                %         epspca = -1.1611;
            else
                epspca = -0.9937;
            end
        elseif(epsSensF==2)
            if(t>tStart && t<tEnd)
                epspca = flin(tStart,tEnd, -0.5773, -0.7447, t);
            elseif(t>=tEnd)
                epspca = flin(tEnd,runtime, -0.7447, -0.5773, t);
            else
                epspca = -0.5773;
            end
        end
    end
    epspcaV (kt) = epspca;
    %======  Calcium-44
    alpca    = epspca/1e3+1;
    % Corg
    %EPLvCAC = alp*Rb(kL).*EPLv;
    %EPHCAC  = alp*Rb(10) *EPH;
    if(CAvflag>0)
        EALvCAC =   1*  alpca*Rbca(kL).*EALv;
        EAHCAC  =     1*alpca*Rbca(10) *ECAH;
        % Ca total: CaCO3
        ECLvCAC = EALvCAC/2;%-alpca*Rbca(kL).*ENLv.*0.01/2;
        ECHCAC  = EAHCAC /2;
    end;
    % calcite export
    EAA   = EALv(1)+EAH*gp(7); % +EAH (no H seds, see below)
    EAI   = EALv(2)+EAH*gp(8);
    EAP   = EALv(3)+EAH*gp(9);


    if(CAvflag>0)
        EAACA  = EALv(1)+ECAH*gp(7); % +EAH (no H seds, see below)
        EAICA   = EALv(2)+ECAH*gp(8);
        EAPCA   = EALv(3)+ECAH*gp(9);
    end;

    % Fpr is rain to sediments, not export
    % Fpr = export - water column diss



    FprA  = (1-nu)*EAA/2/A(1); % mol C/m2/y
    FprI  = (1-nu)*EAI/2/A(2); % mol C/m2/y
    FprP  = (1-nu)*EAP/2/A(3); % mol C/m2/y
    if(CAvflag>0)
        FprAca  = (1-nu)*EAACA/2/A(1); % mol C/m2/y
        FprIca  = (1-nu)*EAICA/2/A(2); % mol C/m2/y
        FprPca  = (1-nu)*EAPCA/2/A(3); % mol C/m2/y
    end;
    %====== Carbon-13
    EAACC   = EALvCC(1)+EAHCC*gp(7); % +EAH (no H seds, above)
    EAICC   = EALvCC(2)+EAHCC*gp(8);
    EAPCC   = EALvCC(3)+EAHCC*gp(9);
    Fpr13A  = (1-nu)*EAACC/2/A(1); % mol C/m2/y
    Fpr13I  = (1-nu)*EAICC/2/A(2); % mol C/m2/y
    Fpr13P  = (1-nu)*EAPCC/2/A(3); % mol C/m2/y

    %====== Calcium-44
    if(CAvflag>0)
        EAACAC   = EALvCAC(1)+EAHCAC*gp(7); % +EAH (no H seds, above)
        EAICAC   = EALvCAC(2)+EAHCAC*gp(8);
        EAPCAC   = EALvCAC(3)+EAHCAC*gp(9);
        Fpr44A  = (1-nu)*EAACAC/2/A(1); % mol C/m2/y
        Fpr44I  = (1-nu)*EAICAC/2/A(2); % mol C/m2/y
        Fpr44P  = (1-nu)*EAPCAC/2/A(3); % mol C/m2/y
    end;


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

            frrfAca = frrfca*ones(1,Ns)*frem(1); %  kg  /m2/y
            frrfIca = frrfca*ones(1,Ns)*frem(2); %  kg  /m2/y
            frrfPca = frrfca*ones(1,Ns)*frem(3); %  kg  /m2/y

            % carbonate rain to sediments AIP
            FprvA = FprA*ones(1,Ns); % mol C/m2/y
            FprvI = FprI*ones(1,Ns); % mol C/m2/y
            FprvP = FprP*ones(1,Ns); % mol C/m2/y

            % FprvA(3:Ns)=0.000001;
            % FprvA(3:Ns)=0.000001;
            % FprvA(3:Ns)=0.000001;

            if(CAvflag>0)
                FprvAca = FprAca*ones(1,Ns); % mol C/m2/y
                FprvIca = FprIca*ones(1,Ns); % mol C/m2/y
                FprvPca = FprPca*ones(1,Ns); % mol C/m2/y
            end;

            Fpr13vA = Fpr13A*ones(1,Ns);
            Fpr13vI = Fpr13I*ones(1,Ns);
            Fpr13vP = Fpr13P*ones(1,Ns);

            if(CAvflag>0)
                Fpr44vA = Fpr44A*ones(1,Ns);
                Fpr44vI = Fpr44I*ones(1,Ns);
                Fpr44vP = Fpr44P*ones(1,Ns);
            end;

            % total C rain AIP
            FAA    = FprA*A(1);
            FII    = FprI*A(2);
            FPP    = FprP*A(3);
            jj     = 2;
            fdpv(1)= ( FAA ...
                -fsh*sum(FprvA(1:jj   ).*asvA(1:jj   )*A(1)) ) ...
                ./    sum(FprvA(jj+1:Ns).*asvA(jj+1:Ns)*A(1))  ;
            fdpv(2)= ( FII ...
                -fshI*sum(FprvI(1:jj   ).*asvI(1:jj   )*A(2)) ) ...
                ./    sum(FprvI(jj+1:Ns).*asvI(jj+1:Ns)*A(2))  ;
            fdpv(3)= ( FPP ...
                -fshP*sum(FprvP(1:jj   ).*asvP(1:jj   )*A(3)) ) ...
                ./    sum(FprvP(jj+1:Ns).*asvP(jj+1:Ns)*A(3))  ;

            fshvA = [fsh*ones(1,jj) fdpv(1)*ones(1,Ns-jj)];
            fshvI = [fshI*ones(1,jj) fdpv(2)*ones(1,Ns-jj)];
            fshvP = [fshP*ones(1,jj) fdpv(3)*ones(1,Ns-jj)];

            % total CA rain AIP
            if(CAvflag>0)
                FAAca    = FprAca*A(1);
                FIIca    = FprIca*A(2);
                FPPca    = FprPca*A(3);
                jj     = 2;
                fdpvca(1)= ( FAAca ...
                    -fsh*sum(FprvAca(1:jj   ).*asvA(1:jj   )*A(1)) ) ...
                    ./    sum(FprvAca(jj+1:Ns).*asvA(jj+1:Ns)*A(1))  ;
                fdpvca(2)= ( FIIca ...
                    -fshI*sum(FprvIca(1:jj   ).*asvI(1:jj   )*A(2)) ) ...
                    ./    sum(FprvIca(jj+1:Ns).*asvI(jj+1:Ns)*A(2))  ;
                fdpvca(3)= ( FPPca ...
                    -fshP*sum(FprvPca(1:jj   ).*asvP(1:jj   )*A(3)) ) ...
                    ./    sum(FprvPca(jj+1:Ns).*asvP(jj+1:Ns)*A(3))  ;

                fshvAca = [fsh*ones(1,jj) fdpvca(1)*ones(1,Ns-jj)];
                fshvIca = [fshI*ones(1,jj) fdpvca(2)*ones(1,Ns-jj)];
                fshvPca = [fshP*ones(1,jj) fdpvca(3)*ones(1,Ns-jj)];
            end;
            frrfA = fshvA.*frrfA;
            frrfI = fshvI.*frrfI;
            frrfP = fshvP.*frrfP;


            AfcA = zeros(1,Ns);
            AfcI = zeros(1,Ns);
            AfcP = zeros(1,Ns);

            AfcA(1) = (A(1)/(asvA(1)*A(1)));
            AfcA(2) = (A(1)/(asvA(2)*A(1)));
            AfcI(1) = (A(2)/(asvI(1)*A(2)));
            AfcI(2) = (A(2)/(asvI(2)*A(2)));
            AfcP(1) = (A(3)/(asvP(1)*A(3)));
            AfcP(2) = (A(3)/(asvP(2)*A(3)));

            % AfcA(jj+1:Ns) = ones(jj+1:Ns).*1;
            % AfcI(jj+1:Ns) = ones(jj+1:Ns).*1;
            % AfcP(jj+1:Ns) = ones(jj+1:Ns).*1;

            if(inorgF)
                FprvA = fshvA.*FprvA.*AfcA;
                FprvI = fshvI.*FprvI.*AfcI;
                FprvP = fshvP.*FprvP.*AfcP;
            else
                FprvA = fshvA.*FprvA;
                FprvI = fshvI.*FprvI;
                FprvP = fshvP.*FprvP;
            end

            if(CAvflag>0)
                frrfAca = fshvAca.*frrfAca;
                frrfIca = fshvIca.*frrfIca;
                frrfPca = fshvPca.*frrfPca;

                if(inorgF)
                    FprvAca = fshvAca.*FprvAca.*AfcA;
                    FprvIca = fshvIca.*FprvIca.*AfcI;
                    FprvPca = fshvPca.*FprvPca.*AfcP;
                else
                    FprvAca = fshvAca.*FprvAca;
                    FprvIca = fshvIca.*FprvIca;
                    FprvPca = fshvPca.*FprvPca;
                end
            end;

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

            if(inorgF)
                Fpr13vA = fshvA.*Fpr13vA.*AfcA;
                Fpr13vI = fshvI.*Fpr13vI.*AfcI;
                Fpr13vP = fshvP.*Fpr13vP.*AfcP;
            else
                Fpr13vA = fshvA.*Fpr13vA;
                Fpr13vI = fshvI.*Fpr13vI;
                Fpr13vP = fshvP.*Fpr13vP;
            end
            if(CAvflag>0)
                if(inorgF)
                    Fpr44vA = fshvAca.*Fpr44vA.*AfcA;
                    Fpr44vI = fshvIca.*Fpr44vI.*AfcI;
                    Fpr44vP = fshvPca.*Fpr44vP.*AfcP;
                else
                    Fpr44vA = fshvAca.*Fpr44vA;
                    Fpr44vI = fshvIca.*Fpr44vI;
                    Fpr44vP = fshvPca.*Fpr44vP;
                end
            end;
            % sum(FprvA.*asvA*A(1))
            % FAA
            % sum(FprvI.*asvI*A(2))
            % FII
            % sum(FprvP.*asvP*A(3))
            % FPP
            % return
            % test: FprvA; sum(FprvA.*asvA*A(1)); FAA;
        end;

        if(ftys)
            EAT    = EALv  (4);

            EATCC  = EALvCC(4);

            FprT   = (1-nu)*EAT  /2/A(11); % mol C/m2/y
            Fpr13T = (1-nu)*EATCC/2/A(11); % mol C/m2/y
            

            if(CAvflag>0)
                ECAT=EALv(4);
                FprTca   = (1-nu)*ECAT /2 /A(11); % mol C/m2/y
                EATCAC  = EALvCAC(4);
                FprvTca  = FprTca*ones(1,Ns);      % mol C/m2/y
                Fpr44T = (1-nu)*EATCAC/2/A(11); % mol C/m2/y
                FprvT  = FprT*ones(1,Ns);      % mol C/m2/y
            end;
            if(shlf)
                fremT  = 6.0;                  % 4. !#!#!#!#!#! clay
                frrfT  = frrf*ones(1,Ns)*fremT;%  kg  /m2/y
                FprvT  = FprT*ones(1,Ns);      % mol C/m2/y
                Fpr13vT= Fpr13T*ones(1,Ns);
                if(CAvflag>0)
                    frrfTca  = frrfca*ones(1,Ns)*fremT;%  kg  /m2/y
                    FprvTca  = FprTca*ones(1,Ns);      % mol C/m2/y
                    Fpr44vT= Fpr44T*ones(1,Ns);
                end;
                fshT   = (fsh)^nshT;           % fshT < fsh !!! 0.6
                FTT    = FprT*A(11);
                fdpv(4)= ( FTT ...
                    -fshT*sum(FprvT(1:jj   ).*asvT(1:jj   )*A(11)) ) ...
                    ./    sum(FprvT(jj+1:Ns).*asvT(jj+1:Ns)*A(11))  ;

                fshvT  = [fshT*ones(1,jj) fdpv(4)*ones(1,Ns-jj)];
                frrfT  = fshvT.*frrfT;
                AfcT = zeros(1,Ns);
                AfcT(1) = (A(11)/(asvT(1)*A(11)));
                AfcT(2) = (A(11)/(asvT(2)*A(11)));
                % AfcT(jj+1:Ns) = ones(jj+1:Ns).*1;
                if(inorgF)
                    FprvT  = fshvT.*FprvT.*AfcT;
                else
                    FprvT  = fshvT.*FprvT;
                end
                % sum(FprvT.*asvT*A(11))
                % FTT
                % return
                if(inorgF)
                    Fpr13vT= fshvT.*Fpr13vT.*AfcT;
                else
                    Fpr13vT= fshvT.*Fpr13vT;
                end
                if(CAvflag>0)
                    if(inorgF)
                        Fpr44vT= fshvT.*Fpr44vT.*AfcT;
                    else
                        Fpr44vT= fshvT.*Fpr44vT;
                    end
                    FTTca    = FprTca*A(11);
                    fdpvca(4)= ( FTTca ...
                        -fshT*sum(FprvTca(1:jj   ).*asvT(1:jj   )*A(11)) ) ...
                        ./    sum(FprvTca(jj+1:Ns).*asvT(jj+1:Ns)*A(11))  ;

                    fshvTca  = [fshT*ones(1,jj) fdpvca(4)*ones(1,Ns-jj)];
                    frrfTca  = fshvTca.*frrfTca;
                    if(inorgF)
                        FprvTca  = fshvTca.*FprvTca.*AfcT;
                    else
                        FprvTca  = fshvTca.*FprvTca;
                    end
                end;
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
            TdvA = TX(klid);
            TdvI = TX(klid+1);
            TdvP = TX(klid+2);
            SdvA = Sv(klid);
            SdvI = Sv(klid+1);
            SdvP = Sv(klid+2);

            % CAdvA = CA(klid);
            % CAdvI = CA(klid+1);
            % CAdvP = CA(klid+2);
            if(ftys)
                TdvT = TX(kliT);
                SdvT = Sv(kliT);
                % CAdvT = CA(kliT);
            end
            zsatv = dsv;

            %==== saturation at sediment levels
            if(CAvflag == 2)
                for k=1:Ns
                    for l=1:Nb
                        CALC(l)=sum(CA(l)*V/Voc);
                    end
                    avCALC=sum(CALC.*V/Voc);

                    if(ftys)
                        [kspcA(k),tmp] = ...
                            kspfun(TdvA(k),SdvA(k),zsatv(k)/10.,CA(k)*1e-3,Mg);
                        [kspcI(k),tmp] = ...
                            kspfun(TdvI(k),SdvI(k),zsatv(k)/10.,CA(k)*1e-3,Mg);
                        [kspcP(k),tmp] = ...
                            kspfun(TdvP(k),SdvP(k),zsatv(k)/10.,CA(k)*1e-3,Mg);
                        [kspcT(k),tmp] = ...
                            kspfun(TdvT(k),SdvT(k),zsatv(k)/10.,CA(k)*1e-3,Mg);
                    else
                        [kspcA(k),tmp] = ...
                            kspfun(TdvA(k),SdvA(k),zsatv(k)/10.,CA(1)*1e-3,Mg);
                        [kspcI(k),tmp] = ...
                            kspfun(TdvI(k),SdvI(k),zsatv(k)/10.,CA(1)*1e-3,Mg);
                        [kspcP(k),tmp] = ...
                            kspfun(TdvP(k),SdvP(k),zsatv(k)/10.,CA(1)*1e-3,Mg);
                    end
                end;
                %     co3satvA  = kspcA/Ca;
                %     co3satvI  = kspcI/Ca;
                %     co3satvP  = kspcP/Ca;
                %     if(ftys)
                %         co3satvT  = kspcT/Ca;
                %     end
                co3satvA  = kspcA./(avCALC*1e-3);
                co3satvI  = kspcI./(avCALC*1e-3);
                co3satvP  = kspcP./(avCALC*1e-3);
                if(ftys)
                    co3satvT  = kspcT./(avCALC*1e-3);
                end
                %     for k=1:Ns
                %         for l=1:Nb
                %             CALC(l)=sum(CA(l)*V/Voc);
                %         end
                %         avCALC=sum(CALC.*V/Voc);
                %         [kspcA(k),tmp] = ...
                %             kspfun(TdvA(k),SdvA(k),zsatv(k)/10.,avCALC*1e-3,Mg);
                %         [kspcI(k),tmp] = ...
                %             kspfun(TdvI(k),SdvI(k),zsatv(k)/10.,avCALC*1e-3,Mg);
                %         [kspcP(k),tmp] = ...
                %             kspfun(TdvP(k),SdvP(k),zsatv(k)/10.,avCALC*1e-3,Mg);
                %         if(ftys)
                %             [kspcT(k),tmp] = ...
                %                 kspfun(TdvT(k),SdvT(k),zsatv(k)/10.,avCALC*1e-3,Mg);
                %         end
                %     end;
                % % co3satv = kspcSed./(avCALC*1e-3);
                % co3satvA  = kspcA./(avCALC*1e-3);
                % co3satvI  = kspcI./(avCALC*1e-3);
                % co3satvP  = kspcP./(avCALC*1e-3);
                %     if(ftys)
                %         co3satvT  = kspcT./(avCALC*1e-3);
                %     end
            elseif(CAvflag == 1)
                for k=1:Ns
                    [kspcA(k),tmp] = ...
                        kspfun(TdvA(k),SdvA(k),zsatv(k)/10.,Ca,Mg);
                    [kspcI(k),tmp] = ...
                        kspfun(TdvI(k),SdvI(k),zsatv(k)/10.,Ca,Mg);
                    [kspcP(k),tmp] = ...
                        kspfun(TdvP(k),SdvP(k),zsatv(k)/10.,Ca,Mg);
                    if(ftys)
                        [kspcT(k),tmp] = ...
                            kspfun(TdvT(k),SdvT(k),zsatv(k)/10.,Ca,Mg);
                    end
                end;
                %     co3satv = kspcSed/Ca;
                co3satvA  = kspcA/Ca;
                co3satvI  = kspcI/Ca;
                co3satvP  = kspcP/Ca;
                if(ftys)
                    co3satvT  = kspcT/Ca;
                end
            elseif(CAvflag == 0)
                for k=1:Ns
                    [kspcA(k),tmp] = ...
                        kspfun(TdvA(k),SdvA(k),zsatv(k)/10.,Ca,Mg);
                    [kspcI(k),tmp] = ...
                        kspfun(TdvI(k),SdvI(k),zsatv(k)/10.,Ca,Mg);
                    [kspcP(k),tmp] = ...
                        kspfun(TdvP(k),SdvP(k),zsatv(k)/10.,Ca,Mg);
                    if(ftys)
                        [kspcT(k),tmp] = ...
                            kspfun(TdvT(k),SdvT(k),zsatv(k)/10.,Ca,Mg);
                    end
                end;
                %     co3satv = kspcSed/Ca;
                co3satvA  = kspcA/Ca;
                co3satvI  = kspcI/Ca;
                co3satvP  = kspcP/Ca;
                if(ftys)
                    co3satvT  = kspcT/Ca;
                end
            end;


            %==== saturation of surface ocean boxes, needed???
            if(1)
                if(CAvflag == 2)
                    lk = 1;
                    for k=kkv
                        CALC(k)=sum(CA(k)*V/Voc);
                        [kspcS(k),kspaS(k)] = ...
                            kspfunCA(TX(k),Sv(k),Pv(k),CA(k)*1e-3,Mg);
                        %         avCALC=sum(CALC.*V/Voc);
                        omegCSv(lk) = co3(k)*(CA(k)*1e-3)/kspcS(k);
                        omegASv = co3(k)*(CA(k)*1e-3)/kspaS(k);
                        lk = lk + 1;
                    end;
                elseif(CAvflag == 1)
                    for k=kkv
                        [kspcS(k),kspaS(k)] = ...
                            kspfunCA(TX(k),Sv(k),Pv(k),Ca,Mg);
                    end;
                    omegCSv = co3(kkv)*Ca./kspcS(kkv);
                    omegASv = co3(kkv)*Ca./kspaS(kkv);

                elseif(CAvflag == 0)
                    for k=kkv
                        [kspcS(k),kspaS(k)] = ...
                            kspfun(TX(k),Sv(k),Pv(k),Ca,Mg);
                    end;
                    omegCSv = co3(kkv)*Ca./kspcS(kkv);
                    omegASv = co3(kkv)*Ca./kspaS(kkv);


                end
            end;
        end;

        % set omega and find supersat indices
        if(CAvflag==0)

            omvA = co3(klid  )./co3satvA; % Atl
            omvI = co3(klid+1)./co3satvI; % Ind
            omvP = co3(klid+2)./co3satvP; % Pac
        else
            omvA = co3(klid  )./co3satvA; % Atl
            omvI = co3(klid+1)./co3satvI; % Ind
            omvP = co3(klid+2)./co3satvP; % Pac

        end

        kdA  = find(omvA >= 1.);
        kdI  = find(omvI >= 1.);
        kdP  = find(omvP >= 1.);


        %====== Carbon-13
        RsA  = f13cvA./fcvA;
        RsI  = f13cvI./fcvI;
        RsP  = f13cvP./fcvP;
        % RsA  = f44cavA./fcvA;
        % RsI  = f44cavI./fcvI;
        % RsP  = f44cavP./fcvP;
        %====== Calcium-44

        if(CAvflag >0)
            RcasA  = f44cavA./fcavA;
            RcasI  = f44cavI./fcavI;
            RcasP  = f44cavP./fcavP;
            RcasAcheck(kt+1,:)=RcasA;
            f44cavAchck (kt+1,:)=f44cavA;
        end;
        fcvAchck(kt+1,:)=fcvA;
        RsAcheck(kt+1,:)=RsA;
        f13cvAchck (kt+1,:)=f13cvA;

        % Ca, Mg/Ca correction for Sigman dissolution
        alpha = 0.0833;
        xm    = Mgm/Cam;
        if(CAvflag == 2)
            for k=1:Ns
                xtv(k)    = Mg /(avCALC*1e-3) ;
                rcak(k)  = ((avCALC*1e-3)/Cam)*(1/(1-alpha*(xm-xtv(k))));
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

        if(CAvflag >0)
            FFca   = (phi1-phi0)/(1-phi1);
            phiAca = (phi0+FFca*fcavA)./(1+FFca*fcavA);
            phiIca = (phi0+FFca*fcavI)./(1+FFca*fcavI);
            phiPca = (phi0+FFca*fcavP)./(1+FFca*fcavP);
        end;

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

        if(CAvflag >0)
            rscavA = FprvAca*m2kg/rhos/(1-phi1);
            rscavI = FprvIca*m2kg/rhos/(1-phi1);
            rscavP = FprvPca*m2kg/rhos/(1-phi1);
            rsrvAca = frrfAca     /rhos/(1-phi0);
            rsrvIca = frrfIca     /rhos/(1-phi0);
            rsrvPca = frrfPca     /rhos/(1-phi0);
            rsvAca  = rscavA+rsrvAca;
            rsvIca  = rscavI+rsrvIca;
            rsvPca  = rscavP+rsrvPca;
        end;

        %====== Carbon-13
        r13scvA = Fpr13vA*m2kg/rhos/(1-phi1);
        r13scvI = Fpr13vI*m2kg/rhos/(1-phi1);
        r13scvP = Fpr13vP*m2kg/rhos/(1-phi1);

        %====== Calcium-44
        if(CAvflag >0)
            r44scvA = Fpr44vA*m2kg/rhos/(1-phi1);
            r44scvI = Fpr44vI*m2kg/rhos/(1-phi1);
            r44scvP = Fpr44vP*m2kg/rhos/(1-phi1);
        end;

        % dissolution
        if(CAvflag == 0);
            dKA   = KS*((co3satvA-co3(klid  ))*rcak).^nc; % mol/m2/y
            dKI   = KS*((co3satvI-co3(klid+1))*rcak).^nc; % mol/m2/y
            dKP   = KS*((co3satvP-co3(klid+2))*rcak).^nc; % mol/m2/y

        elseif(CAvflag == 1)
            dKA   = KS*((co3satvA-co3(klid  ))*rcak).^nc; % mol/m2/y
            dKI   = KS*((co3satvI-co3(klid+1))*rcak).^nc; % mol/m2/y
            dKP   = KS*((co3satvP-co3(klid+2))*rcak).^nc; % mol/m2/y
        elseif(CAvflag == 2)
            dKA   = KS*((co3satvA-co3(klid  )).*rcak(klid)).^nc; % mol/m2/y
            dKI   = KS*((co3satvI-co3(klid+1)).*rcak(klid+1)).^nc; % mol/m2/y
            dKP   = KS*((co3satvP-co3(klid+2)).*rcak(klid+2)).^nc; % mol/m2/y


        end;
        % fc^0.5
        dissA = fcvA.^0.5.*dKA'; % mol/m2/y
        dissI = fcvI.^0.5.*dKI'; % mol/m2/y
        dissP = fcvP.^0.5.*dKP'; % mol/m2/y

        if(CAvflag >0)
            dissAca = fcavA.^0.5.*dKA'; % mol/m2/y
            dissIca = fcavI.^0.5.*dKI'; % mol/m2/y
            dissPca = fcavP.^0.5.*dKP'; % mol/m2/y
        end;

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

            if(CAvflag >0)
                jAca = find(fcavA < 1./Bf^2);
                jIca = find(fcavI < 1./Bf^2);
                jPca = find(fcavP < 1./Bf^2);
            end;

            if(ftys)
                jT = find(fcvT < 1./Bf^2);
                if(CAvflag >0)
                    jTca = find(fcavT < 1./Bf^2);
                end;
            end;
        end;

        dissA(jA) = fcvA(jA).^nf.*dKA(jA)'*Bf; % mol/m2/y
        dissI(jI) = fcvI(jI).^nf.*dKI(jI)'*Bf; % mol/m2/y
        dissP(jP) = fcvP(jP).^nf.*dKP(jP)'*Bf; % mol/m2/y

        if(CAvflag >0)
            dissAca(jAca) = fcavA(jAca).^nf.*dKA(jAca)'*Bf; % mol/m2/y
            dissIca(jIca) = fcavI(jIca).^nf.*dKI(jIca)'*Bf; % mol/m2/y
            dissPca(jPca) = fcavP(jPca).^nf.*dKP(jPca)'*Bf; % mol/m2/y
        end;

        % diss = 0, for omega > 1
        dissA(kdA) = 0.;
        dissI(kdI) = 0.;
        dissP(kdP) = 0.;

        if(CAvflag >0)
            dissAca(kdA) = 0.;
            dissIca(kdI) = 0.;
            dissPca(kdP) = 0.;
        end;

        % diss rate, m/y [(mol/m2/y)*kg/mol / kg*m3 = m/y]
        % pure calcite/(1-phi1) = Delta h1
        rdvA  = dissA'*m2kg/rhos/(1-phi1);     % m/y
        rdvI  = dissI'*m2kg/rhos/(1-phi1);     % m/y
        rdvP  = dissP'*m2kg/rhos/(1-phi1);     % m/y

        if(CAvflag >0)
            rdvAca  = dissAca'*m2kg/rhos/(1-phi1);     % m/y
            rdvIca  = dissIca'*m2kg/rhos/(1-phi1);     % m/y
            rdvPca  = dissPca'*m2kg/rhos/(1-phi1);     % m/y
        end;

        %====== Carbon-13
        r13dvA = RsA.*rdvA';
        r13dvI = RsI.*rdvI';
        r13dvP = RsP.*rdvP';

        %====== Calcium-44
        if(CAvflag >0)
            r44dvA = RcasA.*rdvAca';
            r44dvI = RcasI.*rdvIca';
            r44dvP = RcasP.*rdvPca';
        end;

        % burial rate, m/y
        wvA   = rsvA-rdvA;
        wvI   = rsvI-rdvI;
        wvP   = rsvP-rdvP;

        if(CAvflag >0)
            wvAca   = rsvAca-rdvAca;
            wvIca   = rsvIca-rdvIca;
            wvPca   = rsvPca-rdvPca;
        end;
        % find erosion indices
        lA = find(wvA < 0.);
        lI = find(wvI < 0.);
        lP = find(wvP < 0.);

        if(CAvflag >0)
            lAca = find(wvAca < 0.);
            lIca = find(wvIca < 0.);
            lPca = find(wvPca < 0.);
        end;

        % calcite burial (w>0)
        wcvA  = fcvA.*wvA'.*(1-phiA)/(1-phi1);
        wcvI  = fcvI.*wvI'.*(1-phiI)/(1-phi1);
        wcvP  = fcvP.*wvP'.*(1-phiP)/(1-phi1);
        fbA   = 1*ones(size(rdvA));
        fbI   = 1*ones(size(rdvI));
        fbP   = 1*ones(size(rdvP));

        if(CAvflag >0)
            wcavA  = fcavA.*wvAca'.*(1-phiA)/(1-phi1);
            wcavI  = fcavI.*wvIca'.*(1-phiI)/(1-phi1);
            wcavP  = fcavP.*wvPca'.*(1-phiP)/(1-phi1);
            fbAca   = 1*ones(size(rdvAca));
            fbIca   = 1*ones(size(rdvIca));
            fbPca   = 1*ones(size(rdvPca));
        end;

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

        if(CAvflag >0)
            wcavA(lAca) = -(1-fca0A(lAca)).*wvAca(lAca)'.*(1-phiiAca(lAca))/(1-phi0);
            wcavI(lIca) = -(1-fca0I(lIca)).*wvIca(lIca)'.*(1-phiiIca(lIca))/(1-phi0);
            wcavP(lPca) = -(1-fca0P(lPca)).*wvPca(lPca)'.*(1-phiiPca(lPca))/(1-phi0);
            wcavA(lAca) =   wcavA(lAca) + rsrvAca(lAca)';
            wcavI(lIca) =   wcavI(lIca) + rsrvIca(lIca)';
            wcavP(lPca) =   wcavP(lPca) + rsrvPca(lPca)';
            fbAca (lAca) = 0.;
            fbIca (lIca) = 0.;
            fbPca (lPca) = 0.;
        end;

        %====== Carbon-13
        w13cvA = RsA.*wcvA;
        w13cvI = RsI.*wcvI;
        w13cvP = RsP.*wcvP;

        %====== Calcium-44
        if(CAvflag >0)
            w44cvA = RcasA.*wcavA;
            w44cvI = RcasI.*wcavI;
            w44cvP = RcasP.*wcavP;
        end;

        % dphi/dfc
        dphiA = FF*(1-phi0)./(1+FF*fcvA).^2;
        dphiI = FF*(1-phi0)./(1+FF*fcvI).^2;
        dphiP = FF*(1-phi0)./(1+FF*fcvP).^2;

        if(CAvflag >0)
            dphiAca = FFca*(1-phi0)./(1+FFca*fcavA).^2;
            dphiIca = FFca*(1-phi0)./(1+FFca*fcavI).^2;
            dphiPca = FFca*(1-phi0)./(1+FFca*fcavP).^2;
        end;

        % G's
        GA    = hs*(1-phiA-fcvA.*dphiA)/(1-phi1);
        GI    = hs*(1-phiI-fcvI.*dphiI)/(1-phi1);
        GP    = hs*(1-phiP-fcvP.*dphiP)/(1-phi1);

        if(CAvflag >0)
            GAca    = hs*(1-phiAca-fcavA.*dphiAca)/(1-phi1);
            GIca    = hs*(1-phiIca-fcavI.*dphiIca)/(1-phi1);
            GPca    = hs*(1-phiPca-fcavP.*dphiPca)/(1-phi1);
        end;

        % dissolution (w>0) in mol/y, see above

        % dissolution (w<0) in mol/y
        dissA(lA) = (1*rsvA(lA)-wvA(lA))*(1-phi1)*rhos/m2kg; % mol/m2/y
        dissI(lI) = (1*rsvI(lI)-wvI(lI))*(1-phi1)*rhos/m2kg; % mol/m2/y
        dissP(lP) = (1*rsvP(lP)-wvP(lP))*(1-phi1)*rhos/m2kg; % mol/m2/y

        if(CAvflag >0)
            dissAca(lAca) = (1*rsvAca(lAca)-wvAca(lAca))*(1-phi1)*rhos/m2kg; % mol/m2/y
            dissIca(lIca) = (1*rsvIca(lIca)-wvIca(lIca))*(1-phi1)*rhos/m2kg; % mol/m2/y
            dissPca(lPca) = (1*rsvPca(lPca)-wvPca(lPca))*(1-phi1)*rhos/m2kg; % mol/m2/y
        end;

        %-------------- Tethys --------------------------%
        if(ftys)
            if(CAvflag==0)
                omvT = co3(kliT  )./co3satvT; % Tys
            else
                omvT = co3(kliT  )./co3satvT; % Tys
            end;
            kdT  = find(omvT >= 1.);
            RsT  = f13cvT./fcvT;
            phiT = (phi0+FF*fcvT)./(1+FF*fcvT);
            rscvT = FprvT*m2kg/rhos/(1-phi1);
            rsrvT = frrfT     /rhos/(1-phi0);
            rsvT  = rscvT+rsrvT;
            r13scvT = Fpr13vT*m2kg/rhos/(1-phi1);

            if(CAvflag >0)
                RcasT  = f44cavT./fcavT;
                phiTca = (phi0+FFca*fcavT)./(1+FFca*fcavT);
                rscavT = FprvTca*m2kg/rhos/(1-phi1);
                rsrvTca = frrfTca     /rhos/(1-phi0);
                rsvTca  = rscavT+rsrvTca;
                r44scvT = Fpr44vT*m2kg/rhos/(1-phi1);
            end;
            if(CAvflag == 0)
                dKT   = KS*((co3satvT-co3(kliT))*rcak).^nc; % mol/m2/y
            elseif(CAvflag == 1)
                dKT   = KS*((co3satvT-co3(kliT))*rcak).^nc; % mol/m2/y
            elseif(CAvflag == 2)
                dKT   = KS*((co3satvT-co3(kliT)).*rcak(kliT)).^nc; % mol/m2/y
            end
            dissT = fcvT.^0.5.*dKT';                   % mol/m2/y
            dissT(jT) = fcvT(jT).^nf.*dKT(jT)'*Bf;     % mol/m2/y
            dissT(kdT) = 0.;
            rdvT  = dissT'*m2kg/rhos/(1-phi1);     % m/y

            if(CAvflag >0)
                dissTca = fcavT.^0.5.*dKT';                   % mol/m2/y
                dissTca(jTca) = fcavT(jTca).^nf.*dKT(jTca)'*Bf;     % mol/m2/y
                dissTca(kdT) = 0.;


                rdvTca  = dissTca'*m2kg/rhos/(1-phi1);     % m/y
                r44dvT = RcasT.*rdvTca';
                wvTca   = rsvTca-rdvTca;
                lTca = find(wvTca < 0.);
                wcavT  = fcavT.*wvTca'.*(1-phiTca)/(1-phi1);
                fbTca   = 1*ones(size(rdvTca));
                wcavT(lTca) = -(1-fca0T(lTca)).*wvTca(lTca)'.*(1-phiiTca(lTca))/(1-phi0);
                w44cvT = RcasT.*wcavT;
                fbTca (lTca) = 0.;
            end;

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

            if(CAvflag >0)
                dphiTca = FFca*(1-phi0)./(1+FFca*fcavT).^2;
                GTca    = hs*(1-phiTca-fcavT.*dphiTca)/(1-phi1);
                dissTca(lTca) = (1*rsvTca(lTca)-wvT(lTca))*(1-phi1)*rhos/m2kg; % mol/m2/y



                dissCav(:,:) = [dissAca'.*asvA*A(01);  ...
                    dissIca'.*asvI*A(02);  ...
                    dissPca'.*asvP*A(03);  ...
                    dissTca'.*asvT*A(11)]; % -> mol C/y



                dissC44v(:,:) = ...
                    [RcasA'.*dissAca'.*asvA*A(01);  ...
                    RcasI'.*dissIca'.*asvI*A(02);  ...
                    RcasP'.*dissPca'.*asvP*A(03);  ...
                    RcasT'.*dissTca'.*asvT*A(11)];
            end;
            %      dissC44v(:,:) = ...
            %         [(-1.3/1e3+1)*(RcasA').*dissA'.*asvA*A(01);  ...
            %          (-1.3/1e3+1)*(RcasI').*dissI'.*asvI*A(02);  ...
            %          (-1.3/1e3+1)*(RcasP').*dissP'.*asvP*A(03);  ...
            %          (-1.3/1e3+1)*(RcasT').*dissT'.*asvT*A(11)];
        else
            % dissolution in mol/y
            dissCv(:,:) = [dissA'.*asvA*A(01);  ...
                dissI'.*asvI*A(02);  ...
                dissP'.*asvP*A(03)]; % -> mol C/y


            %====== Carbon-13
            dissC13v(:,:) = [RsA'.*dissA'.*asvA*A(01);  ...
                RsI'.*dissI'.*asvI*A(02);  ...
                RsP'.*dissP'.*asvP*A(03)]; % -> mol C/y
            if(CAvflag>0)
                dissCav(:,:) = [dissAca'.*asvA*A(01);  ...
                    dissIca'.*asvI*A(02);  ...
                    dissPca'.*asvP*A(03)];
                dissC44v(:,:) = [RcasA'.*dissA'.*asvA*A(01);  ...
                    RcasI'.*dissI'.*asvI*A(02);  ...
                    RcasP'.*dissP'.*asvP*A(03)]; % -> mol C/y
            end
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

            if(CAvflag >0)
                rbvAca = wvAca;
                rbvIca = wvIca;
                rbvPca = wvPca;
                dissvAca = dissAca'.*asvA*A(01)*m2kg; % kg/y
                dissvIca = dissIca'.*asvI*A(02)*m2kg; % kg/y
                dissvPca = dissPca'.*asvP*A(03)*m2kg; % kg/y
            end;

            dissv13A = RsA'.*dissvA;
            dissv13I = RsI'.*dissvI;
            dissv13P = RsP'.*dissvP;

            if(CAvflag >0)
                dissv44A = RcasA'.*dissvAca;
                dissv44I = RcasI'.*dissvIca;
                dissv44P = RcasP'.*dissvPca;
            end;

            if(ftys)
                rbvT = wvT;
                dissvT = dissT'.*asvT*A(11)*m2kg; % kg/y
                dissv13T = RsT'.*dissvT;

                if(CAvflag >0)
                    rbvTca = wvTca;
                    dissvTca = dissTca'.*asvT*A(11)*m2kg; % kg/y
                    dissv44T = RcasT'.*dissvTca;
                end;
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
    if(LTflag)
        tstgc=1e3;
        tg=t-kgc*tstgc;
        %if(pflag==0)
        if(tg>=tstgc)
            geoflag=1;
            kgc=kgc+DUM;
            %kgct=t/1e6;
        else
            geoflag=1;
        end
    else
        geoflag = 0;
    end

    if(geoflag)
        if(LTflag)
            
%             [fbgv,fbcv,fwgv,fmcv,fmgv,fwcv,fGGi,fgkc,frkc,fekc,fdkc,rco2,pco2gca]=gcfun12((interp1(1:length(dbc),dbc,(tgc-(kgct)))),tgc-kgct);
            FVC = fmcv*1e12/Aoc;
            Focb0 = fbgv*1e12/Aoc;
            Focw0 = (fwgv+fmgv)*1e12/Aoc;
            FiNgc = fwcv*1e12/Aoc;
            fGG = fGGi;
            frkci = frkc;
            fekci = fekc;
            fdkci = fdkc;
            flakci=flakc;
            aa=exp(-ACT*ws*(tgc-1)/570)*exp(ACT*fGG);
            bb=1-RUN(tgc)*ws*(tgc-1)/570;
            ccc=1-0.087*ws*(tgc-1)/570;
            fb=aa*(pco2a/pco20)^(ACT*gamma(tgc))*(bb+(gamma(tgc)*RUN(tgc))*log(pco2a/pco20)+RUN(tgc)*fGG)^0.65*(2*(pco2a/pco20)/(1+(pco2a/pco20)))^FERT;
            FSigc   = FSi0*(fb+0.0000)*frkci*fekci*(fdkci^0.65);
            
            fbb=(ccc+gamma(tgc)*0.087*log(pco2a/pco20)+0.087*fGG)*(2*(pco2a/pco20)/(1+(pco2a/pco20)))^FERT;
            FiNgc = fbb*flakci*fdkci*fekci*FiN0;
%             FVC=((interp1(1:length(fmcv),fmcv,tgc-kgct))*1.e12)/Aoc;
%             Focb0   =   ((interp1(1:length(fbgv),fbgv,tgc-kgct))*1.e12)/Aoc;        % mol C    /m2/y
%             Focw0  =   (((interp1(1:length(fwgv),fwgv,tgc-kgct))+(interp1(1:length(fmgv),fmgv,tgc-kgct)))*1.e12)/Aoc;
%             FiNgc   =   ((interp1(1:length(fwcv),fwcv,tgc-kgct))*1.e12)/Aoc;
%             aa=exp(-ACT*ws*(tgc-kgct)/570)*exp(ACT*(interp1(1:length(fGG),fGG,tgc-kgct)));
%             bb=1-RUN(tgc)*ws*(tgc-kgct)/570;
%             ccc=1-0.087*ws*(tgc-kgct)/570;
%             fb=aa*(pco2a/pco20)^(ACT*gamma(tgc))*(bb+(gamma(tgc)*RUN(tgc))*log(pco2a/pco20)+RUN(tgc)*(interp1(1:length(fGG),fGG,tgc-kgct)))^0.65*(2*(pco2a/pco20)/(1+(pco2a/pco20)))^FERT;
%             %       fbb=((interp1(1:length(ccv),ccv,tgc-kgct))+gamma(tgc)*0.087*log(pco2a/pco20)+0.087*(interp1(1:length(fGG),fGG,tgc-kgct)))*(2*(pco2a/pco20)/(1+(pco2a/pco20)))^FERT;
%             FSigc   = 5*1e12/Aoc*(fb+0.0000)*(interp1(1:length(frkc),frkc,(tgc-(kgct))))*(interp1(1:length(fekc),fekc,(tgc-(kgct))))*((interp1(1:length(fdkc),fdkc,(tgc-(kgct))))^0.65);
        end

    end
    % volcanic degassing
    Fvc   = FVC;

    if(LTflag)
        FSi   = FSigc;
        Fin   = FiNgc;
    else
        % CaSiO3 weathering
        % if(t<1e4)
        % FSi   = 2*FVC*(pco2a/pCSi)^nSi; % CaSiO3
        % FSichck(kt) = FSi;
        % else
        FSi   = FVC*(pco2a/pCSi)^nSi; % CaSiO3
        % if(t<1e5&&t>100)
        %     FSi   = FVC*3;
        % else
        %     FSi   = FVC*(pco2a/pCSi)^nSi;
        % end

        % CaCO3  weathering
        % if(t<1e4)
        % Fin   = 2*FiN*(pco2a/pCSi)^nCC; % CaC O3
        % Finchck (kt) = Fin;
        % else

        Fin   = FiN*(pco2a/pCSi)^nCC;
        % if(t<1e5&&t>100)
        %     Fin   = FiN*2.4; % CaC O3
        % else
        %     Fin   = FiN*(pco2a/pCSi)^nCC; % CaC O3
        % end
        % if(t<1e5)
        %     CafactCarb   = 1.8714;
        %     CafactSi   = 1.3690;
        % else
        %     CafactCarb   = 1;
        %     CafactSi   = 1;
        % end
    end

    % phoshpate weathering and organic carbon burial
    if(kt > chck1st) % first call with new kt
        first = 1;
    else
        first = 0;
    end

    %  if(first)
    %      display('Called first time around to get initial concentration');
    %      po40 = p(3); % this will need to be area weighted average of surf boxes
    %      chck1st = kt;
    %  end

    if(Pfeed)
        %     if(Pscenario==0)
        %     O0 = 0.123389861337360; % mean dox of all deep boxes at time 0
        %     end
        %     if(Pscenario==1)
        %      O0 = 0.118807766915575;
        %     end
        %     O0 = 0.162748304894886; % dox of deep atlantic at time 0
        %% Lenton and Watson 2000
        %     Focb = Focb0*(p(3)/3.428657287157348e-04)^2;
        %     Focw = Focw0*(4/12*FSi/FVC+8/12*Fin/FiN);
        %
        %     anoxf=1-oxicf0*(3.428657287157348e-04/p(3));
        %     Fpw = Fpw0*(4/12*FSi/FVC+8/12*Fin/FiN);
        %     Fopb = Focb/250;
        %     Ffep = Ffep0/oxicf0*(1-anoxf);
        %     Fcap = Fcap0*(p(3)/3.428657287157348e-04)^2;
        %% Tsandev, Slomp, Van Cappellen

        %% initial [O2] deep Pac = 0.1084 mol/m3

        %sum(EPLv)=3.953740401592986e+014
        %sum(PPLv)=3.041338770456144e+012

        if(ftys)
            doxMeanI =  (dox(7).*V(7)+dox(8).*V(8)+dox(9)*V(9)+dox(13)*V(13))./((V(7)+V(8)+V(9)+V(13)));
            %doxMeanDeep = (dox(4).*V(4)+dox(5).*V(5)+dox(6)*V(6)+ dox(7).*V(7)+dox(8).*V(8)+dox(9)*V(9)+dox(12)*V(12)+dox(13)*V(13))./((V(4)+V(5)+V(6)+V(7)+V(8)+V(9)+V(12)+V(13)));    
%             doxMeanI =  (dox(4).*V(4)+dox(5).*V(5)+dox(6)*V(6)+dox(12)*V(12))./((V(4)+V(5)+V(6)+V(12)));
        else
            doxMeanI = (dox(7).*V(7)+dox(8).*V(8)+dox(9)*V(9))./((V(7)+V(8)+V(9)));
%             doxMeanI = (dox(4).*V(4)+dox(5).*V(5)+dox(6)*V(6))./((V(4)+V(5)+V(6)));
        end

       
%           Ffep = Ffep0*doxMeanI/O0;
        %     Fopb =4.991716e-005*sum(EPLv)/Aoc*(0.25+0.75*doxMeanI/O0);%-3.953740401592986e12/Aoc+0.01*sum(EPLv)/Aoc;
        %     Fcap = 0.005025125628177*(sum(PPLv)/Aoc-Fopb)*(0.1+0.9*doxMeanI/O0);
        %     Fpw = Fpw0*(4/12*FSi/FVC+8/12*Fin/FiN);
        %     Focb = 0;%Fopb*(CPox*CPanox)/(doxMeanI/O0*CPanox+(1-doxMeanI/O0)*CPox)-3.95374e+012/Aoc;
        %     Focw = Focw0*1;%(4/12*FSi/FVC+8/12*Fin/FiN);
        if(~LTflag)
            if(Pscenario==0)
                po4bf0=0.005;
                ocbf0 = 0.01;
            end
            if(Pscenario==1)
                po4bf0=0.01;
                ocbf0 = 0.01;
            end
            if(Pscenario==2)
                po4bf0=0.0025;
                ocbf0 = 0.01;
            end
            if(Pscenario==3)
                po4bf0=0.005;
                ocbf0 = 0.005;
            end
            if(Pscenario==4)
                po4bf0=0.01;
                ocbf0 = 0.005;
            end
            if(Pscenario==5)
                po4bf0=0.0025;
                ocbf0 = 0.005;

            end
            if(Pscenario==6)           
                po4bf0=0.005;
                ocbf0 = 0.02;
            end
            if(Pscenario==8)
                po4bf0=0.0025;
                ocbf0 = 0.02;
            end        
            if(Pscenario==7)
                po4bf0=0.01;
                ocbf0 = 0.02; %*0.817353763464955
            end 
            if(Pscenario==9)
                po4bf0=0.005;
                ocbf0 = 0.04;
            end
            if(Pscenario==10)
                po4bf0=0.01;
                ocbf0 = 0.04;
            end
            if(Pscenario==11)
                po4bf0=0.0025;
                ocbf0 = 0.04;
            end
            if(Pscenario==12)
                po4bf0=0.02;
                ocbf0 = 0.02;
            end
            if(Pscenario==13)
                po4bf0=0.01;
                ocbf0 = 0.02;
            end
        end
        if(Floegel)
            if(Floegel == 1)
                po4bf0=0.02;
                ocbf0 = 0.01;
            elseif(Floegel == 2)
                po4bf0=0.02;
                ocbf0 = 0.02;
            elseif(Floegel == 3)
                po4bf0=0.02;
                ocbf0 = 0.005;
            elseif(Floegel == 4)
                po4bf0=0.04;
                ocbf0 = 0.01;
            elseif(Floegel == 5)
                po4bf0=0.01;
                ocbf0 = 0.01;
            end
            po4bf=po4bf0;
            ocbf = po4bf*ocbf0/po4bf0; %Organic C burial factor

            %%% Floegel and Wallmann
            Yf=123; % +-24
            Af=-112; % +-24 
            rf=32;   % +-19
            rREG = Yf+Af*exp(-doxMeanI*1e3/rf);
            if(counter==0)
                rREG0 = Yf+Af*exp(-doxMeanI*1e3/rf);
                totCexp0 = (sum(EPLv)+EPH)/Aoc;
                Focw0 = (sum(EPLv)+EPH)*ocbf0/Aoc;
                Fopb0 = (sum(PPLv)+PPH)*po4bf0/Aoc;
                Ffep0 = 0;
                Fcap0 = 0;
                Fpw0 = Fopb0+Ffep0+Fcap0;
                capk= (Fcap0)/(sum(PPLv)/Aoc-sum(PPLv)/Aoc*po4bf); % rate const. for Ca sorbed P burial
                REDPC0 = REDPC;
            end
            totCexp=(sum(EPLv)+EPH)/Aoc;
            Focb = 0;
            Focw = Focw0;%*(2/12*FSi/FVC+5/12*Fin/FiN+5/12*1); %Focw0
            Fpw = Fpw0*((FSi+Fin)/(FiN+FVC)); %Floegel, wallmann 2011
            Fopb = 0;%Fopb0;
            Ffep = Ffep0*doxMeanI/O0; %Ffep0

            Fcap = capk*(sum(PPLv)/Aoc-sum(PPLv)/Aoc*po4bf*(0.1+0.9*doxMeanI/O0)); %Fcap0
            Ffepv(kt) = Ffep*Aoc;
            Fcapv(kt) = Fcap*Aoc;
            counter = counter+1;
            ocbf =ocbf0;%*(totCexp/totCexp0)^1.05;
            oI = 1-eI-ocbf;
           
            oIp = 1-eI-ocbf*po4bf0/ocbf0*(rREG/rREG0)*((REDPC0)/REDPC);
            rREGv(kt) = rREG;
        else
           
            Qfac=(1/Q10)^((TCv0(9)-2.0)/10);
%             Qfac=(Q10)^((TCv0(9)-2.0)/10);
            
            % Can also make Q10 dependent on temperature
            %%%% 0.78.*(2.*2.^(-(TCvt(:,9)-2)./14)).^((TCvt(:,9)-2.0)/10)
%                Qfac=1;
%             if(tgc<53)
%                Qfac = (1/interp1([1 59],[1.5 3],53,'pchip'))^((TCv0(9)-2.0)/10);
%             end
            % Q10 effect *1/3^((TCv(9)-TCv0(9))/10)
            po4bf=po4bf0*(0.25+0.75*doxMeanI/O0)*ocbf/ocbf0; %*ocbf/ocbf0 fraction of PO4 rain that gets burried
            % org PO4 burial needs to be
            % 260 times smaller than org C
            % burial. It is already 130
            % times smaller due to
            % redfield in the surface ocean.
            % And additional two
            % times smaller due two the
            % initial factor of 0.005
            % compared to 0.01 for org C
            
%             ocbf = po4bf*ocbf0/po4bf0; %Organic C burial factor
%             ocbf = ocbft;
            CPox=1/(REDPC)*(ocbf0/po4bf0); %260
            CPanox=1100;  %4000
            if(counter==0)
%                 Focw0 = (sum(EPLv)+EPH)*ocbf0/Aoc;
%                 fbgv
%                 Focb0*Aoc
%                 (sum(EPLv)+EPH)*0.016346752281632

%             (sum(EPLv)+EPH)*0.02*0.817353763464955
%             return
%                 Focb0*Aoc/(sum(EPLv)+EPH)
%                 return
                if(~LTflag)
                    Fopb0 = (sum(PPLv)+PPH)*po4bf0/Aoc;
                    Ffep0 = Fopb0;
                    Fcap0 = Ffep0*2;
                    Fpw0 = Fopb0+Ffep0+Fcap0;
                    capk= (Fcap0)/(sum(PPLv)/Aoc-sum(PPLv)/Aoc*po4bf); % rate const. for Ca sorbed P burial
                end
            end
            if(~LTflag)
                Focb = 0;
                Focw = Focw0;%*(2/12*FSi/FVC+5/12*Fin/FiN+5/12*1); %Focw0
                %     Fpw = Fpw0*(2/12*FSi/FVC+5/12*Fin/FiN+5/12*1); %Fpw0
                Fpw = Fpw0*((FSi+Fin)/(FiN0+FSi0)); %Floegel, wallmann 2011
                Fopb = 0;%Fopb0;
                Ffep = Ffep0*doxMeanI/O0; %Ffep0
                Fcap = capk*(sum(PPLv)/Aoc-sum(PPLv)/Aoc*po4bf*(0.1+0.9*doxMeanI/O0)); %Fcap0
            else
                Focb = 0; % zero because it is calculated
                Focw = Focw0;%*(2/12*FSi/FVC+5/12*Fin/FiN+5/12*1); %Focw0
                %     Fpw = Fpw0*(2/12*FSi/FVC+5/12*Fin/FiN+5/12*1); %Fpw0
                Fpw = Fpw0*((FSi*beta+Fin*(beta-1))/(FSi0*beta+FiN0*(beta-1))); %Floegel, wallmann 2011
                Fopb = 0;%zero because it is calculated
                Ffep = Ffep0*doxMeanI/O0; %Ffep0
                Fcap = capk0*((sum(PPLv)+PPH)/Aoc-((sum(PPLv)+PPH)/Aoc)*po4bf*(0.1+0.9*doxMeanI/O0)); %Fcap0 

            end
            Ffepv(kt) = Ffep*Aoc;
            Fcapv(kt) = Fcap*Aoc;
            Fpwv(kt) = Fpw*Aoc;
            counter = counter+1;
            
%             eIi = 0.78;%*Qfac;
            % the condition below should never be met as long as Q10 is
            % less than 1.20 because temp increase is no more than 11
            % degrees across the Cenozoic
%             if(eI>=0.99 || (eI+((1-eI)/0.22)*ocbf*((*CPanox)/(doxMeanI/O0*CPanox+(1-doxMeanI/O0)*))/)>=1.0)
%                 eI=0.99;
%                 ocbf = 0.009;
%             end
%             if(eI>=0.99 || ((1-eI)/0.22)*(eI+po4bf)>=1.0)
%                 eI=0.99;
%                 po4bf = 0.009;
%             end
%             oI = 1-eI-((1-eI)/0.22)*ocbf*((*CPanox)/(doxMeanI/O0*CPanox+(1-doxMeanI/O0)*))/;
%             
%             oIp = 1-eI-((1-eI)/0.22)*po4bf;
%             oIo = oI;
            fbd = 1.0; %fraction of total fraction buried in the deep
            fbi = 1-fbd; %fraction of total fraction buried in intermediate
            %total C burial fraction
            ocbft = ocbf*((CPox*CPanox)/(doxMeanI/O0*CPanox+(1-doxMeanI/O0)*CPox))/CPox*po4bf/po4bf0;%*po4bf/po4bf0;
            % C and P burial fraction in deep
            ocbfd = ocbft*fbd;  
            po4bfd = po4bf*fbd;
            % C and P burial fraction in intermediate
            ocbfi=ocbft*fbi;
            po4bfi = po4bf*fbi;
            eI = eIi-ocbfi;
            oI = 1-eIi-ocbfd; 
            %oIpold = 1-eIi-po4bf;
            eIp=eIi-po4bfi;
            oIp = 1-eIi-po4bfd;
            
            if(oI<0)
                oI=0;
            end
            if(oIp<0)
                oIp=0;
            end
            
            oIo = oI;
            
%             orgCb = ocbf.*(sum(EPLv)+EPH)
%             oI
%             eI
%             ocbf
%             FSi*Aoc
%             FSi0*Aoc
%             ((FSi+Fin)/(FiN0+FSi0))
%             (sum(EPLv)+EPH)
%             Focb0*Aoc/(sum(EPLv)+EPH)
%             return
        end
    else

        Focb = Focb0;
        Focw = Focw0;

        anoxf=1-oxicf0;
        Fpw = Fpw0;
        Fopb = Fopb0;
        Ffep = Ffep0/oxicf0*(1-anoxf);
        Fcap = Fcap0;
        O0 = dox(9);
        % fraction EPL, remineralized in I boxes
        oI = 1-eI;
        oIp = 1-eI;
    end

    Focbv(kt) = Focb*Aoc;
    FSichck(kt) = FSi*Aoc;  %mol C/y
    Finchck (kt) = Fin*Aoc; %mol C/y

    % Used for testing when I compared simple Payne and Loscar Ca model
    CafactCarb   = 1;
    CafactSi   = 1;

    %end
    % All C13 Fluxes
    FSi13 = FSi*Rin;
    Fin13 = Fin*Rin;
    Fvc13 = Fvc*Rvc;
    Fkg13  = Rkg*Focw;


    % All Ca44 Fluxes
    FSi44 = FSi*RinCA*CafactSi;
    Fin44 = Fin*RinCA*CafactCarb;
    % Fvc44 = Fvc*Rvc;

    % if(t<50e3)
    % %     Focb0 = flin(0,50e3,1.0*05.e12/Aoc,0.0*05.e12/Aoc,t) ;
    %     Focb0 = 1*02.5e12/Aoc;
    % else
    %     Focb0 = 1*05.e12/Aoc;
    % end
    %

    if(~fsed)
        Fin   = 0*FiN;
        Fin13 = 0*FiN13;
        Fin44 = 0*Fin*RinCA;
        FSi   = 0.;
        FSi13 = 0.;
        FSi44 = 0.;
        Fvc   = 0.;
        Fvc13 = 0.;
        Focb   = 0.;
        Focw  = 0.;
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
            wcvtA(it,:) = wcvA;
            wcvtI(it,:) = wcvI;
            wcvtP(it,:) = wcvP;
            if(ftys)
                rsedvtT(it,:) = rbvT;
                wcvtT(it,:) = wcvT;
            end;

            if(CAvflag>0)
                rburvtAca(it,:) = rbvAca;
                rburvtIca(it,:) = rbvIca;
                rburvtPca(it,:) = rbvPca;
            end;
            if(ftys)
                rsedvtT(it,:) = rbvT;
                if(CAvflag>0)
                    rsedvtTca(it,:) = rbvTca;
                end;
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

            if(CAvflag>0)
                dissvtAca(it,:) = dissvAca;
                dissvtIca(it,:) = dissvIca;
                dissvtPca(it,:) = dissvPca;

                dissv44tA(it,:) = dissv44A;
                dissv44tI(it,:) = dissv44I;
                dissv44tP(it,:) = dissv44P;

                phivtAca (it,:) = phiAca';
                phivtIca (it,:) = phiIca';
                phivtPca (it,:) = phiPca';
            end;

            FprtA  (it  ) = FprA;
            FprtI  (it  ) = FprI;
            FprtP  (it  ) = FprP;
            Fpr13tA(it  ) = Fpr13A;
            Fpr13tI(it  ) = Fpr13I;
            Fpr13tP(it  ) = Fpr13P;

            if(CAvflag>0)
                FprtAca  (it  ) = FprAca;
                FprtIca  (it  ) = FprIca;
                FprtPca  (it  ) = FprPca;
                Fpr44tA(it  ) = Fpr44A;
                Fpr44tI(it  ) = Fpr44I;
                Fpr44tP(it  ) = Fpr44P;
            end;
            if(ftys)
                dissvtT(it,:) = dissvT;
                dissv13tT(it,:) = dissv13T;
                phivtT (it,:) = phiT';
                FprtT  (it  ) = FprT;
                Fpr13tT(it  ) = Fpr13T;

                if(CAvflag>0)
                    dissvtTca(it,:) = dissvTca;
                    dissv44tT(it,:) = dissv44T;
                    phivtTca (it,:) = phiTca';
                    FprtTca  (it  ) = FprTca;
                    Fpr44tT(it  ) = Fpr44T;
                end;

            end;
            Fint   (it  ) = Fin;
            Fkgst  (it  ) = Focw;
            Fin13t (it  ) = Fin13;
            Fin44t (it  ) = Fin44;
            FSit   (it  ) = FSi;
            FSi13t (it  ) = FSi13;
            FSi44t (it  ) = FSi44;
        end; % sediments
        %%%    it = it+1; % increased below (blflag)
    end; % dYflag

    % set all derivs to zero
    %yp(1:3*Nb+1) = 0;
    %yp(1:3*Nb+1+Ns) = 0;
    %yp(1:3*Nb+1+3*Ns) = 0;
    %yp(1:4*Nb+2+3*Ns) = 0;
    if(CAvflag == 0)
        if(fsed)
            if(ftys)
                yp(1:4*Nb+2+8*Ns) = 0;
            else
                yp(1:4*Nb+2+6*Ns) = 0;
            end;
            EXLv   = EPLv;
            EXLvCC = EPLvCC;
            EXH    = EPH;
            EXHCC  = EPHCC;
            x      = 0;
        else % fsed
            yp(1:4*Nb+2     ) = 0;
            EXLv   = ECLv;
            EXLvCC = ECLvCC;
            EXH    = ECH;
            EXHCC  = ECHCC;
            x      = 1;
        end;
    else
        if(fsed)
            if(ftys)
                yp(1:7*Nb+2+16*Ns) = 0;
            else
                yp(1:7*Nb+2+12*Ns) = 0;
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
    end;


    % omegCSv(1)=2.1031;
    % omegCSv(2)=2.1141;
    % omegCSv(3)=2.1035;
    % omegCSv(4)=1.1923;
    % omegCSv(5)=1.6344;
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
    % sum(EPLv)
    % sum(ECLv)
    % sum(EXLv)
    % sum(EALv)
    % return
    cp = THmfun(c,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
    % air-sea
    for k=kkv
        cp(k ) = cp(k ) + ...                           %
            kasv(k)*(pco2a-pco2(k))/V(k);        %
    end;
    % bio pump Corg
    for k=1:3
        cp(k  ) = cp(k )  -    ECLv(k)/V(k  );    % L
        cp(k+3) = cp(k+3) + eI*EXLv(k)/V(k+3);    % I #!  EC or EP
        cp(k+6) = cp(k+6) + oI*EXLv(k)/V(k+6) ... % D #!(tot or Corg)
            + nu*EALv(k)/V(k+6)/2;  % D ClmnDiss
    end;
    cp(10)  = cp(10)  -    ECH/V10;           % H
    for k=7:9                                 % DA,DI,DP
        cp(k )  = cp(k )  +    (eI+oI)*EXH/V(k)*gp(k)...  % #!  EC or EP
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
            if(omegCSv(k)<1)
                kcarb=kcarb0;
            end
            cp(k  ) = cp(k  ) + 2*Fin*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                + 2*FSi*Aoc/V(k)/nOC         ... % #! Si
                + 1*Faom*Aoc/V(k)/nOC...
                - 1*Fmeth*Aoc/V(k)/nOC...
                - 1*Focb*Aoc/V(k)/nOC         ... % #! krgn
                - 0*kcarb*(omegCSv(k)-1)^2*Aoc/V(k)/nOC;
            cp(k  ) = cp(k  ) + sum(dissCv(k,1   :n1))/V(k  ); % diss L
            cp(k+3) = cp(k+3) + sum(dissCv(k,n1+1:n2))/V(k+3); % diss I
            cp(k+6) = cp(k+6) + sum(dissCv(k,n2+1:Ns))/V(k+6); % diss D
        end;
        if(ftys)
            k=11;
            if(omegCSv(5)<1)
                kcarb=kcarb0;
            end
            cp(k  ) = cp(k  ) + 2*Fin*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                + 2*FSi*Aoc/V(k)/nOC         ... % #! Si
                + 1*Faom*Aoc/V(k)/nOC...
                - 1*Fmeth*Aoc/V(k)/nOC...
                - 1*Focb*Aoc/V(k)/nOC         ... % #! krgn
                - 0*kcarb*(omegCSv(5)-1)^2*Aoc/V(k)/nOC;
            cp(k  ) = cp(k  ) + sum(dissCv(4,1   :n1))/V(k  ); % diss L
            cp(k+1) = cp(k+1) + sum(dissCv(4,n1+1:n2))/V(k+1); % diss I
            cp(k+2) = cp(k+2) + sum(dissCv(4,n2+1:Ns))/V(k+2); % diss D
        end;
    end;
    EPLvv(kt)  = sum(EPLv);
    oIv (kt)   = oI;
    oIpv (kt)  = oIp;
    eIv(kt)    = eI;
    eIpv(kt)    = eIp;
    PPLvv (kt) = sum(PPLv);
    % if(t<10e6)
    add=0;
    % else
    %     add=0;
    % end
    aalk=1;
    if(tgc>52)
       aalk =0; 
    end
    %=================== TA   ======================%
    % TH & mixing
    ap = THmfun(a,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
    % bio pump CaCO3, Aorg
    for k=1:3
        ap(k  ) = ap(k ) +add*2 -      EALv(k)/V(k  )+ENLv(k)/V(k  ) ;
        ap(k+3) = ap(k+3)+eIp*(x*EALv(k)/V(k+3)-ENLv(k)/V(k+3))...
            +fbi*Pfeed*Ffep*REDNC/REDPC*Aoc/V(k+3)/nOC...% Alkalinity source from iron-sorbed burial
            +fbi*Pfeed*Fcap*REDNC/REDPC*Aoc/V(k+3)/nOC; % Alkalinity source from fluorapatite burial;  %;   % 1#!
        ap(k+6) = ap(k+6)+oIp*(x*EALv(k)/V(k+6)-ENLv(k)/V(k+6))... % 1#!
            +   nu*EALv(k)/V(k+6);                   % D ClmnDiss
    end;
    ap(10)  = ap(10)  - EAH/V10        + ENH/V10;
    for k=7:9                                              % DA,DI,DP
        ap(k )  = ap(k )  +(eIp+oIp)*(x*EAH/V(k)     - ENH/V(k))*gp(k)...% 1#!
            +nu*EAH/V(k)                *gp(k)...
            +fbd*Pfeed*Ffep*REDNC/REDPC*Aoc/V(k)/nOC...% Alkalinity source from iron-sorbed burial
            +fbd*Pfeed*Fcap*REDNC/REDPC*Aoc/V(k)/nOC; % Alkalinity source from fluorapatite burial;  %
    end;
    if(ftys)
        k = 11;
        ap(k  ) = ap(k  )-      EALv(4)/V(k  )+ENLv(4)/V(k  ) ;
        ap(k+1) = ap(k+1)+eIp*(x*EALv(4)/V(k+1)-ENLv(4)/V(k+1))...
            +fbi*Pfeed*Ffep*REDNC/REDPC*Aoc/V(k+1)/nOC...% Alkalinity source from iron-sorbed burial
            +fbi*Pfeed*Fcap*REDNC/REDPC*Aoc/V(k+1)/nOC; % Alkalinity source from fluorapatite burial;;
        ap(k+2) = ap(k+2)+oIp*(x*EALv(4)/V(k+2)-ENLv(4)/V(k+2))... %
            +   nu*EALv(4)/V(k+2)...
            +fbd*Pfeed*Ffep*REDNC/REDPC*Aoc/V(k+2)/nOC...% Alkalinity source from iron-sorbed burial
            +fbd*Pfeed*Fcap*REDNC/REDPC*Aoc/V(k+2)/nOC; % Alkalinity source from fluorapatite burial;
    end;
    % riverine & sediment fluxes
    if(fsed)                                           % #!
        for k=1:3
            if(omegCSv(k)<1)
                kcarb=kcarb0;
            end
            ap(k  ) = ap(k  )+2*Fin*Aoc/V(k)/nOC          ...  % Fin  L
                +2*FSi*Aoc/V(k)/nOC          ...   % #! Si
                + aalk*Faom*Aoc/V(k)/nOC...
                - 1*Fmeth*Aoc/V(k)/nOC...
                -0*kcarb*(omegCSv(k)-1)^2*Aoc/V(k)/nOC...
                -Pfeed*Fpw*REDNC/REDPC*Aoc/V(k)/nOC;%...%Alkalinity sink from phosphate weathering
            %                  +Pfeed*Ffep*REDNC/REDPC*Aoc/V(k)/nOC...% Alkalinity source from iron-sorbed burial
            %                  +Pfeed*Fcap*REDNC/REDPC*Aoc/V(k)/nOC; % Alkalinity source from fluorapatite burial
            ap(k  ) = ap(k  )+2*sum(dissCv(k,1   :n1))/V(k  ); % diss L
            ap(k+3) = ap(k+3)+2*sum(dissCv(k,n1+1:n2))/V(k+3); % diss I
            ap(k+6) = ap(k+6)+2*sum(dissCv(k,n2+1:Ns))/V(k+6); % diss D
        end;
        if(ftys)
            k=11;
            if(omegCSv(5)<1)
                kcarb=kcarb0;
            end
            ap(k  ) = ap(k  )+2*Fin*Aoc/V(k)/nOC          ...  % Fin  L
                +2*FSi*Aoc/V(k)/nOC          ...   % #! Si
                + aalk*Faom*Aoc/V(k)/nOC...
                - 1*Fmeth*Aoc/V(k)/nOC...
                -0*kcarb*(omegCSv(5)-1)^2*Aoc/V(k)/nOC...
                -Pfeed*Fpw*REDNC/REDPC*Aoc/V(k)/nOC;%...%Alkalinity sink from phosphate weathering
            %                  +Pfeed*Ffep*REDNC/REDPC*Aoc/V(k)/nOC...% Alkalinity source from iron-sorbed burial
            %                  +Pfeed*Fcap*REDNC/REDPC*Aoc/V(k)/nOC; % Alkalinity source from iron-sorbed burial
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
        pp(k+3) = pp(k+3) + eIp*PPLv(k)/V(k+3)...
            - fbi*Ffep*Aoc/V(k+3)/nOC...
            - fbi*Fcap*Aoc/V(k+3)/nOC...
            - fbi*Fopb*Aoc/V(k+3)/nOC;             % DA,DI,DP;
        pp(k+6) = pp(k+6) + oIp*PPLv(k)/V(k+6);
    end;
    pp(10)  = pp(10)  - PPH/V10;
    for k=7:9
        pp(k )  = pp(k )  + (eIp+oIp)*PPH/V(k)*gp(k)...
            - fbd*Ffep*Aoc/V(k)/nOC...
            - fbd*Fcap*Aoc/V(k)/nOC...
            - fbd*Fopb*Aoc/V(k)/nOC;             % DA,DI,DP
    end;
    if(ftys)
        k = 11;
        pp(k  ) = pp(k  ) -    PPLv(4)/V(k  );
        pp(k+1) = pp(k+1) + eIp*PPLv(4)/V(k+1)...
            - fbi*Ffep*Aoc/V(k+1)/nOC...
            - fbi*Fcap*Aoc/V(k+1)/nOC...
            - fbi*Fopb*Aoc/V(k+1)/nOC;
        pp(k+2) = pp(k+2) + oIp*PPLv(4)/V(k+2)...
            - fbd*Ffep*Aoc/V(k+2)/nOC...
            - fbd*Fcap*Aoc/V(k+2)/nOC...
            - fbd*Fopb*Aoc/V(k+2)/nOC;
    end;
    % riverine PO4 fluxes and PO4 burial
    for k=1:3
        pp(k) = pp(k) + Fpw*Aoc/V(k)/nOC...
            ;%...
        %                  - Ffep*Aoc/V(k)/nOC...
        %                  - Fcap*Aoc/V(k)/nOC;
    end
    if(ftys)
        k=11;
        pp(k) = pp (k) + Fpw*Aoc/V(k)/nOC...
            ;%...
        %                  - Ffep*Aoc/V(k)/nOC...
        %                  - Fcap*Aoc/V(k)/nOC;
    end

    % EALvv ECALvv ECAHv EAHv dissCavv
    %=================== CALCIUM  ======================%
    % TH & mixing
%     if(CAvflag == 1 || CAvflag == 2)
%         CAp = THmfun(CA,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
%         % Checking balances mol/yr
%         % sum(ECALv)-nu*sum(EALv)/2 - Total Available for deposition = 3.8856e+013
%         % sum(dissCav(:)) - total diss/bur = 1.7278e+013
%         % (Fin+FSi)*Aoc - Total riverine fluxes = 2.1578e+013
%         % Total available for deposition = total diss + total riverine; but should
%         % it be? ot tot riv = total avail for dep - total diss?
%         % bio pump CaCO3, Aorg
%         for k=1:3
%             CAp(k  ) = CAp(k ) -      ECALv(k)/V(k  );
%             CAp(k+3) = CAp(k+3)+eI*(x*ECALv(k)/V(k+3)); % 1#!
%             CAp(k+6) = CAp(k+6)+oI*(x*ECALv(k)/V(k+6))...% 1#!
%                 +   nu*EALv(k)/V(k+6)/2;%+ 0.01*ENLv(k)/V(k+6)/2;            % D ClmnDiss
%         end;
%         CAp(10)  = CAp(10)  - ECAH/V10;
%         for k=7:9                                              % DA,DI,DP
%             CAp(k )  = CAp(k )  +(eI+oI)*(x*ECAH/V(k))...%   1#!
%                 +nu*EAH/V(k)   *gp(k) / 2;  %
%         end;
%         if(ftys)
%             k = 11;
%             CAp(k  ) = CAp(k  )-      ECALv(4)/V(k  );
%             CAp(k+1) = CAp(k+1)+eI*(x*ECALv(4)/V(k+1));
%             CAp(k+2) = CAp(k+2)+oI*(x*ECALv(4)/V(k+2))...
%                 +   nu*EALv(4)/V(k+2)/2;% +0.01*ENLv(4)/V(k+2)/2;
%         end;
%         % riverine & sediment fluxes
%         if(fsed)                                           % #!
%             for k=1:3
%                 if(omegCSv(k)<1)
%                     kcarb=kcarb0;
%                 end
%                 CAp(k  ) = CAp(k  )+1*Fin*CafactCarb*Aoc/V(k)/nOC          ...  % Fin  L
%                     +1*FSi*CafactSi*Aoc/V(k)/nOC            ...    % #! Si
%                     -0*kcarb*(omegCSv(k)-1)^2*Aoc/V(k)/nOC...
%                     +4/V(k)/nOC...
%                     -4/V(k)/nOC;
% 
%                 CAp(k  ) = CAp(k  )+1*sum(dissCav(k,1   :n1))/V(k  ); % diss L
%                 % dissCavv(kt) = sum(dissCav(k,1   :n1))/V(k  );
%                 CAp(k+3) = CAp(k+3)+1*sum(dissCav(k,n1+1:n2))/V(k+3); % diss I
%                 % dissCavv(kt) = dissCavv(kt)+sum(dissCav(k,n1+1:n2))/V(k+3);
%                 CAp(k+6) = CAp(k+6)+1*sum(dissCav(k,n2+1:Ns))/V(k+6); % diss D
%                 % dissCavv(kt) = dissCavv(kt)+sum(dissCav(k,n2+1:Ns))/V(k+6);
%             end;
%             if(ftys)
%                 k=11;
%                 if(omegCSv(5)<1)
%                     kcarb=kcarb0;
%                 end
%                 CAp(k  ) = CAp(k  )+1*Fin*CafactCarb*Aoc/V(k)/nOC          ...  % Fin  L
%                     +1*FSi*CafactSi*Aoc/V(k)/nOC       ...% #! Si
%                     -0*kcarb*(omegCSv(5)-1)^2*Aoc/V(k)/nOC...
%                     +4/V(k)/nOC...
%                     -4/V(k)/nOC;
% 
%                 CAp(k  ) = CAp(k  )+1*sum(dissCav(4,1   :n1))/V(k  ); % diss L
%                 % dissCavv(kt) = dissCavv(kt)+sum(dissCav(4,1   :n1))/V(k  );
%                 CAp(k+1) = CAp(k+1)+1*sum(dissCav(4,n1+1:n2))/V(k+1); % diss I
%                 % dissCavv(kt) = dissCavv(kt)+sum(dissCav(4,n1+1:n2))/V(k+1);
%                 CAp(k+2) = CAp(k+2)+1*sum(dissCav(4,n2+1:Ns))/V(k+2); % diss D
%                 % dissCavv(kt) = dissCavv(kt)+sum(dissCav(4,n2+1:Ns))/V(k+2);
%             end;
%         end;
%     end;
    
    if(CAvflag == 1 || CAvflag == 2)
        CAp = THmfun(CA,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
        % Checking balances mol/yr
        % sum(ECALv)-nu*sum(EALv)/2 - Total Available for deposition = 3.8856e+013
        % sum(dissCav(:)) - total diss/bur = 1.7278e+013
        % (Fin+FSi)*Aoc - Total riverine fluxes = 2.1578e+013
        % Total available for deposition = total diss + total riverine; but should
        % it be? ot tot riv = total avail for dep - total diss?
        % bio pump CaCO3, Aorg
        for k=1:3
            CAp(k  ) = CAp(k );% -      ECALv(k)/V(k  );
            CAp(k+3) = CAp(k+3);%+eI*(x*ECALv(k)/V(k+3)); % 1#!
            CAp(k+6) = CAp(k+6);%+oI*(x*ECALv(k)/V(k+6))...% 1#!
                %+   nu*EALv(k)/V(k+6)/2;%+ 0.01*ENLv(k)/V(k+6)/2;            % D ClmnDiss
        end;
        CAp(10)  = CAp(10);%  - ECAH/V10;
        for k=7:9                                              % DA,DI,DP
            CAp(k )  = CAp(k );%  +(eI+oI)*(x*ECAH/V(k))...%   1#!
                %+nu*EAH/V(k)   *gp(k) / 2;  %
        end;
        if(ftys)
            k = 11;
            CAp(k  ) = CAp(k  );%-      ECALv(4)/V(k  );
            CAp(k+1) = CAp(k+1);%+eI*(x*ECALv(4)/V(k+1));
            CAp(k+2) = CAp(k+2);%+oI*(x*ECALv(4)/V(k+2))...
                %+   nu*EALv(4)/V(k+2)/2;% +0.01*ENLv(4)/V(k+2)/2;
        end;
        % riverine & sediment fluxes
        if(fsed)                                           % #!
            for k=1:3
                if(omegCSv(k)<1)
                    kcarb=kcarb0;
                end
                CAp(k  ) = CAp(k  );%+1*Fin*CafactCarb*Aoc/V(k)/nOC          ...  % Fin  L
                    %+1*FSi*CafactSi*Aoc/V(k)/nOC            ...    % #! Si
                    %-0*kcarb*(omegCSv(k)-1)^2*Aoc/V(k)/nOC...
                   % +4/V(k)/nOC...
                   % -4/V(k)/nOC;

                CAp(k  ) = CAp(k  );%+1*sum(dissCav(k,1   :n1))/V(k  ); % diss L
                % dissCavv(kt) = sum(dissCav(k,1   :n1))/V(k  );
                CAp(k+3) = CAp(k+3);%+1*sum(dissCav(k,n1+1:n2))/V(k+3); % diss I
                % dissCavv(kt) = dissCavv(kt)+sum(dissCav(k,n1+1:n2))/V(k+3);
                CAp(k+6) = CAp(k+6);%+1*sum(dissCav(k,n2+1:Ns))/V(k+6); % diss D
                % dissCavv(kt) = dissCavv(kt)+sum(dissCav(k,n2+1:Ns))/V(k+6);
            end;
            if(ftys)
                k=11;
                if(omegCSv(5)<1)
                    kcarb=kcarb0;
                end
                CAp(k  ) = CAp(k  );%+1*Fin*CafactCarb*Aoc/V(k)/nOC          ...  % Fin  L
                    %+1*FSi*CafactSi*Aoc/V(k)/nOC       ...% #! Si
                    %-0*kcarb*(omegCSv(5)-1)^2*Aoc/V(k)/nOC...
                    %+4/V(k)/nOC...
                    %-4/V(k)/nOC;

                CAp(k  ) = CAp(k  );%+1*sum(dissCav(4,1   :n1))/V(k  ); % diss L
                % dissCavv(kt) = dissCavv(kt)+sum(dissCav(4,1   :n1))/V(k  );
                CAp(k+1) = CAp(k+1);%+1*sum(dissCav(4,n1+1:n2))/V(k+1); % diss I
                % dissCavv(kt) = dissCavv(kt)+sum(dissCav(4,n1+1:n2))/V(k+1);
                CAp(k+2) = CAp(k+2);%+1*sum(dissCav(4,n2+1:Ns))/V(k+2); % diss D
                % dissCavv(kt) = dissCavv(kt)+sum(dissCav(4,n2+1:Ns))/V(k+2);
            end;
        end;
    end;
    if(CAvflag>0)
        dissCavv(kt) = sum(dissCav(:));
    end
    EALvv(kt) = sum (EALv);
    ECALvv(kt) = sum(ECALv);
    ECAHv(kt) = sum (ECAH);
    EAHv(kt) = sum (EAH);
    EPHv(kt) = EPH;
    PPHv(kt) = PPH;
 

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
            doxp(k+6) = doxp(k+6) - oIo*EXLv(k)/V(k+6)*REDO2C*mmoxD; % D #!(tot or Corg)
        end;
        for k=7:9 % DA,DI,DP
            mmoxD     = dox(k  )/(dox(k  )+KMMOX)*sign(sign(dox(k  ))+1);
            doxp(k )  = doxp(k  ) -    EXH/V(k)*gp(k)*REDO2C*mmoxD; % #!  EC or EP
        end;
        doxp(10)  = doxp(10 ) +    (eI+oIo)*ECH/V10       *REDO2C;       % H
        if(ftys)
            k = 11;
            mmoxI     = dox(k+1)/(dox(k+1)+KMMOX)*sign(sign(dox(k+1))+1);
            mmoxD     = dox(k+2)/(dox(k+2)+KMMOX)*sign(sign(dox(k+2))+1);
            doxp(k  ) = doxp(k  ) +    ECLv(4)/V(k  )*REDO2C;
            doxp(k+1) = doxp(k+1) - eI*EXLv(4)/V(k+1)*REDO2C*mmoxI;
            doxp(k+2) = doxp(k+2) - oIo*EXLv(4)/V(k+2)*REDO2C*mmoxD;
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
        ccp(k )  = ccp(k )  +    (eI+oI)*EXHCC/V(k)*gp(k)...    % #!  EC or EP
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
            if(omegCSv(k)<1)
                kcarb=kcarb0;
            end
            ccp(k  ) = ccp(k  ) + 2*Fin13*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                + 2*FSi13*Aoc/V(k)/nOC         ... % #! Si
                + 1*Faom13*Aoc/V(k)/nOC...
                - 1*Fmeth13*Aoc/V(k)/nOC...
                - 1*Focb  *Aoc/V(k)/nOC*alp*Rb(k)...  % #! krgn
                - 0*kcarb*Rin*(omegCSv(k)-1)^2*Aoc/V(k)/nOC;

            ccp(k  ) = ccp(k  ) + sum(dissC13v(k,1   :n1))/V(k  ); % diss L
            ccp(k+3) = ccp(k+3) + sum(dissC13v(k,n1+1:n2))/V(k+3); % diss I
            ccp(k+6) = ccp(k+6) + sum(dissC13v(k,n2+1:Ns))/V(k+6); % diss D
        end;
        if(ftys)
            k=11;
            if(omegCSv(5)<1)
                kcarb=kcarb0;
            end
            ccp(k  ) = ccp(k  ) + 2*Fin13*Aoc/V(k)/nOC         ... % Fin  L wthr:2
                + 2*FSi13*Aoc/V(k)/nOC         ... % #! Si
                + 1*Faom13*Aoc/V(k)/nOC...
                - 1*Fmeth13*Aoc/V(k)/nOC...
                - 1*Focb  *Aoc/V(k)/nOC*alp*Rb(k)...  % #! krgn
                - 0*kcarb*Rin*(omegCSv(5)-1)^2*Aoc/V(k)/nOC;
            ccp(k  ) = ccp(k  ) + sum(dissC13v(4,1   :n1))/V(k  ); % diss L
            ccp(k+1) = ccp(k+1) + sum(dissC13v(4,n1+1:n2))/V(k+1); % diss I
            ccp(k+2) = ccp(k+2) + sum(dissC13v(4,n2+1:Ns))/V(k+2); % diss D
        end;
    end;

    %=================== 44CALCIUM  ======================%
    if(CAvflag == 1 || CAvflag == 2)
        cacp = THmfun(cac,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI);
        % bio pump CaCO3, Aorg
        for k=1:3
            cacp(k  ) = cacp(k ) -      ECLvCAC(k)/V(k  );%+ENLv(k)/V(k  ) ;
            cacp(k+3) = cacp(k+3)+eI*(x*ECLvCAC(k)/V(k+3));%-ENLv(k)/V(k+3));   % 1#!
            cacp(k+6) = cacp(k+6)+oI*(x*ECLvCAC(k)/V(k+6))...%-ENLv(k)/V(k+6))... % 1#!
                +   nu*ECLvCAC(k)/V(k+6);                   % D ClmnDiss
        end;
        cacp(10)  = cacp(10)  - ECHCAC/V10;%        + ENH/V10;
        for k=7:9                                              % DA,DI,DP
            cacp(k )  = cacp(k )  +(eI+oI)*(x*ECHCAC/V(k))...%     - ENH/V(k))*gp(k)...% 1#!
                +nu*ECHCAC/V(k)                *gp(k);  %
        end;
        if(ftys)
            k = 11;
            cacp(k  ) = cacp(k  )-      ECLvCAC(4)/V(k  );%+ENLv(4)/V(k  ) ;
            cacp(k+1) = cacp(k+1)+eI*(x*ECLvCAC(4)/V(k+1));%-ENLv(4)/V(k+1));
            cacp(k+2) = cacp(k+2)+oI*(x*ECLvCAC(4)/V(k+2))...%-ENLv(4)/V(k+2))... %
                +   nu*ECLvCAC(4)/V(k+2);
        end;
        % riverine & sediment fluxes
        if(fsed)                                           % #!
            for k=1:3
                if(omegCSv(k)<1)
                    kcarb=kcarb0;
                end
                cacp(k  ) = cacp(k  )+1*Fin44*Aoc/V(k)/nOC          ...  % Fin  L
                    +1*FSi44*Aoc/V(k)/nOC              ...  % #! Si
                    - 0*kcarb*RinCA*(omegCSv(k)-1)^2*Aoc/V(k)/nOC...
                    +4*RinCA/V(k)/nOC...
                    -4*RinCA/V(k)/nOC;

                cacp(k  ) = cacp(k  )+1*sum(dissC44v(k,1   :n1))/V(k  ); % diss L
                cacp(k+3) = cacp(k+3)+1*sum(dissC44v(k,n1+1:n2))/V(k+3); % diss I
                cacp(k+6) = cacp(k+6)+1*sum(dissC44v(k,n2+1:Ns))/V(k+6); % diss D
            end;
            if(ftys)
                k=11;
                if(omegCSv(5)<1)
                    kcarb=kcarb0;
                end
                cacp(k  ) = cacp(k  )+1*Fin44*Aoc/V(k)/nOC          ...  % Fin  L
                    +1*FSi44*Aoc/V(k)/nOC       ...% #! Si
                    - 0*kcarb*RinCA*(omegCSv(5)-1)^2*Aoc/V(k)/nOC...
                    +4*RinCA/V(k)/nOC...
                    -4*RinCA/V(k)/nOC;

                cacp(k  ) = cacp(k  )+1*sum(dissC44v(4,1   :n1))/V(k  ); % diss L
                cacp(k+1) = cacp(k+1)+1*sum(dissC44v(4,n1+1:n2))/V(k+1); % diss I
                cacp(k+2) = cacp(k+2)+1*sum(dissC44v(4,n2+1:Ns))/V(k+2); % diss D
            end;
        end;
    end;


    %=================== C  atm ====================%
    Cp    = sum(                                 ...%
        -kasv(kkv).*(pco2a-pco2(kkv))          ...% mol/m2/y
        )/Aoc                                ...%
        -1*Fin   + Fvc   - 2*FSi   + Focw;         % wthr #!
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
    %=================== O atm =====================%

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

            if(CAvflag>0)
                fcavAp(l)   = ( fbAca(l)*rscavA(l) - fbAca(l)*rdvAca(l) - wcavA(l) )/GAca(l);
                fcavIp(l)   = ( fbIca(l)*rscavI(l) - fbIca(l)*rdvIca(l) - wcavI(l) )/GIca(l);
                fcavPp(l)   = ( fbPca(l)*rscavP(l) - fbPca(l)*rdvPca(l) - wcavP(l) )/GPca(l);

                f44cavAp(l) = ( fbAca(l)*r44scvA(l) - fbAca(l)*r44dvA(l) - w44cvA(l) )/GAca(l);
                f44cavIp(l) = ( fbIca(l)*r44scvI(l) - fbIca(l)*r44dvI(l) - w44cvI(l) )/GIca(l);
                f44cavPp(l) = ( fbPca(l)*r44scvP(l) - fbPca(l)*r44dvP(l) - w44cvP(l) )/GPca(l);
            end;

            if(ftys)
                fcvTp(l)   = ( fbT(l)*rscvT(l) - fbT(l)*rdvT(l) - wcvT(l) )/GT(l);
                f13cvTp(l) = ( fbT(l)*r13scvT(l) - fbT(l)*r13dvT(l) - w13cvT(l) )/GT(l);

                if(CAvflag>0)
                    fcavTp(l)   = ( fbTca(l)*rscavT(l) - fbTca(l)*rdvTca(l) - wcavT(l) )/GTca(l);
                    f44cavTp(l) = ( fbTca(l)*r44scvT(l) - fbTca(l)*r44dvT(l) - w44cvT(l) )/GTca(l);
                end;
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
            if(CAvflag>0)
                f44cavAp(l)= 0;
                f44cavIp(l)= 0;
                f44cavPp(l)= 0;

                m44cAp(l) = 0;
                m44cIp(l) = 0;
                m44cPp(l) = 0;
            end;
            if(ftys)
                fcvTp(l)  = 0;
                mcTp(l)   = 0;
                f13cvTp(l)= 0;
                m13cTp(l) = 0;
                if(CAvflag>0)
                    f44cvTp(l)= 0;
                    m44cTp(l) = 0;
                end;
            end;
        end;
    end; % sediments


    tmpCO(1) = cp(kb);
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
        % DTS = 06.e3;   % n06 [10 1] e3 years 04
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
            %CAp(kb)=0;%%% THIS ADDED
            %CAp(kT)=0; %%%THIS ADDED
            Cinpc(kb)=(CBl/12/DTS)*(12/1e15);
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
            fB  = .734*(  fx);  %.4 .75 Atl fraction .9 .7 .47 .6 n.5 .734
            fM  = .734*(1-fx);  %.4 .75 Atm fraction .9 .7 .47 .6 n.5 .734
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
        %Cinpc(kt)=cp(kb);
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
        cacp(11:13) =          0;         % T
    end;

    % all in one
    if(CAvflag == 0)
        if(0)
            cp(11:13) =          0;         % T
            ap(11:13) =          0;         % T
            pp(11:13) =          0;         % T
            ccp(11:13) =          0;         % T
        end;

        % all in one

        if(fdox)
            yp(     1:  Nb     ) = cp;
            yp(  Nb+1:2*Nb     ) = ap;
            yp(2*Nb+1:3*Nb     ) = pp;
            yp(3*Nb+1:4*Nb     ) = doxp;
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
        else % fdox
            yp(     1:  Nb     ) = cp;
            yp(  Nb+1:2*Nb     ) = ap;
            yp(2*Nb+1:3*Nb     ) = pp;
            yp(3*Nb+1:4*Nb     ) = ccp;
            yp(4*Nb+1          ) = Cp;
            yp(4*Nb+2          ) = CCp;
            if(fsed)
                yp(4*Nb+3     :4*Nb+2+1*Ns) = fcvAp;
                yp(4*Nb+3+1*Ns:4*Nb+2+2*Ns) = fcvIp;
                yp(4*Nb+3+2*Ns:4*Nb+2+3*Ns) = fcvPp;
                yp(4*Nb+3+3*Ns:4*Nb+2+4*Ns) = f13cvAp;
                yp(4*Nb+3+4*Ns:4*Nb+2+5*Ns) = f13cvIp;
                yp(4*Nb+3+5*Ns:4*Nb+2+6*Ns) = f13cvPp;
                if(ftys)
                    yp(4*Nb+3+6*Ns:4*Nb+2+7*Ns) = fcvTp;
                    yp(4*Nb+3+7*Ns:4*Nb+2+8*Ns) = f13cvTp;
                end;
            end;
        end; % fdox
    else
        if(fdox)
            yp(     1:  Nb     ) = cp;
            yp(  Nb+1:2*Nb     ) = ap;
            yp(2*Nb+1:3*Nb     ) = pp;
            yp(3*Nb+1:4*Nb     ) = CAp ;
            yp(4*Nb+1:5*Nb     ) = doxp;
            yp(5*Nb+1:6*Nb     ) = cacp;
            yp(6*Nb+1:7*Nb     ) = ccp;
            yp(7*Nb+1          ) = Cp;
            yp(7*Nb+2          ) = CCp;
            if(fsed)
                yp(7*Nb+3     :7*Nb+2+1*Ns) = fcvAp;
                yp(7*Nb+3+1*Ns:7*Nb+2+2*Ns) = fcvIp;
                yp(7*Nb+3+2*Ns:7*Nb+2+3*Ns) = fcvPp;
                yp(7*Nb+3+3*Ns:7*Nb+2+4*Ns) = f13cvAp;
                yp(7*Nb+3+4*Ns:7*Nb+2+5*Ns) = f13cvIp;
                yp(7*Nb+3+5*Ns:7*Nb+2+6*Ns) = f13cvPp;
                yp(7*Nb+3+6*Ns:7*Nb+2+7*Ns) = fcavAp;
                yp(7*Nb+3+7*Ns:7*Nb+2+8*Ns) = fcavIp;
                yp(7*Nb+3+8*Ns:7*Nb+2+9*Ns) = fcavPp;
                yp(7*Nb+3+9*Ns:7*Nb+2+10*Ns) = f44cavAp;
                yp(7*Nb+3+10*Ns:7*Nb+2+11*Ns) = f44cavIp;
                yp(7*Nb+3+11*Ns:7*Nb+2+12*Ns) = f44cavPp;
                if(ftys)
                    yp(7*Nb+3+12*Ns:7*Nb+2+13*Ns) = fcvTp;
                    yp(7*Nb+3+13*Ns:7*Nb+2+14*Ns) = f13cvTp;
                    yp(7*Nb+3+14*Ns:7*Nb+2+15*Ns) = fcavTp;
                    yp(7*Nb+3+15*Ns:7*Nb+2+16*Ns) = f44cavTp;
                end;
            end;
        else % fdox
            yp(     1:  Nb     ) = cp;
            yp(  Nb+1:2*Nb     ) = ap;
            yp(2*Nb+1:3*Nb     ) = pp;
            yp(3*Nb+1:4*Nb     ) = CAp;
            yp(4*Nb+1:5*Nb     ) = cacp;
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
        end; % fdox
    end; %CAvflag

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
