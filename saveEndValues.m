function []=saveEndValues(caMean, epspcaV)
        dsvp = dsv';
        save                    caMean.DAT  caMean -ASCII -DOUBLE -TABS;
        save                    epspcaV.DAT  epspcaV -ASCII -DOUBLE -TABS;
        save                    d44cabulkV.DAT  d44cabulkV -ASCII -DOUBLE -TABS;
        save vars.tex  Nb Ns d13CBl -ASCII -DOUBLE -TABS;
        save dsv.DAT   dsvp  -ASCII -DOUBLE -TABS;
        save tv.DAT    tv    -ASCII -DOUBLE -TABS;
        save tv11.DAT    tv11    -ASCII -DOUBLE -TABS;
        save c.DAT     c     -ASCII -DOUBLE -TABS;
        save a.DAT     a     -ASCII -DOUBLE -TABS;
        save p.DAT     p     -ASCII -DOUBLE -TABS;
        save kkv.DAT   kkv   -ASCII -DOUBLE -TABS;
        save pco2t.DAT pco2t -ASCII -DOUBLE -TABS;
        save pco2v.DAT pco2v -ASCII -DOUBLE -TABS;
        save co3tv.DAT co3tv -ASCII -DOUBLE -TABS;
        save phtv.DAT  phtv  -ASCII -DOUBLE -TABS;
        save d13c.DAT  d13c  -ASCII -DOUBLE -TABS;
        save d13CA.DAT d13C  -ASCII -DOUBLE -TABS;
        save temp.DAT TCvt  -ASCII -DOUBLE -TABS;

        if(Pfeed)
            save oIv.DAT oIv -ASCII -DOUBLE -TABS;
            save oIpv.DAT  oIpv  -ASCII -DOUBLE -TABS;
            save EPLvv.DAT  EPLvv  -ASCII -DOUBLE -TABS;
            save PPLvv.DAT PPLvv  -ASCII -DOUBLE -TABS;
            save Fcapv.DAT  Fcapv  -ASCII -DOUBLE -TABS;
            save Ffepv.DAT Ffepv  -ASCII -DOUBLE -TABS;
        end
        save V.DAT V  -ASCII -DOUBLE -TABS;

        save TCvt.DAT  TCvt  -ASCII -DOUBLE -TABS;

        if(0)
            save THt.DAT   THt   -ASCII -DOUBLE -TABS;
        end
        
        if(Floegel)
          save REDPC.DAT REDPCv  -ASCII -DOUBLE -TABS;  
          save meanSpco2.DAT meanSpco2v  -ASCII -DOUBLE -TABS; 
          save rREG.DAT rREGv  -ASCII -DOUBLE -TABS; 
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

        save EPH.DAT  EPH  -ASCII -DOUBLE -TABS;
        save PPH.DAT  PPH  -ASCII -DOUBLE -TABS;
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
return