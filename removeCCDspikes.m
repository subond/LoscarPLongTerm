function [ccdA,ccdI,ccdP] = removeCCDspikes(fcA_ ,fcI_, fcP_)
global gctl
%REMOVECCDSPIKES Summary of this function goes here
%   Detailed explanation goes here

zv   = [000.:10:6000.]; % z, continuous (1 or 10 m)
  %------------ CCD --------------------%
        fccd  = 0.05; % 0.10 0.05
        dd    = 0.01;
        nd    = 0;
        fccdv = [fccd-nd*dd:dd:fccd+nd*dd];
        ld    = length(fccdv);
        clear jccdA jccdI jccdP;
        dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000;
        for i=1:gctl
            yA(i,:)         = interp1(dsv,fcA_(i,:),zv,'PCHIP');
            yI(i,:)         = interp1(dsv,fcI_(i,:),zv,'PCHIP');
            yP(i,:)         = interp1(dsv,fcP_(i,:),zv,'PCHIP');
%             yT(i,:)         = interp1(dsv,fcT(i,:),zv,'PCHIP');
            for k=1:ld
                [tmp, jccdA(i,k)] = min(abs(yA(i,150:end)-fccdv(k)));
                [tmp, jccdI(i,k)] = min(abs(yI(i,150:end)-fccdv(k)));
                [tmp, jccdP(i,k)] = min(abs(yP(i,150:end)-fccdv(k)));
%                 [tmp, jccdT(i,k)] = min(abs(yT(i,:)-fccdv(k)));
            end;
        end;

        if(nd == 0)
            ccdA = zv(jccdA)+zv(150);
            ccdI = zv(jccdI)+zv(150);
            ccdP = zv(jccdP)+zv(150);
%             ccdT = zv(jccdT);
        else
            ccdA = sum(zv(jccdA),2)/ld;
            ccdI = sum(zv(jccdI),2)/ld;
            ccdP = sum(zv(jccdP),2)/ld;
%             ccdT = sum(zv(jccdT),2)/ld;
        end
end

