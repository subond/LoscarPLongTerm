function ydt = payneD44CaDiff(t,y)

global  Friv Fhyd Fpw driv dhyd dpw epsc kcarb MCa0 Friv0 Fcarbv kt Frivv...
    Moc ksp omgCv flag epscv kcarb2 drivv;
MCa=y(1);
dsw = y(2);
co3 = y(3);
TC = y(4);
TA = y(5);
omgC = (co3 * MCa/Moc)/ksp; % MCa/Moc converting from total moles Ca to moles/kg

tStart =1e6;
tEnd= tStart+1e5;
if(t>=tStart && t<=tEnd)
    Friv = Friv0*3.3; %3.3  2.0 1.7143
    if(flag==0)
        Fcarb = kcarb*(0.823);  % 0.823 Payne '10, 17.7% reduction
        epsc=-1.4;
%         epsc=flin(tStart, tEnd, -1.4, -0.7, t);
%         driv=flin(tStart, tEnd, -0.6, -0.1, t);
    else
        Fcarb = kcarb*(1);
        epsc=-1.4;
    end;
    if(flag==2)
        Fcarb = kcarb2*(omgC-1)^2;
%          epsc=flin(1e6,1e6+10e4,-1.4,-1.8,t);
            epsc=-1.4;
    end
elseif (t>tEnd)
    Friv = Friv0;
    if(flag==2)
        Fcarb = kcarb2*(omgC-1)^2;
         epsc=-1.4;
    else
        Fcarb = kcarb*(MCa/MCa0)^2;
         epsc=-1.4;
%         driv=-0.1;
    end;  
else
    Friv = Friv0;
    if(flag==2)
        Fcarb = kcarb2*(omgC-1)^2;
         epsc=-1.4;
    else
        Fcarb = kcarb*(MCa/MCa0)^2;
        epsc=-1.4;
%         driv=-0.6;
    end; 
end
if(t>1.1e6 && t<=1.2e6)
    Friv = Friv0*1;  %0.4
end
% if(t>=1e6 && t<=1.15e6)
%    epsc=flin(1e6, 1.15e6, -1.4, -1.2, t); 
% end
% if(t>1.15e6 && t<=1.3e6)
%    epsc=flin(1.15e6, 1.3e6, -1.2, -1.7, t); 
% end
% if(t>1.3e6)
%    epsc=flin(1.3e6, 1.5e6, -1.7, -1.4, t); 
% end

Fcarbv(kt+1)=Fcarb;
Frivv(kt+1) = Friv;
omgCv(kt+1) = omgC;
epscv(kt+1) = epsc;
drivv(kt+1) = driv;

ydt = [Friv + Fhyd + Fpw - Fcarb % Total Moles Ca per year
    (Friv*(driv-dsw) + Fhyd*(dhyd-dsw) + Fpw*(dpw-dsw) - Fcarb*(epsc))/MCa
%     (Friv*(driv) + Fhyd*(dhyd) + Fpw*(dpw) - Fcarb*(dsw+epsc))/MCa
    ((Friv + Fhyd + Fpw - Fcarb)/(2*Moc)) % Moles/kg CO32=
    (Friv + Fhyd + Fpw - Fcarb)/Moc   % TC - Not used
    (2*(Friv + Fhyd + Fpw) - 2*Fcarb)/Moc]; %TA not used


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