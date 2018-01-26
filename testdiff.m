function ydt=testdiff(t,y)
global C
Fin0=1e12;
Fin=Fin0;
k=Fin0/C;

if(t>=60e3 & t<=100e3)
%     Fin=flin(60e3,100e3,1e12,0,t);
    Fin=interp1([60e3 100e3],[1e12 0], t);
else
    Fin=1e12;
end
ydt = Fin - k*y(1);

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