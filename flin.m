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