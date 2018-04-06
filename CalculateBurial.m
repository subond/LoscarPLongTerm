%% Calculates Organic Carbon burial between two ocean depths when supplied with Organic Carbon export and "b" term in Martin Curve
function [Fbct, Fraint,Fbci,Fraini,Fbcd,Fraind,Fremi,Fremd]=CalculateBurial(Fp,b,m)
global fbetta

Fbct=0;
Fraint=0;
Fbci=0;
Fbcd=0;
Fraini=0;
Fraind=0;
for i=1:2
    if(i==1)
        lb=-1.0;
        ub=-0.1;
    else
        lb=-4.5;
        ub=-1.0;
    end       
if(lb>=-1.0)
    ps1 = 0.307; ps2 = 0.624; ps3 = 0.430;ps4 = 0.991;
    A=ps1*lb.^3 + ps2*lb.^2 + ps3*lb + ps4;
    
else
    ps1 = 0.02; ps2 = 0.103; ps3 = 0.219;ps4 = 1.025;
end
% A=ps1*zs.^3 + ps2*zs.^2 + ps3*zs + ps4;
b1 = 0.411; b2 = 0.153; b3 = 1.0; % default: b1 = 0.411; b2 = 0.153; b3 = 1.0;
sig=(1/600)*1000;
% Martin curve
if (m)
    bur = @(x) ((fbetta).*(b1*exp(-b3*(-x)) + b2).*(3*ps1*x.^2 + 2*ps2*x + ps3).*((x/(-0.1)).^(-b)));
    rain = @(x) ((x/(-0.1)).^(-b));
else
% Shaffer curve (see Bjerrum et al., 2006)
    bur = @(x) ((b1*exp(-b3*(-x)) + b2).*(3*ps1*x.^2 + 2*ps2*x + ps3).*exp(-sig*(-x)));
    rain = @(y) (3*ps1*y.^2 + 2*ps2*y + ps3).*exp(-sig*(-y));
end
if(i==1)
    Fbci  = Fp*integral(bur,lb,ub);
    Fraini= Fp*integral(rain,lb,ub);
    if(m)
        Fremi = Fp*(1-(-1.0/(-0.1)).^(-b))-Fbci;
    else
        Fremi = Fp*(1-A*exp(sig*(-1.0)))-Fbci;
    end
else
    Fbcd  = Fp*integral(bur,lb,ub);
    Fraind= Fp*integral(rain,lb,ub);
    Fremd= Fp-(Fremi+Fbci+Fbcd);
end
% Fbct  = Fbct+Fp*integral(fun,lb,ub);
% Fraint= Fraint+Fp*integral(fun2,lb,ub);
end
Fbct  = Fbci+Fbcd;
Fraint= Fraini+Fraind;