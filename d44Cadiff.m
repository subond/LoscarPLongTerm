function ydt=d44Cadiff(t,y)
global fFw X
 Nca=y(1);
 dsw=y(2);

Fw0=2.15e13; % 1.5e13 mol/year (Fwc + Fsi)
dw0=-1.3;    % d44Ca of weathering in permil -1.3 griffith et al 2008
Dsed= -1.3;
Fsed0=2.15e13; % 1.5e13 mol/year
if(1)
if(t<50.e3)
    Fw=Fw0*fFw;  % incresing weathering by a factor
    Fsed=Fsed0;
    dw=dw0+X;
else
    Fw=Fw0;
    Fsed=Fsed0;
    dw=dw0;
end
end
%     Fw=Fw0;
%     Fsed=Fsed0;

ydt=[(Fw-Fsed)
(Fw*(dw-dsw)-(Fsed*Dsed))/Nca]; % Evalute Ca concentration at time t