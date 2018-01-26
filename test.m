clear all
global C
C = 5000e15; %Pg C

y0=[C];

t0=0;
tfinal=1e5;
tspan=[tfinal - t0];

options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
    'InitialStep',1.,'Maxstep',tspan/1000.);
tic
[T,Y] = myode15s(@testdiff,[t0 tfinal],y0,options);
toc
figure
hold on
plot(T,Y)
clear all

C0 = 5000e15; %Pg C
Fin0=1e12;
Fin=Fin0;
k=Fin0/C0;
 t=[1:1000:1e5];
% Fin=ones(1,100)*Fin0;
counter = 1;
for i=t(1):1000:t(end)
%    Fin(counter)=flin(60e3,100e3,1e12,0,i);
    if(i>=60e3 & i<=100e3)
        Fin(counter)=0;
    else
        Fin(counter)=Fin0;
    end
%    mytest(counter)=interp1([60e3 100e3],[1e12 0], i);
   counter= counter+1;
end
% C=Fin/k+exp(-k*t); %Wolfram-alpha solution
Ct = (C0-Fin/k).*exp(-k*t)+Fin/k;
% figure
plot(t,Ct,'r')
hold off
