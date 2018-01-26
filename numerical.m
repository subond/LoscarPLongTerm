%tspan=[0:10]*10^5; %is this the way to define time?
% %[T,Y] = ode23s(@vdp1000,tspan,10); %solve ODE
% t0=0.;
% tfinal=10e6;
% Ca=0
% tspan=[tfinal - t0];
% format compact
% options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
%     'InitialStep',1.,'Maxstep',tspan/1000.);
% tic
% [T,Y] = myode15s(@vdp1000,[t0 tfinal],Ca,options);
% toc
% plot(T,Y);
% title('numerical solution, k=Fin/10.3*10^-3')
% xlabel('years');
% ylabel('Ca concentration [mol/kg]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global N
N=1.71; % N=1.71 used to get Nca increase of 2 percent
t0=0.;
tfinal=10e4;

Moc=1.3243e21;

y0=[2.6486e19 0.00];
tspan=[tfinal - t0];
format compact
options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
    'InitialStep',1.,'Maxstep',tspan/1000.);
tic
[T,Y] = myode15s(@vdp1000,[t0 tfinal],y0,options);
toc
(max(Y(:,1))/y0(1,1))*100-100  %inventory increase in percent
DdCa=y0(1,2)-((min(Y(:,2))+y0(1,2))/2)        % Delta delta seawater calcium change in permil
subplot(211)
plot(T,Y(:,1),'-m');
title('numerical solution, \d_w')
xlabel('years');
ylabel('Ca concentration [mol/kg]');
subplot(212)
plot(T,Y(:,2),'g');