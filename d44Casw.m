global fFw X
fFw=2.0; % Factor used to increase weathering flux
X=-0.0;  % Changing d44Ca of the the input (d44w+X)
t0=0.;
tfinal=10e5;

Moc=1.3243e21;

y0=[2.6486e19 0.00];
tspan=[tfinal - t0];
format compact
options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
    'InitialStep',1.,'Maxstep',tspan/1000.);
tic
[T,Y] = myode15s(@d44Cadiff,[t0 tfinal],y0,options);
toc
(max(Y(:,1))/y0(1,1))*100-100  %inventory increase in percent
DdCa=((min(Y(:,2))+y0(1,2)))        % Delta delta seawater calcium change in permil
DdCaM=((max(Y(:,2))+y0(1,2)))
subplot(211)
plot(T,Y(:,1)/Moc*1000,'-m');
title('numerical solution, \delta_w')
xlabel('years');
ylabel('Ca concentration [mmol/kg]');
subplot(212)
plot(T,Y(:,2),'g');
xlabel('years');
ylabel('\delta^{44}_{sw}');