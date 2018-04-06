
%% martin curve
clear all
close all
z=[100:100:4500]; %,
Aoc = 3.49e+014; %m2
b=0.51; %attenuation factor default: 0.858; upper limit 1.33 buesseler et al 2007
b0=0.858; %attenuation factor default: 0.858; upper limit 1.33 buesseler et al 2007
b2=1.33; %attenuation factor default: 0.858; upper limit 1.33 buesseler et al 2007
% LOSCAR low latitude export
Cexp = 395.19e12; %mol C/y  %460 534 | 524 611
F100 = 100;%Cexp/Aoc; % mol C/m2/y
F= F100.*(z./100).^-b;
F2= F100.*(z./100).^-b2;

% Hypsometry polyn. coefficients
ps1 = 0.307; ps2 = 0.624; ps3 = 0.430;ps4 = 0.991;
% Burial efficiency polyn. coefficients
b1 = 0.411; b2a = 0.153; b3 = 1.0;

zkm= [-100:-100:-4500]./1000; % depth (z) in km. negative down
zp = zkm(1); % export layer depth
Dz = zkm(2); % The change in depth, assuming that the beginning depth is zp 
A=ps1*zkm.^3 + ps2*zkm.^2 + ps3*zkm + ps4; % Hypsometry curve for 0 to 1000 m
dA = 3*ps1*zkm.^2 + 2*ps2.*zkm + ps3; % derivative of the Hypso curve
Bz = F100*((b1*exp(-b3*-Dz) + b2a)*dA(2)*((Dz/zp)^(-b))); % Burial at Dz
RMz = F100*(1-(Dz/zp)^(-b))-Bz; % Remineralization at depth Dz
dRNz = F100-F100*(Dz/zp)^(-b); % Change in rain from depth zp do Dz
dRNz
(Bz + RMz)

% dRA = F100*(-(b*(z/100)^(-b)))/z;

% RA(1:10)+diff(RA(1:11))
% za=[100:100:4500];
% syms f(za)
% f(za) = F100.*(za./100).^-b;
% dfs = diff(f,za);
% return
% Berger '87. formulation
% NOTE TO MYSELF 5/12/2017:
% Berger curve at 100 meters is 16% of primary production. I need
% to normalize the berger flux at this depth so that both Martin 
% and Berger are exactly the same at 100m in mol C/m2/y
Fb = 9.*F100/0.16./z+0.7*F100/0.16./z.^0.5;

F4000 = F(40)*Aoc;

% Shaffer 1996
sig=1/600;
SL=0;
% *1./(exp(-sig*(100))) this part added in order to normalize the curve so
% there is 100% rain at 100 meters depth
Fs = F100.*exp(-sig*(SL--z)).*1./(exp(-sig*(F100)));

sig2=(1/600)*1000
z2=[-0:-0.1:-4.5];
Fs2 = F100.*exp(-sig2*(SL-z2));
Fm2 = F100.*(z2./-0.1).^-b;


% My formulation 
Q10 = 2.6;
k = log(Q10)/10; % per degree C. so 0.693 (or log(2)/2) per 10 degrees C 
T = [20:-1:0];
T0=20;

Dt=10;
kDt = k*Dt;
% exp(k*Dt)
% if(kDt==0)
%     kDt=1;
% end
Fy = F100.*(exp(-sig*exp(kDt)*(SL--z))).*1./(exp(-sig*exp(kDt)*(F100))); 


figure
% plot(F/F100*100,z,'r',Fb/F100*100,z, 'b')
plot(F,z,'r',Fb,z, 'b',Fs,z, 'g')
legend('Martin', 'Berger','Shaffer')
set(gca, 'YDir','reverse')
% xlabel('% of Primary Production')
xlabel('Global org. C rain [mol C m^{-2} y^{-1}]')
ylabel('Depth (m)')
title ('Global ocean')

figure
% plot(F/F100*100,z,'r',F2/F100*100,z, 'b')
plot(F,z./-1000,'r--',F2,z./-1000, 'bd-')
legend(['Martin, b=',num2str(b)], ['Martin, b=',num2str(b2)])
% set(gca, 'YDir','reverse')
 xlabel('% of Carbon Export at 100 m reaching the sea floor')
% xlabel('Global org. C rain [mol C m^{-2} y^{-1}]')
ylabel('Depth (10^3 m)')
title ('Global ocean')

figure
plot(Fs,z,'r-',Fy,z,'+b')
set(gca, 'YDir','reverse')

legend('Shaffer','My dependent on T')

%%% Calculating the difference between the two curves for two different b's
% 
% x0(1) = fzero(F-F2,1);
% x0(2) = fzero(F-F2,1.5);

% x = linspace(x0(1),x0(2),100);
% y1 = x.^2;
% y2 = 2*x;
% figure
% fill([x x(end:-1:1)],[y1 y2(end:-1:1)],'r')
% hold on
% plot(x,y1,x,y2,'LineWidth',2)
% grid on
% a = trapz(x,y2)-trapz(x,y1);

Xmax=16;
Xmin=1.054;
mya=1;
mya=2.981;
myb=100;
mya+(8-Xmin)*(myb-mya)/(Xmax-Xmin)


%Buesseler et al., 2007 temperature profiles
zb = [150 300 500 1000];
t1 = [21.93 13.55 7.62 3.94];
t2 = [2.17 3.37 3.17 2.57];

figure
plot(t1,zb,'r', t2,zb,'g')
legend('Station Aloha', 'K2')
set(gca, 'YDir','reverse')
xlabel('Temperature')
ylabel('Depth')

% [a,b]=find(EPLvva{8}>=max(EPLvva{8}))
return
% Organic export in Atlantic
CexpA = EPLv(1)+EPH*gp(7); %Organic export in Atlantic
A1=A(1); %Atlantic area
F100A = CexpA/A1
FA= F100A.*(z./100).^-b;
% fshA = 0.5602;
% Berger '87. formulation
FAb = 9.*F100A/0.16./z+0.7*F100A/0.16./z.^0.5;

figure
plot(FA,z,'r',FAb,z, 'b')
legend('Martin', 'Berger')
set(gca, 'YDir','reverse')
xlabel('Org C. rain in the Atlantic Ocean [mol C m^{-2} y^{-1}]')
ylabel('Depth (m)')
title ('Atlantic')

%%%%%% GLOBAL %%%%%%%%

% 395.19e12 organic C export in low latitude
% FAA = 8.8928e+012 total inorganic C rain in the deep atlantic
% FAA + FII +... 4.0698e+013 total inorg C for all oceans
% above divided by Aoc is 0.1166 mol C/m2/y
% Martin curve at 4000 meters for the enitre ocean gives us 0.0478 mol
% C/m2/y, which is 4.22% of surface export. So inorg C to org C on the
% surface of the sediments is 0.1166/0.0478 = 2.4393 times. This ratio is
% identical to the simple calculations:

% Org rain: 100 * 0.0422 = 4.22
% CaCO3 rain: 100 * (1-0.31)/6.7 = 10.2985
% 10.2985/4.22 = 2.44

%%%%%% Atlantic %%%%%%%

% Corg: 0.08482 mol C/m2/y but 0.06963 if no addition from high latitude
% which is the EPH*gp(7) part in the equation
% CaCO3: 8.8928e+012/A(1) = 0.1699
% CaCO3 twice the size 

% Initial Atlantic org c burial fraction 7.6667e+011,standard run
% 4.6e12*V(1)/(V(1)+V(2)+V(3)+V(11))
% Initial Atlantic org rain 7.6335e+13:
% 4.5801e+014*V(1)/(V(1)+V(2)+V(3)+V(11))
% rain at 4km = *0.0422 = 3.2213e+12
%7.6667e+011/3.2213e+012 = 0.238 => 23.8 percent is burial

% The same but at max export (tv=125):
% orgCb(125): 5.5166e+012
%orgCb(125)*V(1)/(V(1)+V(2)+V(3)+V(11)) = 9.1943e+11
% organin rain (EPLvv(d)+EPH)): 5.3396e+014
% (EPLvv(d)+EPH)*V(1)/(V(1)+V(2)+V(3)+V(11)) = 8.8993e+013
% at 4km = 3.7555e+012
% 9.1943e+11/3.7555e+012 => 0.2448 => 24.5 percent





