%% martin curve
z=[100:100:6000]; %,
Aoc = 3.49e+014; %m2
b=0.858; %attenuation factor default: 0.858; upper limit 1.33 buesseler et al 2007
% LOSCAR low latitude export
Cexp = 395.19e12; %mol C/y  %460 534 | 524 611
F100 = Cexp/Aoc % mol C/m2/y
F= F100.*(z./100).^-b;

% Berger '87. formulation
% NOTE TO MYSELF 5/12/2017:
% Berger curve at 100 meters is 16% of primary production. I need
% to normalize the berger flux at this depth so that both Martin 
% and Berger are exactly the same at 100m in mol C/m2/y
Fb = 9.*F100/0.16./z+0.7*F100/0.16./z.^0.5;

F4000 = F(40)*Aoc;
figure
% plot(F/F100*100,z,'r',Fb/F100*100,z, 'b')
plot(F,z,'r',Fb,z, 'b')
legend('Martin', 'Berger')
set(gca, 'YDir','reverse')
% xlabel('% of Primary Production')
xlabel('Global org. C rain [mol C m^{-2} y^{-1}]')
ylabel('Depth (m)')
title ('Global ocean')


Xmax=16;
Xmin=1.054;
mya=1;
mya=2.981;
myb=100;
mya+(8-Xmin)*(myb-mya)/(Xmax-Xmin)

% [a,b]=find(EPLvva{8}>=max(EPLvva{8}))

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





