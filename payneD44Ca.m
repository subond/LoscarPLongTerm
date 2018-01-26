%%%% Payne d44Ca
%%%%
clear all
% close all
global Friv Fhyd Fpw driv dhyd dpw epsc kcarb MCa0 Friv0 Fcarbv kt Frivv...
    Moc phflag k1k2flag Mgm Cam ksp omgCv flag epscv kcarb2 drivv;
addpath('myode');
kt =1;

%%%Simulation Flags
flag = 1; %% 0: standard Payne et al sim, 1: weathering only (no acid.)  
          %% 2: weathering only but with carbonate saturation feedback



%%% Fluxes in mol/yr
Friv = 14e12;
Friv0=Friv;
Fhyd = 4e12;
Fpw  = 6.3e12;
kcarb = 24.3e12;


% d44Ca in permil
driv = -0.6;
dhyd = -0.25;
dpw = -0.45;
dsw = 0.8965;

epsc = -1.4;

t0=0.;
tfinal=4e6;
Moc=1.3243e21; % mass ocean [kg]

carbc = 180*1e-6; %carbonate ion conc. [moles/kg]
Ca0=10; % Initial calcium conc. moles/kg
MCa0=Moc*Ca0/1000; % Total ocean Ca inventory in moles



%%%%% for Ksp calculation
Cam=10.0e-3;  %modern
Mgm=53e-3;      %modern
Cap=10.0e-3;    %paleo
Mgp=53e-3;      %paleo
TC=20.;                     % temperature
S=35;                  % salinity
P=15;                 % pressure
phflag = 0;
k1k2flag = 1;
[mykspc,x] = ...
     kspfun(TC,S,P,Cap,Mgp);

ksp = mykspc;

omgC = carbc*Cam/ksp(1);

kcarb2 = 24.3e12/(omgC-1)^2;

%%% SOLVER
y0=[MCa0 dsw carbc 2.2828e-3 2.3896e-3];
tspan=[tfinal - t0];
format compact
options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
    'InitialStep',1.,'Maxstep',tspan/1000.);
tic
[T,Y] = myode15s(@payneD44CaDiff,[t0 tfinal],y0,options);
toc



Ca=Y(:,1)/Moc; % Calcium concentration in moles/kg
d44Ca = Y(:,2); % d44 calcium in sea water
myco3 = Y(:,3)*1e6; % carbonate ion concentration [microl mol/kg]using simplified eq. from Zeebe and Westbroek
TCAR=Y(:,4); %Dissolved inorganic carbon in moles/kg
TALK=Y(:,5); % Total alkalinity



% dic1=TCAR;
% for i=1:length(dic1) 
%      alk(i)=TALK(i); 
%      [co2(i),pco2(i),co31(i),ph(i),kh(i),o2(i),kspc(i),kspa(i)]=...
%          dafunPE(dic1(i),alk(i),TC,S,P,Cam,Mgm); 
% end;
% 
% 
% omgC = co31*Cam/kspc(1);

%  T=253.6e6-T;
%%% Data from Payne
pd =  csvread('PayneData.csv');
stel = pd(:,1); % stratigraphic elevation in meters
d44ca = pd(:,2); 
stdp = pd(:,3); % standard deviation of calcium data
depth = [0 160]; % used for calculating time
for i=1:length(stel)
time(i) = interp1(depth, [0 1] , stel(i));
end

% Plots
lw = 4; %2 specifies line width
ms = 20;%10 marker size
fs = 18;%12 font size
figurehndl=1;    
FigHandle = figure(figurehndl);
  set(FigHandle, 'Position', [0, 0, 600, 800]);
  
subplot(411)

epscv(1)=-1.4;
plot(T,Ca*1000,'-m','LineWidth',lw);
% plot(T,d44Ca,'g');
ylabel({'[Ca^{2+}]';' [mmol kg^{-1}]'},'fontsize',fs);
% title('numerical solution, \delta_w')
set(gca,'YLim',[8 13],'XTickLabel',[])
set(gca,'FontSize',fs)
text(1.00,12.55,' a)','FontSize',12,'FontWeight','bold');
% title('Simple, 1-box model, no acidification')

subplot(412)
plot(T,d44Ca+epscv','g','LineWidth',lw);
hold on
plot((time+1-0.27)*1e6,d44ca,'.','LineWidth',lw,'markersize',ms,'Color','r');
h=errorbar((time+1-0.27)*1e6,d44ca,stdp,'r'); 
set(h,'linestyle','none','LineWidth',0.5)
hold off
ylabel('\delta^{44} Ca_{sed} (‰)','fontsize',fs);
set(gca,'YLim',[-0.9 -0.3],'XTickLabel',[])

set(gca,'FontSize',fs)
if(flag == 0)
lgnd =legend('Acidification Scenario','\delta^{44}Ca - data','Location','SouthEast');
set(lgnd, 'fontsize',12)
end
if(flag == 1)
legend('Weathering Scenario','\delta^{44}Ca - data','Location','SouthEast')
end
if(flag == 2)
legend('Weathering + [CO_3^{2-} feedback]','\delta^{44}Ca - data','Location','SouthEast')
end
text(1.00,-0.350,' b)','FontSize',12,'FontWeight','bold');

% plot(T,d44Ca--[epsc epscv(2:end)]','g');
% % plot(T,d44Ca,'g');
% ylabel('\delta^{44}Ca_{sw}');
% % set(gca,'YLim',[0.6 1.2],'XTickLabel',[])
% text(1.00,1.150,' b)','FontSize',12,'FontWeight','bold');



Fcarbv(1) = kcarb/1e12;
Frivv1(1) = Friv0/1e12;
subplot(413)
plot(T,[kcarb/1e12 Fcarbv(2:end)/1e12],'g','LineWidth',lw);
set(gca,'YLim',[15 40],'XTickLabel',[])
ylabel({'Fcarb';' (10^{12} mol yr^{-1})'},'fontsize',fs);
text(1.00,36,' c)','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',fs)

subplot(414)
plot(T/1e6,[Friv0/1e12 Frivv(2:end)/1e12],'g','LineWidth',lw);
% xlabel('years');
%  set(gca,'XTickLabel',{'253.6','253.1','252.6','252.1','251.6','251.1','250.6','250.1','249.6'})
 set(gca,'XTickLabel',{'253.6','252.6','251.6','250.6','249.6'})
ylabel({'Friv ';'(10^{12} mol yr^{-1})'},'fontsize',fs);
xlabel('Time (Ma)')
text(0.00,46,' d)','FontSize',12,'FontWeight','bold');
xlhand = get(gca,'xlabel');
set(xlhand,'string','Time (Ma)','fontsize',fs)
set(gca,'FontSize',fs)
% return
figurehndl=2;    
FigHandle = figure(figurehndl);
  set(FigHandle, 'Position', [0, 0, 600, 1000]);
subplot(611)

plot(T,Ca*1000,'-m','LineWidth',lw);
ylabel({'[Ca^{2+}]';' (mmol kg^{-1})'});
% title('numerical solution, \delta_w')
set(gca,'YLim',[9.9 10.2],'XTickLabel',[])
text(1.00,10.15,' a)','FontSize',12,'FontWeight','bold');
% title('Simple, 1-box model, [CO_3^{2-}] included')
subplot(612)
plot(T,d44Ca+epsc,'g','LineWidth',lw);
hold on
plot((time+1-0.27)*1e6,d44ca,'.','LineWidth',lw,'markersize',10,'Color','r');
h=errorbar((time+1-0.27)*1e6,d44ca,stdp,'r'); 
set(h,'linestyle','none','LineWidth',1)
hold off
ylabel('\delta^{44} Ca_{sed} (â€°)');
set(gca,'YLim',[-0.9 -0.3],'XTickLabel',[])
if(flag == 0)
legend('Acidification Scenario','\delta^{44}Ca - data','Location','SouthEast')
end
if(flag == 1)
legend('Weathering Scenario','\delta^{44}Ca - data','Location','SouthEast')
end
if(flag == 2)
legend('Weathering Scenario + [CO_3^{2-}] feedback','\delta^{44}Ca - data','Location','SouthEast')
end

set(gca,'XTickLabel',[])
text(1.00,-0.350,' b)','FontSize',12,'FontWeight','bold');



Fcarbv(1) = kcarb/1e12;
Frivv1(1) = Friv0/1e12;
subplot(613)
plot(T,[kcarb/1e12 Fcarbv(2:end)/1e12],'g','LineWidth',lw);
set(gca,'XTickLabel',[])
ylabel({'Fcarb ';'(10^{12} mol yr^{-1})'});
text(1.00,55,' c)','FontSize',12,'FontWeight','bold');

subplot(614)
plot(T,[Friv0/1e12 Frivv(2:end)/1e12],'g','LineWidth',lw);
% xlabel('years');
set(gca,'XTickLabel',[],'YLim',[10 50])
ylabel({'Friv';' (10^{12} mol yr^{-1})'})
text(0.00,46,' d)','FontSize',12,'FontWeight','bold');

subplot(615)
plot(T,myco3,'g','LineWidth',lw);
% xlabel('years');
set(gca,'XTickLabel',[])
ylabel({'[CO_3^{2-}]';' [\mu mol kg^{-1}]'});
text(0.00,230,' e)','FontSize',12,'FontWeight','bold');


subplot(616)
plot(T/1e6,[omgC(1) omgCv(2:end)],'g','LineWidth',lw);
xlabel('Time (Ma)');
ylabel('\Omega');
set(gca,'XTickLabel',{'253.6','253.1','252.6','252.1','251.6','251.1','250.6','250.1','249.6'})
text(0.00,5.5,' f)','FontSize',12,'FontWeight','bold');

% subplot(716)
% plot(T,TCAR*1e3,'g');
% % xlabel('years');
% set(gca,'XTickLabel',[])
% ylabel('TC [mmol kg^{-1}]');
% 
% subplot(717)
% plot(T,TALK*1e3,'g');
% xlabel('years');
% ylabel('TA [mmol kg^{-1}]');

% figure(3)
% 
% plot(T,myco3,'g', T,co31*1e6, 'k');
% legend('my CO3=', 'csys CO3=');
% 
% 
% 
figure(3)
subplot(411)
plot(T/1e6,[Friv0/1e12 Frivv(2:end)/1e12],'g','LineWidth',lw);
% xlabel('years');
% set(gca,'XTickLabel',{'253.6','253.1','252.6','252.1','251.6','251.1','250.6','250.1','249.6'})
ylabel('Friv (10^{12} mol yr^{-1})');
set(gca,'XTickLabel',[])
text(0.00,27,' a)','FontSize',12,'FontWeight','bold');

subplot(412)
plot(T/1e6,d44Ca,'g','LineWidth',lw);
ylabel('\delta^{44}Ca_{sw}');
% xlabel('Time (Ma)')
% plot(T/1e6,myco3,'g','LineWidth',lw);
% ylabel('[CO_3^{2-}] [\mu mol kg^{-1}]');
set(gca,'XTickLabel',[])
text(0.00,210,' b)','FontSize',12,'FontWeight','bold');

subplot(413)
plot(T/1e6,epscv','g','LineWidth',lw);
ylabel('\epsilon_{carb}');
set(gca,'XTickLabel',[])
text(0.00,-1.20,' c)','FontSize',12,'FontWeight','bold');

subplot(414)

plot(T/1e6,d44Ca+epscv','g','LineWidth',lw);
hold on
plot((time+1),d44ca,'.','LineWidth',lw);
plot([1.23 1.23],[-1 -0.2],'-k');
hold off
legend('See figure caption','\delta^{44}Ca - data','extinction horizon','Location','SouthEast')
ylabel('\delta^{44}Ca_{carb}');
xlabel('Time (Ma)')
set(gca,'XTickLabel',{'253.9','253.4','252.9','252.4','251.9','251.4','250.9','250.4','249.9'})
text(0.00,-0.25,' d)','FontSize',12,'FontWeight','bold');


lw = 4; 
%%% For EGU presentation
figurehndl=5;    
fs = 18;
FigHandle = figure(figurehndl);
  set(FigHandle, 'Position', [0, 0, 600, 800]);
  
subplot(211)

epscv(1)=-1.4;
plot(T,Ca*1000,'-m','LineWidth',lw);
% plot(T,d44Ca,'g');
ylabel('[Ca^{2+}] [mmol kg^{-1}]','fontsize',fs);
% title('numerical solution, \delta_w')
set(gca,'YLim',[9 13],'XTickLabel',[])
text(1.00,12.55,' a)','FontSize',12,'FontWeight','bold');
% title('Simple, 1-box model, no acidification')
set(gca,'FontSize',fs)

subplot(212)
p1=plot(T,d44Ca+epscv','g','LineWidth',lw);
hold on
p2=plot((time+1-0.27)*1e6,d44ca,'.','LineWidth',lw,'markersize',20,'Color','r');
h=errorbar((time+1-0.27)*1e6,d44ca,stdp,'r'); 
set(h,'linestyle','none','LineWidth',0.5)
hold off
ylabel('\delta^{44} Ca_{sed} (‰)','fontsize',fs);
set(gca,'YLim',[-0.9 -0.3])
xlhand = get(gca,'xlabel');
set(xlhand,'string','Time (Ma)','fontsize',fs)
set(gca,'FontSize',fs)
% set(gca,'XTickLabel',{'253.9','253.4','252.9','252.4','251.9','251.4','25
% 0.9','250.4','249.9'})
set(gca,'XTickLabel',{'253.6','252.6','251.6','250.6','249.6'})
if(flag == 0)
legend('Acidification Scenario','\delta^{44}Ca - data','Location','SouthEast')
end
if(flag == 1)
legend([p2],{'\delta^{44}Ca - data'},'Location','SouthEast')
end
if(flag == 2)
legend([p2],'\delta^{44}Ca - data','Location','SouthEast')
end
text(1.00,-0.350,' b)','FontSize',12,'FontWeight','bold');















