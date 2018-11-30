clear all
% close all
global fbetta gctl
fbetta = 0.3987;

%% Chapter3 results and comparisons
%%% Not smoothed kurtz data
kurtz=load('dat/Cenozoicd13c/kurtz03/d13c.csv');
t=kurtz(:,1);
d13c=kurtz(:,2);

%%% smoothed kurtz data
c13kurtz=load('dat/LPEEkurtz/Sim2/c13kurtz2.DAT');
dbc=c13kurtz(2,:);

%%% Zachos Benthic 08
zachos=csvread('dat/Cenozoicd13c/zac08.csv');
tzachos=zachos(:,1); %time zachos
c13zac=zachos(:,2); %c13 zachos
o18zac=zachos(:,3); %o18 zachos
rz=ksrlin(tzachos,c13zac,1); %Local linear kernel smoothing regression zachos
rzo=ksrlin(tzachos,o18zac,1); %Local linear kernel smoothing regression zachos


%%% pH Data Foster et al. 2012, Raitzsch 2013, Pearson 2009
foster=csvread('dat/Cenozoicd13c/CenozoicpH/FosterpH.csv');
% Raitzsch is deep ocean pH
raitP=csvread('dat/Cenozoicd13c/CenozoicpH/RaitzschpHPacific.csv');
raitA=csvread('dat/Cenozoicd13c/CenozoicpH/RaitzschpHAtlantic.csv');
pear=csvread('dat/Cenozoicd13c/CenozoicpH/Pearson09.csv');
stap=csvread('dat/Cenozoicd13c/CenozoicpH/Stap16.csv');
seki=csvread('dat/Cenozoicd13c/CenozoicpH/Seki10.csv');

tpHf=foster(:,1);
pHf = foster(:,2);

tpHrP = raitP(:,1);
pHrP = raitP(:,2);
tpHrA = raitA(:,1);
pHrA = raitA(:,2);

tpHp = pear(:,1);
pHp = pear(:,2);

tpHs = stap(:,1);
pHs = stap(:,2);

tpHsk=seki(:,1);
pHsk=seki(:,2);

%###################### CCD Palike and Van Andel data ######################
palikeCCD=csvread('dat/Cenozoicd13c/CenozoicCCD/palike12ccd.csv');
palikeCCDoff=csvread('dat/Cenozoicd13c/CenozoicCCD/palike12ccdoffeq.csv');
% Van Andel North Atlantic CCD
VANA=csvread('dat/Cenozoicd13c/CenozoicCCD/VanAndelNAtlantic.csv');
VASA=csvread('dat/Cenozoicd13c/CenozoicCCD/VanAndelSAntlantic.csv');
VAI=csvread('dat/Cenozoicd13c/CenozoicCCD/VanAndelndian.csv');
VAP=csvread('dat/Cenozoicd13c/CenozoicCCD/VanAndelPacific.csv');
% Slotnick et al., 2015 Indian ocean data excluding hyperthermals
SL=csvread('dat/Cenozoicd13c/CenozoicCCD/ccdISlotNoH.csv');

timeP = palikeCCD(:,1); % million years
timePo = palikeCCDoff(:,1); % million years
timeVNA = VANA(:,1);
timeVSA = VASA(:,1);
timeVI = VAI(:,1);
timeVP = VAP(:,1);
timeSL = SL(:,1);

%CCDs meters 
PccdP = palikeCCD(:,2);
PccdPo=palikeCCDoff(:,2);
VAccdNA = VANA(:,2)*1000;
VAccdSA = VASA(:,2)*1000;
VAccdI = VAI(:,2)*1000;
VAccdP = VAP(:,2)*1000;
SLccdI = SL(:,2);

%Campbell 2018 CCD range pacific
[timeC, CCDrp]=campbell18CCD();
% ####################### pCO2 Berling & Royer 2011 compilation ###########

BCa=csvread('dat/Cenozoicd13c/CenozoicpCO2/BCa.csv');
Boron=csvread('dat/Cenozoicd13c/CenozoicpCO2/Boron.csv');
L=csvread('dat/Cenozoicd13c/CenozoicpCO2/Liverworts.csv');
N=csvread('dat/Cenozoicd13c/CenozoicpCO2/Nacholite.csv');
Pal=csvread('dat/Cenozoicd13c/CenozoicpCO2/Paleosols.csv');
Phy=csvread('dat/Cenozoicd13c/CenozoicpCO2/Phytoplankton.csv');
St=csvread('dat/Cenozoicd13c/CenozoicpCO2/Stomata.csv');
Alk=csvread('dat/Cenozoicd13c/CenozoicpCO2/Alkenone.csv');

% time data
timeB  = (BCa(:,1));
timeBo = (Boron(:,1));
timeL = (L(:,1));
timeN = (N(:,1));
timePal = (Pal(:,1));
timePhy = (Phy(:,1));
timeSt = (St(:,1));
timeAlk = (Alk(:,1));

%time error bars
timeBn = (BCa(:,2));
timeBp = (BCa(:,3));

% pco2 data
co2B  = (BCa(:,4));
co2Bo  = (Boron(:,4));
co2L= (L(:,4));
co2N= (N(:,4));
co2Pal= (Pal(:,4));
co2Phy= (Phy(:,4));
co2St= (St(:,4));
co2Alk= (Alk(:,4));

% pco2 error bars
co2Bn = (BCa(:,5));
co2Bp = (BCa(:,6));
co2Bon = (Boron(:,5));
co2Bop = (Boron(:,6));

co2Ln = (L(:,5));
co2Lp = (L(:,6));

co2Nn = (N(:,5));
co2Np = (N(:,6));

co2Paln = (Pal(:,5));
co2Palp = (Pal(:,6));

co2Phyn = (Phy(:,5));
co2Phyp = (Phy(:,6));

co2Stn = (St(:,5));
co2Stp = (St(:,6));

co2Alkn = (Alk(:,5));
co2SAlkp = (Alk(:,6));

% new pco2&pH data from Richard

X = load('dat/paleopH.dat');
tvec = X(:,1); lt = length(tvec);
[An16,Ed15,Pe09,F761,F926,Ba11,Se10] = fpHAuth(X);


gct = 1;
gctl=59; %last year in Ma
crit = 0.0001; % Newton-Raphson accuracy


dir = 'dat/LoscarLT/LOSCAR/105/';

TCvt =  (importdata([dir 'TCvt.DAT']));
V =  (importdata([dir 'V.DAT']));
FiN =  (importdata([dir 'FiN.DAT']));
FSi =  (importdata([dir 'FSi.DAT']));    
oIpv=  (importdata([dir 'oIpv.DAT']));
oIv=  (importdata([dir 'oIv.DAT']));
PPLvv=  (importdata([dir 'PPLvv.DAT']));
EPLvv=  (importdata([dir 'EPLvv.DAT'])); 
Fcapv=  (importdata([dir 'Fcapv.DAT']));
Ffepv=  (importdata([dir 'Ffepv.DAT']));
EPH =  (importdata([dir 'EPH.DAT']));
PPH =  (importdata([dir 'PPH.DAT']));
Fpwv=  (load([dir 'Fpwv.DAT']));

dox=  (importdata([dir 'dox.DAT']));
p=  (importdata([dir 'p.DAT']));
a =  (importdata([dir 'a.DAT']));
c =  (importdata([dir 'c.DAT']));

% ccdA=  (load([dir 'ccdA.DAT']));
% ccdI=  (load([dir 'ccdI.DAT']));
% ccdP=  (load([dir 'ccdP.DAT']));
% % ccdT=zeros(59,1);
% ccdT(1:24)=NaN;
% ccdT(25:gctl)=  (load([dir 'ccdT.DAT']));
pco2t=  (load([dir 'pco2t.DAT']));
phtv=  (importdata([dir 'phtv.DAT']));
co3tv=  (importdata([dir 'co3tv.DAT']));

d13cL=  (importdata([dir 'd13c.DAT']));
d13fcA =  (importdata([dir 'd13fcA.DAT']));
d13fcI =  (importdata([dir 'd13fcI.DAT']));
d13fcP =  (importdata([dir 'd13fcP.DAT']));
d13fcT =  (importdata([dir 'd13fcT.DAT']));
epsp=  (importdata([dir 'epspv.DAT']));

eI=0.78;
eIv=  (importdata([dir 'eIv.DAT']));
eIpv=  (importdata([dir 'eIpv.DAT']));
orgCb = ((1-(eIv+oIv)).*(EPLvv+EPH));
orgPb=(1-(eIpv+oIpv)).*(PPLvv+PPH);

meanDoxDeep = (dox(:,7).*V(:,7)+dox(:,8).*V(:,8)+dox(:,9).*V(:,9))./((V(:,7)+V(:,8)+V(:,9)));

dsv   = [.1 .6 1 1.5 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6.5]*1000;
fcA=  (load([dir 'fcA.DAT']));
fcI=  (load([dir 'fcI.DAT']));
fcP=  (load([dir 'fcP.DAT']));
fcT(1:gctl,1:13)=NaN;
fcT(1:24,1:13)=NaN;
fcT(25:gctl,:) = (load([dir 'fcT.DAT']));

[ccdA, ccdI, ccdP] = removeCCDspikes(fcA,fcI,fcP);

% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'CaCO3OverTime.gif';
% for n = 1:gctl
%         plot(fcA(n,:)*1e2,dsv,'b-d');
%         hold on;
%         plot(fcI(n,:)*1e2,dsv,'m-s');
%         plot(fcP(n,:)*1e2,dsv,'k-o')
%         plot(fcT(n,:)*1e2,dsv,'g-p');
%         drawnow
%         hold off;
%         set(gca,'YDir','reverse');
%         set(gca,'FontSize',10);
%         set(gca,'YDir','reverse');
%         xlabel('CaCO_3 (wt %)');
%         ylabel('Depth (m)');
%         %axis([0 100 0 5.5]);
%          legend('Atlantic','Indic','Pacific','Tethys');
%          title(['Year: ' num2str(n-1) ' Ma'] ,'FontSize',16)
%          % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
% end
for gct = 1:1:gctl
    if(gct>1)
%         myFbg(gct) = orgCb(end)/1e12;
        d13G(gct)=myNR(3,crit, gct, orgCb(gct)/1e12);
    else
%         myFbg(gct) =5.0;
        d13G(gct) = 0.2470;
    end
    [fbgv(gct),fbcv(gct),fwgv(gct),fmcv(gct),fmgv(gct),fwcv(gct),fGGi(gct),fgkc(gct),frkc(gct),fekc(gct),fdkc(gct),flakc(gct),rco2(gct),pco2gca(gct), acv(gct), G(gct), dcv(gct), fmv(gct), famv(gct)]=gcfun13(d13G(gct), gct);
end


fs=12;
lw=2;
cs = 'gggkkkrrrbgkr';
sstr = '- ---.- ---.- ---. -: : : ';
time = [0:length(p)-1];
lstr0 = 'LALILPIAIIIPDADIDP HLTITDT';
lstrR = {'Anagnostou+ ''16';'Edgar ''15?';'Pearson+ ''09';'Foster+ ''12';...
'Foster+ ''12';'Bartoli+ ''11';'Seki+ ''10'};
for k=1:13
    lstr(k,:) = sprintf('%s',lstr0(2*k-1:2*k));
end;

FigHandle = figure(102);
  set(FigHandle, 'Position', [200, 0, 800, 900]);
% figure
subplot(321)
box on
hold on
% plot(timeB,co2B,'ro')
% plot(timeBo,co2Bo,'bd')
% plot(timeL,co2L,'r+')
% plot(timeN,co2N,'b*')
% plot(timePal,co2Pal,'m.')
% plot(timePhy,co2Phy,'kx')
plot(timeSt,co2St,'g^')
% plot(timeAlk,co2Alk,'bx')


plot(An16.Age,An16.CO2,'s','MarkerEdgeColor','k','MarkerFaceColor','r')
plot(Ed15.Age,Ed15.CO2,'s','MarkerEdgeColor','k','MarkerFaceColor','w')
plot(Pe09.Age,Pe09.CO2,'^','MarkerEdgeColor','k','MarkerFaceColor','g')
plot(F761.Age,F761.CO2,'^','MarkerEdgeColor','k','MarkerFaceColor','b')
plot(F926.Age,F926.CO2,'s','MarkerEdgeColor','w','MarkerFaceColor','b')
plot(Ba11.Age,Ba11.CO2,'o','MarkerEdgeColor','k','MarkerFaceColor','b')
plot(Se10.Age,Se10.CO2,'d','MarkerEdgeColor','w','MarkerFaceColor','b')
plot(time,pco2t,'r','LineWidth',lw)
hold off
set(gca,'XLim',[0 60],'XTickLabel',[])
% legend('B/Ca','Boron', 'Liverworts','Nacholite','Paleosols','Phytoplankton','Stomata','Model')

H = legend(['Stomata';lstrR;'Model'],'Location','NorthEast');
ylabel('pCO_2 (ppmv)')
text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')




subplot(322)
box  on;
hold on
plot(t,d13c,'o')
% plot(time, d13G,'m','LineWidth',lw)
plot(time, d13fcA,'g','LineWidth',lw)
% plot(tzachos,c13zac,'+k')
hold off
text(0.02,0.9,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
% set(gca,'XLim',axx);
% xlabel('Time (Ma)');
ylabel('Bulk carbonate \delta^{13}C (‰)')
legend('Kurtz et al. ''03', 'Model','Location','SouthEast')
set(gca,'XLim',[0 60],'XTickLabel',[])

subplot (323)
tRH=[tpHrP;tpHrA];
phRH = [pHrP;pHrA];
box  on;
hold on;
for k=1:13
    h1=plot(time  ,phtv (:,k),sstr(2*k-1:2*k),'Color',cs(k));
%     key1{k}=lstr(k,1:2);
end;
% h2=plot(tpHf,pHf,'go');
% h3=plot(tpHp,pHp,'gd');
% h4=plot(tpHrP,pHrP,'r*');
% h5=plot(tpHrA,pHrA,'r*');
% h6=plot(tpHs,pHs,'ms');
% h7=plot(tpHsk,pHsk,'k.');

h2=plot(An16.Age,An16.pH,'s','MarkerEdgeColor','k','MarkerFaceColor','r');
h3=plot(Ed15.Age,Ed15.pH,'s','MarkerEdgeColor','k','MarkerFaceColor','w');
h4=plot(Pe09.Age,Pe09.pH,'^','MarkerEdgeColor','k','MarkerFaceColor','g');
h5=plot(F761.Age,F761.pH,'^','MarkerEdgeColor','k','MarkerFaceColor','b');
h6=plot(F926.Age,F926.pH,'s','MarkerEdgeColor','k','MarkerFaceColor','b');
h7=plot(Ba11.Age,Ba11.pH,'o','MarkerEdgeColor','k','MarkerFaceColor','r');
h8=plot(Se10.Age,Se10.pH,'d','MarkerEdgeColor','k','MarkerFaceColor','b');
h9=plot(tRH,phRH,'r*');
hold off;
text(0.02,0.98,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
% set(gca,'FontSize',fs);
set(gca,'XTickLabel',[])
% xlabel('Time (y)');
% legend([h2 h3 h4 h6 h7],{'Foster et al., ''12','Pearson et al., ''09','Raitzsch & Hoenisch, ''13','Stap et al., 16','Seki et al., ''10'},'Location','SouthWest');
legend([h2 h3 h4 h5 h6 h7 h8 h9],[lstrR;'Raitzsch & Hoenisch, ''13'],'Location','SouthWest');

% H1 = legend([lstrR;'Raitzsch & Hoenisch, ''13'],"Location","NorthEast");% gridLegend(hdlY,5,key1,'location','north','Fontsize',8,'Box','off');
ylabel('pH');

subplot (324)
box on
hold on;
for k=1:13
    hdlY(k)=plot(time  ,d13cL (:,k),sstr(2*k-1:2*k),'Color',cs(k));
     key1{k}=lstr(k,1:2);
end;
% plot(tzachos,c13zac,'+k')
hold off;
% set(gca,'FontSize',fs);
text(0.02,0.98,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
% set(gca,'YDir','reverse');
% xlabel('Time (Ma)');
ylabel('\delta^{13}C (‰)');
set(gca,'XTickLabel',[])
gridLegend(hdlY,5,key1,'location','north','Fontsize',8,'Box','off');
% Hl=legend(lstr,1);
% set(Hl,'FontSize',10);

subplot(325)
box on
hold on
fill(timeC,CCDrp,'y')
plot(timeP, PccdP,'r','LineWidth',lw)
plot(timePo, PccdPo,'r--','LineWidth',lw)
% plot(timeSL, SLccdI,'g-','LineWidth',lw)
plot(timeVP, VAccdP,'b','LineWidth',lw)
% plot(timeVNA, VAccdNA,'b--','LineWidth',lw)
% plot(timeVSA, VAccdSA,'b-.','LineWidth',lw)
% plot(timeVI, VAccdI,'b:','LineWidth',lw)


% plot(time, ccdI,'k--','LineWidth',lw)
plot(time, ccdP,'k.-','LineWidth',lw)
% plot(time, ccdA,'k','LineWidth',lw)
% plot(time, ccdT,'k:','LineWidth',lw)
hold off

% legend('Pacific (Palike et al. ''12)','Pacific (Van Andel ''75)','North
% Atlantic (Van Andel ''75)','South Atlantic (Van Andel ''75)','Indic (Van Andel ''75)')
% legend('Pacific','Pacific','North Atlantic','South Atlantic','Indic')
legend('Campbel+ ''18','Pälike+ ''12 Eq','Pälike+ ''12 oE' ,'Van Andel ''75','Model','Location','SouthEast')
ylabel('Pacific CCD (m)')
xlabel('Year (Ma)')
text(0.02,0.98,'e)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'YDir','reverse')
set(gca,'xlim',[0 60])


subplot(326)
box on
hold on

plot(timeVNA, VAccdNA,'b--','LineWidth',lw)
plot(timeVSA, VAccdSA,'b-.','LineWidth',lw)
plot(time, ccdA,'k','LineWidth',lw)
hold off

% legend('Pacific (Palike et al. ''12)','Pacific (Van Andel ''75)','North
% Atlantic (Van Andel ''75)','South Atlantic (Van Andel ''75)','Indic (Van Andel ''75)')
% legend('Pacific','Pacific','North Atlantic','South Atlantic','Indic')
% legend('North Atlantic Van Andel ''75','South Atlantic Van Andel ''75','Atlantic model')
legend('Van Andel ''75 NA','Van Andel ''75 SA','Model')
ylabel('Atlantic CCD (m)')
xlabel('Year (Ma)')
text(0.02,0.98,'f)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'YDir','reverse')
set(gca,'xlim',[0 60])

% figure
% % subplot(414)
% box on
% hold on
% plot(timeVI, VAccdI,'b:','LineWidth',lw)
% plot(time, ccdI,'k--','LineWidth',lw)
% 
% hold off
% 
% % legend('Pacific (Palike et al. ''12)','Pacific (Van Andel ''75)','North
% % Atlantic (Van Andel ''75)','South Atlantic (Van Andel ''75)','Indic (Van Andel ''75)')
% % legend('Pacific','Pacific','North Atlantic','South Atlantic','Indic')
% legend('Indian Van Andel ''75','Indian model')
% ylabel('CCD (m)')
% xlabel('Age (Ma)')
% text(0.02,0.98,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
% set(gca,'YDir','reverse')
% set(gca,'xlim',[0 60])
% 
% figure
% box on
% hold on
% plot(time, flakc,'b-','LineWidth',lw)
% plot(time, fdkc,'g.','LineWidth',lw)
% plot(time, fekc,'r:','LineWidth',lw)
% plot(time, frkc,'k--','LineWidth',lw)
% hold off
% legend('f_{LA}', 'f_{D}','f_E','f_R')
% xlabel('Age (Ma)')
% ylabel('Dimensionless GEOCARB parameters')





FigHandle = figure(202);
  set(FigHandle, 'Position', [1000, 0, 800, 900]);
subplot(321)
box on
hold on
% for k=1:13
%     plot(time  ,TCvt (:,k)   ,sstr(2*k-1:2*k),'Color',cs(k));
% end;
plot(time  ,TCvt (:,9)   ,'r','LineWidth',lw);
hold off
set(gca,'XTickLabel',[],'YLim',[0 15])
text(0.02,0.9,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
ylabel('Deep ocean temperature (^oC)')
% Hl=legend(lstr,1);
% set(Hl,'FontSize',10);
% legend('Deep Ocean')

subplot (322)
box on
hold on
% Org. C burial

plot(time, FiN/1e12,'b','LineWidth',lw);
plot(time, FSi/1e12,'r','LineWidth',lw);
plot(time, Fpwv/1e10,'g','LineWidth',lw);
hold off
% set(gca,'FontSize',fs);
set(gca,'XTickLabel',[],'YLim',[0 25])
%     Hl=legend('F_{OCB}');
%     set(Hl,'FontSize',10);
text(0.02,0.9,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
legend('F_{wc}','F_{Si}','F_{pw}')
% xlabel('Time (Ma)');
% ylabel({'F_{wc}, F_{Si} [10^{12} mol yr^{-1}]';'F_{pw} [10^{10} mol yr^{-1}]'});    
ylabel({'F_{wc}, F_{Si}, F_{pw} (see caption)'});    
    
% figure
subplot (323)
box on
hold on
% Org. C burial

[hAx,hLine1,hLine2] =plotyy(time, orgCb./1e12,time,meanDoxDeep);
% plot(time(1:gctl), fbgv1(1:gctl));
hold off
% set(gca,'FontSize',fs);
text(0.02,0.9,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')

%     Hl=legend('F_{OCB}');
%     set(Hl,'FontSize',10);
% xlabel('Time (Ma)');
% ylabel('Org. C burial [10^{12}mol y^{-1}]');
set(hLine1,'LineWidth',lw,'Color','k');
set(hLine2,'LineWidth',lw,'Color','r','LineStyle','--');
set(hAx,{'ycolor'},{'k';'r'})  % Left color red, right color blue..
set(hAx(1),'YLim',[3 8],'box','off','YTick',[3 4 5 6 7 8],'XTickLabel',[])
set(hAx(2),'XTickLabel',[])
set(gca,'XTickLabel',[])
ylabel(hAx(1),'Organic C burial (10^{12} mol yr^{-1})')
ylabel(hAx(2),'Mean deep ocean (O_2) (mol m^{-3})')


subplot(324)
plot(time, (EPLvv+EPH)/1e12,'g','LineWidth',lw)
text(0.02,0.9,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
% xlabel('Time (Ma)');
set(gca,'XTickLabel',[],'YAxisLocation', 'right','YLim',[200 1200])
ylabel('Organic C export (10^{12} mol yr^{-1})')

% PO4    
% figure
subplot (325)
box  on;
hold on
counter=1;
for k=1:13
    hdlY(k)=plot(time  ,p (:,k)*1e3,sstr(2*k-1:2*k),'Color',cs(k));
%     key{k}=sprintf('trace %d',k);
 key{k}=lstr(k,1:2);
end;
hold off
% set(gca,'FontSize',fs);
text(0.02,0.9,'e)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
% set(gca,'XLim',axx);
xlabel('Time (Ma)');
set(gca,'YLim',[0 8])
ylabel('PO_4 (\mumol kg^{-1})');
% Hl=legend(lstr,4);

gridLegend(hdlY,5,key,'location','northoutside','Fontsize',8,'Box','off');
% set(hdlY,'FontSize',10,'Location','Best');

% P bur
% figure
subplot (326)
box on
hold on
% Org. P burial
plot(time, orgPb./1e10,'g');

plot(time, Fcapv./1e10,'--b');
plot(time, Ffepv./1e10,'-.r');
totPbur = orgPb+Fcapv+Ffepv;
plot(time, totPbur./1e10,'k');
hold off
% set(gca,'FontSize',fs);
text(0.02,0.9,'f)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'YLim',[0 10])
Hl=legend('F_{pb}','F_{CaP}','F_{FeP}','Total P burial');
set(Hl,'FontSize',10);
xlabel('Time (Ma)');
ylabel('P burial fluxes (10^{10}mol y^{-1})');


% figure
% for k=1:13
%     hold on
%     h1=plot(time  ,TCvt (:,k),sstr(2*k-1:2*k),'Color',cs(k));
%     hold off
% %     key1{k}=lstr(k,1:2);
% end;
% 
clear fcA fcI fcP fcT ccdA ccdI ccdP ccdT;
for i=1:1:5

    if(i==1)
        dir = 'dat/LoscarLT/LOSCAR/105/';
    elseif(i==2)
        dir = 'dat/LoscarLT/LOSCAR/104/';
    elseif(i==3)
        dir = 'dat/LoscarLT/LOSCAR/103/';
    elseif(i==4)
        dir = 'dat/LoscarLT/LOSCAR/102/';
    else
        dir = 'dat/LoscarLT/LOSCAR/101/';
    end
    
    fcA {i,:}=  (load([dir 'fcA.DAT']));
    fcI {i,:}=  (load([dir 'fcI.DAT']));
    fcP {i,:}=  (load([dir 'fcP.DAT']));
    fcT{i,:}(1:gctl,1:13)=NaN;
    fcT{i,:}(1:24,1:13)=NaN;
    fcT{i,:}(25:gctl,:) = (load([dir 'fcT.DAT']));
% 
    [Accd, Iccd, Pccd] = removeCCDspikes(fcA{i,:},fcI{i,:},fcP{i,:});
    ccdA{i,:} = Accd;
    ccdI{i,:} = Iccd;
    ccdP{i,:} = Pccd;
    
    mpco2 {i,:}=  (load([dir 'pco2t.DAT']));
    
    my_time{i} = [0:length(Accd)-1];
end
% Campbell is a fill curve so break into upper and lower
% component
[t0, d0] = interpCCDs(timeC(1:length(timeC)/2), CCDrp(1:length(timeC)/2));
[t0b, d0b] = interpCCDs(timeC(length(timeC)/2:length(timeC)), CCDrp(length(timeC)/2:length(timeC)));
[t1, d1] = interpCCDs(timeP, PccdP);
[t2, d2] = interpCCDs(timePo, PccdPo);
[t3, d3] = interpCCDs(timeVP, VAccdP);
d3(1)=VAccdP(end);

ccd_min = min([d0;d0b;d1;d2;d3]);
ccd_max = max([d0;d0b;d1;d2;d3]);

%remove all NaN's from both X and Y
ccd_min = ccd_min(~isnan(ccd_min));
ccd_max = ccd_max(~isnan(ccd_max));

figure
x2 = [t3(1:length(ccd_min)), fliplr(t3(1:length(ccd_min)))];
inBetween = [ccd_min, fliplr(ccd_max)];


figure
hold on
fill(x2,inBetween,[211/255 211/255 211/255],'FaceColor',[211/255 211/255 211/255])
plot(t0,d0,'y')
plot(t0b,d0b,'k')
plot(t1,d1,'r')
plot(t2,d2,'r--')
plot(t3,d3,'b')

% plot(t3,ccd_min,'m')
% plot(t3,ccd_max,'m')
hold off

figure
subplot (211)
box on
hold on
plot(timeSt,co2St,'g^','HandleVisibility','off')
plot(An16.Age,An16.CO2,'s','MarkerEdgeColor','k','MarkerFaceColor','r','HandleVisibility','off')
plot(Ed15.Age,Ed15.CO2,'s','MarkerEdgeColor','k','MarkerFaceColor','w','HandleVisibility','off')
plot(Pe09.Age,Pe09.CO2,'^','MarkerEdgeColor','k','MarkerFaceColor','g','HandleVisibility','off')
plot(F761.Age,F761.CO2,'^','MarkerEdgeColor','k','MarkerFaceColor','b','HandleVisibility','off')
plot(F926.Age,F926.CO2,'s','MarkerEdgeColor','w','MarkerFaceColor','b','HandleVisibility','off')
plot(Ba11.Age,Ba11.CO2,'o','MarkerEdgeColor','k','MarkerFaceColor','b','HandleVisibility','off')
plot(Se10.Age,Se10.CO2,'d','MarkerEdgeColor','w','MarkerFaceColor','b','HandleVisibility','off')
plot(my_time{1},mpco2{1},'k--','LineWidth',lw)
plot(my_time{2},mpco2{2},'bx-','LineWidth',lw)
plot(my_time{3},mpco2{3},'cd-','LineWidth',lw)
plot(my_time{4},mpco2{4},'m.-','LineWidth',lw)
plot(my_time{5},mpco2{5},'k*-','LineWidth',lw)
hold off
legend('Q10 = 1.0, cons fsh','Q10 = 1.5, cons fsh','Q10 = 1.5, fsh(sea-level)','Q10 = 1.5, fsh(sea-level) + CaCO_3 prolif', 'preferred')
text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
ylabel('Pacific CCD (m)')
set(gca,'XLim',[0 60],'XTickLabel',[])
% legend('B/Ca','Boron', 'Liverworts','Nacholite','Paleosols','Phytoplankton','Stomata','Model')

% H = legend(['Stomata';lstrR;'Model'],'Location','NorthEast');
ylabel('pCO_2 (ppmv)')



subplot (212)
box on
hold on
fill(x2,inBetween,[211/255 211/255 211/255],'FaceColor',[211/255 211/255 211/255])
% fill(timeC,CCDrp,'y','HandleVisibility','off')
% plot(timeP, PccdP,'r','LineWidth',lw,'HandleVisibility','off')
% plot(timePo, PccdPo,'r--','LineWidth',lw,'HandleVisibility','off')
% plot(timeVP, VAccdP,'b','LineWidth',lw,'HandleVisibility','off')

% plot(time, ccdP,'k.-','LineWidth',lw)
plot(my_time{1},ccdP{1},'k--','LineWidth',lw)
plot(my_time{2},ccdP{2},'bx-','LineWidth',lw)
plot(my_time{3},ccdP{3},'cd-','LineWidth',lw)
plot(my_time{4},ccdP{4},'m.-','LineWidth',lw)
plot(my_time{5},ccdP{5},'k*-','LineWidth',lw)
% plot(time3,ccdP3,'k.','LineWidth',lw)
% plot(timeP, PccdP,'r','LineWidth',lw)
% plot(timePo, PccdPo,'r--','LineWidth',lw)
% plot(timeSL, SLccdI,'g-','LineWidth',lw)
hold off
legend('Data: Range','Q10 = 1.0, cons fsh','Q10 = 1.5, cons fsh','Q10 = 1.5, fsh(sea-level)','Q10 = 1.5, fsh(sea-level) + CaCO_3 prolif', 'preferred')
text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
ylabel('Pacific CCD (m)')
xlabel('Year (Ma)')
% text(0.02,0.98,'e)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'YDir','reverse')
set(gca,'xlim',[0 60])



% return
figure
box on
hold on
for k=1:13
    plot(time  ,co3tv (:,k)*1e6   ,sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off
Hl=legend(lstr);
set(Hl,'FontSize',10);
xlabel('Time (Ma)');
ylabel('(CO_3^{=})');

figure
box on
hold on
for k=1:13
    plot(time  ,dox (:,k)   ,sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off


return
 figure
subplot(211)
box  on;
hold on;
for k=1:13
    plot(time  ,c  (:,k),sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off;
set(gca,'FontSize',fs);
xlabel('Time (Ma)');
ylabel('TCO_2 (mmol kg^{-1})');
Hl=legend(lstr);
set(Hl,'FontSize',10);

% figure
subplot(212)
box  on;
hold on;
for k=1:13
    plot(time  ,a  (:,k),sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off;
set(gca,'FontSize',fs);
xlabel('Time (Ma)');
ylabel('TA (mmol kg^{-1})');
Hl=legend(lstr);
% set(Hl,'FontSize',10);


