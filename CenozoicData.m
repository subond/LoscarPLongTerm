% Cenozoic data

c13kurtz=load('dat/LPEEkurtz/Sim2/c13kurtz2.DAT');
dbc=c13kurtz(2,:);
%%########## DATA ###############################
%%% Zachos 2008 benthic d13C record

zachos=csvread('dat\Cenozoicd13c\zac08.csv');
tzachos=zachos(:,1); %time zachos
c13zac=zachos(:,2); %time zachos
o18zac=zachos(:,3); %time zachos
rz=ksrlin(tzachos,c13zac,1); %Local linear kernel smoothing regression zachos
rzo=ksrlin(tzachos,o18zac,1); %Local linear kernel smoothing regression zachos

% Linear regression
p = polyfit(tzachos,c13zac,1);
slope=p(1);
intercept=p(2);
linreg=slope*tzachos+intercept;



d18od=inpaint_nans(o18zac);
d18om=runmean(d18od,100);
dosw=zeros(size(d18om));
doswa=zeros(size(d18om));
dosw(find(tzachos<10.07))=-0.2;
dosw(find(tzachos>=10.07))=-1.0;
dosw(find(tzachos>34.04))=-1.2;

doswa(find(tzachos>34.04))=12;
doswa(find(tzachos<=34.04))=-4.25;
k = zeros(size(d18om));
oTemp1= zeros(size(d18om));
k(find(tzachos>34.04))=-4;
k(find(tzachos<=34.04))=-2;

oTemp = 16.9 - 4.8*(d18om-dosw);
oTemp1(find(tzachos>34.04)) = k(find(tzachos>34.04)).*d18om(find(tzachos>34.04))+doswa(find(tzachos>34.04));
oTemp1(find(tzachos<=34.04)) = k(find(tzachos<=34.04)).*(d18om(find(tzachos<=34.04))+doswa(find(tzachos<=34.04)));


% % Oxygen isotope
% figure
% box on
% plot(tzachos,o18zac,'c+',rzo.x,rzo.f,'k-',tzachos,d18om ,'r');
% legend('Zachos ''08','Local linear regression', '100 point run mean')
% ylabel('\delta^{18}O (‰)')
% set(gca,'xlim',[0 60])

% Ice free temp
figure
subplot(411)
box on
plot(tzachos,oTemp1);
ylabel('Ice-free temperature [^oC]')
text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'xlim',[0 60],'XTickLabel',[])

% figure
% box on
% plot(tzachos,c13zac,'c+',rz.x,rz.f,'k-');
% legend('Zachos ''08','Local linear regression')
% ylabel('\delta^{13}C (‰)')
% set(gca,'xlim',[0 60])

CO=zachos;
XX=csvread('dat\Cenozoicd13c\ds01.csv');
d13cd=CO(:,2);
d13cd=inpaint_nans(d13cd);
Y=runmean(d13cd,50);

% figure
subplot(412)
box on
hold on
plot(rz.x,rz.f,'k-','LineWidth',2);
% plot(CO(:,1)+0.53,Y) %Mean deep zachos
plot(XX(:,1),XX(:,2),'r','LineWidth',2)
plot(c13kurtz(1,:),dbc,'g','LineWidth',2)
hold off

text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'xlim',[0 60],'XTickLabel',[])
legend('Benthic','Surface','Bulk')
ylabel('\delta^{13}C (‰)')
% xlabel('Age (Ma)')
% set(gca,'XDir','reverse')


%###################### CCD Palike and Van Andel data ######################
palikeCCD=csvread('dat\Cenozoicd13c\CenozoicCCD\palike12ccd.csv');
% Van Andel North Atlantic CCD
VANA=csvread('dat\Cenozoicd13c\CenozoicCCD\VanAndelNAtlantic.csv');
VASA=csvread('dat\Cenozoicd13c\CenozoicCCD\VanAndelSAntlantic.csv');
VAI=csvread('dat\Cenozoicd13c\CenozoicCCD\VanAndelndian.csv');
VAP=csvread('dat\Cenozoicd13c\CenozoicCCD\VanAndelPacific.csv');

timeP = palikeCCD(:,1)/1000; % million years
timeVNA = VANA(:,1);
timeVSA = VASA(:,1);
timeVI = VAI(:,1);
timeVP = VAP(:,1);

%CCDs meters 
PccdP = palikeCCD(:,2);
VAccdNA = VANA(:,2)*1000;
VAccdSA = VASA(:,2)*1000;
VAccdI = VAI(:,2)*1000;
VAccdP = VAP(:,2)*1000;


% figure
subplot(414)
box on
hold on
plot(timeP, PccdP,'r','LineWidth',lw)
plot(timeVP, VAccdP,'b','LineWidth',lw)
plot(timeVNA, VAccdNA,'b--','LineWidth',lw)
plot(timeVSA, VAccdSA,'b-.','LineWidth',lw)
plot(timeVI, VAccdI,'b:','LineWidth',lw)
hold off

% legend('Pacific (Palike et al. ''12)','Pacific (Van Andel ''75)','North
% Atlantic (Van Andel ''75)','South Atlantic (Van Andel ''75)','Indic (Van Andel ''75)')
legend('Pacific','Pacific','North Atlantic','South Atlantic','Indic')
ylabel('CCD (m)')
xlabel('Age (Ma)')
text(0.02,0.98,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'YDir','reverse')
set(gca,'xlim',[0 60])



% ####################### pCO2 Berling & Royer 2011 compilation ###########

BCa=csvread('dat\Cenozoicd13c\CenozoicpCO2\BCa.csv');
Boron=csvread('dat\Cenozoicd13c\CenozoicpCO2\Boron.csv');
L=csvread('dat\Cenozoicd13c\CenozoicpCO2\Liverworts.csv');
N=csvread('dat\Cenozoicd13c\CenozoicpCO2\Nacholite.csv');
Pal=csvread('dat\Cenozoicd13c\CenozoicpCO2\Paleosols.csv');
Phy=csvread('dat\Cenozoicd13c\CenozoicpCO2\Phytoplankton.csv');
St=csvread('dat\Cenozoicd13c\CenozoicpCO2\Stomata.csv');

% time data
timeB  = (BCa(:,1));
timeBo = (Boron(:,1));
timeL = (L(:,1));
timeN = (N(:,1));
timePal = (Pal(:,1));
timePhy = (Phy(:,1));
timeSt = (St(:,1));

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


% figure
% hold on
% errorbar(timeB,co2B,co2Bn,co2Bp,'ro')
% errorbar(timeBo,co2Bo,co2Bon,co2Bop,'bd')
% errorbar(timeL,co2L,co2Ln,co2Lp,'y+')
% errorbar(timeN,co2N,co2Nn,co2Np,'b*')
% errorbar(timePal,co2Pal,co2Paln,co2Palp,'m.')
% errorbar(timePhy,co2Phy,co2Phyn,co2Phyp,'kx')
% errorbar(timeSt,co2St,co2Stn,co2Stp,'g^')
% hold off

% figure
subplot(413)
box on
hold on
plot(timeB,co2B,'ro')
plot(timeBo,co2Bo,'bd')
plot(timeL,co2L,'r+')
plot(timeN,co2N,'b*')
plot(timePal,co2Pal,'m.')
plot(timePhy,co2Phy,'kx')
plot(timeSt,co2St,'g^')
hold off

legend('B/Ca','Boron', 'Liverworts','Nacholite','Paleosols','Phytoplankton','Stomata')
ylabel('pCO2')
text(0.02,0.98,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
set(gca,'xlim',[0 60],'XTickLabel',[])

