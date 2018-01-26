dir = 'dat/FinalPERMsims/Sens4/a/';

Voc  = 1.2918235e18;    % (m3) volume ocean
rho  = 1.025e3;
Aoc  = 3.49e14;         % (m2) area ocean
%%% load all the files
tv= load([dir 'tv.DAT']);
tv11= load([dir 'tv11.DAT']);
epspcaV= load([dir 'epspcaV.DAT']);
d44ca= load([dir 'd44ca.DAT']);
d44cabulkV= load([dir 'd44cabulkV.DAT']);
caMean = load([dir 'caMean.DAT']);
dincaV = load([dir 'dincaV.DAT']);

FiN= load([dir 'FiN.DAT']);
FVC= load([dir 'FVC.DAT']);
V = load([dir 'V.DAT']);


epspcaV1  = [epspcaV(1)' epspcaV(1)'];
dincaV1  = [dincaV(1)' dincaV(1)'];
dcaMean=(d44ca(:,1).*V(1)+d44ca(:,2).*V(2)+d44ca(:,3)*V(3)+d44ca(:,10)*V(10)+d44ca(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
dcaMean1  = [dcaMean(1)' dcaMean(1)'];
d44cabulkV1 = [d44cabulkV(1)'  d44cabulkV(1)'];


tauCA=(caMean(1)*1e-3*Voc*rho)/((FiN+FVC)*Aoc)
figurehndl=55;    
FigHandle = figure(figurehndl);
box on
  set(FigHandle, 'Position', [800, 0, 600, 500]);  
subplot(311)
plot(tv  ,dincaV ,'r-','LineWidth',2);
hold on;
plot(tv11,dincaV1,'r-','LineWidth',2);
set(gca,'XTickLabel',[])
% ylabel('\epsilon_{carb} (‰)');
ylabel('\delta^{44}Ca_{riv} (‰)');
 text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
% title('Fractionation effect')
title('\delta^{44}Ca_{riv} effect')

subplot(312)
plot(tv  ,dcaMean ,'r-','LineWidth',2);
hold on;
plot(tv11,dcaMean1,'r-','LineWidth',2);
set(gca,'XTickLabel',[])
ylabel({'' ;'\delta^{44} Ca_{sw}(‰) '})
 text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')

subplot(313)
plot(tv  ,d44cabulkV ,'r-','LineWidth',2);
hold on;
plot(tv11,d44cabulkV1,'r-','LineWidth',2);
ylabel({'';'\delta^{44} Ca_{carb}(‰)'});
xlabel('Time (ky)')
set(gca,'YLim',[-0.65 -0.45]);
set(gca,'XTickLabel',[-50 0 50 100 150 200 250 300 350 400])
 text(0.02,0.98,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')