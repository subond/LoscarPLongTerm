figure(111)
% pCO2
subplot(511)
plot(tv  ,pco2t ,'r-','LineWidth',2);
hold on;
plot(tv11,pco2t1,'r-','LineWidth',2);
hold off;
set(gca,'XTickLabel',[])
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
ylabel('pCO_2 (\muatm)');

subplot(512)
box on
hold on
hndl=plot(tv, FSichck/1e12,'m',tv,Finchck/1e12,'g');
set(hndl,'LineWidth',1)
legend('F_S_i','F_C')
hndl1=plot(tv11, FSichck1(9,:)./1e12,'m',tv11,Finchck1(9,:)./1e12,'g');
set(hndl1,'LineWidth',1)
%xlabel('Time (y)');
%ylabel('Silicate flux increase');
%hold on
%plot(tv,Finchck,'g','LineWidth',1)
%plot(tv11,Finchck1(9,:),'g','LineWidth',1)
hold off
ylabel('F_C, F_S_i (10^{12} mol C y^{-1})')
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
set(gca,'XTickLabel',[])

%pH
subplot(513)
box  on;
hold on;
for k=1:Nb
    plot(tv  ,phtv (:,k),sstr(2*k-1:2*k),'Color',cs(k));
end;
for k=1:Nb
    plot(tv11,phtv1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
set(gca,'XTickLabel',[])
xlabel('Time (y)');
ylabel('pH');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);

% omega
subplot(514)
box on
hold on
for k=1:5
    plot(tv  ,omegCSvt (:,k),sstr(2*k-1:2*k),'Color',cs(k));
end;
for k=1:5
    plot(tv11  ,omegCSvt1 (k,:),sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off

set(gca,'XTickLabel',[])
set(gca,'FontSize',fs);
% set(gca,'XLim',axx);
ylabel('Calcite saturation');

subplot(515);
plot(tv/1e3  ,ccdP,'k-','LineWidth',2);
hold on
plot(tv11/1e3,ccdP1,'k-','LineWidth',2);
hold off
set(gca,'XDir','normal','YDir','reverse')
xlabel('ky')
ylabel('CCD (m)')
%     set(gca,'XTickLabel',[])
legend('Pacific')

figure(112)

subplot(211)
plot(tv  ,CA  (:,9),'-b','LineWidth',2);
hold on
plot(tv11,CA1 (9,:),'-b','LineWidth',2);
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
set(gca,'XTickLabel',[])
% set(gca,'XTickLabel',[-50 0 50 100 150 200])
% xlabel('Time (ky)');
ylabel('Ca^2^+ (mmol kg^{-1})');
legend('Ca^2^+ (DP)')
refresh;
% xlabel('Time (ky)')
ylabel('[Ca^2^+] mmol/kg')
subplot(212)
box  on;
hold on;
for k=1:Nb
    plot(tv  ,d44ca (:,k),sstr(2*k-1:2*k),'Color',cs(k));
end;
for k=1:Nb
    plot(tv11,d44ca1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
end;
%     plot(tv  ,d44ca (:,9),sstr(2*9-1:2*9),'Color',cs(9));
%     plot(tv11,d44ca1(9,:),sstr(2*9-1:2*9),'Color',cs(9));
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
%set(gca,'YDir','reverse');
xlabel('Time (y)');
ylabel('\delta^{44}Ca (‰)');
Hl=legend(lstr,1);
set(Hl,'FontSize',10);

