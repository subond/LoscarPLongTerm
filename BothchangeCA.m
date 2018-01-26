%clear all all;
Fin = [1:0.1:1.7];
diss =[1000:1000:6000];
Cachng = csvread('BothchangeCA.csv');
Ca = Cachng;
figure(1)
[Fin,diss]=contourf(Fin,diss,Ca,'w','LineWidth',1);
b=autumn;
z=sortrows(b,-2);
colormap(z)
xlabel('Car and SI weathering Flux (norm.)','FontSize',12);
%set(gca,'XTickLabel','FontSize',14)
ylabel('Dissolved Carbon (Pg)','FontSize',12);
%zlabel('fgds');
clabel(Fin,diss,'FontSize',12)
title('% change in [Ca] (Deep Pacific)')
