%clear all all;
diss = [-1:-0.2:-2];
 Fin=[0.0:1:2];
Cachng = csvread('isotope2.csv');
Ca = Cachng;
figure(2)
[Fin,diss]=contourf(Fin,diss,Ca,'w','LineWidth',1);
b=cool;
z=sortrows(b,2);
colormap(z)
xlabel('[Ca^{2+}] increase (%)','FontSize',12);
%set(gca,'XTickLabel','FontSize',14)
ylabel('\delta^{44}Ca_w (‰)','FontSize',12);
%zlabel('fgds');
clabel(Fin,diss,'FontSize',10)
title('\Delta\delta^{44}Ca_{sw} (Deep Pacific)')