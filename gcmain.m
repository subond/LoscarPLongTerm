clear all;
close all
global time ws gamma RUN myfwc pco20  dbckv acv FERT ACT fckc dbcv
tp=59;    %time point at which we want to see all the values
FERT=0.40;
ACT=0.05;%0.045
pco20=300;
% CALLING fgeo10 to calculate all the fluxes
myT = tp;
tic
% %
% %
[ybbv,ybv,ycv,aav,bbv,ccv,fGG,fgkc,frkc,fekc,fdkc,flakc,fbch,fbbv,fwcv,fwsv1,fwsv2,fbcv,fbgv,rco2,pco2gca,fwgv,fmcv,fmgv,ggc,cgc,dg,dc]=gcfun10(59,0);
% [fbgv,fbcv,fwgv,fmcv,fmgv,fwcv,fGG,fgkc,frkc,fekc,fdkc,rco2,pco2gca] = gcfun11(3.54, 59);
% [ybbv,ybv,ycv,aav,bbv,ccv,fGG,fbbv,fwcv,fwsv1,fwsv2,fbcv,fbgv,rco2,pco2gc
% a,fwgv,fmcv,fmgv,ggc,cgc,dg,dc]=gcfun11(59,3.54);
toc
TP=59;
return;
MASSB=fwcv(TP)+fmcv(TP)+fwgv(TP)+fmgv(TP)-fbcv(TP)-fbgv(TP)
ISOMASSB=fwcv(TP)*dc(TP)+fmcv(TP)*(-4)+(fwgv(TP)+fmgv(TP))*dg(TP)-dbckv(TP)*fbcv(TP)-(dbckv(TP)-acv(TP))*fbgv(TP)
d13Cvc = -4.0;
d13Ckg = -21.2;
d13Cin = +2.0;
epsp57   = -31.8;
epsp52   = -30.97;
DBC=(dc(TP)*fwcv(TP)+d13Cvc*fmcv(TP)+dg(TP)*(fwgv(TP)+fmgv(TP))+acv(TP)*fbgv(TP))/(fbcv(TP)+fbgv(TP))
 dcNK57=(d13Cvc*fmcv(58)+d13Cin*fwcv(58)+d13Ckg*(fwgv(58)+fmgv(58))-0*fbcv(58)-epsp57*fbgv(58))/(fbcv(58)+fbgv(58));
 dcNK52=(d13Cvc*fmcv(53)+d13Cin*fwcv(53)+d13Ckg*(fwgv(53)+fmgv(53))-0*fbcv(53)-epsp52*fbgv(53))/(fbcv(53)+fbgv(53));
 DBC57=(dc(58)*fwcv(58)+d13Cvc*fmcv(58)+dg(58)*(fwgv(58)+fmgv(58))+acv(58)*fbgv(58))/(fbcv(58)+fbgv(58));
 DBC52=(dc(53)*fwcv(53)+d13Cvc*fmcv(53)+dg(53)*(fwgv(53)+fmgv(53))+acv(53)*fbgv(53))/(fbcv(53)+fbgv(53));
 %DBC=(dc(58)*(fwcv(58)+fmcv(58))+dg(58)*(fwgc(58)+fmgc(58))+acv(58)*fbgv(58))/(fbcv(58)+fbgv(58))
year=sprintf('%4.0f Ma',tp-1)
% FLUXES = ...
%     sprintf('  %d    %d    %d     %d      %d     %d    %d     %d',...
%     rco2(tp), fwcc(tp), fwgc(tp)5, fbcc(tp)-fwcc(tp), fmcc(tp), fmgc(tp),...
%     fbcc(tp), fbgc(tp))
FLUXES=[rco2(tp) fwcv(tp) fwgv(tp) fbcv(tp)-fwcv(tp) fmcv(tp) fmgv(tp)...
    fbcv(tp) fbgv(tp)]
disp('    RCO2      Fwc         Fwg       Fws       Fmc       Fmg     Fbc      Fbg');
%return
figure(1)
% plot(time,BIGD,'k');
% xlabel('Millions of years from present');
% ylabel('weighted \delta');
%pause;

%disp([c(p) g(p)]);
%disp([dc(p) dg(p) fbc fbg rco2(p)]);
figure(2)
plot(time,rco2,'k');
xlabel('Millions of years from present');
ylabel('RCO_2');
%print -deps fig711;
title('GEOCARB II with \Deltat = 1');
%disp([time' rco2']);
%pause;
%printsto -deps prob749;
figure(3)
DT=gamma(57)*log(rco2)+ws*(time)/570;
plot(time,DT,'k',time,ws*time/570,'k--');
xlabel('Millions of years from present');
ylabel('T - T_0');
title('GEOCARB II with \Deltat = 1');
%pause;
%printsto -deps prob746;
figure(4)
plot(time,ggc,'k');
axis([-600 0,1000 max(ggc)+10]);
xlabel('Millions of years from present');
ylabel('Organic carbon (10^1^8 mol)');
title('GEOCARB II with \Deltat = 1');
%pause;
%printsto -deps prob743;
figure(5)
plot(time,cgc,'k');
xlabel('Time (Millions of years from present)');
ylabel('Carbonate carbon (10^1^8 mol)');
title('GEOCARB II with \Deltat = 1');
%pause;
%printsto -deps prob742;
figure(6)
plot(time,dg,'k');
xlabel('Time (Millions of years from present)');
%sylabel('\black{{\symbol d}{_\10 g} ({^\12 o}/oo)}');
ylabel('\delta_g (per mil)');
title('GEOCARB II with \Deltat = 1');
%pause;
%printsto -deps prob745;
figure(7)
plot(time,dc,'k');
xlabel('Time (Millions of years from present)');
%sylabel('{\symbol d}{_\12 c} ({^\12 o}/oo)');
ylabel('\delta_C (per mil)');
title('GEOCARB II with \Deltat = 1');
figure(8)
plot(time,fwcv,'b',time,fwsv1,'g',time,fmcv,'k',time,fwgv,'r',time,fmgv,'m',time,fbcv,'--k',time,fbgv,'-.k')
legend('Fwc','Fws','Fmc','Fwg','Fmg','Fbc','Fbg')
figure(100)
 geoc3=[570	11.7
560	16.3
550	18.0
540	17.2
530	25.5
520	26.2
510	22.4
500	18.9
490	17.3
480	17.3
470	17.7
460	15.5
450	15.9
440	16.7
430	17.0
420	13.9
410	11.0
400	11.3
390	13.5
380	15.3
370	8.0
360	6.1
350	4.3
340	2.7
330	1.7
320	1.3
310	1.3
300	1.2
290	1.3
280	1.3
270	1.4
260	1.9
250	6.1
240	7.1
230	5.2
220	5.8
210	4.9
200	5.4
190	4.4
180	4.8
170	8.6
160	9.1
150	7.6
140	8.2
130	6.6
120	6.1
110	5.9
100	5.3
90	4.3
80	4.2
70	3.2
60	2.8
50	3.2
40	2.1
30	1.4
20	1.2
10	1.0
0	1.0
];
L=[0:10:570]';          %
P=[geoc3];             %
X=[0:1:570]';
geoc3=interp1(L,P,X);
plot(geoc3(1:10:end,1)*(-1),geoc3(1:10:end,2),'-ok', 'MarkerFaceColor','g')
hold on
plot(time(1:10:end),rco2(1:10:end),'-ok', 'MarkerFaceColor','k')
%set(gca,'XDir','reverse')
title('RCO2 GEOCARB III')
legend('B&K 2001 ', 'myGEOCARB III')
%printsto -deps prob744;