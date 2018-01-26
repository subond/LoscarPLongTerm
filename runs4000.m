clear all all
%global CBl nCC

for slj=1:3
  % clear all all
    nCC=0.3+slj/10
    for po=1:1
CBl=1000.e15*po
Loscar10
dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
save                    dCalcium.DAT  dCalcium -ASCII -DOUBLE -TABS;
load dCalcium.dat
 bDP(:,:)=dCalcium
% ovco=bDP
% dCADP=ovco(slj,1)'
       bDPstr = [ 'myDataFile' num2str(slj) '.dat' ];
        %save(dCalcium);
save(bDPstr,'bDP', '-ASCII')
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
% save                    dCalcium.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load dCalcium.dat
% slj
%  bDP(slj,:)=dCalcium
%  ovco=bDP
%  dCADP=ovco(:,1)'
% save                    dCADP.DAT  dCADP -ASCII -DOUBLE -TABS;
 keep CBl  slj nCC bDP ovco dCADP
%%return
end
end
load myDataFile1.dat
load myDataFile2.dat
load myDataFile3.dat
load myDataFile4.dat
load myDataFile5.dat
load myDataFile6.dat
load myDataFile7.dat
load myDataFile8.dat
load myDataFile9.dat
A=myDataFile1;
B=myDataFile2;
C=myDataFile3;
D=myDataFile4;
E=myDataFile5;
F=myDataFile6;
G=myDataFile7;
H=myDataFile8;
I=myDataFile9;
CAallruns=[A 
    B
    C
    D
    E
    F
    G
    H
    I];
dCAallbasins=CAallruns(:,:);
save                    dCAallbasins.DAT  dCAallbasins -ASCII -DOUBLE -TABS;
csvwrite('allbasinsdCADP', dCAallbasins)
xlswrite('allbasinsdCADP', dCAallbasins)
dCADP=CAallruns(:,9);
save                    dCADP.DAT  dCADP -ASCII -DOUBLE -TABS;
csvwrite('allrunsdCADP', dCADP')
xlswrite('allrunsdCADP', dCADP')
%clear all all

% % end
% clear all all
% nSi=0.2;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium02.DAT  dCalcium -ASCII -DOUBLE -TABS;
%  load Pg5000dCalcium02.dat
%  beg=Pg5000dCalcium02;
%   xlswrite('Pg5000dCalcium02', beg) 
% clear all all
% 
% nSi=0.3;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium03.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium03.dat
%  beg=Pg5000dCalcium03;
%   xlswrite('Pg5000dCalcium03', beg)
%  clear all all
% 
% nSi=0.4;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium04.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium04.dat
%  beg=Pg5000dCalcium04;
%   xlswrite('Pg5000dCalcium04', beg)
%  clear all all
% 
% nSi=0.5;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium05.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium05.dat
%  beg=Pg5000dCalcium05;
%   xlswrite('Pg5000dCalcium05', beg)
%  clear all all
% 
% nSi=0.6;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium06.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium06.dat
%  beg=Pg5000dCalcium06;
%   xlswrite('Pg5000dCalcium06', beg)
%  clear all all
% 
% nSi=0.7;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium07.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium07.dat
%  beg=Pg5000dCalcium07;
%   xlswrite('Pg5000dCalcium07', beg)
%  clear all all
% 
% nSi=0.8;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium08.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium08.dat
%  beg=Pg5000dCalcium08;
%   xlswrite('Pg5000dCalcium08', beg)
%  clear all all
% 
% nSi=0.9;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium09.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium09.dat
%  beg=Pg5000dCalcium09;
%   xlswrite('Pg5000dCalcium09', beg)
%  clear all all
% 
% nSi=1.0;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium10.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium10.dat
%  beg=Pg5000dCalcium10;
%   xlswrite('Pg5000dCalcium10', beg)
%  clear all all
% 
% nSi=1.1;
% Loscar10
% dCalcium = ((max(CA)./(CA(1,:)))-1)*100;
%  save                    Pg5000dCalcium11.DAT  dCalcium -ASCII -DOUBLE -TABS;
% load Pg5000dCalcium11.dat
%  beg=Pg5000dCalcium11;
%   xlswrite('Pg5000dCalcium11', beg)
%  clear all all
%  load Pg5000dCalcium02.dat
%  b(1,:)=Pg5000dCalcium02;
%  load Pg5000dCalcium03.dat
%  b(2,:)=Pg5000dCalcium03;
%  load Pg5000dCalcium04.dat
%  b(3,:)=Pg5000dCalcium04;
%  load Pg5000dCalcium05.dat
%  b(4,:)=Pg5000dCalcium05;
%  load Pg5000dCalcium06.dat
%  b(5,:)=Pg5000dCalcium06;
%  load Pg5000dCalcium07.dat
%  b(6,:)=Pg5000dCalcium07;
%  load Pg5000dCalcium08.dat
%  b(7,:)=Pg5000dCalcium08;
%  load Pg5000dCalcium09.dat
%  b(8,:)=Pg5000dCalcium09;
%  load Pg5000dCalcium10.dat
%  b(9,:)=Pg5000dCalcium10;
%  load Pg5000dCalcium11.dat
%  b(10,:)=Pg5000dCalcium11;
%  dCADP=b(:,9)';