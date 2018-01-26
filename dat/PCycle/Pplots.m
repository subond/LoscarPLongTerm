clear all
close all

fs = 10;
eI=0.78;

cs = 'gggkkkrrrbgkr';
sstr = '- ---.- ---.- ---. -: : : ';
lstr0 = 'LALILPIAIIIPDADIDP HLTITDT';
Nb=13;

% to look at the individual runs change the path to the
% folder you are interested in e.g. PonW, Poff, P0, Leak1...

dir = 'rain1/';
tv= load([dir 'tv.DAT']);
tv11= load([dir 'tv11.DAT']);
pco2t= load([dir 'pco2t.DAT']);

oIpv= load([dir 'oIpv.DAT']);
oIv= load([dir 'oIv.DAT']);
PPLvv= load([dir 'PPLvv.DAT']);
EPLvv= load([dir 'EPLvv.DAT']);
Fcapv= load([dir 'Fcapv.DAT']);
Ffepv= load([dir 'Ffepv.DAT']);
EPH = load([dir 'EPH.DAT']);
PPH = load([dir 'PPH.DAT']);

dox= load([dir 'dox.DAT']);
p= load([dir 'p.DAT']);
T= load([dir 'temp.DAT']);









% return;
%------------------ d13C PETM data
%== Zachos data
ZchF3a = load('../Zachos/ZchF3aNEW.DAT');
age3  = ZchF3a(:,3);
c13c3 = ZchF3a(:,04);
%== Roehl data
% new age model, G^3, 2007
NEW690     = load('../Zachos/NEW690.txt');
d690n      = NEW690(:,2);
a690n      = NEW690(:,3);
% c13 data
Roehl690   = load('../Zachos/Roehl690.txt');
d690       = Roehl690(:,1);
% interpolate
age690blk  = interp1(d690n,a690n,d690);
c13c690blk = Roehl690(:,5);


myRoehl690   = load('../Zachos/690roehl07.csv');
myd690       = myRoehl690(:,3);
mya690       = myRoehl690(:,2);
% interpolate
myage690blk  = interp1(d690n,a690n,d690);
myc13c690blk = myRoehl690(:,4);


%Farley 2003 3He age model Hole 690
Fdata= csvread('../Zachos/Farley690Age.csv');
Fdepth = Fdata(:,1);
Fage = Fdata(:,3);
Fage690 = interp1(Fdepth,Fage, d690);
myFage690 = interp1(Fdepth,Fage, myd690);

% Murphy, Farley Zachos 2010
Mdata= csvread('../Zachos/Murphy1266.csv');
depth = Mdata(:,1);
d13c1266 = Mdata(:,2);

% Murph, Farley, Zachos age model
Murphy1266 = csvread('../Zachos/ageMurphy.csv');
depth1266 = Murphy1266(:,1);
age1266 = Murphy1266(:,3);
age1266blk = interp1(depth1266,age1266, depth);

% Roehl 2007, site 1266 age model
Roehl1266     = load('../Zachos/ageRoehl1266.csv');
Rd690n      = Roehl1266(:,1);
Ra690n      = Roehl1266(:,4);
% interpolate
Rage1266blk  = interp1(Rd690n,Ra690n,depth);
% c13c1266blk = Roehl690(:,5);


cramer690     = load('../Zachos/Cramer690BIsotopes.csv');
cd690n      = cramer690(:,1);
C13c690      = cramer690(:,5);
% interpolate
cage690blk  = interp1(d690n,a690n,cd690n)*1e3;

penData1409 = load('../Zachos/1409CaCO3.csv');
penAge1409 = load('../Zachos/1409age.csv');
penD = penData1409(:,1);
penCaCO31409 = penData1409(:,2);
penAgeD = penAge1409(:,1);
penAgeAge = penAge1409(:,2);
pAge1409 = interp1(penAgeD, penAgeAge, penD);

penData1403 = load('../Zachos/1403CaCO3.csv');
penAge1403 = load('../Zachos/1403age.csv');
penD = penData1403(:,1);
penCaCO31403 = penData1403(:,2);
penAgeD = penAge1403(:,1);
penAgeAge = penAge1403(:,2);
pAge1403 = interp1(penAgeD, penAgeAge, penD);

% figure
% plot(pAge1409, penCaCO31409);
% figure
% plot(pAge1403, penCaCO31403);

pen1409 = load('../Zachos/pen1409.csv');
pen1403 = load('../Zachos/pen1403.csv');

% figure
% hold on
% plot(pen1409(:,1),pen1409(:,2),'r')
% plot(pAge1409, penCaCO31409);
% hold off
%
% figure
% hold on
% plot(pen1403(:,1),pen1403(:,2),'r')
% plot(pAge1403, penCaCO31403);
% hold off

MS  = 'MarkerSize';
ms  = 5; %4
MFC = 'MarkerFaceColor';
LW ='LineWidth';
Dam = 9.; %7

FS='FontSize';
fs=18; %12

FW = 'FontWeight';
fw = 'Bold'; %noramal
%
% % plot(cage690blk,C13c690)
% % return
% FigHandle=figure;
% set(FigHandle, 'Position', [100, 100, 950, 800]);
%
% subplot(211)
% box on
% hold on
% % plot(Fage690,c13c690blk,'k- d',LW,2,MFC,'k',MS,ms);
% % plot(age690blk,c13c690blk,'g- o',LW,2,MFC,'g',MS,ms);
% plot(myFage690,myc13c690blk,'k- d',LW,2,MFC,'k',MS,ms);
% plot(mya690,myc13c690blk,'g- o',LW,2,MFC,'g',MS,ms);
% hold off
% ylabel('\delta^{13}C (‰; ODP Site 690)',FS,fs, FW, fw);
% text(0.02,0.95,'a)','Units', 'Normalized',FS,fs, FW, fw);
% legend('^3He-based age model','Cycle-based age model','Location','Southeast');
% set(gca,'XTickLabel',[],'XLim',[-50 300], FS,fs, FW, fw);
%
% subplot(212)
% box on
% hold on
% plot(age1266blk,d13c1266,'k- d',LW,2,MFC,'k',MS,ms);
% plot(Rage1266blk,d13c1266,'g- o',LW,2,MFC,'g',MS,ms);
% hold off
% ylabel('\delta^{13}C (‰; ODP Site 1266)',FS,fs, FW, fw);
% text(0.02,0.95,'b)','Units', 'Normalized',FS,fs, FW, fw);
% legend('^3He-based age model','Cycle-based age model','Location','Southeast');
% set(gca,'XLim',[-50 300], FS, fs, FW, fw);
% return

p1     = [p(1,:)' p(1,:)'];
dox1   = [dox(1,:)' dox(1,:)'];

orgPb = (1-(eI+oIpv)).*(PPLvv+PPH);
totPbur = orgPb+Fcapv+Ffepv;
orgCb = (1-(eI+oIv)).*(EPLvv+EPH);

orgPb1 = [orgPb(1) orgPb(1)];
totPbur1 = [totPbur(1) totPbur(1)];
orgCb1 = [orgCb(1) orgCb(1)];
Fcapv1 = [Fcapv(1) Fcapv(1)];
Ffepv1 = [Ffepv(1) Ffepv(1)];

CPbur = orgCb./orgPb;
CPbur1 = [CPbur(1) CPbur(1)];

% ind = find(tv<=4.0e5);
ind = find(tv>=0.61e5 & tv<=1.0e5 );
integ = trapz(tv(ind),orgCb(ind));
excessC = integ - (orgCb(1)*(tv(ind(end))-tv(ind(1)))); % mol C
excessCg = excessC*12 %g C




%%% array to plot results from multiple runs (3 in this case)
sensFlag =5; % 0: Comparing model with P cy cle and without
% 1: Sensitivity studies of Standard P
% 2: Oxygen sensitivity
% 3: Floegel Oxygen fraction sensitivity
% 4: Floegel C and P fractions buried sensitivity
% 5: Leak (Default) PETM scenario, original Loscar vs. Loscar with P burial
% 6: Original original LOSCAR version, 0: default leak; 1: doubling Corg
% butial over 32ky during the onset to get ~2000 PgC burial increase
% 7: Original Loscar with leak, LOSCAR-P with the same leak, LOSCAR-P with
% the leak (control run, not preferred run) but CaCO3 rain constant
% 8: Sensitivity studies but with leak included
% 9: Preferred Leak scenario with weathering sensitivity runs
% 10: Comparing LOSCAR-P variable CaCO3 and constant CaCO3 rain
axx = [-50 200];
if(sensFlag == 1)
    ind = 8;
    %     clrs ='krbgcymkr' ;
    %     symb= '.--x---+s---.- ---. -: : : ';
    clrs ='kkkgggrrr' ;
    symb= '-- --.-- --.-- --.';
    Legend=cell(9,1);%  five positions
    lgnd = {'f_{OP} = 0.50%; f_{OC} = 1.00% ','f_{OP} = 1.00%; f_{OC} = 1.00%','f_{OP} = 0.25%; f_{OC} = 1.00% ',...
        'f_{OP} = 0.50%; f_{OC} = 0.50% ','f_{OP} = 1.00%; f_{OC} = 0.50% ','f_{OP} = 0.25%; f_{OC} = 0.50% '...
        'f_{OP} = 0.50%; f_{OC} = 2.00% ','f_{OP} = 1.00%; f_{OC} = 2.00% ','f_{OP} = 0.25%; f_{OC} = 2.00% ' };
elseif(sensFlag == 0)
    ind = 2;
    clrs ='brkg' ;
    symb = '- -- .-';
    Legend=cell(3,1);%  three positions
    lgnd = {'Original LOSCAR', 'With P cycle','LOSCAR-P'};
elseif(sensFlag == 2)
    ind = 2;
    clrs ='krb' ;
    symb = '- -- -';
    Legend=cell(3,1);%  three positions
    lgnd = {'oxA = 0.4', 'oxA = 0.0','oxA = 1.0'};
elseif(sensFlag == 3)
    ind = 2;
    clrs ='rbkgcmyr' ;
    symb = '- -- - - - - - -';
    Legend=cell(3,1);%  five positions
    lgnd = {'oxA = 0.4', 'oxA = 0.0','oxA = 1.0'};
    %     lgnd = {'f_{OP} = 0.50%; f_{OC} = 1.00% ','f_{OP} = 1.00%; f_{OC} = 1.00% ','f_{OP} = 0.25%; f_{OC} = 1.00% ','f_{OP} = 0.50%; f_{OC} = 2.00% ','f_{OP} = 0.50%; f_{OC} = 0.50% ','f_{OP} = 0.50%; f_{OC} = 0.50% ' };
elseif(sensFlag == 4)
    ind = 4;
    clrs ='rbkgcmyr' ;
    symb = '- -- - - - - - -';
    Legend=cell(5,1);%  five positions
    lgnd = {'f_{OP} = 2.00%; f_{OC} = 1.00% ','f_{OP} = 2.00%; f_{OC} = 2.00% ','f_{OP} = 2.00%; f_{OC} = 0.5% ','f_{OP} = 4.00%; f_{OC} = 1.00% ','f_{OP} = 1.00%; f_{OC} = 1.00% '};
elseif(sensFlag == 5)
    ind = 3;
    clrs ='bkrG' ;
    symb = '- -- -- ';
    Legend=cell(4,1);%  three positions
    lgnd = {'Original LOSCAR', 'LOSCAR-P, f_{OP} = 0.50%; f_{OC} = 1.00% ','LOSCAR-P, f_{OP} = 1.00%; f_{OC} = 2.00% ','LOSCAR-P, C leak = 2,500 PgC'};
elseif(sensFlag == 6)
    ind = 1;
    clrs ='bg' ;
    symb = '- - ';
    Legend=cell(2,1);%  three positions
    lgnd = {'Original LOSCAR', 'Enhanced F_{OCB}'};
elseif(sensFlag == 7)
    ind = 2;
    clrs ='bkrg' ;
    symb = '- -- -- ';
    Legend=cell(3,1);%  three positions
    lgnd = {' Original LOSCAR', ' LOSCAR-P control',' LOSCAR-P const. CaCO_3 rain'};
elseif(sensFlag==8)
    ind = 8;
    %     clrs ='krbgcymkr' ;
    %     symb= '.--x---+s---.- ---. -: : : ';
    clrs ='kkkgggrrr' ;
    symb= '-- --.-- --.-- --.';
    Legend=cell(9,1);%  five positions
    lgnd = {'f_{OP} = 0.50%; f_{OC} = 1.00% ','f_{OP} = 1.00%; f_{OC} = 1.00%','f_{OP} = 0.25%; f_{OC} = 1.00% ',...
        'f_{OP} = 0.50%; f_{OC} = 0.50% ','f_{OP} = 1.00%; f_{OC} = 0.50% ','f_{OP} = 0.25%; f_{OC} = 0.50% '...
        'f_{OP} = 0.50%; f_{OC} = 2.00% ','f_{OP} = 1.00%; f_{OC} = 2.00% ','f_{OP} = 0.25%; f_{OC} = 2.00% ' };
elseif(sensFlag==9)
    ind = 2;
    clrs ='ggg' ;
    symb = '- -- :';
    Legend=cell(3,1);%  three positions
    lgnd = {'nC=0.4; nS=0.2', 'nC=1.0; nS=0.8','nC=0.1; nS=0.1'};
elseif(sensFlag == 10)
    ind = 1;
    clrs ='gk' ;
    symb = '- -- -- ';
    Legend=cell(2,1);%  three positions
    lgnd = {'LOSCAR-P, C leak = 2,500 PgC','CaCO_3 rain = const.'};

end

for i=1:1:ind+1
    V{i,:}= load([dir 'V.DAT']);
    if(sensFlag == 1)
        dira = ['P' num2str(i-1) '/'];
    elseif(sensFlag == 0)
        if (i==1)
            dira = 'Poff/';
        end
        if(i==2)
            dira = 'Pon/';
        end
        if (i==3)
            dira = 'PonW/';
        end
    elseif(sensFlag == 2)
        if (i==1)
            dira = ['P' num2str(i-1) '/'];
        else
            dira = ['POx' num2str(i-2) '/'];
        end
    elseif(sensFlag == 3)
        dira = ['Pfloegel' num2str(i-1) '/'];
    elseif(sensFlag == 4)
        k=3;
        dira = ['Pfloegel' num2str(i-1+k) '/'];
    elseif(sensFlag == 5)
        dira = ['Leak' num2str(i-1) '/'];
    elseif(sensFlag == 6)
        dira = ['OR' num2str(i-1) '/'];
    elseif(sensFlag == 7)
        if (i==1 || i==2)
            dira = ['Leak' num2str(i-1) '/'];
        else
            dira = 'ConstInorgRain/';
        end
    elseif(sensFlag==8)
        dira = ['LeakSens' num2str(i-1) '/'];
    elseif(sensFlag==9)
        dira = ['WLeak' num2str(i-1) '/'];
    elseif(sensFlag==10)
        dira = ['rain' num2str(i-1) '/'];
    end
    tva{i,:}= load([dira 'tv.DAT'])./1e3;
    tv11a{i,:}= load([dira 'tv11.DAT'])./1e3;
    pco2ta{i,:}= load([dira 'pco2t.DAT']);
    pco2ta1{i,:}= [pco2ta{i}(1) pco2ta{i}(1)];
    d13c{i,:}= load([dira 'd13c.DAT']);
    d13c1{i,:}= [d13c{i}(1) d13c{i}(1)];
    d13cmean{i,:} = (d13c{i,:}(:,1).*V{i,:}(1)+d13c{i,:}(:,2).*V{i,:}(2)+d13c{i,:}(:,3)*V{i,:}(3)+d13c{i,:}(:,11)*V{i,:}(11))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(11)));
    d13cmean1{i,:} = [d13cmean{i,:}(1,1)'  d13cmean{i,:}(1,1)'];
    fcAa{i,:}= load([dira 'fcA.DAT']);
    fcAa1{i,:}= [fcAa{i}(1,:);fcAa{i}(1,:)]';





    CCDa{i,:}= load([dira 'ccdA.DAT'])./1e3;
    CCDa1{i,:}= [CCDa{i}(1) CCDa{i}(1)];
    if(sensFlag==0 && i==1)
        CCDa{i} = CCDa{i}+0.24;
        CCDa1{i} = CCDa1{i}+0.24;
    end
    CCDi{i,:}= load([dira 'ccdI.DAT'])./1e3;
    CCDi1{i,:}= [CCDi{i}(1) CCDi{i}(1)];

    CCDp{i,:}= load([dira 'ccdP.DAT'])./1e3;
    CCDp1{i,:}= [CCDp{i}(1) CCDp{i}(1)];

    CCDt{i,:}= load([dira 'ccdT.DAT'])./1e3;
    CCDt1{i,:}= [CCDt{i}(1) CCDt{i}(1)];

    doxa{i,:}= load([dira 'dox.DAT']);
    dox1a{i,:}= [doxa{i,:}(1,:)' doxa{i,:}(1,:)'];

    omegCa{i,:} = load([dira 'OmegCS.DAT']);
    omegC1a{i,:}= [omegCa{i}(1,:);omegCa{i}(1,:)]';


    if(sensFlag ~= 6)



        pa{i,:}= load([dira 'p.DAT']).*1e3;
        pa1{i,:}= [pa{i,:}(1,:)' pa{i,:}(1,:)'];
        PO4mean{i,:} = (pa{i,:}(:,1).*V{i,:}(1)+pa{i,:}(:,2).*V{i,:}(2)+pa{i,:}(:,3)*V{i,:}(3)+pa{i,:}(:,4).*V{i,:}(4)+pa{i,:}(:,5).*V{i,:}(5)+pa{i,:}(:,6)*V{i,:}(6)+pa{i,:}(:,7).*V{i,:}(7)+pa{i,:}(:,8).*V{i,:}(8)+pa{i,:}(:,9)*V{i,:}(9)+pa{i,:}(:,10)*V{i,:}(10)+pa{i,:}(:,11)*V{i,:}(11)+pa{i,:}(:,12).*V{i,:}(12)+pa{i,:}(:,13).*V{i,:}(13))./(sum(V{i,:}));
        Ta= load([dira 'temp.DAT']);
        if(sensFlag>0 && (sensFlag~=5 )&& sensFlag~=7)
            oIpva{i,:}= load([dira 'oIpv.DAT']);
            oIva{i,:}= load([dira 'oIv.DAT']);
            PPLvva{i,:}= load([dira 'PPLvv.DAT']);
            EPLvva{i,:}= load([dira 'EPLvv.DAT']);
            Fcapva{i,:}= load([dira 'Fcapv.DAT']);
            Ffepva{i,:}= load([dira 'Ffepv.DAT']);
            EPHa{i,:} = load([dira 'EPH.DAT']);
            PPHa{i,:} = load([dira 'PPH.DAT']);
        elseif(sensFlag == 0 || sensFlag == 5 || sensFlag == 7)
            if(i==1)
                oIpva{i,:}= zeros(1,size(tva{i}));
                oIva{i,:}= zeros(1,size(tva{i}));
                PPLvva{i,:}= zeros(1,size(tva{i}));
                EPLvva{i,:}= zeros(1,size(tva{i}));
                Fcapva{i,:}= zeros(1,size(tva{i}));
                Ffepva{i,:}= zeros(1,size(tva{i}));
                EPHa{i,:} = zeros(1,size(tva{i}));
                PPHa{i,:} = zeros(1,size(tva{i}));
            else
                oIpva{i,:}= load([dira 'oIpv.DAT']);
                oIva{i,:}= load([dira 'oIv.DAT']);
                PPLvva{i,:}= load([dira 'PPLvv.DAT']);
                EPLvva{i,:}= load([dira 'EPLvv.DAT']);
                Fcapva{i,:}= load([dira 'Fcapv.DAT']);
                Ffepva{i,:}= load([dira 'Ffepv.DAT']);
                EPHa{i,:} = load([dira 'EPH.DAT']);
                PPHa{i,:} = load([dira 'PPH.DAT']);

                EPLvva1{i,:} = [EPLvva{i,:}(1,1)' EPLvva{i,:}(1,1)'];
                EPHa1{i,:} = [EPHa{i,:}(1,1)' EPHa{i,:}(1,1)'];
            end
            if(sensFlag == 5)
                RlsCtva {i,:} = load([dira 'RlsCtv.DAT']);
                DTSva {i,:} = load([dira 'fDTS.DAT'])/1.e3;
                Rma{i,:}  =  (RlsCtva{i,:}(1:end-1)+RlsCtva{i,:}(2:end  ))/2./1.e15; % Gt C y
                dtva{i,:} =  (    tva{i,:}(2:end  )-    tva{i,:}(1:end-1));
                tma{i,:}  =  (    tva{i,:}(1:end-1)+    tva{i,:}(2:end  ))/2.;
                trla{i,:} = [0.; tma{i,:}];
                rla{i,:} = [0.  Rma{i,:}];
                % total release
                TRa{i,:} = sum(Rma{i,:}*diff(tva{i,:}));  % Gt C
                % split release

            end
        end
        orgPba{i,:} = (1-(eI+oIpva{i,:})).*(PPLvva{i,:}+PPHa{i,:});
        totPbura{i,:} = orgPba{i,:}+Fcapva{i,:}+Ffepva{i,:};
        orgCba{i,:} = (1-(eI+oIva{i,:})).*(EPLvva{i,:}+EPHa{i,:});

        orgPb1a{i,:} = [orgPba{i,:}(1,1)' orgPba{i,:}(1,1)'];
        totPbur1a{i,:} = [totPbura{i,:}(1,1)' totPbura{i,:}(1,1)'];
        orgCb1a{i,:} = [orgCba{i,:}(1,1)'  orgCba{i,:}(1,1)'];
        Fcapv1a{i,:} = [Fcapva{i,:}(1,1)' Fcapva{i,:}(1,1)'];
        Ffepv1a{i,:} = [Ffepva{i,:}(1,1)' Ffepva{i,:}(1,1)'];

        CPbura{i,:} = orgCba{i,:}./orgPba{i,:};
        CPbur1a{i,:} = [CPbura{i,:}(1,1)' CPbura{i,:}(1,1)'];

        inda{i,:} = find(tva{i,:}<=400);
        intega{i,:} = trapz(tva{i,:}(inda{i,:}),orgCba{i,:}(inda{i,:}));
        excessCa{i,:} = intega{i,:} - (orgCba{i,:}(1)*tva{i,:}(inda{i,:}(end))); % mol C
        excessCga{i,:} = excessCa{i,:}*12*1e3/1e15 %Pg C

        if(~sensFlag)
            orgCba{1}=ones(1,size(tva{1})).*4.58e12;
            orgCb1a{1} = [orgCba{1,:}(1,1)'  orgCba{1,:}(1,1)'];
        end
    end
end


labelFS=16;
legendFS=12;
textFS=12;
tickFS=12;
linew=2.5;
% symb='- - - - - ';
% clrs ='bgg' ;
% %%%% GOLDSCHMIDT PLOTS
% % figures 1 & 2, sensFlag = 0;
% % figure 3, sensFlag = 5;
% % figure 4, sensFlag = 10;
% FigHandle=figure;
% set(FigHandle, 'Position', [400, 100, 500, 950]);
% subplot(211)
% box on
% hold on
% mark = 'o x d o x d o x d ';
% for i=1:1:ind+1
%     if(i==1 || i==3)
%     h1=plot(tva{i}, pco2ta{i},[symb(2*i-1:2*i)],'Color',clrs(i),LW,linew);
%     Legend{i}=lgnd{i} ;
%     end
% end
% for i=1:1:ind+1
%     if(i==1 || i==3)
%     h2=plot(tv11a{i}, pco2ta1{i},symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
%     end
% end
% hold off
% 
% ylabel('pCO_2 (\muatm)','FontSize',labelFS);
% set(gca,'XTickLabel',[],'XLim',axx,'fontWeight','bold')
% legend(lgnd{1},lgnd{3});
% if(sensFlag==0)
% %     title('No P burial vs P')
% elseif(sensFlag==1)
%     title('org P and C fraction buried Sensitivity runs')
% elseif(sensFlag==2)
%     title('Oxygen sensitivity runs')
% elseif(sensFlag==2)
%     title('Floegel runs')
% end
% text(0.02,0.90,'a)','Units', 'Normalized','fontWeight','bold');
% subplot(212)
% box on
% hold on
% for i=1:1:ind+1
%     % scaling d13C data so they have the same initial value
%     d13cmean{i}=d13cmean{i}-d13cmean{i}(1)+d13cmean{1}(1);
%     if(i==1 || i==3)
%     plot(tva{i}, d13cmean{i},symb(2*i-1:2*i),'Color',clrs(i),LW,linew)
%     end
% end
% for i=1:1:ind+1
%     % scaling d13C data so they have the same initial value
%     d13cmean1{i}=d13cmean1{i}-d13cmean1{i}(1)+d13cmean{1}(1);
%     if(i==1 || i==3)
%     plot(tv11a{i}, d13cmean1{i},symb(2*i-1:2*i),'Color',clrs(i),LW,linew)
%     end
% end
% hold off
% ylabel('Mean surface \delta^{13}C (‰)','FontSize',labelFS);
% xlabel('Time (ky)','FontSize',labelFS)
% text(0.02,0.90,'b)','Units', 'Normalized','fontWeight','bold');
% set(gca,'XLim',axx,'fontWeight','bold');
% 
% 
% 
% FigHandle=figure;
% set(FigHandle, 'Position', [400, 100, 500, 950]);
% subplot(211)
% box on
% hold on
% for i=1:1:ind+1
%     if(i==1 || i==3)
%     plot(tva{i}, orgCba{i}/1e12,symb(2*i-1:2*i),'Color',clrs(i),LW,linew)
%     end
% end
% for i=1:1:ind+1
%     if(i==1 || i==3)
%     plot(tv11a{i}, orgCb1a{i}/1e12,symb(2*i-1:2*i),'Color',clrs(i),LW,linew)
%     end
% end
% hold off
% legend(lgnd{1},lgnd{3});
% ylabel('F_{OCB} (10^{12} mol yr^{-1})','FontSize',labelFS)
% set(gca,'XTickLabel',[],'XLim',axx,'fontWeight','bold')
% 
% text(0.02,0.90,'c)','Units', 'Normalized','fontWeight','bold');
% subplot(212)
% box on
% hold on
% for i=1:1:ind+1
%     if(i==1 || i==3)
%     plot(tva{i}, totPbura{i}/1e10,symb(2*i-1:2*i),'Color',clrs(i),LW,linew)
%     end
% end
% for i=1:1:ind+1
%     if(i==1 || i==3)
%     plot(tv11a{i}, totPbur1a{i}/1e10,symb(2*i-1:2*i),'Color',clrs(i),LW,linew)
%     end
% end
% hold off
% text(0.02,0.90,'d)','Units', 'Normalized','fontWeight','bold');
% set(gca,'XLim',axx,'fontWeight','bold');
% ylabel({'Total P burial';' (10^{10} mol yr^{-1})'},'FontSize',labelFS)
% xlabel('Time (ky)','FontSize',labelFS)

if(sensFlag==5)
    %%% Calculating mass accumulation
    myInd = 4; % array index of simulation we want to calculate this for
    tvi=[1:1:400];
    Ci=interp1(tva{myInd},orgCba{myInd},tvi);
    integi = (trapz(tvi(1:400),Ci(1:400))*1e3-(orgCba{myInd}(1)*(tvi(400)-tvi(1))*1e3))*12;

    lti=length(tvi);
    mycum(1)=0;
    for i=1:lti-1
        mycum(i+1)=mycum(i)+(trapz(tvi(i:i+1),Ci(i:i+1))*1e3-(orgCba{myInd}(1)*(tvi(i+1)-tvi(i))*1e3))*12;
    end



    john08b=load('../Zachos/johnExcessC1.csv'); %Bass river data
    john08l=load('../Zachos/johnExcessC2.csv'); %Lodo Gulch data

    timeb=john08b(:,1);
    cumCb=john08b(:,2);
    timel=john08l(:,1);
    cumCl=john08l(:,2);
    figure
    box on
    hold on
    plot(timeb,cumCb,'k- d',LW,2,MFC,'k',MS,ms)
    plot(timel,cumCl,'r- o',LW,2,MFC,'k',MS,ms)
    % plot(tvi,mycum/(52e6*1e10),'g',LW,2) %mass accum. in g/cm^2
    [AX,H1,H2] = plotyy(tvi,mycum/(52e6*1e10),tvi,mycum/1e15,'plot');
    hold off
    xlabel('Time (ky)', 'FontSize',labelFS)
    ylabel('Global organic carbon sequestration (g/cm^2)','FontSize',labelFS)
    legend('Bass River','Lodo Gulch','LOSCAR-P','Location','Southeast','FontSize',14)
    set(gca,'fontWeight','bold');
    set(get(AX(1),'Ylabel'),'String','Global organic carbon sequestration (g/cm^2)')
    set(get(AX(2),'Ylabel'),'String','Pg C')
    set(H1,'LineStyle','-','Color','g',LW,linew)
    set(H2,'LineStyle','none','Color','r',LW,linew)
    set(AX,{'ycolor'},{'k';'k'},'XLim',[0 250])      % ...and to adjust the axis color
end
% cumCburi(1)=0;
% for i=1:lti-1
%     cumCburi(i+1)=cumCburi(i)+(((Ci(i+1)+Ci(i))*12))*1e3-Ci(1)*(1e3)*12;
% end
%
%
% figure
% plot(tvi,cumCburi,'r')
%
%
%
% lt=length(tva{myInd});
% for i=1:lt-1
%    dtva(i)     = (tva{myInd}(i+1)-tva{myInd}(i))*1e3;
% end
%
%
% cumCbur(1)=0;
% for i=1:lt-1
%     cumCbur(i+1)=cumCbur(i)+(((orgCba{myInd}(i+1)+orgCba{myInd}(i))*12)*dtva(i)/2-orgCba{myInd}(1)*12);
% end
% figure
% plot(tva{myInd},cumCbur) %/(52e6*1e10)


%---------------------- Calc Release
RlsCtv = load([ 'Leak1/' 'RlsCtv.DAT']);
RlsCtv2 = load([ 'Leak3/' 'RlsCtv.DAT']);
newtv= load(['Leak1/' 'tv.DAT']);
newtv2= load(['Leak3/' 'tv.DAT']);
DTSv = load([ 'Leak3/' 'fDTS.DAT'])/1.e3;

Rm  =  (RlsCtv(1:end-1)+RlsCtv(2:end  ))/2./1.e15; % Gt C y
dtv =  (    newtv(2:end  )-    newtv(1:end-1));
tm  =  (    newtv(1:end-1)+    newtv(2:end  ))/2.;
trl = [0.; tm];
rl = [0.  Rm];
% total release
TR = sum(Rm*diff(newtv));  % Gt C
% split release

Rm2  =  (RlsCtv2(1:end-1)+RlsCtv2(2:end  ))/2./1.e15; % Gt C y
dtv2 =  (    newtv2(2:end  )-    newtv2(1:end-1));
tm2  =  (    newtv2(1:end-1)+    newtv2(2:end  ))/2.;
trl2 = [0.; tm2];
rl2 = [0.  Rm2];
% total release
TR2 = sum(Rm2*diff(newtv2));  % Gt C




% % split release
%
% DTS  = DTSv(1)+0.1; % ky duration 1st bump
% DTS2 = DTSv(2)+0.1; % ky duration 1st cont release
% ts3  = DTSv(3)    ; % ky onset    2nd bump
% DTS3 = DTSv(4)+0.1; % ky duration 2nd bump
% DTS4 = DTSv(5)+0.1; % ky duration 2nd cont release
%
% % indices
% [a,kt1] = min(abs(tv-DTS     )); % end   1st bump
% [a,kt3] = min(abs(tv-ts3     )); % start 2nd bump
% [a,kt4] = min(abs(tv-ts3-DTS3)); % end   2nd bump
%
%
% R1 = sum(Rm(    1:kt1).*dtv(    1:kt1)'*1e3); % Gt C
% R2 = sum(Rm(kt1+1:kt3).*dtv(kt1+1:kt3)'*1e3); % Gt C
% R3 = sum(Rm(kt3+1:kt4).*dtv(kt3+1:kt4)'*1e3); % Gt C
% R4 = sum(Rm(kt4+1:end).*dtv(kt4+1:end)'*1e3); % Gt C
%
% % bumps
% str1 = sprintf('%4.0f Pg C',10*round(R1/10));
% str3 = sprintf('%4.0f Pg C',10*round(R3/10));
% % sum cont release
% str2 = sprintf('%4.0f Pg C',10*(round(R2/10)+round(R4/10)));
% % sum cont release + bump2
% str4 = sprintf('%4.0f Pg C',10*round(R3/10)+10*(round(R2/10)+round(R4/10)));
% % sumtotal
% str  = sprintf('%4.0f Pg C',10*(round(R1/10)+...
%                                 round(R2/10)+...
%                                 round(R3/10)+...
%                                 round(R4/10)));
% cstr = sprintf('%+3.0f',[-50]);
%
%


if(sensFlag==5)
    fcA = load('orfcA.DAT');
    fcA1     = [fcA(1,:)' fcA(1,:)'];
    ortv= load('ortv.DAT');
    ortv11= load('ortv11.DAT');

    figure
    %----------------------- Release --------------------------%
    subplot(311)
    % axxn=[-50 200]*1e3;

    %  plot(trla{1},rla{1},'-k');
    %  plot(trl,rl,'-k');
    box on;
    hold on;
    % for i=1:1:ind+1
    %     fill(trla{1},rla{1},'r');
    % end
    fill(trl2,rl2,'g ');
    fill(trl,rl,'b ');

    hold off;
    ylabel({'Carbon Input';'(Pg C y^{-1})'},'FontSize',labelFS);
    % set(gca,'FontSize',labelFS);
    set(gca,'XLim',axx.*1e3);
    set(gca,'YLim',[0 .7],'FontWeight','bold');
    set(gca,'XTickLabel',[]);
    % set(gca,'XTickL',[]);

    % text(020,0.400,['      ' str1],'FontSize',fs);
    % text(070,0.100,['      ' str4],'FontSize',fs);
    % annotation('arrow',[.37 .32],[.905 .905]);
    % annotation('arrow',[.54 .49],[.86 .86]);
    text(0.02,0.90,'a)','Units', 'Normalized','FontSize',textFS);
    % text(005,0.400,['  \leftarrow' str1],'FontSize',fs);
    % text(070,0.100,['  \leftarrow' str4],'FontSize',fs);
    % text(0.7,0.800,['Total = '     str ],'FontSize',fs,'Units','n');
    % text(1.0,0.800,['\delta^{13}C_{^{inp}}  = ' cstr '‰ '],'FontSize',...
    %     fs,'Units','n','Hor','r');

    % text(xyl(1),xyl(2),Flb(1),'FontS',ffs,'FontW','b','Units','n','Col','k');
    myH=legend( 'LOSCAR-P','Original LOSCAR');
    set(myH,'FontSize',legendFS);

    subplot(312)
    box on
    hold on
    %     Legend=cell(3,1);%  three positions
    for i=1:1:ind+1
        if(i==1 || i==4)
            h1=plot(tva{i}, CCDa{i},symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
            Legend{i}=lgnd{i} ;
            CCDovershoot = max(CCDa{i})-CCDa{i}(1)
        end
    end
    for i=1:1:ind+1
        if(i==1 || i==4)
            h2=plot(tv11a{i}, CCDa1{i},symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
        end
    end
    hold off
    ylabel('Atlantic CCD (km)', 'FontSize',labelFS);

    text(0.02,0.90,'b)','Units', 'Normalized','FontSize',textFS);
    set(gca, 'YDir','reverse','XLim',axx,'XTickLabel',[], 'fontWeight','bold')
    myH=legend(Legend{1}, Legend{4});
    set(myH,'FontSize',legendFS);
    if(sensFlag == 6)

        subplot(313)
        box on
        hold on
        css = 'bgrmckbbgbmckybgrmcky';
        Ns=13;
        % for l=1:Ns
        %     plot(tva{1}  ,fcAa{1} (:,7)*1e2,'-','Color','b','LineWidth',2);
        plot(tva{1}  ,fcAa{1} (:,10)*1e2,'--','Color','b',LW,linew);
        %     plot(tva{2}  ,fcAa{2} (:,7)*1e2,'-','Color','g','LineWidth',2);
        plot(tva{2}  ,fcAa{2} (:,10)*1e2,'--','Color','g',LW,linew);
        % end;
        % for l=1:Ns
        %     plot(tv11a{1}  ,fcAa1{1} (7,:)*1e2,'-','Color','b','LineWidth',2);
        plot(tv11a{1},fcAa1{1}(10,:)*1e2,'--','Color','b',LW,linew);
        %         plot(tv11a{2}  ,fcAa1{2} (7,:)*1e2,'-','Color','g','LineWidth',2);
        plot(tv11a{2},fcAa1{2}(10,:)*1e2,'--','Color','g',LW,linew);
        % end;
        for l=1:Ns
        end;
        for l=1:Ns
        end;
        hold off
        ylabel('Atlantic CaCO_3 (wt%)','FontSize',labelFS);
        xlabel('Time (ky)','FontSize',labelFS)
        text(0.02,0.90,'c)','Units', 'Normalized','FontSize',textFS);
        dsv =[100 600 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6500];
        dstr(1,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
        dstr(2,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
        dstr(3,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
        dstr(4,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
        Hl=legend(dstr,4);

        set(Hl,'FontSize',legendFS);
        set(gca,'XLim',axx, 'fontWeight','bold','YLim',[0 60])
    elseif(sensFlag==9)

        subplot(313)
        box on
        hold on
        for i=1:1:ind+1
            if(i==1 || i==4)
                h1=plot(tva{i}, fcAa{i}(:,10)*1e2,symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
                Legend{i}=lgnd{i} ;
            end
        end
        for i=1:1:ind+1
            if(i==1 || i==4)
                h2=plot(tv11a{i}, fcAa1{i}(10,:)*1e2,symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
            end
        end
        hold off
        ylabel('Atlantic CaCO_3 at 4500m (wt%)','FontSize',labelFS);
        xlabel('Time (ky)','FontSize',labelFS)
        text(0.02,0.90,'c)','Units', 'Normalized','FontSize',textFS);
        set(gca ,'XLim',axx, 'fontWeight','bold','YLim',[0 60])
    end


    if(sensFlag == 5 || sensFlag==6|| sensFlag==10)
        if(sensFlag == 5 || sensFlag==6)
            dInd = 4;
            clr='g';
        else
            dInd=2;
            clr='k';
        end
        % figure
        subplot(313)
        box on
        hold on
        css = 'bgrmckbbgbmckybgrmcky';
        Ns=13;
        % for l=1:Ns
        plot(ortv./1e3  ,fcA (:,10)*1e2,'-','Color','b','LineWidth',2);
        plot(tva{dInd}  ,fcAa{dInd} (:,10)*1e2,'-','Color',clr,LW,linew);
        plot(pen1403(:,1),pen1403(:,2),'r--',LW,linew);
        % end;
        % for l=1:Ns
        plot(ortv11./1e3  ,fcA1 (10,:)*1e2,'-','Color','b',LW,linew);
        plot(tv11a{dInd},fcAa1{dInd}(10,:)*1e2,'-','Color',clr,LW,linew);

        % end;

        hold off
        ylabel('Atlantic CaCO_3 (wt%)','FontSize',labelFS);
        xlabel('Time (ky)','FontSize',labelFS)
        text(0.02,0.90,'c)','Units', 'Normalized','FontSize',textFS);

        dsv =[100 600 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6500];

        dstr(1,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
        dstr(2,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];

        Hl=legend('Original LOSCAR 4500m','LOSCAR-P 4500m','Data, Site 1403: 4374m');
        set(Hl,'FontSize',legendFS);
        set(gca,'XLim',axx, 'fontWeight','bold','YLim',[0 60])
    end


    figure
    subplot(313)
    box on
    hold on
    for i=1:1:ind+1
        h1=plot(tva{i}, fcAa{i}(:,10)*1e2,symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
        Legend{i}=lgnd{i} ;
    end
    for i=1:1:ind+1
        h2=plot(tv11a{i}, fcAa1{i}(10,:)*1e2,symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
    end
    hold off
    ylabel('Atlantic CaCO_3 at 4500m (wt%)','FontSize',labelFS);
    xlabel('Time (ky)','FontSize',labelFS)
    text(0.02,0.90,'c)','Units', 'Normalized','FontSize',textFS);
    set(gca ,'XLim',axx,'YLim',[0 60])
    set(gca, 'YDir','reverse','XLim',axx,'XTickLabel',[], 'fontWeight','bold')
    myH=legend(Legend);
    set(myH,'FontSize',legendFS);
end

if(sensFlag==10)
    figure
    subplot(211)
    box on
    hold on
    %     Legend=cell(3,1);%  three positions
    for i=1:1:ind+1
        h1=plot(tva{i}, CCDa{i},symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
        Legend{i}=lgnd{i} ;
        CCDovershoot = max(CCDa{i})-CCDa{i}(1)
    end
    for i=1:1:ind+1
        h2=plot(tv11a{i}, CCDa1{i},symb(2*i-1:2*i),'Color',clrs(i));
    end
    hold off
    ylabel('Atlantic CCD (km)','FontSize',labelFS);
    % xlabel('Time (ky)')
    text(0.02,0.90,'a)','Units', 'Normalized','FontSize',textFS);
    set(gca, 'YDir','reverse','XLim',axx,'XTickLabel',[], 'fontWeight','bold')
    legend(Legend,'FontSize',legendFS);
    if(sensFlag == 6)

        subplot(212)
        box on
        hold on
        %     Legend=cell(3,1);%  three positions
        % css = 'bgrmckybgrmckybgrmcky';
        css = 'bgrmckbbgbmckybgrmcky';
        Ns=13;
        for l=1:Ns
            %     if(l==3 || l==4 || l==5 || l==6 || l==7)
            plot(tva{1}  ,fcAa{1} (:,7)*1e2,'-','Color','b',LW,linew);
            plot(tva{1}  ,fcAa{1} (:,10)*1e2,'--','Color','b',LW,linew);
            plot(tva{2}  ,fcAa{2} (:,7)*1e2,'-','Color','g',LW,linew);
            plot(tva{2}  ,fcAa{2} (:,10)*1e2,'--','Color','g',LW,linew);
            %     end
        end;
        for l=1:Ns
            %     if(l==3 || l==4 || l==5 || l==6 || l==7)
            plot(tv11a{1}  ,fcAa1{1} (7,:)*1e2,'-','Color','b',LW,linew);
            plot(tv11a{1},fcAa1{1}(10,:)*1e2,'--','Color','b',LW,linew);
            plot(tv11a{2}  ,fcAa1{2} (7,:)*1e2,'-','Color','g',LW,linew);
            plot(tv11a{2},fcAa1{2}(10,:)*1e2,'--','Color','g',LW,linew);
            %     end
        end;
        for l=1:Ns
            %     if(l==3 || l==4 || l==5 || l==6 || l==7)


            %     end
        end;
        for l=1:Ns
            %     if(l==3 || l==4 || l==5 || l==6 || l==7)

            %     end
        end;
        % plot(pen1409(:,1),pen1409(:,2),'r');
        hold off
        ylabel('Atlantic CaCO_3 (wt%)','FontSize',labelFS);
        xlabel('Time (ky)')
        text(0.02,0.90,'c)','Units', 'Normalized','FontSize',textFS);
        % set(gca, 'YDir','reverse')
        dsv =[100 600 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6500];
        % for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)
        %     dstr(l,:) = sprintf('%4.0f m',dsv(10));
        dstr(1,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
        dstr(2,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
        dstr(3,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
        dstr(4,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
        Hl=legend(dstr,4);

        set(Hl,'FontSize',legendFS);
        set(gca,'XLim',axx, 'fontWeight','bold','YLim',[0 60])
        %    return
    elseif(sensFlag==9)

        subplot(212)
        box on
        %  hold on
        %  for i=1:1:ind+1
        %     h1=plot(tva{i}, omegCa{i}(:,10),symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2);
        %     Legend{i}=lgnd{i} ;
        % end
        % for i=1:1:ind+1
        %     h2=plot(tv11a{i}, omegC1a{i}(10,:),symb(2*i-1:2*i),'Color',clrs(i));
        % end
        % hold off
        % ylabel('Carbonate Saturation (4500m)');
        % xlabel('Time (ky)')
        % text(0.02,0.90,'b)','Units', 'Normalized');
        % set(gca ,'XLim',axx)
        % legend(Legend);
        hold on
        for i=1:1:ind+1
            h1=plot(tva{i}, fcAa{i}(:,10)*1e2,symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
            Legend{i}=lgnd{i} ;
        end
        for i=1:1:ind+1
            h2=plot(tv11a{i}, fcAa1{i}(10,:)*1e2,symb(2*i-1:2*i),'Color',clrs(i),LW,linew);
        end
        hold off
        ylabel('Atlantic CaCO_3 at 4500m (wt%)','FontSize',labelFS);
        xlabel('Time (ky)','FontSize',labelFS)
        text(0.02,0.90,'c)','Units', 'Normalized','FontSize',textFS);
        set(gca ,'XLim',axx, 'fontWeight','bold','YLim',[0 60])
    end


    if(sensFlag == 5 || sensFlag==6|| sensFlag==10)
        if(sensFlag == 5 || sensFlag==6)
            dInd = 4;
            clr='g';
        else
            dInd=2;
            clr='k';
        end
        % figure
        subplot(212)
        box on
        hold on
        %     Legend=cell(3,1);%  three positions
        % css = 'bgrmckybgrmckybgrmcky';
        css = 'bgrmckbbgbmckybgrmcky';
        Ns=13;
        for l=1:Ns
            %     if(l==3 || l==4 || l==5 || l==6 || l==7)
            %     plot(tva{dInd}  ,fcAa{dInd} (:,7)*1e2,'-','Color',clr,'LineWidth',2);
            %     plot(tva{dInd}  ,fcAa{dInd} (:,8)*1e2,':','Color',clr,'LineWidth',2);
            %     plot(tva{dInd}  ,fcAa{dInd} (:,9)*1e2,'-.','Color',clr,'LineWidth',2);
            plot(tva{dInd}  ,fcAa{dInd} (:,10)*1e2,'--','Color',clr,LW,linew);
            %     plot(pen1409(:,1),pen1409(:,2),'r');
            plot(pen1403(:,1),pen1403(:,2),'r--',LW,linew);
            %     end
        end;
        for l=1:Ns
            %     if(l==3 || l==4 || l==5 || l==6 || l==7)
            %     plot(tv11a{dInd}  ,fcAa1{dInd} (7,:)*1e2,'-','Color',clr,'LineWidth',2);
            %     plot(tv11a{dInd},fcAa1{dInd}(8,:)*1e2,':','Color',clr,'LineWidth',2);
            %     plot(tv11a{dInd},fcAa1{dInd}(9,:)*1e2,'-.','Color',clr,'LineWidth',2);
            plot(tv11a{dInd},fcAa1{dInd}(10,:)*1e2,'--','Color',clr,LW,linew);
            %     end
        end;

        hold off
        ylabel('Atlantic CaCO_3 (wt%)','FontSize',labelFS);
        xlabel('Time (ky)','FontSize',labelFS)
        text(0.02,0.90,'b)','Units', 'Normalized','FontSize',textFS);
        % set(gca, 'YDir','reverse')
        dsv =[100 600 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6500];
        % for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)
        %     dstr(l,:) = sprintf('%4.0f m',dsv(10));
        dstr(1,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
        dstr(2,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
        % dstr(1,:) = [dstr(1,:) 'sdfsdf'];
        %     end
        % end;
        set(gca,'XLim',axx, 'fontWeight','bold','YLim',[0 60])
        Hl=legend('LOSCAR-P 4500m','Data, Site 1403: 4374m');
        set(Hl,'FontSize',legendFS);

    end

end

% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOTTING
%CCD Atlantic only figure
figure
subplot(211)
box on
hold on
%     Legend=cell(3,1);%  three positions
for i=1:1:ind+1
    h1=plot(tva{i}, CCDa{i},symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2);
    Legend{i}=lgnd{i} ;
    CCDovershoot = max(CCDa{i})-CCDa{i}(1)
end
for i=1:1:ind+1
    h2=plot(tv11a{i}, CCDa1{i},symb(2*i-1:2*i),'Color',clrs(i));
end
hold off
ylabel('Atlantic CCD (km)');
% xlabel('Time (ky)')
text(0.02,0.90,'a)','Units', 'Normalized');
set(gca, 'YDir','reverse','XLim',axx,'XTickLabel',[])
legend(Legend);
if(sensFlag == 6)

    subplot(212)
    box on
    hold on
    %     Legend=cell(3,1);%  three positions
    % css = 'bgrmckybgrmckybgrmcky';
    css = 'bgrmckbbgbmckybgrmcky';
    Ns=13;
    for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)
        plot(tva{1}  ,fcAa{1} (:,7)*1e2,'-','Color','b','LineWidth',2);
        plot(tva{1}  ,fcAa{1} (:,10)*1e2,'--','Color','b','LineWidth',2);
        plot(tva{2}  ,fcAa{2} (:,7)*1e2,'-','Color','g','LineWidth',2);
        plot(tva{2}  ,fcAa{2} (:,10)*1e2,'--','Color','g','LineWidth',2);
        %     end
    end;
    for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)
        plot(tv11a{1}  ,fcAa1{1} (7,:)*1e2,'-','Color','b','LineWidth',2);
        plot(tv11a{1},fcAa1{1}(10,:)*1e2,'--','Color','b','LineWidth',2);
        plot(tv11a{2}  ,fcAa1{2} (7,:)*1e2,'-','Color','g','LineWidth',2);
        plot(tv11a{2},fcAa1{2}(10,:)*1e2,'--','Color','g','LineWidth',2);
        %     end
    end;
    for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)


        %     end
    end;
    for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)

        %     end
    end;
    % plot(pen1409(:,1),pen1409(:,2),'r');
    hold off
    ylabel('Atlantic CaCO_3 (wt%)');
    xlabel('Time (ky)')
    text(0.02,0.90,'b)','Units', 'Normalized');
    % set(gca, 'YDir','reverse')
    dsv =[100 600 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6500];
    % for l=1:Ns
    %     if(l==3 || l==4 || l==5 || l==6 || l==7)
    %     dstr(l,:) = sprintf('%4.0f m',dsv(10));
    dstr(1,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
    dstr(2,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
    dstr(3,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
    dstr(4,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
    Hl=legend(dstr,4);

    set(Hl,'FontSize',10);
    set(gca,'XLim',axx)
    %    return
elseif(sensFlag==9)

    subplot(212)
    box on
    %  hold on
    %  for i=1:1:ind+1
    %     h1=plot(tva{i}, omegCa{i}(:,10),symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2);
    %     Legend{i}=lgnd{i} ;
    % end
    % for i=1:1:ind+1
    %     h2=plot(tv11a{i}, omegC1a{i}(10,:),symb(2*i-1:2*i),'Color',clrs(i));
    % end
    % hold off
    % ylabel('Carbonate Saturation (4500m)');
    % xlabel('Time (ky)')
    % text(0.02,0.90,'b)','Units', 'Normalized');
    % set(gca ,'XLim',axx)
    % legend(Legend);
    
    hold on
    for i=1:1:ind+1
        h1=plot(tva{i}, fcAa{i}(:,10)*1e2,symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2);
        Legend{i}=lgnd{i} ;
    end
    for i=1:1:ind+1
        h2=plot(tv11a{i}, fcAa1{i}(10,:)*1e2,symb(2*i-1:2*i),'Color',clrs(i));
    end
    hold off
    ylabel('Atlantic CaCO_3 at 4500m (wt%)');
    xlabel('Time (ky)')
    text(0.02,0.90,'b)','Units', 'Normalized');
    set(gca ,'XLim',axx)
end


if(sensFlag == 5 || sensFlag==6|| sensFlag==10)
    if(sensFlag == 5 || sensFlag==6)
        dInd = 4;
        clr='g';
    else
        dInd=2;
        clr='k';
    end
    
    % figure
    subplot(212)
    box on
    hold on
    %     Legend=cell(3,1);%  three positions
    % css = 'bgrmckybgrmckybgrmcky';
    css = 'bgrmckbbgbmckybgrmcky';
    
    Ns=13;
    for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)
        plot(tva{dInd}  ,fcAa{dInd} (:,7)*1e2,'-','Color',clr,'LineWidth',2);
        plot(tva{dInd}  ,fcAa{dInd} (:,8)*1e2,'--','Color',clr,'LineWidth',2);
        plot(tva{dInd}  ,fcAa{dInd} (:,9)*1e2,'-.','Color',clr,'LineWidth',2);
        plot(tva{dInd}  ,fcAa{dInd} (:,10)*1e2,':','Color',clr,'LineWidth',2);
        plot(pen1409(:,1),pen1409(:,2),'r');
        plot(pen1403(:,1),pen1403(:,2),'r--');
        %     end
    end;
    for l=1:Ns
        %     if(l==3 || l==4 || l==5 || l==6 || l==7)
        plot(tv11a{dInd}  ,fcAa1{dInd} (7,:)*1e2,'-','Color',clr,'LineWidth',2);
        plot(tv11a{dInd},fcAa1{dInd}(8,:)*1e2,'--','Color',clr,'LineWidth',2);
        plot(tv11a{dInd},fcAa1{dInd}(9,:)*1e2,'-.','Color',clr,'LineWidth',2);
        plot(tv11a{dInd},fcAa1{dInd}(10,:)*1e2,':','Color',clr,'LineWidth',2);
        %     end
    end;

    hold off
    ylabel('Atlantic CaCO_3 (wt%)');
    xlabel('Time (ky)')
    text(0.02,0.90,'b)','Units', 'Normalized');
    % set(gca, 'YDir','reverse')
    dsv =[100 600 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6500];
    % for l=1:Ns
    %     if(l==3 || l==4 || l==5 || l==6 || l==7)
    %     dstr(l,:) = sprintf('%4.0f m',dsv(10));
    dstr(1,:) = [sprintf('%4.0f m',dsv(7)) ' Site U1409'];
    dstr(2,:) = [sprintf('%4.0f m',dsv(10)) ' Site U1403'];
    % dstr(1,:) = [dstr(1,:) 'sdfsdf'];
    %     end
    % end;
    Hl=legend('LOSCAR-P 3000m', 'LOSCAR-P 3500m', 'LOSCAR-P 4000m', 'LOSCAR-P 4500m','Data, Site 1409: 2913m','Data, Site 1403: 4374m');
    set(Hl,'FontSize',10);
    set(gca,'XLim',axx)
end


















FigHandle=figure;
set(FigHandle, 'Position', [400, 100, 500, 950]);
% subplot(511)
% box on;
% hold on;
% fill(trl2,rl2,'g ');
% fill(trl,rl,'r ');
%
% hold off;
% set(gca,'FontSize',fs);
% set(gca,'XLim',axx.*1e3);
% set(gca,'XTickLabel',[]);
% set(gca,'YLim',[0 .7]);
% % set(gca,'XTickLabel',[]);
% % set(gca,'XTickL',[]);
% ylabel({'Carbon Input';' (Pg C y^{-1})'});
% text(0.02,0.90,'a)','Units', 'Normalized');
% annotation('arrow',[.37 .32],[.905 .905]);
% annotation('arrow',[.54 .49],[.86 .86]);
subplot(411)
box on
hold on
mark = 'o x d o x d o x d ';
for i=1:1:ind+1
    h1=plot(tva{i}, pco2ta{i},[symb(2*i-1:2*i)],'Color',clrs(i),'LineWidth',2);
    Legend{i}=lgnd{i} ;
end
for i=1:1:ind+1
    h2=plot(tv11a{i}, pco2ta1{i},symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2);
end
hold off

ylabel('pCO_2 (\muatm)');
set(gca,'XTickLabel',[],'XLim',axx)
% columnlegend(2,Legend);
legend(Legend);
if(sensFlag==0)
%     title('No P burial vs P')
elseif(sensFlag==1)
    title('org P and C fraction buried Sensitivity runs')
elseif(sensFlag==2)
    title('Oxygen sensitivity runs')
elseif(sensFlag==2)
    title('Floegel runs')
end
text(0.02,0.90,'a)','Units', 'Normalized');
subplot(412)
box on
hold on
for i=1:1:ind+1
    % scaling d13C data so they have the same initial value
    d13cmean{i}=d13cmean{i}-d13cmean{i}(1)+d13cmean{1}(1);
    plot(tva{i}, d13cmean{i},symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2)
end
for i=1:1:ind+1
    % scaling d13C data so they have the same initial value
    d13cmean1{i}=d13cmean1{i}-d13cmean1{i}(1)+d13cmean{1}(1);
    plot(tv11a{i}, d13cmean1{i},symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2)
end
% plot(age3,     c13c3     ,'b--s',MFC,'b',MS,ms);
% plot(age690blk,c13c690blk,'r- d',MFC,'r',MS,ms);
% plot(age1266blk,d13c1266,'g- d',MFC,'r',MS,ms);
% plot(Fage690,c13c690blk,'k- d',MFC,'k',MS,ms);
hold off
ylabel('Mean surface \delta^{13}C (‰)');
text(0.02,0.90,'b)','Units', 'Normalized');
set(gca,'XTickLabel',[],'XLim',axx);

subplot(413)
box on
hold on
for i=1:1:ind+1
    plot(tva{i}, orgCba{i}/1e12,symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2)
end
for i=1:1:ind+1
    plot(tv11a{i}, orgCb1a{i}/1e12,symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2)
end
hold off

ylabel('F_{OCB} (10^{12} mol yr^{-1})')
set(gca,'XTickLabel',[],'XLim',axx)
%     text(0.8,1,'Total Pg C over 200ky','Units','normalized');
%     text(0.8,1-(i*0.2), num2str( excessCga{i}),'Color',clrs(i),'Units','normalized')
text(0.02,0.90,'c)','Units', 'Normalized');
subplot(414)
box on
hold on
for i=1:1:ind+1
    plot(tva{i}, totPbura{i}/1e10,symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2)
end
for i=1:1:ind+1
    plot(tv11a{i}, totPbur1a{i}/1e10,symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2)
end
hold off
text(0.02,0.90,'d)','Units', 'Normalized');
set(gca,'XLim',axx);
ylabel({'Total P burial';' (10^{10} mol yr^{-1})'})
xlabel('Time (ky)')

return
if(sensFlag==1)
    figure
    subplot(311)
    hold on
    box on
    for i=1:1:ind+1
        plot(orgCba{i,1}(1)./1e12, excessCga{i},'MarkerFaceColor',clrs(i),'MarkerEdgeColor',clrs(i),...
            'MarkerSize',20,...
            'Marker',mark(2*i-1:2*i),...
            'LineStyle','none')
    end
    ylabel('Excess C buried (Pg C)')
    xlabel('Initial F_{OCB} (10^{12} mol yr^{-1})')
    text(0.02,0.90,'a)','Units', 'Normalized');
    hold off
    subplot(312)
    % figure
    box on
    hold on
    for i=1:1:ind+1
        plot(totPbura{i,1}(1)./1e10, excessCga{i},'MarkerFaceColor',clrs(i),'MarkerEdgeColor',clrs(i),...
            'MarkerSize',20,...
            'Marker',mark(2*i-1:2*i),...
            'LineStyle','none')
    end
    ylabel('Excess C buried (Pg C)')
    xlabel('Initial Total P burial (10^{10} mol yr^{-1})')
    text(0.02,0.90,'b)','Units', 'Normalized');
    hold off

    Ndecimals = 2;
    f = 10.^Ndecimals;
    count = 1;
    for j=1:3
        X(j)=round(((1-(eI(1)+oIva{count}(1)))*100)*f)/f;
        for k=1:3

            Y(k)=round(((1-(eI(1)+oIpva{count}(1)))*100)*f)/f;
            Z(j,k) = excessCga{count}(1);
            count=count+1;
        end
    end
    % Sort Org C and org P fractions and their respective excess org C burial
    % so that we can contour plot it
    sortedX=sort(X);
    sortedY=sort(Y);
    for j=1:3
        for k=1:3
            sortedZ(k,j) = Z(find(X==sortedX(j)), find(Y==sortedY(k)));
        end
    end
    % figure
    subplot(313)
    [C,h] = contourf(sortedX,sortedY,sortedZ);
    set(h,'ShowText','on','Color',[1 1 1],'LineStyle','--')
    colormap cool
    xlabel('f_{OC} (%)')
    ylabel('f_{OP} (%)')
    text(0.02,0.90,'c)','Units', 'Normalized');
    % Excess C burial as a function of initial P and C fluxes
    % figure
    % X1=[0.5 1.0 2.0];
    % Y1 = [0.25 0.5 1.0];
    % Z1 = [861.7 1361 2129;1227 1839 2815; 1499 2200 2840];
    %
    % [C,h] = contourf(X1,Y1,Z1);
    % set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    % colormap cool
    % xlabel('f_{OC} (%)')
    % ylabel('f_{OP} (%)')
end
if(sensFlag == 0)
    %PCO2 only figure
    % figure(2)
    figure
    box on
    hold on
    Legend=cell(3,1);%  three positions
    for i=1:1:ind+1
        h1=plot(tva{i}, pco2ta{i},symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2);
        Legend{i}=lgnd{i} ;
    end
    for i=1:1:ind+1
        h2=plot(tv11a{i}, pco2ta1{i},symb(2*i-1:2*i),'Color',clrs(i));
    end
    hold off
    ylabel('pCO_2 (\muatm)');
    xlabel('Time (ky)')
    legend(Legend);

    %CCD Atlantic only figure
    figure
    box on
    hold on
    Legend=cell(3,1);%  three positions
    for i=1:1:ind+1
        h1=plot(tva{i}, CCDa{i},symb(2*i-1:2*i),'Color',clrs(i),'LineWidth',2);
        Legend{i}=lgnd{i} ;
    end
    for i=1:1:ind+1
        h2=plot(tv11a{i}, CCDa1{i},symb(2*i-1:2*i),'Color',clrs(i));
    end
    hold off
    ylabel('Atlantic (km)');
    xlabel('Time (ky)')
    set(gca, 'YDir','reverse')
    legend(Legend);

    % Oxygen only for Standard P model
    for k=1:Nb
        lstr(k,:) = sprintf('%s',lstr0(2*k-1:2*k));
    end;



    % Phosphate only for Standard P model
    figure
    subplot(211)
    %     clf;
    box  on;
    hold on;
    for k=1:Nb
        plot(tva{3}  ,pa {3}(:,k),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',2);
    end;
    for k=1:Nb
        plot(tv11a{3},pa1{3}(k,:),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',2);
    end;
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %     xlabel('Time (ky)');
    ylabel('PO_4 (\mumol kg^{-1})');
    Hl=legend(lstr,4);
    set(Hl,'FontSize',10);
    set(gca,'XTickLabel',[]);
    text(0.02,0.90,'a)','Units', 'Normalized');
    % Export
    subplot(212)
    %     figure
    %     clf;
    box  on;
    hold on;
    plot(tva{3}  ,(EPLvva {3} +EPHa{3})./1e12,'g','LineWidth',2);
    plot(tv11a{3}  ,(EPLvva1 {3}+EPHa1{3})./1e12,'g','LineWidth',2);

    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    xlabel('Time (ky)');
    ylabel('Global export production (10^{12} mol C)');
    Hl=legend('Global ocean production');
    set(Hl,'FontSize',10);
    text(0.02,0.90,'b)','Units', 'Normalized');



    figure
    subplot(211)
    %     clf;
    box  on;
    hold on;
    for k=1:Nb
        plot(tva{3}  ,doxa {3}(:,k),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',2);
    end;
    for k=1:Nb
        plot(tv11a{3},dox1a{3}(k,:),sstr(2*k-1:2*k),'Color',cs(k),'LineWidth',2);
    end;
    hold off;
    text(0.02,0.90,'a)','Units', 'Normalized');
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %     xlabel('Time (ky)');
    ylabel('O_2 (mol m^{-3})');
    Hl=legend(lstr,4);
    set(Hl,'FontSize',10);
    set(gca,'XTickLabel',[])

    %     figure
    subplot(212)
    box on
    hold on
    % Org. P burial
    plot(tv./1e3, totPbur./1e10,'.k','LineWidth',2);
    plot(tv./1e3, Fcapv./1e10,'--b','LineWidth',2);
    plot(tv./1e3, orgPb./1e10,'g','LineWidth',2);
    plot(tv./1e3, Ffepv./1e10,'-.r','LineWidth',2);

    plot(tv11./1e3, totPbur1./1e10,'-.k','LineWidth',2);
    plot(tv11./1e3, Fcapv1./1e10,'--b','LineWidth',2);
    plot(tv11./1e3, orgPb1./1e10,'g','LineWidth',2);
    plot(tv11./1e3, Ffepv1./1e10,'-.r','LineWidth',2);

    hold off
    set(gca,'FontSize',fs,'XLim',axx);
    text(0.02,0.90,'b)','Units', 'Normalized');
    Hl=legend('Total P burial','F_{CaP}','F_{OPB}','F_{FeP}',4);
    set(Hl,'FontSize',10);
    xlabel('Time (ky)');
    ylabel('P burial fluxes (10^{10} mol P y^{-1})');

end


return
figure(1)
subplot(311)
hold on
phndl = plot(tva{1}, pco2ta{1},'r',tva{2}, pco2ta{2},'b--',tva{3}, pco2ta{3},'k.-');
legend('W/out cycle', 'With P cycle','With P cycle + P Weath.')
phndl2=plot(tv11a{1}, pco2ta1{1},'r',tv11a{2}, pco2ta1{2},'b--',tv11a{3}, pco2ta1{3},'k.-');
hold off

ylabel('pCO_2 (\muatm)');
set(gca,'XTickLabel',[])
subplot(312)
hold on
phnd3 = plot(tva{1}, d13cmean{1},'r',tva{2}, d13cmean{2},'b--',tva{3}, d13cmean{3},'k.-');
phnd4 = plot(tv11a{1}, d13cmean1{1},'r',tv11a{2}, d13cmean1{2},'b--',tv11a{3}, d13cmean1{3},'k.-');
hold off
ylabel('Mean surface \delta^{13}C');
set(gca,'XTickLabel',[]);
subplot(313)
hold on
plot(tva{1}, orgCba{1},'r', tva{2}, orgCba{2},'b--', tva{3}, orgCba{3},'k.-')
plot(tv11a{1}, orgCb1a{1},'r',tv11a{2}, orgCb1a{2},'b--',tv11a{3}, orgCb1a{3},'k.-');
text(200,4.5e12, num2str( excessCga{1}),'Color',[1 0 0])
text(150,4.9e12, num2str( excessCga{2}),'Color',[0 0 1])
text(100,5.2e12, num2str( excessCga{3}),'Color',[0 0 0])
% set(gca,'Color',[0 1 1]);
hold off
xlabel('Time (ky)')
ylabel('Corgb (mol yr^{-1})')

return

for k=1:Nb
    lstr(k,:) = sprintf('%s',lstr0(2*k-1:2*k));
end;
axx = [-0.5e5 400000];
figure
clf;
box  on;
hold on;
% dir = 'Poff/';
% tv= load([dir 'tv.DAT']);
% tv11= load([dir 'tv11.DAT']);
% p= load([dir 'p.DAT']);
% p1     = [p(1,:)' p(1,:)'];
for k=1:Nb
    plot(tv  ,p (:,k)*1e3,sstr(2*k-1:2*k),'Color',cs(k));
end;
for k=1:Nb
    plot(tv11,p1(k,:)*1e3,sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (ky)');
ylabel('PO_4 (\mumol kg^{-1})');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);


figure
clf;
box  on;
hold on;
% dir = 'Poff/';
% tv= load([dir 'tv.DAT']);
% tv11= load([dir 'tv11.DAT']);
% p= load([dir 'p.DAT']);
% p1     = [p(1,:)' p(1,:)'];
for k=1:Nb
    plot(tv  ,T (:,k),sstr(2*k-1:2*k),'Color',cs(k));
end;
% for k=1:Nb
%     plot(tv11,p1(k,:)*1e3,sstr(2*k-1:2*k),'Color',cs(k));
% end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('Temperature (C^o)');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);


figure(311)
hold on
% Org. P burial
plot(tv, orgPb,'g');
plot(tv, Fcapv,'--b');
plot(tv, Ffepv,'-.r');
plot(tv, totPbur,'k');

plot(tv11, orgPb1,'g');
plot(tv11, Fcapv1,'--b');
plot(tv11, Ffepv1,'-.r');
plot(tv11, totPbur1,'k');
hold off
set(gca,'FontSize',fs);
Hl=legend('F_{OPB}','F_{CaP}','F_{FeP}','Total P burial',4);
set(Hl,'FontSize',10);
xlabel('Time (y)');
ylabel('P burial fluxes (mol P y^{-1})');

figure(32)
hold on
% Org. C burial
plot(tv, orgCb,'r');
plot(tv11, orgCb1,'r');
hold off
set(gca,'FontSize',fs);
%     Hl=legend('F_{OCB}');
%     set(Hl,'FontSize',10);
xlabel('Time (y)');
ylabel('Org. C burial (mol C y^{-1})');

figure(33)
hold on
% Org. C burial
plot(tv, CPbur,'r');
plot(tv11, CPbur1,'r');
hold off
set(gca,'FontSize',fs);
%     Hl=legend('F_{OCB}');
%     set(Hl,'FontSize',10);
xlabel('Time (y)');
ylabel('C:P ratio in deep sea sediments');

figure(35)
clf;
box  on;
hold on;
for k=1:Nb
    plot(tv  ,dox (:,k),sstr(2*k-1:2*k),'Color',cs(k));
end;
for k=1:Nb
    plot(tv11,dox1(k,:),sstr(2*k-1:2*k),'Color',cs(k));
end;
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
xlabel('Time (y)');
ylabel('O_2 (mol m^{-3})');
Hl=legend(lstr,4);
set(Hl,'FontSize',10);
