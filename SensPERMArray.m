clear all
close all
pd =  csvread('PayneData.csv');
        stel = pd(:,1); % stratigraphic elevation in meters
        d44cap = pd(:,2);
        stdp = pd(:,3); % standard deviation of calcium data
        d13cp = pd(:,4);
        
        depth = [0 160]; % used for calculating time
        for i=1:length(stel)
            time(i) = interp1(depth, [0 1] , stel(i));
        end
%%%
axx = [-0.5e5 4e5];
cs = 'gggkkkrrrbgkr';
sstr = '- ---.- ---.- ---. -: : : ';
fs =10;
ftys =1;
CAvflag = 2;
inorgF = 1;

% dir = 'dat/FinalPERMsims/Sens1/a/'; %Sim1/, Sim2/, Sim3/, BiopSim/
currentRun = 1; % whether to plot current run or saved files, 0 or 1
figurehndl=42;
if(currentRun)
    figurehndl=52;
%%% load all the files
for i=1:1:3

    if(i==1)
        dir = 'dat/FinalPERMsims/Sens1/a/';
    elseif(i==2)
        dir = 'dat/FinalPERMsims/Sens1/b/';
    elseif(i==3)
        dir = 'dat/FinalPERMsims/Sens1/c/';
    end
tv{i,:}= load([dir 'tv.DAT']);
tv11{i,:}= load([dir 'tv11.DAT']);
pco2t{i,:}= load([dir 'pco2t.DAT']);

FSichck{i,:}= load([dir 'FSi.DAT']);
Finchck{i,:}= load([dir 'FiN.DAT']);
FSichck1{i,:} = load([dir 'FSi1.DAT']);
Finchck1{i,:}= load([dir 'FiN1.DAT']);

a{i,:} = load([dir 'a.DAT']);
c{i,:}= load([dir 'c.DAT']);
V{i,:}= load([dir 'V.DAT']);
co3tv{i,:} = load([dir 'co3tv.DAT']);

d13fcA{i,:} = load([dir 'd13fcA.DAT']);
d13fcI{i,:}= load([dir 'd13fcI.DAT']);
d13fcP{i,:}= load([dir 'd13fcP.DAT']);
d13fcT{i,:}= load([dir 'd13fcT.DAT']);

BurialCA{i,:}= load([dir 'BurialCA.DAT']);
BurialCI{i,:}= load([dir 'BurialCI.DAT']);
BurialCP{i,:}= load([dir 'BurialCP.DAT']);
BurialCT{i,:}= load([dir 'BurialCT.DAT']);

d44fcA{i,:}= load([dir 'd44fcA.DAT']);
d44fcI{i,:}= load([dir 'd44fcI.DAT']);
d44fcP{i,:}= load([dir 'd44fcP.DAT']);
d44fcT{i,:}= load([dir 'd44fcT.DAT']);

d44ca{i,:} = load([dir 'd44ca.DAT']);

ccdA{i,:}= load([dir 'ccdA.DAT']);
ccdI{i,:}= load([dir 'ccdI.DAT']);
ccdP{i,:}= load([dir 'ccdP.DAT']);
ccdT{i,:}= load([dir 'ccdT.DAT']);

omegCSvt{i,:}= load([dir 'OmegCS.DAT']);
CA{i,:}= load([dir 'CA.DAT']);
epspcaV{i,:}= load([dir 'epspcaV.DAT']);
pco2t1{i,:}  = [pco2t{i,:}(1,1)' pco2t{i,:}(1,1)'];
epspcaV1{i,:}  = [epspcaV{i,:}(1,1)' epspcaV{i,:}(1,1)'];
% MEAN SURFACE OCEAN carbonate saturation
oMean{i,:}=(omegCSvt{i,:}(:,1).*V{i,:}(1)+omegCSvt{i,:}(:,2).*V{i,:}(2)+...
    omegCSvt{i,:}(:,3)*V{i,:}(3)+omegCSvt{i,:}(:,11)*V{i,:}(11))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(11)));
    % MEAN SURFACE OCEAN alkalinity and total carbon over time
      % MEAN SURFACE OCEAN carbonate ion conc. over time
       % MEAN SURFACE OCEAN calcium ion conc. over time
    if(ftys)
        aMean{i,:}=(a{i,:}(:,1).*V{i,:}(1)+a{i,:}(:,2).*V{i,:}(2)+a{i,:}(:,3)*V{i,:}(3)+a{i,:}(:,10)*V{i,:}(10)+a{i,:}(:,11)*V{i,:}(11))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(10)+V{i,:}(11)));
        cMean{i,:}=(c{i,:}(:,1).*V{i,:}(1)+c{i,:}(:,2).*V{i,:}(2)+c{i,:}(:,3)*V{i,:}(3)+c{i,:}(:,10)*V{i,:}(10)+c{i,:}(:,11)*V{i,:}(11))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(10)+V{i,:}(11)));
        co3Mean{i,:}=(co3tv{i,:}(:,1).*V{i,:}(1)+co3tv{i,:}(:,2).*V{i,:}(2)+co3tv{i,:}(:,3)*V{i,:}(3)+co3tv{i,:}(:,10)*V{i,:}(10)+co3tv{i,:}(:,11)*V{i,:}(11))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(10)+V{i,:}(11)));
         d13cbulkV{i,:} = ((d13fcA{i,:}(:,1)*BurialCA{i,:}+d13fcP{i,:}(:,1)*BurialCP{i,:}+d13fcI{i,:}(:,1)*BurialCI{i,:}+d13fcT{i,:}(:,1)*BurialCT{i,:})./(BurialCT{i,:}+BurialCP{i,:}+BurialCI{i,:}+BurialCA{i,:})+...
         (d13fcA{i,:}(:,2)*BurialCA{i,:}+d13fcP{i,:}(:,2)*BurialCP{i,:}+d13fcI{i,:}(:,2)*BurialCI{i,:}+d13fcT{i,:}(:,2)*BurialCT{i,:})./(BurialCT{i,:}+BurialCP{i,:}+BurialCI{i,:}+BurialCA{i,:}))/2;
            d44cabulkV{i,:} = ((d44fcA{i,:}(:,1)*BurialCA{i,:}+d44fcP{i,:}(:,1)*BurialCP{i,:}+d44fcI{i,:}(:,1)*BurialCI{i,:}+d44fcT{i,:}(:,1)*BurialCT{i,:})./(BurialCT{i,:}+BurialCP{i,:}+BurialCI{i,:}+BurialCA{i,:})+...
         (d44fcA{i,:}(:,2)*BurialCA{i,:}+d44fcP{i,:}(:,2)*BurialCP{i,:}+d44fcI{i,:}(:,2)*BurialCI{i,:}+d44fcT{i,:}(:,2)*BurialCT{i,:})./(BurialCT{i,:}+BurialCP{i,:}+BurialCI{i,:}+BurialCA{i,:}))/2;  
        caMean{i,:}=(CA{i,:}(:,1).*V{i,:}(1)+CA{i,:}(:,2).*V{i,:}(2)+CA{i,:}(:,3)*V{i,:}(3)+CA{i,:}(:,10)*V{i,:}(10)+CA{i,:}(:,11)*V{i,:}(11))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(10)+V{i,:}(11)));
        d44caMean{i,:}=(d44ca{i,:}(:,1).*V{i,:}(1)+d44ca{i,:}(:,2).*V{i,:}(2)+d44ca{i,:}(:,3)*V{i,:}(3)+d44ca{i,:}(:,10)*V{i,:}(10)+d44ca{i,:}(:,11)*V{i,:}(11))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(10)+V{i,:}(11)));
    else
       aMean{i,:}=(a{i,:}(:,1).*V{i,:}(1)+a{i,:}(:,2).*V{i,:}(2)+a{i,:}(:,3)*V{i,:}(3)+a{i,:}(:,10)*V{i,:}(10))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(10)));
       cMean{i,:}=(c{i,:}(:,1).*V{i,:}(1)+c{i,:}(:,2).*V{i,:}(2)+c{i,:}(:,3)*V{i,:}(3)+c{i,:}(:,10)*V{i,:}(10))./((V{i,:}(1)+V{i,:}(2)+V{i,:}(3)+V{i,:}(10)));     
    end
    
    aMean1{i,:} = [aMean{i,:}(1,1)'  aMean{i,:}(1,1)'];
    cMean1{i,:} = [cMean{i,:}(1,1)'  cMean{i,:}(1,1)'];
    co3Mean1{i,:} = [co3Mean{i,:}(1,1)'  co3Mean{i,:}(1,1)'];
    co3Mean{i,:} = [co3Mean1{i,:} co3Mean{i,:}']';
    oMean1{i,:} = [oMean{i,:}(1,1)'  oMean{i,:}(1,1)'];
    oMean{i,:} = [oMean1{i,:} oMean{i,:}']';
    ctv{i,:}=[tv11{i,:} tv{i,:}']';
    epspcaVf{i,:} = [epspcaV1{i,:} epspcaV{i,:}]';
    d13cbulkV1{i,:} = [d13cbulkV{i,:}(1,1)'  d13cbulkV{i,:}(1,1)'];
    d44cabulkV1{i,:} = [d44cabulkV{i,:}(1,1)'  d44cabulkV{i,:}(1,1)'];

    d44cabulkVf{i,:} = [d44cabulkV1{i,:}  d44cabulkV{i,:}'];
    
     caMean1{i,:} = [caMean{i,:}(1:1)'  caMean{i,:}(1:1)'];
     d44caMean1{i,:} = [d44caMean{i,:}(1:1)'  d44caMean{i,:}(1:1)'];

            

      

end
end







%%% Plot
% figure(42)
FigHandle = figure(figurehndl);
  set(FigHandle, 'Position', [200, 0, 600, 1000]);
    clf;
%     X=[tv{1},fliplr(tv{1})];
%     Y=[pco2t{1},fliplr(pco2t{2})];
    subplot(711);
    plot(tv{1}  ,pco2t{1} ,'r-','LineWidth',2);
    hold on;
%     plot(tv{2}  ,pco2t{2} ,'r--','LineWidth',1);
%     plot(tv{3}  ,pco2t{3} ,'r.','LineWidth',1);

    plot(tv11{1},pco2t1{1},'r-','LineWidth',2);
    hold off;
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('pCO_2 (\muatm)');

    text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    title('LOSCAR Simulation: 12,000 Pg C, Strangelove')
    subplot(712)
    box on
    hold on
    hndl=plot(tv{1}, FSichck{1}/.1e13,'m',tv{1},Finchck{1}/.1e13,'g');
%     hndl2=plot(tv{2}, FSichck{2}/.1e13,'m--',tv{2},Finchck{2}/.1e13,'g--');
%     hndl3=plot(tv{3}, FSichck{3}/.1e13,'m.',tv{3},Finchck{3}/.1e13,'g.');
    set(hndl,'LineWidth',1)
%      set(hndl3,'LineWidth',1,'markersize',3)
    legend('F_S_i','F_C')
    hndl1=plot(tv11{1}, FSichck1{1}(9,:)/.1e13,'m',tv11{1},Finchck1{1}/.1e13,'g');
    set(hndl1,'LineWidth',1)

    hold off
    ylabel({'F_C, F_S_i';' (10^{12} mol C y^{-1})'})
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'YLim',[3 15]);
    set(gca,'XTickLabel',[])
    text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')


    subplot(713)
    box  on;
    hold on;
%     hndl1=plot(tv  ,a  (:,3),'-b',tv  ,c  (:,3),'-r');
    hndl1=plot(tv{1}  ,aMean{1},'-b',tv{1}  ,cMean{1},'-r');
%     hndl2=plot(tv{2}  ,aMean{2},'--b',tv{2}  ,cMean{2},'--r');
%     hndl4=plot(tv{3}  ,aMean{3},'.b',tv{3}  ,cMean{3},'.r');
    
    legend('TA','TC ')
%     hndl=plot(tv11,a1 (3,:),'-b',tv11,c1 (3,:),'-r');
    hndl=plot(tv11{1},aMean1{1},'-b',tv11{1},cMean1{1},'-r');
%     hndl3=plot(tv11{2}  ,aMean1{2},'-b',tv11{2}  ,cMean1{2},'-r');
%     hndl5=plot(tv11{3}  ,aMean1{3},'-b',tv11{3}  ,cMean1{3},'-r');
    set(hndl,'LineWidth',1)
    %legend('TA')
    hold off;
    set(hndl1,'LineWidth',1)
%     set(hndl2,'LineWidth',1)
%     set(hndl3,'LineWidth',1)
%     set(hndl5,'LineWidth',1)
%     set(hndl4,'LineWidth',1,'markersize',3)
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx,'YLim',[1 5]);
    set(gca,'XTickLabel',[])
    % set(gca,'XTickLabel',[])
    % plot(tv,TCAR*1e3,'k',tv,TALK*1e3,'r');
    ylabel({'TA, TC';' (mmol kg^{-1})'})
    text(0.02,0.98,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
  
    subplot(714);
   
    box  on;
    hold on;

    [ax,p1,p2] = plotyy(ctv{1},co3Mean{1}.*1e6,ctv{1},oMean{1},'plot','plot');
%       [ax2,p3,p4] = plotyy(ctv{2},co3Mean{2}.*1e6,ctv{2},oMean{2},'plot','plot');
%     plot(tv  ,co3Mean.*1e6,sstr(2*3-1:2*3),'Color',cs(3));
%     plot(tv11,co3Mean1*1e6,sstr(2*3-1:2*3),'Color',cs(3));
    %end;
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
%     ylabel('[CO_3^{2-}] (\mumol kg^{-1})');
    ylabel(ax(1),{'[CO_3^{2-}]';' (\mumol kg^{-1})'}) % label left y-axis
    ylabel(ax(2),'\Omega') % label right y-axis
    p1.LineWidth = 2;
    p2.LineWidth = 2; 
%     p3.LineWidth = 1;
%     p4.LineWidth = 1; 
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    legend('mean surface ocean','mean surface ocean')
    set(gca,'XTickLabel',[])
    text(0.02,0.98,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
   
 
    subplot(715);
    box on;
    hold on;
   myh1= plot(tv{1}, d13cbulkV{1},'-b','LineWidth',2 );
   myh2= plot(tv11{1}, d13cbulkV1{1},'-b','LineWidth',2 );
   myh3= plot((time-0.26)*1e6,d13cp,'o','markersize', 3);
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('\delta^{13}C_{carb} (‰)');
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    legend([myh1 myh3], 'Mean shallow', 'data - Payne et al. 2010')
    set(gca,'XTickLabel',[])
    
 text(0.02,0.98,'e)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    if(CAvflag~=0)
        subplot(716);
        box on;
%         plot(tv  ,d44ca (:,9),sstr(2*9-1:2*9),'Color',cs(9));
        hold on
        [ax,p1,p2] = plotyy(ctv{1},d44cabulkVf{1},ctv{1},epspcaVf{1},'plot','plot');
%         [ax2,p3,p4] = plotyy(tv11{1},d44cabulkV1{1},tv11{1},epspcaV1{1},'plot','plot');
        ylabel(ax(1),{'[CO_3^{2-}]';' (\mumol kg^{-1})'}) % label left y-axis
    ylabel(ax(2),'\epsilon_{carb}') % label right y-axis
    p1.LineWidth = 2;
    %Line 2 color
    color2= [204/255 51/255 255/255];
    p2.LineWidth = 2;
    p2.Color = color2;
    set(ax(2), 'YColor',color2 );
%         myh1 = plot(tv{1}, d44cabulkV{1},'-b','LineWidth',2);
%         myh2 = plot(tv11{1},d44cabulkV1{1},'-b', 'LineWidth', 2);
        myh3 = plot((time-0.26)*1e6,d44cap,'o','Color','r','markersize', 3);
        h=errorbar((time-0.26)*1e6,d44cap,stdp,'b'); 
        set(h,'linestyle','none','LineWidth',0.5, 'Color',[1 0 0])
%         plot(tv11,d44ca1(9,:),sstr(2*9-1:2*9),'Color',cs(9));
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx,'YLim',[-1 -0.4]);
    ylabel('Ca^2^+ (mmol kg^{-1})');
    legend([p1 p2 myh3],'Mean shallow','Calculated \epsilon_{carb}','data - Payne et al. 2010')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    set(gca,'XTickLabel',[])
    ylabel('\delta^{44} Ca_{carb} (‰)')
    text(0.02,0.98,'f)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
    

    if(CAvflag~=0)
        
        
        
        subplot(717);
%         plot(tv  ,CA  (:,3),'-b','LineWidth',2);
        plot(tv{1}  ,caMean{1},'-b','LineWidth',2);
        hold on
%         plot(tv11,CA1 (3,:),'-b','LineWidth',2);
        plot(tv11{1},caMean1{1},'-b','LineWidth',2);
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
     set(gca,'XTickLabel',[-50 0 50 100 150 200 250 300 350 400])
    xlabel('Time (ky)');
    ylabel('Ca^2^+ (mmol kg^{-1})');
    legend('mean surface ocean')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    xlabel('Time (ky)')
    ylabel({'[Ca^2^+]';'(mmol kg^{-1})'})
    
    text(0.02,0.98,'g)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
% epsilon and d44in sensitivityies, respectively  

figurehndl=53;    
FigHandle = figure(figurehndl);
  set(FigHandle, 'Position', [800, 0, 500, 500]);
    clf;    
    subplot(211)

        box on;
        hold on
        myh1 = plot(tv{1}, d44cabulkV{1},'-b','LineWidth',2);
        myh2 = plot(tv11{1},d44cabulkV1{1},'-b', 'LineWidth', 2);
        
        myh3 = plot(tv{2}, d44cabulkV{2},'-k','LineWidth',2);
        myh4 = plot(tv11{2},d44cabulkV1{2},'-b', 'LineWidth', 1);
        
        myh5 = plot(tv{3}, d44cabulkV{3},'-g','LineWidth',2);
        myh6 = plot(tv11{3},d44cabulkV1{3},'-b', 'LineWidth', 1);
        
        hold off;

    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    ylabel('Ca^2^+ (mmol kg^{-1})');
   legend([myh1 myh3 myh5],'\tau_{Ca} \approx 1 Ma (standard)','\tau_{Ca} \approx 2 Ma','\tau_{Ca} \approx 0.5 Ma')
    refresh;
    set(gca,'XTickLabel',[])
%     set(gca,'YLim',[-0.9 -0.5]);
    ylabel('\delta^{44} Ca_{carb} (‰)')
    text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
    
    
    
    for i=1:1:3

    if(i==1)
        dir = 'dat/FinalPERMsims/Sens1/a/';
    elseif(i==2)
        dir = 'dat/FinalPERMsims/Sens2/b/';
    elseif(i==3)
        dir = 'dat/FinalPERMsims/Sens2/c/';
    end
    BurialCA{i,:}= load([dir 'BurialCA.DAT']);
BurialCI{i,:}= load([dir 'BurialCI.DAT']);
BurialCP{i,:}= load([dir 'BurialCP.DAT']);
BurialCT{i,:}= load([dir 'BurialCT.DAT']);

d44fcA{i,:}= load([dir 'd44fcA.DAT']);
d44fcI{i,:}= load([dir 'd44fcI.DAT']);
d44fcP{i,:}= load([dir 'd44fcP.DAT']);
d44fcT{i,:}= load([dir 'd44fcT.DAT']);
    tv{i,:}= load([dir 'tv.DAT']);
tv11{i,:}= load([dir 'tv11.DAT']);    
 d44cabulkV{i,:} = ((d44fcA{i,:}(:,1)*BurialCA{i,:}+d44fcP{i,:}(:,1)*BurialCP{i,:}+d44fcI{i,:}(:,1)*BurialCI{i,:}+d44fcT{i,:}(:,1)*BurialCT{i,:})./(BurialCT{i,:}+BurialCP{i,:}+BurialCI{i,:}+BurialCA{i,:})+...
         (d44fcA{i,:}(:,2)*BurialCA{i,:}+d44fcP{i,:}(:,2)*BurialCP{i,:}+d44fcI{i,:}(:,2)*BurialCI{i,:}+d44fcT{i,:}(:,2)*BurialCT{i,:})./(BurialCT{i,:}+BurialCP{i,:}+BurialCI{i,:}+BurialCA{i,:}))/2;      
 d44cabulkV1{i,:} = [d44cabulkV{i,:}(1,1)'  d44cabulkV{i,:}(1,1)'];   
    end

    
    subplot(212)
    box on;
    hold on;
       myh1 = plot(tv{1}, d44cabulkV{1},'-b','LineWidth',2);
        myh2 = plot(tv11{1},d44cabulkV1{1},'-b', 'LineWidth', 2);
        
        myh3 = plot(tv{2}, d44cabulkV{2},'-k','LineWidth',2);
        myh4 = plot(tv11{2},d44cabulkV1{2},'-b', 'LineWidth', 2);
        
        myh5 = plot(tv{3}, d44cabulkV{3},'-g','LineWidth',2);
        myh6 = plot(tv11{3},d44cabulkV1{3},'-b', 'LineWidth', 2);

        hold off;

    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
     set(gca,'YLim',[-0.9 -0.5]);
      set(gca,'XTickLabel',[-50 0 50 100 150 200 250 300 350 400])
    xlabel('Time (ky)');
    ylabel('\delta^{44} Ca_{carb} (‰)');
    legend([myh1 myh3 myh5],'\delta^{44}Ca_{riv} = -0.6 (standard)','\delta^{44}Ca_{riv} = -0.3','\delta^{44}Ca_{riv} = -0.9')
    refresh;
    
    text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
    
 figurehndl=54;    
FigHandle = figure(figurehndl);
box on
  set(FigHandle, 'Position', [1600, 0, 500, 500]);   
  %Hinojosa data
hd =  csvread('hinojosadata.csv');
stelH = hd(:,1); % stratigraphic elevation in meters
stelH1 = hd(7:end,1); % stratigraphic elevation in meters
d44caH = hd(:,2); 
stdpH = hd(:,3); % standard deviation of calcium data
depthH = [-0.88 0.5];
%calculating time from known time points
for i=1:length(stelH)
timeH(i) = interp1(depthH, [252.87 252.36] , stelH(i),'linear','extrap');
end
  
meand44Pre=mean(d44caH(1:13));
meand44Post=mean(d44caH(14:17));
meanStd=mean(stdpH(1:13));
semPre=(meanStd)/sqrt(length(d44caH(1:13)));  

meanStdpost=mean(stdpH(14:17));
semPost=(meanStdpost)/sqrt(length(d44caH(14:17)));

hold on
% plot(252.6-tv{1}/1e6,d44caMean{1})
h1=plot([-20 0.0008]*1e4,d44caMean1{1},'-b');
h2 = plot(tv{1},d44caMean{1},'-b');
  set(h1,'LineWidth',2)
  set(h2,'LineWidth',2)
% plot([252.6 252.6],[-0.2 1],'-k')
% plot(timeH,d44caH(7:end)+1.1,'o','Color','g')
h3=plot((252.6-timeH)*1e6,d44caH+1.05,'o','Color','g');
h = errorbar((252.6-timeH)*1e6,d44caH+1.05,stdpH);
plot([0 0],[-0.1 0.7],'k')
set(h,'linestyle','none','Color','g')
hold off
set(gca,'XLim',[-200000 400000])
set(gca,'XTickLabel',[-200 -100  0 100 200  300 400])
% set(hndl2,'XDir','reverse')
legend([h2 h3],'LOSCAR mean ocean','data - Hinojosa et al. 2012 (scaled)')
xlabel('Time (ky)')
ylabel('\delta^{44} Ca_{sw} (‰)');
% title('Hinojosa') 

figure
box on
mt=[252.9254 252.5734];
hold on
plot(mt(1),meand44Pre+1.05,'.','markersize',10,'Color','g')
eh1=errorbar(mt(1),meand44Pre+1.05,semPre,'Color','g');    
set(eh1,'XData',[252.9254-0.5 252.9254+0.5]+0.5)
plot(mt(2),meand44Post+1.05,'.','markersize',10,'Color','r')
eh2=errorbar(mt(2),meand44Post+1.05,semPost,'Color','r') 
set(eh2,'XData',[252.5734-0.5 252.5734+0.5]+0.5)
ylabel(['scaled \delta^{44} Ca_{conod}' '(' char(8240) ')']);
set(gca,'XDir','reverse','XLim',[251.8 253])%,'XLim',[251.8 253]
 legend('pre-extinction mean','standard error of the mean','excursion mean','standard error of the mean')
hold off    
    
% % Fractionation and din sensitivities
% dir = 'dat/FinalPERMsims/Sens3/a/';   
% figure
% subplot(421)
% plot(tv,eps)
% subplot(422)
% subplot(423)
% subplot(424)


% hf = figure;
% X = 0:pi/10:pi;
% Y = sin(X);
% E = std(Y)*ones(size(X));
% ha = errorbar(X,Y,E);
% hb = get(ha,'children');  
% Xdata = get(hb(2),'Xdata');
% temp = 4:3:length(Xdata);
% temp(3:3:end) = [];
% % xleft and xright contain the indices of the left and right
% %  endpoints of the horizontal lines
% xleft = temp; xright = temp+1; 
% % Increase line length by 0.2 units
% Xdata(xleft) = Xdata(xleft) - .1;
% Xdata(xright) = Xdata(xright) + .1;
% set(hb(2),'Xdata',Xdata)











%#########################################################################
%#########################################################################
 fs = 18;
 ms = 20;
 lw = 4;
FigHandle = figure(6);
  set(FigHandle, 'Position', [200, 0, 600, 1000]);
    clf;
%     X=[tv{1},fliplr(tv{1})];
%     Y=[pco2t{1},fliplr(pco2t{2})];
    subplot(411);
    plot(tv{1}  ,pco2t{1} ,'r-','LineWidth',lw);
    hold on;
%     plot(tv{2}  ,pco2t{2} ,'r--','LineWidth',1);
%     plot(tv{3}  ,pco2t{3} ,'r.','LineWidth',1);

    plot(tv11{1},pco2t1{1},'r-','LineWidth',lw);
    hold off;
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('pCO_2 (\muatm)','fontsize',fs);

    text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    title('LOSCAR Simulation: 12,000 Pg C, Strangelove')
    
  
    subplot(412);
   
    box  on;
    hold on;

    [ax,p1,p2] = plotyy(ctv{1},co3Mean{1}.*1e6,ctv{1},oMean{1},'plot','plot');
%       [ax2,p3,p4] = plotyy(ctv{2},co3Mean{2}.*1e6,ctv{2},oMean{2},'plot','plot');
%     plot(tv  ,co3Mean.*1e6,sstr(2*3-1:2*3),'Color',cs(3));
%     plot(tv11,co3Mean1*1e6,sstr(2*3-1:2*3),'Color',cs(3));
    %end;
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
%     ylabel('[CO_3^{2-}] (\mumol kg^{-1})');
    ylabel(ax(1),{'[CO_3^{2-}]';' (\mumol kg^{-1})'}) % label left y-axis
    ylabel(ax(2),'\Omega','fontsize',fs) % label right y-axis
    set(ax(2),'fontsize',fs );
    p1.LineWidth = lw;
    p2.LineWidth = lw; 
%     p3.LineWidth = 1;
%     p4.LineWidth = 1; 
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    lgn=legend('mean surface ocean','mean surface ocean');
    set(lgn,'Location','southeast')
    set(gca,'XTickLabel',[])
    text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
   
 
    subplot(413);
    box on;
    hold on;
   myh1= plot(tv{1}, d13cbulkV{1},'-b','LineWidth',lw );
   myh2= plot(tv11{1}, d13cbulkV1{1},'-b','LineWidth',lw );
   myh3= plot((time-0.26)*1e6,d13cp,'.','markersize', ms);
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('\delta^{13}C_{carb} (�)');
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    legend([myh1 myh3], 'Mean shallow', 'data - Payne et al. 2010')
    set(gca,'XTickLabel',[])
    
 text(0.02,0.98,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    if(CAvflag~=0)
        subplot(414);
        box on;
%         plot(tv  ,d44ca (:,9),sstr(2*9-1:2*9),'Color',cs(9));
        hold on
        [ax,p1,p2] = plotyy(ctv{1},d44cabulkVf{1},ctv{1},epspcaVf{1},'plot','plot');
%         [ax2,p3,p4] = plotyy(tv11{1},d44cabulkV1{1},tv11{1},epspcaV1{1},'plot','plot');
        ylabel(ax(1),{'[CO_3^{2-}]';' (\mumol kg^{-1})'}) % label left y-axis
    ylabel(ax(2),'\epsilon_{carb}') % label right y-axis
    p1.LineWidth = lw;
    %Line 2 color
    color2= [204/255 51/255 255/255];
    p2.LineWidth = lw;
    p2.Color = color2;
    set(ax(2), 'YColor',color2,'fontsize',fs );
%         myh1 = plot(tv{1}, d44cabulkV{1},'-b','LineWidth',2);
%         myh2 = plot(tv11{1},d44cabulkV1{1},'-b', 'LineWidth', 2);
        myh3 = plot((time-0.26)*1e6,d44cap,'.','Color','r','markersize', ms);
        h=errorbar((time-0.26)*1e6,d44cap,stdp,'b'); 
        set(h,'linestyle','none','LineWidth',0.5, 'Color',[1 0 0])
%         plot(tv11,d44ca1(9,:),sstr(2*9-1:2*9),'Color',cs(9));
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx,'YLim',[-1 -0.4]);
    ylabel('Ca^2^+ (mmol kg^{-1})');
    legend([p1 p2 myh3],'Mean shallow','Calculated \epsilon_{carb}','data - Payne et al. 2010')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    set(gca,'XTickLabel',[-50 0 50 100 150 200 250 300 350 400])
    xlabel('Time (ky)');
    ylabel('\delta^{44} Ca_{carb} (�)')
    text(0.02,0.98,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
    

    