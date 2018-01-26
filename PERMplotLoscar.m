%%% Permian Big Plot
%%%%Payne data set for d44Ca and d13C end-Permian
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

dir = 'dat/FinalPERMsims/Sim1/'; %Sim2/, Sim2/, Sim3/, BiopSim/
currentRun = 1; % whether to plot current run or saved files, 0 or 1
figurehndl=42;
if(currentRun)
    figurehndl=52;
%%% load all the files
tv= load([dir 'tv.DAT']);
tv11= load([dir 'tv11.DAT']);
pco2t= load([dir 'pco2t.DAT']);

FSichck= load([dir 'FSi.DAT']);
Finchck= load([dir 'FiN.DAT']);
FSichck1 = load([dir 'FSi1.DAT']);
Finchck1= load([dir 'FiN1.DAT']);

a = load([dir 'a.DAT']);
c= load([dir 'c.DAT']);
V= load([dir 'V.DAT']);
co3tv = load([dir 'co3tv.DAT']);

d13fcA = load([dir 'd13fcA.DAT']);
d13fcI= load([dir 'd13fcI.DAT']);
d13fcP= load([dir 'd13fcP.DAT']);
d13fcT= load([dir 'd13fcT.DAT']);

BurialCA= load([dir 'BurialCA.DAT']);
BurialCI= load([dir 'BurialCI.DAT']);
BurialCP= load([dir 'BurialCP.DAT']);
BurialCT= load([dir 'BurialCT.DAT']);

d44fcA= load([dir 'd44fcA.DAT']);
d44fcI= load([dir 'd44fcI.DAT']);
d44fcP= load([dir 'd44fcP.DAT']);
d44fcT= load([dir 'd44fcT.DAT']);

d44ca= load([dir 'd44ca.DAT']);

ccdA= load([dir 'ccdA.DAT']);
ccdI= load([dir 'ccdI.DAT']);
ccdP= load([dir 'ccdP.DAT']);
ccdT= load([dir 'ccdT.DAT']);

omegCSvt= load([dir 'OmegCS.DAT']);
CA= load([dir 'CA.DAT']);
% epspcaV= load([dir 'epspcaV.DAT']);
end
% MEAN SURFACE OCEAN carbonate saturation
oMean=(omegCSvt(:,1).*V(1)+omegCSvt(:,2).*V(2)+omegCSvt(:,3)*V(3)+omegCSvt(:,11)*V(11))./((V(1)+V(2)+V(3)+V(11)));

pco2t1  = [pco2t(1)' pco2t(1)'];

%%% Plot
% figure(42)
FigHandle = figure(figurehndl);
  set(FigHandle, 'Position', [200, 0, 600, 1000]);
    clf;
    
    subplot(711);
    plot(tv  ,pco2t ,'r-','LineWidth',2);
    hold on;
    %for k=kkv
    %plot(tv(lt),pco2v(k),systr(2*k-1:2*k));
    %end;
    plot(tv11,pco2t1,'r-','LineWidth',2);
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
    hndl=plot(tv, FSichck/.1e13,'m',tv,Finchck/.1e13,'g');
    set(hndl,'LineWidth',1)
    legend('F_S_i','F_C')
    hndl1=plot(tv11, FSichck1(9,:)/.1e13,'m',tv11,Finchck1/.1e13,'g');
    set(hndl1,'LineWidth',1)
    %xlabel('Time (y)');
    %ylabel('Silicate flux increase');
    %hold on
    %plot(tv,Finchck,'g','LineWidth',1)
    %plot(tv11,Finchck1(9,:),'g','LineWidth',1)
    hold off
    ylabel({'F_C, F_S_i';' (10^{12} mol C y^{-1})'})
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[])
    text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    % MEAN SURFACE OCEAN alkalinity and total carbon over time
    if(ftys)
        aMean=(a(:,1).*V(1)+a(:,2).*V(2)+a(:,3)*V(3)+a(:,10)*V(10)+a(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
        cMean=(c(:,1).*V(1)+c(:,2).*V(2)+c(:,3)*V(3)+c(:,10)*V(10)+c(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
    else
       aMean=(a(:,1).*V(1)+a(:,2).*V(2)+a(:,3)*V(3)+a(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));
       cMean=(c(:,1).*V(1)+c(:,2).*V(2)+c(:,3)*V(3)+c(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));     
    end
    
    aMean1 = [aMean(1)'  aMean(1)'];
    cMean1 = [cMean(1)'  cMean(1)'];

    subplot(713)
    box  on;
    hold on;
%     hndl1=plot(tv  ,a  (:,3),'-b',tv  ,c  (:,3),'-r');
    hndl1=plot(tv  ,aMean,'-b',tv  ,cMean,'-r');
    set(hndl1,'LineWidth',1)
    legend('TA','TC ')
%     hndl=plot(tv11,a1 (3,:),'-b',tv11,c1 (3,:),'-r');
    hndl=plot(tv11,aMean1,'-b',tv11,cMean1,'-r');
    set(hndl,'LineWidth',1)
    %legend('TA')
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[])
    % set(gca,'XTickLabel',[])
    % plot(tv,TCAR*1e3,'k',tv,TALK*1e3,'r');
    ylabel({'TA, TC';' (mmol kg^{-1})'})
    text(0.02,0.98,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
    % MEAN SURFACE OCEAN carbonate ion conc. over time
    if(ftys)
        co3Mean=(co3tv(:,1).*V(1)+co3tv(:,2).*V(2)+co3tv(:,3)*V(3)+co3tv(:,10)*V(10)+co3tv(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
    else
        co3Mean=(co3tv(:,1).*V(1)+co3tv(:,2).*V(2)+co3tv(:,3)*V(3)+co3tv(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));
    end
    co3Mean1 = [co3Mean(1)'  co3Mean(1)'];
    co3Mean = [co3Mean1 co3Mean']';
    oMean1 = [oMean(1)'  oMean(1)'];
    oMean = [oMean1 oMean']';
    ctv=[tv11 tv']';

    subplot(714);
   
    box  on;
    hold on;

    [ax,p1,p2] = plotyy(ctv,co3Mean.*1e6,ctv,oMean,'plot','plot');
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
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    legend('mean surface ocean','mean surface ocean')
    set(gca,'XTickLabel',[])
    text(0.02,0.98,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    if(ftys)
    d13cbulkV = ((d13fcA(:,1)*BurialCA+d13fcP(:,1)*BurialCP+d13fcI(:,1)*BurialCI+d13fcT(:,1)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA)+...
    (d13fcA(:,2)*BurialCA+d13fcP(:,2)*BurialCP+d13fcI(:,2)*BurialCI+d13fcT(:,2)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA))/2;
    d44cabulkV = ((d44fcA(:,1)*BurialCA+d44fcP(:,1)*BurialCP+d44fcI(:,1)*BurialCI+d44fcT(:,1)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA)+...
    (d44fcA(:,2)*BurialCA+d44fcP(:,2)*BurialCP+d44fcI(:,2)*BurialCI+d44fcT(:,2)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA))/2;  
    else
       d13cbulkV = ((d13fcA(:,1)*BurialCA+d13fcP(:,1)*BurialCP+d13fcI(:,1)*BurialCI)./(BurialCP+BurialCI+BurialCA)+...
    (d13fcA(:,2)*BurialCA+d13fcP(:,2)*BurialCP+d13fcI(:,2)*BurialCI)./(BurialCP+BurialCI+BurialCA))/2;
    d44cabulkV = ((d44fcA(:,1)*BurialCA+d44fcP(:,1)*BurialCP+d44fcI(:,1)*BurialCI)./(BurialCP+BurialCI+BurialCA)+...
    (d44fcA(:,2)*BurialCA+d44fcP(:,2)*BurialCP+d44fcI(:,2)*BurialCI)./(BurialCP+BurialCI+BurialCA))/2;   
        
    end

    d13cbulkV1 = [d13cbulkV(1)'  d13cbulkV(1)'];
    d44cabulkV1 = [d44cabulkV(1)'  d44cabulkV(1)'];
 
    subplot(715);
    box on;
    hold on;
   myh1= plot(tv, d13cbulkV,'-b','LineWidth',2 );
   myh2= plot(tv11, d13cbulkV1,'-b','LineWidth',2 );
   myh3= plot((time-0.26)*1e6,d13cp,'o','markersize', 3);
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('\delta^{13}C_{carb} (â€°) ');
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
        myh1 = plot(tv, d44cabulkV,'-b','LineWidth',2);
        myh2 = plot(tv11,d44cabulkV1,'-b', 'LineWidth', 2);
        myh3 = plot((time-0.26)*1e6,d44cap,'o','Color','r','markersize', 3);
        h=errorbar((time-0.26)*1e6,d44cap,stdp,'b'); 
        set(h,'linestyle','none','LineWidth',0.5, 'Color',[1 0 0])
%         plot(tv11,d44ca1(9,:),sstr(2*9-1:2*9),'Color',cs(9));
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx,'YLim',[-1 -0.4]);
    ylabel('Ca^2^+ (mmol kg^{-1})');
    legend([myh1 myh3],'Mean shallow','data - Payne et al. 2010')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    set(gca,'XTickLabel',[])
    ylabel('\delta^{44} Ca_{carb} (â€°)')
    text(0.02,0.98,'f)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
    

    if(CAvflag~=0)
        
        % MEAN SURFACE OCEAN calcium ion conc. over time
        if(ftys)
            caMean=(CA(:,1).*V(1)+CA(:,2).*V(2)+CA(:,3)*V(3)+CA(:,10)*V(10)+CA(:,11)*V(11))./((V(1)+V(2)+V(3)+V(10)+V(11)));
        else
            caMean=(CA(:,1).*V(1)+CA(:,2).*V(2)+CA(:,3)*V(3)+CA(:,10)*V(10))./((V(1)+V(2)+V(3)+V(10)));
        end
        caMean1 = [caMean(1)'  caMean(1)'];
        
        subplot(717);
%         plot(tv  ,CA  (:,3),'-b','LineWidth',2);
        plot(tv  ,caMean,'-b','LineWidth',2);
        hold on
%         plot(tv11,CA1 (3,:),'-b','LineWidth',2);
        plot(tv11,caMean1,'-b','LineWidth',2);
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
    
    %%% saturation plots
    
    
%#########################################################################    
%#########################################################################    
%#########################################################################    
    
    fs=18;
    ms=20;
    lw = 4;
    FigHandle = figure(2);
  set(FigHandle, 'Position', [200, 0, 600, 1000]);
    clf;
    
    subplot(411);
    plot(tv  ,pco2t ,'r-','LineWidth',lw);
    hold on;
    %for k=kkv
    %plot(tv(lt),pco2v(k),systr(2*k-1:2*k));
    %end;
    plot(tv11,pco2t1,'r-','LineWidth',lw);
    hold off;
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('pCO_2 (\muatm)');

    
    text(0.02,0.98,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
%     title('LOSCAR Simulation: 13,200 Pg C')
    title('LOSCAR Simulation: 43,200 Pg C')
%     title({'LOSCAR Simulation: 12,000 Pg ';' Strangelove'})


    subplot(412);
   
    box  on;
    hold on;

    [ax,p1,p2] = plotyy(ctv,co3Mean.*1e6,ctv,oMean,'plot','plot');
%     plot(tv  ,co3Mean.*1e6,sstr(2*3-1:2*3),'Color',cs(3));
%     plot(tv11,co3Mean1*1e6,sstr(2*3-1:2*3),'Color',cs(3));
    %end;
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    set(gca,'XTickLabel',[])
    set(ax(2),'fontsize',fs );
    %xlabel('Time (y)');
%     ylabel('[CO_3^{2-}] (\mumol kg^{-1})');
    ylabel(ax(1),{'[CO_3^{2-}]';' (\mumol kg^{-1})'}) % label left y-axis
    ylabel(ax(2),'\Omega') % label right y-axis
    p1.LineWidth = lw;
    p2.LineWidth = lw; 
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    lgn =legend('mean surface ocean','mean surface ocean');
    set(lgn,'fontsize',12,'Location','southeast')
    text(0.02,0.98,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    if(ftys)
    d13cbulkV = ((d13fcA(:,1)*BurialCA+d13fcP(:,1)*BurialCP+d13fcI(:,1)*BurialCI+d13fcT(:,1)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA)+...
    (d13fcA(:,2)*BurialCA+d13fcP(:,2)*BurialCP+d13fcI(:,2)*BurialCI+d13fcT(:,2)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA))/2;
    d44cabulkV = ((d44fcA(:,1)*BurialCA+d44fcP(:,1)*BurialCP+d44fcI(:,1)*BurialCI+d44fcT(:,1)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA)+...
    (d44fcA(:,2)*BurialCA+d44fcP(:,2)*BurialCP+d44fcI(:,2)*BurialCI+d44fcT(:,2)*BurialCT)./(BurialCT+BurialCP+BurialCI+BurialCA))/2;  
    else
       d13cbulkV = ((d13fcA(:,1)*BurialCA+d13fcP(:,1)*BurialCP+d13fcI(:,1)*BurialCI)./(BurialCP+BurialCI+BurialCA)+...
    (d13fcA(:,2)*BurialCA+d13fcP(:,2)*BurialCP+d13fcI(:,2)*BurialCI)./(BurialCP+BurialCI+BurialCA))/2;
    d44cabulkV = ((d44fcA(:,1)*BurialCA+d44fcP(:,1)*BurialCP+d44fcI(:,1)*BurialCI)./(BurialCP+BurialCI+BurialCA)+...
    (d44fcA(:,2)*BurialCA+d44fcP(:,2)*BurialCP+d44fcI(:,2)*BurialCI)./(BurialCP+BurialCI+BurialCA))/2;   
        
    end

    d13cbulkV1 = [d13cbulkV(1)'  d13cbulkV(1)'];
    d44cabulkV1 = [d44cabulkV(1)'  d44cabulkV(1)'];
    subplot(413);
    box on;
    hold on;
   myh1= plot(tv, d13cbulkV,'-b','LineWidth',lw );
   myh2= plot(tv11, d13cbulkV1,'-b','LineWidth',lw );
   myh3= plot((time-0.26)*1e6,d13cp,'.','markersize', ms);
    hold off;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
    %xlabel('Time (y)');
    ylabel('\delta^{13}C_{carb} (‰) ');
    %Hl=legend(lstr,4);
    %set(Hl,'FontSize',10);
    lgn=legend([myh1 myh3], 'Mean shallow', 'data - Payne et al. 2010');
    set(lgn,'fontsize',12,'Location','southeast')
    set(gca,'XTickLabel',[])
    
 text(0.02,0.98,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')


    if(CAvflag~=0)
        subplot(414);
        box on;
%         plot(tv  ,d44ca (:,9),sstr(2*9-1:2*9),'Color',cs(9));
        hold on
        myh1 = plot(tv, d44cabulkV,'-b','LineWidth',lw);
        myh2 = plot(tv11,d44cabulkV1,'-b', 'LineWidth', lw);
        myh3 = plot((time-0.26)*1e6,d44cap,'.','Color','r','markersize', ms);
        h=errorbar((time-0.26)*1e6,d44cap,stdp,'b'); 
        set(h,'linestyle','none','LineWidth',0.5, 'Color',[1 0 0])
%         plot(tv11,d44ca1(9,:),sstr(2*9-1:2*9),'Color',cs(9));
        hold off;
    end;
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx,'YLim',[-1 -0.4]);
    ylabel('Ca^2^+ (mmol kg^{-1})');
    lgnd = legend([myh1 myh3],'Mean shallow','data - Payne et al. 2010');
    set(lgnd,'color','none','fontsize',12,'Location','southeast')
    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    ylabel('\delta^{44} Ca_{carb} (‰)')
    text(0.02,0.98,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontw','b')
    
    

    
    set(gca,'FontSize',fs);
    set(gca,'XLim',axx);
     set(gca,'XTickLabel',[0  100  200  300  400])
%      set(gca,'XTickLabel',[-50 0 50 100 150 200 250 300 350 400])


    %text(0.7,0.600,['\delta^{13}C  = ' dstr '?'],'FontSize',fs,'Units','n');
    refresh;
    xlabel('Time (ky)')

    
    
    
    

