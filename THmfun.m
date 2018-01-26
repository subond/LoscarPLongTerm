%-------------------------------------
% function [yp] = THmfun(y,TH,TT,mv,mhd,V,gA,tA,gI,tI)
%
% returns yp for given conveyor configuration
% and mixing. same for all ocean tracers
%
% updates:
%
% 05/03/06 Pac CCD shallowed (x in THmfun)
% 03/01/06 not changed SO circ (fcon = 3)
% 12/30/05 new file
%-------------------------------------
function [yp] = THmfun(y,TH,TS,TT,mv,mhd,V,gA,tA,gI,tI)


global Nb fcon ftys;


% BOXES
% LA LI LP IA II IP DA DI DP  H LT IT DT
%  1  2  3  4  5  6  7  8  9 10 11 12 13

%THX      = TH;
%TTX      = TT;
yp(1:Nb) = 0;


% Fraction SO Source A I P
%x= [1 1 1]./3; % .3 .3 .4
x = [.4 .3 .3]; % .4 .3 .3 (05/03/06)

% TH branches
% tA        upwelled into intermdt Atl 0.27 0.15
% tI        upwelled into intermdt Ind 0.29 0.30
% gA = (1-tA)        export into     Deep Ind
% gI = (1-tA-tI)        export into     Deep Pac

% TH  
if     (fcon == 1) %%% NADW
yp(04) =     ( gA*TH* y(05)+ ...                % IA
               tA*TH* y(07)- ...                %
                  TH* y(04)      )/V(04);       % 
yp(05) =     ( gI*TH* y(06)+ ...                % II
               tI*TH* y(08)- ...                %
               gA*TH* y(05)      )/V(05);       %  
yp(06) =       gI*TH*(y(09)-y(06))/V(06);       % IP
yp(07) =          TH*(y(10)-y(07))/V(07);       % DA
yp(08) =       gA*TH*(y(07)-y(08))/V(08);       % DI
yp(09) =       gI*TH*(y(08)-y(09))/V(09);       % DP
yp(10) =          TH*(y(04)-y(10))/V(10);       % H

    if(1) % add SO ???
    fa = x(1); fi = x(2);fp = x(3); 
    yp(04) =yp(04)+fa*TS*(y(07)-y(04)) /V(04);      % IA
    yp(05) =yp(05)+fi*TS*(y(08)-y(05)) /V(05);      % II
    yp(06) =yp(06)+fp*TS*(y(09)-y(06)) /V(06);      % IP
    yp(07) =yp(07)+fa*TS*(y(10)-y(07)) /V(07);      % DA
    yp(08) =yp(08)+fi*TS*(y(10)-y(08)) /V(08);      % DI
    yp(09) =yp(09)+fp*TS*(y(10)-y(09)) /V(09);      % DP
    yp(10) =yp(10)+(fa*TS*(y(04)-y(10))     ...     % H
                  +fi*TS*(y(05)-y(10))      ...     % H
                  +fp*TS*(y(06)-y(10)))/V(10);      % H
    end;

elseif(fcon == 2) %%% NPDW
if(1)    
%TH = 0.4*THX; % 0.3(1,1) 0.4(2,3)
%TS = 0.6*THX;
%TT = 1.0*TTX;

% BOXES
% LA LI LP IA II IP DA DI DP  H LT IT DT
%  1  2  3  4  5  6  7  8  9 10 11 12 13
yp(06) =     ( gA*TH* y(05)+ ...                % IP
               tA*TH* y(09)- ...                %
                  TH* y(06)      )/V(06);       % 
yp(05) =     ( gI*TH* y(04)+ ...                % II
               tI*TH* y(08)- ...                %
               gA*TH* y(05)      )/V(05);       %  
yp(04) =       gI*TH*(y(07)-y(04))/V(04);       % IA
yp(09) =          TH*(y(10)-y(09))/V(09);       % DP
yp(08) =       gA*TH*(y(09)-y(08))/V(08);       % DI
yp(07) =       gI*TH*(y(08)-y(07))/V(07);       % DA
yp(10) =          TH*(y(06)-y(10))/V(10);       % H
end;
if(1) % add SO ???
fa = x(1); fi = x(2);fp = x(3); 
yp(04) =yp(04)+fa*TS*(y(07)-y(04)) /V(04);      % IA
yp(05) =yp(05)+fi*TS*(y(08)-y(05)) /V(05);      % II
yp(06) =yp(06)+fp*TS*(y(09)-y(06)) /V(06);      % IP
yp(07) =yp(07)+fa*TS*(y(10)-y(07)) /V(07);      % DA
yp(08) =yp(08)+fi*TS*(y(10)-y(08)) /V(08);      % DI
yp(09) =yp(09)+fp*TS*(y(10)-y(09)) /V(09);      % DP
yp(10) =yp(10)+(fa*TS*(y(04)-y(10))     ...     % H
              +fi*TS*(y(05)-y(10))      ...     % H
              +fp*TS*(y(06)-y(10)))/V(10);      % H
end;
% BOXES
% LA LI LP IA II IP DA DI DP  H LT IT DT
%  1  2  3  4  5  6  7  8  9 10 11 12 13
elseif(fcon == 3) %%% SO
fa = x(1); fi = x(2);fp = x(3); 
yp(04) =       fa*TH*(y(07)-y(04)) /V(04);      % IA
yp(05) =       fi*TH*(y(08)-y(05)) /V(05);      % II
yp(06) =       fp*TH*(y(09)-y(06)) /V(06);      % IP
yp(07) =       fa*TH*(y(10)-y(07)) /V(07);      % DA
yp(08) =       fi*TH*(y(10)-y(08)) /V(08);      % DI
yp(09) =       fp*TH*(y(10)-y(09)) /V(09);      % DP
yp(10) =      (fa*TH*(y(04)-y(10))      ...     % H
              +fi*TH*(y(05)-y(10))      ...     % H
              +fp*TH*(y(06)-y(10)))/V(10);      % H
end; % fcon

% Tethys
% LI > LT > DT > DI > II > LI ... etc
if(ftys) %TT
yp(11) =          TT*(y(02)-y(11))/V(11);       % LT
yp(13) =          TT*(y(11)-y(13))/V(13);       % DT
yp(08) = yp(08)+  TT*(y(13)-y(08))/V(08);       % DI
yp(05) = yp(05)+  TT*(y(08)-y(05))/V(05);       % II
yp(02) = yp(02)+  TT*(y(05)-y(02))/V(02);       % LI
end;  

% mixing
for k=1:3 % L-I
yp(k  ) = yp(k  )+ mv(k)*(y(k+3)-y(k  ))/V(k  );% L-I
yp(k+3) = yp(k+3)+ mv(k)*(y(k  )-y(k+3))/V(k+3);% I-L
end;
for k=1:3 % H-D
yp(k+6) = yp(k+6)+mhd(k)*(y(010)-y(k+6))/V(k+6);% D-H
yp(010) = yp(010)+mhd(k)*(y(k+6)-y(010))/V(010);% H-D
end;
% Tethys
if(ftys)
yp(012) =           mv(4)*(y(011)-y(012))/V(012);% IT-LT
yp(011) = yp(011)+  mv(4)*(y(012)-y(011))/V(011);% LT-IT
yp(012) = yp(012)+  mv(5)*(y(005)-y(012))/V(012);% IT-II
yp(005) = yp(005)+  mv(5)*(y(012)-y(005))/V(005);% II-IT
yp(011) = yp(011)+ mhd(4)*(y(013)-y(011))/V(011);% LT-DT
yp(013) = yp(013)+ mhd(4)*(y(011)-y(013))/V(013);% DT-LT
end;


return;



