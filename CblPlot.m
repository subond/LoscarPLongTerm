figure(777)
RlsCtv = load(['dat/PCycle/Leak0/' 'RlsCtv.DAT']);

%---------------------- Calc Release 

Rm  =  (RlsCtv(1:end-1)+RlsCtv(2:end  ))/2./1.e15; % Gt C y
dtv =  (    tv(2:end  )-    tv(1:end-1));
tm  =  (    tv(1:end-1)+    tv(2:end  ))/2.;
    
% total release
TR = sum(Rm*diff(tv)*1e3);  % Gt C
% split release
 DTSv = load(['dat/PCycle/Leak0/' 'fDTS.DAT'])/1.e3;
DTS  = DTSv(1)+0.1; % ky duration 1st bump
DTS2 = DTSv(2)+0.1; % ky duration 1st cont release
ts3  = DTSv(3)    ; % ky onset    2nd bump
DTS3 = DTSv(4)+0.1; % ky duration 2nd bump
DTS4 = DTSv(5)+0.1; % ky duration 2nd cont release

% indices
[a,kt1] = min(abs(tv-DTS     )); % end   1st bump
[a,kt3] = min(abs(tv-ts3     )); % start 2nd bump
[a,kt4] = min(abs(tv-ts3-DTS3)); % end   2nd bump


R1 = sum(Rm(    1:kt1).*dtv(    1:kt1)'*1e3); % Gt C
R2 = sum(Rm(kt1+1:kt3).*dtv(kt1+1:kt3)'*1e3); % Gt C
R3 = sum(Rm(kt3+1:kt4).*dtv(kt3+1:kt4)'*1e3); % Gt C
R4 = sum(Rm(kt4+1:end).*dtv(kt4+1:end)'*1e3); % Gt C

% bumps
str1 = sprintf('%4.0f Pg C',10*round(R1/10));
str3 = sprintf('%4.0f Pg C',10*round(R3/10));
% sum cont release
str2 = sprintf('%4.0f Pg C',10*(round(R2/10)+round(R4/10)));
% sum cont release + bump2
str4 = sprintf('%4.0f Pg C',10*round(R3/10)+10*(round(R2/10)+round(R4/10)));
% sumtotal 
str  = sprintf('%4.0f Pg C',10*(round(R1/10)+...
                                round(R2/10)+...
                                round(R3/10)+...
                                round(R4/10)));
dstr = sprintf('%+3.0f',[-50]);

%----------------------- Release --------------------------% 
subplot(111)
axx=[-50 200]*1e3;
trl = [0.; tm];
 rl = [0.  Rm];
plot(trl,rl,'-k');
hold on;
fill(trl,rl,'r ');
hold off;
set(gca,'FontSize',fs);
set(gca,'XLim',axx);
set(gca,'YLim',[0 .55]);
% set(gca,'XTickLabel',[]);
% set(gca,'XTickL',[]);
ylabel('(Pg C y^{-1})');
text(020,0.400,['      ' str1],'FontSize',fs);
text(070,0.100,['      ' str4],'FontSize',fs);
annotation('arrow',[.37 .32],[.905 .905]);
annotation('arrow',[.54 .49],[.86 .86]);

%text(005,0.400,['  \leftarrow' str1],'FontSize',fs);
%text(070,0.100,['  \leftarrow' str4],'FontSize',fs);
%text(0.7,0.800,['Total = '     str ],'FontSize',fs,'Units','n');
text(1.0,0.800,['\delta^{13}C_{^{inp}}  = ' dstr '‰ '],'FontSize',...
    fs,'Units','n','Hor','r');

% text(xyl(1),xyl(2),Flb(1),'FontS',ffs,'FontW','b','Units','n','Col','k');