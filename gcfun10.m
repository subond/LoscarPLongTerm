function [ybbv,ybv,ycv,aav,bbv,ccv,fGG,fgkc,frkc,fekc,fdkc,flakc,fbch,fbbv,fwcv,fwsv1,fwsv2,fbcv,fbgv,rco2,pco2gca,fwgv,fmcv,fmgv,ggc,cgc,dg,dc]=gcfun10(t,nFBG)
global time ws gamma RUN ACT myfwc pco20 FERT acv dbckv fckc crit tnr dbcv

XX=load('dat/LPEEkurtz/Sim2/GeoCarb3Dat.dat');
% GG=load('GEOGt.txt');
%  GC=load('GC3intdata.dat');     %%% 1-million year step GEOCARB data -d13c
%(linearly interpolated between 57 and 52 Ma)and eps constant after 52Ma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USED FOR MASTER THESIS
%GC=load('GC3intdataOR.dat');    %%% 1 million year step ORIGINAL GEOCARB III
GC=load('dat/LPEEkurtz/Sim2/GC3intdataHil.dat');    %%% 1-million year step GEOCARB data -d13c
%                                      from Hilting '08 and eps
%                                      recalculated
% data
% Kurtz data + dickens13 for period between 50 and 62 Ma
%c13kurtz=load('C:\Users\Nemanja\Dropbox\geocarb\c13kurtz.dat');
c13kurtz=load('dat/LPEEkurtz/Sim2/c13kurtz2.DAT');
% dbc=[1.0000    1.7000    2.0000    2.2000    2.2000    1.3000    5.3000    2.4000    2.4000    2.5000 2.6000    2.7000    2.5000    2.0000    1.2000    1.7000    1.5400    1.3800    1.2200    1.0600 0.9000    0.7600    1.0000    2.0000    2.5000    3.0000    4.0000    4.5000    4.8000    5.0000 5.0000    4.8000    4.5000    4.0000    3.5000    3.0000    2.0000    1.5000    1.0000    1.0000 1.5000    2.0000    2.0000    1.5000    0.4000         0   -0.4000   -0.8000   -0.8000   -0.8000 -0.8000   -0.8000   -0.8000   -0.8000   -0.8000   -0.8000   -0.8000   -0.8000];
% fla=[1.0000    1.0000    0.9900    1.2900    1.0010    1.0    1.10    0.7656 0.6952 0.7128 0.8528 0.8840 0.9152 0.8944 0.8736 1.0292 1.0625 1.0875 1.1125 1.2376 1.2376 1.1070 1.1070 1.2232 1.2104 1.2298 1.1004 1.0611 1.0611 1.1600 1.1600 1.1600 1.1139 1.0857 1.0716 1.0293 0.9800 1.0143 1.0143 1.0584 1.2062 1.2062 1.1396 0.8866 0.9009 0.7705 0.6634 0.6634 0.6138 0.6435 0.6435 0.6402 0.6720 0.6930 0.5355 0.5670 0.5985 0.6300];
% fa=XX(:,11)';
% fd=XX(:,3)';
% fd=[fd.*fa];
% fg=XX(:,6)';
% alpha=[22    25    27    29    29    30.472    32.132   28    29    29    28    29    30    29    29 29    30    32    32    31    31    31    31    30    31    31    32    32    32    32 31    31    31    31    31    31    31    30    30    29    28    28    29    31    31 30    30    29    29    29    29    29    29    29    28    28    28    28];
% myFBGV=[5.0000    4.9243    5.0554    5.1703    5.2641    5.3883    5.4987    5.4957    5.3363    5.1974    5.1363 5.1394    5.1720    5.2954    5.4644    5.5581    5.5874    5.5220    5.3945    5.3004    5.2778    5.2328 5.1947    5.0973    4.8754    4.6745    4.5711    4.4540    4.4440    4.4584    4.4590    4.4269    4.3984 4.3974    4.4347    4.3709    4.3631    4.3479    4.3523    4.3394    4.2826    4.2878    4.4076    4.5178 4.5609    4.5640    4.5244    4.5125    4.5733    4.6654    4.6622    4.2529    3.7638    3.7620    4.0784 4.8746    5.5573    5.9508    6.0300    5.7816    5.4080    5.0751    4.8874];
myFBGV=[5.0000    4.7606    4.8889    5.0030    5.0984    5.2217    5.3327    5.3399    5.2038    5.0844    5.0350  5.0461    5.0835    5.2018    5.3619    5.4558    5.4912    5.4402    5.3333    5.2567    5.2442    5.2270  5.2162    5.1546    4.9893    4.8467    4.7879    4.7165    4.7306    4.7634    4.7844    4.7507    4.7280  4.7277    4.7577    4.7051    4.6948    4.6767    4.6740    4.6563    4.6010    4.6529    4.7694    4.8779  4.9303    4.9499    4.9345    4.9418    5.0090    5.1020    5.1164    4.7911    4.4083    4.4309    4.7054  5.3541    5.9103    6.2350    6.3038    6.1044    5.8098    5.5252    5.3527];
nC=0.23;
nS=0.26;
geogt=XX(:,10)';
% fr=XX(:,4)';
% fc=XX(:,7)';
% fe=XX(:,5)';
dbc=GC(:,1)';
dbc=c13kurtz(2,:);
fla=GC(:,2)';
fd=GC(:,3)';
fg=GC(:,4)';
fr=GC(:,5)';
fe=GC(:,6)';
fc=GC(:,7)';
% alpha=GC(:,8)';
% See komar et al. 2013 to see how alpha is calculated
% to calculate we need d13c of total organic carbon at time t (dtoc)
% I am not sure where I saved this data but I know I calculated alpha 
% correctly from it. That means I can use alpha to solve eq. 3 from komar
% et al. for dtoc. Even the dtoc is calculated it is the actual data. Now
% if we ever have a new dbc data we can alpha based on dtoc 
% dtoc=((c13kurtz(2,:)+1000)./((alpha/1000)+1))-1000
alpha=c13kurtz(3,:); %new correct ep data based on kurtz d13c




mend=10;

% bf = 1 if you want to integrate from 570 BP forward and -1 if
% you want to integrate from present time backwards

bf=-1;
dk=57*mend*(1+bf)/2;

% Initial values and constants

% i is initial entry for variables
% j is initial entry for matrices

i=1+dk;
j=1+dk/mend;
ac=alpha(j);
kwc=0.0024;   % default 0.00267; 0.0024 to get 12x10^12 mol/y
kwg=0.003;     % default 0.003
% ACT=0.06;      % default 0.09
% FERT=0.15;     % default 0.4

ws=7.4;        %GEOCARB II = 12.9, GEOCARB III = 7.4
gamma=ones(1,570)*4.33; % standard 3.3, 4.33=3 deg C
gamma(1,1:41)=4.0;
gamma(1,261:341)=4.0;
RUN=ones(1,570)*0.025;
RUN(1,1:41)=0.045;
RUN(1,261:341)=0.045;

% For testing, MAKE SURE TO REMOVE THIS AFTER!
% and make to remove "570" in 4 places in replace it
% with appropriate variables (t and 3 p's)
% dbc=ones(1,length(dbc)).*dbc(58);
% dbc=3+[1:71]./15;
% geogt=ones(1,length(geogt)).*geogt(58);
% fla=ones(1,length(fla)).*fla(58);
% fd=ones(1,length(fd)).*fd(58);
% fg=ones(1,length(fg)).*fg(58);
% fr=ones(1,length(fr)).*fr(58);
% fe=ones(1,length(fe)).*fe(58);
% fc=ones(1,length(fc)).*fc(58);
% alpha=ones(1,length(alpha)).*alpha(58);
% gamma=ones(1,570)*4.00;
% RUN=ones(1,570)*0.045;
% dbc(58)=6;

aav(1)=1;
bbv(1)=1;
ccv(1)=1;
ybbv(1)=1;
ybv(1)=1;
ycv(1)=1;


%ws=0;

% Use this section if you want C + G = 6250

x= 3.0;         %default 2.0   2.06      1.9311
x1= 24.1875;       %default 27.7  27.01   24.1875 24.755
% x= -.0264;         %default 2.0   2.06      1.9311 -.0264
% x1= 21.2;       %default 27.7  27.01   27.1875    21.2
dcv(1:80)=-4.0;
cgc(i)=5000*(1-bf)/2+5005*(1+bf)/2;
ggc(i)=1250*(1-bf)/2+1245*(1+bf)/2;
dc(i)=x*(1-bf)/2+1.4454*(1+bf)/2;
dg(i)=-x1*(1-bf)/2-23.3787*(1+bf)/2;
fbc=17*(1-bf)/2+14.7133*(1+bf)/2;
fbg=5*(1-bf)/2+3.4767*(1+bf)/2;
r=(1-bf)/2+16.74*(1+bf)/2;
BIGD(i)=(fbg*(dbc(j)-alpha(j))+fbc*dbc(j))/(fbg+fbc);
fbc0=fbc;
fbg0=fbg;


rco2(i)=r;
fwg=fd(j)*1*kwg*ggc(i);
time(i)=-570*(1+bf)/2;
t=-time(i);
fwg0=fwg;

cc=1-0.087*ws*t/570;
fbb=(cc+gamma(t+1)*0.087*log(r))*sqrt(r);
fwc=fbb*fla(j)*fd(j)*fe(j)*kwc*cgc(i);
dbck=dbc(j);
fmg0=((dbc(1)*fbc+(dbc(1)-alpha(1))*fbg-dc(1)*fwc+dcv(1)*fwc-dcv(1)*(fbc+fbg))/(dg(1)-dcv(1)))-kwg*1250;
fmg=fg(j)*fmg0;
fmc0=fbc+fbg-kwc*5000-kwg*1250-fmg0;
fwc0=fwc;

fmc=fg(j)*fc(j)*fmc0;
fwsi0=fbc-fwc;
fwsv2(1)=fwsi0;
fwsv1(1)=fwsi0;
dbcv(1)=dbc(1);
p=i;
tnr=0.0001;
L=[0:10:570]';          %
X=[0:1:570]';



for k=1:62;

    cgc(p-bf)=cgc(p);%-bf*(fwc+fmc-fbc)*10/mend;
    ggc(p-bf)=ggc(p);%-bf*(fwg+fmg-fbg)*10/mend;
    dg(p-bf)=dg(p);%-bf*(fbg*(ac+dg(p)-dbck)/ggc(p))*10/mend;
    dc(p-bf)=dc(p);%-bf*(fbc*(dc(p)-dbck)/cgc(p))*10/mend;

    p=p-bf;
    time(p)=-10*(p-1)/mend;
    t=-time(p);
    if(p>=0)
        dg(p)=-21.2;
        dc(p)=+2.0;
        dcv(p)=-4.0;
    end

    % dbck=interp1(L,dbc,k);
    % flak=interp1(L,fla,k);
    % fdk=interp1(L,fd,k);
    % fgk=interp1(L,fg,k);
    % frk=interp1(L,fr,k);
    % fek=interp1(L,fe,k);
    % fck=interp1(L,fc,k);
    % ac=interp1(L,alpha,k);

    dbck=dbc(p);
    flak=fla(p);
    fdk=fd(p);
    fgk=fg(p);
    frk=fr(p);
    fek=fe(p);
    fck=fc(p);
    ac=alpha(p);

    P=[geogt];             %

    fGG=interp1(L,P,X);

    flakc(p)=flak;
    fdkc(p)=fdk;
    fckc(p)=fck;
    fekc(p)=fek;
    frkc(p)=frk;
    fgkc(p)=fgk;
    acv(p)=ac;
    dbckv(p)=dbck;

    % Solve for fmc, fmg, and fwg using equations 4-6

    fmc=fgk*fck*fmc0;%8.7660;
    fmg=fgk*fmg0;
    fwg=fdk*frk*kwg*ggc(p);

    fmcv(p)=fmc;
    fmgv(p)=fmg;
    fwgv(p)=fwg;
    % Set up constants needed to solve for r, fB, and fBB

    aa=exp(-ACT*ws*p/570)*exp(ACT*fGG(p));
    bb=1-RUN(p)*ws*p/570;
    cc=1-0.087*ws*p/570;
    ybb=flak*fdk*fek*kwc*cgc(p)*(dbck-dc(p));
    yb=ac*frk*fek*fwsi0*fdk^0.65;
    yc=(fmc+fwg+fmg)*(dbck-ac)-dcv(p)*fmc-dg(p)*(fwg+fmg);

    aav(p)=aa;
    bbv(p)=bb;
    ccv(p)=cc;
    ybbv(p)=ybb;
    ybv(p)=yb;
    ycv(p)=yc;

    % Iterate to solve for r, fB and fBB using equations 3 and 11
    % Crit must be less than tnr to terminate the Newton-Raphson

    crit=1;

    if t <= 350

        % This is during time of vascular plants
%         f=@(r) aa*r^(ACT*gamma(t))*(bb+(gamma(t)*RUN(t))*log(r)+RUN(t)*fGG(t))^0.65*(2*r/(1+r))^0.4
%         
%         return
        while crit >= tnr;
            fb=aa*r^(ACT*gamma(p))*(bb+(gamma(p)*RUN(p))*log(r)+RUN(p)*fGG(p))^0.65*(2*r/(1+r))^FERT;
            fbb=(cc+gamma(p)*0.087*log(r)+0.087*fGG(p))*(2*r/(1+r))^FERT;
            % f=@(r) aa*r^(ACT*gamma(t))*(bb+(gamma(t)*RUN(t))*log(r)+RUN(t)*fGG(t))^0.65*(2*r/(1+r))^0.4;
            % dfb=deriv(f,r,1e-3,'complex');
            % f2=@(r) (cc+gamma(t)*0.087*log(r)+0.087*fGG(t))*(2*r/(1+r))^0.4;
            %  dfbb=deriv(f2,r,1e-3,'complex');
            % f=@(r) (r/1)^nS;
            % f2=@(r) (r/1)^nC;
            % dfb=deriv(f,r,1e-3,'complex')
            % dfbb=deriv(f2,r,1e-3,'complex');
            % dfb=nS*(r/1)^(nS-1);
            % dfbb=nC*(r/1)^(nC-1);
            dfb=(fb/r)*((ACT*gamma(p))+(0.65*RUN(p)*gamma(p))/(bb+(gamma(p)*RUN(p))*log(r)+RUN(p)*fGG(p))+FERT/(1+r));
            dfbb=(fbb/r)*((gamma(p)*0.087)/(cc+gamma(p)*0.087*log(r)+0.087*fGG(p))+FERT/(1+r));

            y=fbb*ybb+fb*yb+yc;
            dy=dfbb*ybb+dfb*yb;
            crit=abs(y/dy);
            r=r-y/dy;
        end;

    end
    % This is before time of vascular plants

    pco2gca(1)=rco2(1)*pco20;
    rco2(p)=r;
    pco2gca(p)=pco20*rco2(p);
    fbbv(p)=fbb;
    fbch(p)=fb;
    % Solve for Fwc using equation 3

    fwc=fbb*flak*fdk*fek*kwc*cgc(p);
    fws=fb*frk*fek*fdk^0.65*fwsi0;
    fwcv(p)=fwc;
    myfwc(p)= (rco2(p)/rco2(1))^nC*kwc*cgc(p)*flak*fdk*fek;
    % Solve for Fbc and Fbg using equations 1-2
    fbcv(1)=fbc0;
    fbgv(1)=fbg0;
    fwcv(1)=fwc0;
    fmcv(1)=fmc0;
    fwgv(1)=fwg0;
    fmgv(1)=fmg0;
    flakc(1)=1;
    fdkc(1)=1;
    fckc(1)=1;
    fekc(1)=1;
    frkc(1)=1;
    fgkc(1)=1;
    acv(1)=c13kurtz(3,1);
%     fbc=(dc(p)*fwc+dcv(p)*fmc+dg(p)*(fwg+fmg)-(fwc+fmc+fwg+fmg)*(dbck-ac))/ac;
%     fbg=fwc+fmc+fwg+fmg-fbc;
    f = @(X)fwc+fmc+fwg+fmg-myFBGV(p)-((dc(p)*fwc+dcv(p)*fmc+dg(p)*(fwg+fmg)-(fwc+fmc+fwg+fmg)*(X-ac))/ac);
    z = fzero(f,3);
    fbc=(dc(p)*fwc+dcv(p)*fmc+dg(p)*(fwg+fmg)-(fwc+fmc+fwg+fmg)*(z-ac))/ac;
%     BalCheck=fwc+fmc+fwg+fmg-fbc-myFBGV(p)
    BIGD(p)=(fbg*(dbc(j)-alpha(j))+fbc*dbc(j))/(fbg+fbc);
    fbcv(p)=fbc;
    fbgv(p)=myFBGV(p); %myFBGV(p)
    fwsv1(p)=fbc-fwc;
    fwsv2(p)=fws;%fwsi0*((pco2gca(p)/pco20)^nS)*frk*fek*fdk^0.65;
    myfws(p)=(rco2(p)/rco2(1))^nS*6.65;
    dbcv(p)=z;
end
DBC=(dc(58)*fwcv(58)+dcv(58)*fmcv(58)+dg(58)*(fwgv(58)+fmgv(58))+acv(58)*fbgv(58))/(fbcv(58)+fbgv(58))
return;
% Integrate over time with dt of 10/mend my









