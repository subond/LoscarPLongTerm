
function ydt = CappellenDiff(t,y)
global  koa vmix k54 CPred k36 kup CPoxic CPanoxic k1011...
      k58 k59 k25 kt Fcbur
 

M14 = y(1);
M1=y(2);
M2=y(3);
M3=y(4);
M4=y(5);
M5=y(6);
M6=y(7);
M7=y(8);
M8=y(9);
M9=y(10);
M10=y(11);
M12=y(12);
M13=y(13);


F54t =  k54*vmix *M5;
F13t = CPred*F54t;
F36t = k36*F13t^2.5;

Fcbur(kt+1) = F36t;

kupt = kup;
if(t>0)
    kupt =1.5*kup;
end
F1210 = kupt*M12;
F1310 = kupt*M13;
F72=kupt*M7;
F82 = kupt*M8;
F92 = kupt*M9;

DOA = 1-koa*(vmix*M14/F13t);
x = 0.44+0.56*(DOA);
CPbur = (CPoxic*CPanoxic)/((1-DOA)*CPanoxic+DOA*CPoxic);

%Iron
F1011 = k1011*M10;
F1113 = x*F1011;
F1112 = (1-x)*F1011;


F47 = F36t/CPbur;
F45 = F54t - F47;

F58 = k58*F45^2.5;
F59 = k59*F1112;

% if(t>0)
%   k25=  2.57*1e-9*1.5;
% end

F25 = k25*M2;


% Carbon

F31 = F13t-F36t;

F61 = kupt*M6;

ydt = [F13t-(F31-(7.5/2)*F1113)-(7.5/2)*F1310-F61
    F61+F31-F13t
    F92+F82+F72-F25
    F13t-F36t-F31
    F54t-F45-F47
    F45-F54t-F58-F59+F25
    F36t-F61
    F47-F72
    F58-F82
    F59-F92
    F1310+F1210-F1011
    F1112-F1210
    F1113-F1310
    ];