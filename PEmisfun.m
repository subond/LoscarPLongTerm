%-------------------------------------
% function y = PEmisfun(t)
%
% returns  Fem (CO2 emissions at time t)
%      
%-------------------------------------
function [y] = PEmisfun(t)

global tem em;   

t0 = tem(1);
te = tem(end);

if    (t > t0  & t < te)
  y = interp1(tem,em,t,'linear');
else
  y = 0.;
end;   


return;

