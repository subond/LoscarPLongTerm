function out = odemxv(t,y,Mfun,Margs,v)
%ODEMXV  Helper function -- evaluates M(t,y)*v
%   Used to get d(M(t,y)*v)/dy when the property MStateDependence is 'strong'  
%
%   See also IC3DAE, ODE15S, ODE23T, ODE23TB.

%   Jacek Kierzenka, Lawrence Shampine
%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2000/06/02 00:11:22 $

out = feval(Mfun,t,y,Margs{:})*v;

