function [Jconstant,Janalytic,JacFcn,JacArgs,JacPattern,jthresh,vectorized] = ...
    odejacobian(FcnHandlesUsed,ode,t0,y0,options,extras);
%ODEJACOBIAN  Helper function for the Jacobian function in ODE solvers
%    ODEJACOBIAN determines whether the Jacobian is constant and if so,
%    returns its value as JacFcn. If an analytical Jacoban is available from
%    a function, ODEJACOBIAN initiallizes JacFcn and creates a cell array of
%    additional input arguments. For numerical Jacobian, ODEJACOBIAN tries to
%    extract JacPattern and initializes jthresh and vectorized.
%
%   See also ODE15S, ODE23S, ODE23T, ODE23TB.

%   Jacek Kierzenka
%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 1.4 $  $Date: 2000/06/02 00:11:22 $


true = 1;
false = ~true;

Jconstant = strcmp(odeget(options,'JConstant','off','fast'),'on');
Janalytic = false;   % use numjac
JacFcn = [];
JacArgs = {};
JacPattern = [];
atol = odeget(options,'AbsTol',1e-6,'fast');
jthresh = zeros(size(y0))+ atol(:);
vectorized = false;   % ode is not vectorized

if FcnHandlesUsed
  JacFcn = odeget(options,'Jacobian',[],'fast');
  if ~isempty(JacFcn)
    if isnumeric(JacFcn)
      Jconstant = true;
    else
      Janalytic = true;
      JacArgs = extras;
    end
  end  
else  % ode-file used  
  Joption = odeget(options,'Jacobian','off','fast');  
  switch lower(Joption)
    case 'on'    % ode(t,y,'jacobian',p1,p2...)
      Janalytic = true;
      JacFcn = ode;
      JacArgs = [{'jacobian'} extras];  
    case 'off'   % use numjac
    otherwise
      error(['Unrecognized option ''Jacobian'': ' Joption]);
  end    
end  

if ~Janalytic   % numjac will be used
  vectorized = strcmp(odeget(options,'Vectorized','off','fast'),'on');
  if FcnHandlesUsed  
    JacPattern = odeget(options,'JPattern',[],'fast'); 
  else  % ode-file used
    JPoption = odeget(options,'JPattern','off','fast'); 
    switch lower(JPoption)
      case 'on'
        JacPattern = feval(ode,[],[],'jpattern',extras{:}); 
      case 'off'  % no pattern provided
      otherwise
        error(['Unrecognized option ''JPattern'': ' JPoption]);        
    end          
  end  
end
    
