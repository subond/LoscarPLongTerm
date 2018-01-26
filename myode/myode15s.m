function [tout,yout,varargout] = myode15s(ode,tspan,y0,options,varargin)


%REZ
global kt tst Dtst yst dYst stflag;

%n = 100;
%tst    = zeros(1,n);
%yst    = zeros(1,n);

 tst(1)   = tspan(1);
% yst(1,:) = y0;
%dYst(:,1) = PEboxdif(tspan(1),y0');

DtPr = tspan/10.;
kPr  = 1;
%kt = 1;

%ODE15S Solve stiff differential equations and DAEs, variable order method.
%   [T,Y] = ODE15S(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates the
%   system of differential equations y' = f(t,y) from time T0 to TFINAL with
%   initial conditions Y0. Function ODEFUN(T,Y) must return a column vector
%   corresponding to f(t,y). Each row in the solution array Y corresponds to
%   a time returned in the column vector T. To obtain solutions at specific
%   times T0,T1,...,TFINAL (all increasing or all decreasing), use 
%   TSPAN = [T0 T1  ... TFINAL].     
%   
%   [T,Y] = ODE15S(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the ODESET function. See ODESET for details. Commonly used options
%   are scalar relative error tolerance 'RelTol' (1e-3 by default) and vector
%   of absolute error tolerances 'AbsTol' (all components 1e-6 by default).  
%   
%   [T,Y] = ODE15S(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2...) passes the additional
%   parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
%   all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
%   no options are set.    
%   
%   The Jacobian matrix df/dy is critical to reliability and efficiency. Use
%   ODESET to set 'Jacobian' to a function FJAC if FJAC(T,Y) returns the
%   Jacobian df/dy or to the matrix df/dy if the Jacobian is constant. If the
%   'Jacobian' option is not set (the default), df/dy is approximated by
%   finite differences. Set 'Vectorized' 'on' if the ODE function is coded so
%   that ODEFUN(T,[Y1 Y2 ...]) returns [ODEFUN(T,Y1) ODEFUN(T,Y2) ...]. If
%   df/dy is a sparse matrix, set 'JPattern' to the sparsity pattern of
%   df/dy, i.e., a sparse matrix S with S(i,j) = 1 if component i of f(t,y)
%   depends on component j of y, and 0 otherwise.    
%
%   ODE15S can solve problems M(t,y)*y' = f(t,y) with mass matrix M(t,y). Use
%   ODESET to set the 'Mass' property to a function MASS if MASS(T,Y) returns
%   the value of the mass matrix. If the mass matrix is constant, the matrix
%   can be used as the value of the 'Mass' option. Problems with
%   state-dependent mass matrices are more difficult. If the mass matrix does
%   not depend on the state variable Y and the function MASS is to be called
%   with one input argument T, set 'MStateDependence' to 'none'. If the mass
%   matrix depends weakly on Y, set 'MStateDependence' to 'weak' (the
%   default) and otherwise, to 'strong'. In either case the function MASS is
%   to be called with the two arguments (T,Y). If there are many differential
%   equations, it is important to exploit sparsity: Return a sparse
%   M(t,y). Either supply the sparsity pattern of df/dy using the 'JPattern'
%   property or a sparse df/dy using the Jacobian property. For strongly
%   state-dependent M(t,y), set 'MvPattern' to a sparse matrix S with S(i,j)
%   = 1 if for any k, the (i,k) component of M(t,y) depends on component j of
%   y, and 0 otherwise.    
%
%   If the mass matrix is non-singular, the solution of the problem is
%   straightforward. See examples FEM1ODE, FEM2ODE, BATONODE, or
%   BURGERSODE. If M(t0,y0) is singular, the problem is a differential-
%   algebraic equation (DAE). ODE15S solves DAEs of index 1. DAEs have
%   solutions only when y0 is consistent, i.e., there is a yp0 such that
%   M(t0,y0)*yp0 = f(t0,y0). Use ODESET to set 'MassSingular' to 'yes', 'no',
%   or 'maybe'. The default of 'maybe' causes ODE15S to test whether M(t0,y0)
%   is singular. You can provide yp0 as the value of the 'InitialSlope'
%   property. The default is the zero vector. If y0 and yp0 are not
%   consistent, ODE15S treats them as guesses, tries to compute consistent
%   values close to the guesses, and then goes on to solve the problem. See
%   examples HB1DAE or AMP1DAE.  
%
%   [T,Y,TE,YE,IE] = ODE15S(ODEFUN,TSPAN,Y0,OPTIONS) with the 'Events'
%   property in OPTIONS set to a function EVENTS, solves as above while also
%   finding where functions of (T,Y), called event functions, are zero. For
%   each function you specify whether the integration is to terminate at a
%   zero and whether the direction of the zero crossing matters. These are
%   the three vectors returned by EVENTS: [VALUE,ISTERMINAL,DIRECTION] =
%   EVENTS(T,Y). For the I-th event function: VALUE(I) is the value of the
%   function, ISTERMINAL(I)=1 if the integration is to terminate at a zero of
%   this event function and 0 otherwise. DIRECTION(I)=0 if all zeros are to
%   be computed (the default), +1 if only zeros where the event function is
%   increasing, and -1 if only zeros where the event function is
%   decreasing. Output TE is a column vector of times at which events occur.
%   Rows of YE are the corresponding solutions, and indices in vector IE
%   specify which event occurred.    
%   
%   Example
%         [t,y]=ode15s(@vdp1000,[0 3000],[2 0]);   
%         plot(t,y(:,1));
%     solves the system y' = vdp1000(t,y), using the default relative error
%     tolerance 1e-3 and the default absolute tolerance of 1e-6 for each
%     component, and plots the first component of the solution.
%
%   See also 
%     other ODE solvers:   ODE23S, ODE23T, ODE23TB, ODE45, ODE23, ODE113
%     options handling:    ODESET, ODEGET
%     output functions:    ODEPLOT, ODEPHAS2, ODEPHAS3, ODEPRINT
%     ODE examples:        VDPODE, FEM1ODE, BRUSSODE, HB1DAE
%
%   NOTE: 
%     The interpretation of the first input argument of the ODE solvers and
%     some properties available through ODESET have changed in this version
%     of MATLAB. Although we still support the v5 syntax, any new
%     functionality is available only with the new syntax. To see the v5 help
%     type in the command line  
%         more on, type ode15s, more off

%   NOTE:
%     This portion describes the v5 syntax of ODE15S.
%
%   [T,Y] = ODE15S('F',TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates the
%   system of differential equations y' = F(t,y) from time T0 to TFINAL with
%   initial conditions Y0.  'F' is a string containing the name of an ODE
%   file.  Function F(T,Y) must return a column vector.  Each row in
%   solution array Y corresponds to a time returned in column vector T.  To
%   obtain solutions at specific times T0, T1, ..., TFINAL (all increasing
%   or all decreasing), use TSPAN = [T0 T1 ... TFINAL].   
%
%   [T,Y] = ODE15S('F',TSPAN,Y0,OPTIONS) solves as above with default
%   integration parameters replaced by values in OPTIONS, an argument
%   created with the ODESET function.  See ODESET for details.  Commonly
%   used options are scalar relative error tolerance 'RelTol' (1e-3 by
%   default) and vector of absolute error tolerances 'AbsTol' (all
%   components 1e-6 by default).   
%
%   [T,Y] = ODE15S('F',TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
%   parameters P1,P2,... to the ODE file as F(T,Y,FLAG,P1,P2,...) (see
%   ODEFILE).  Use OPTIONS = [] as a place holder if no options are set.   
%
%   It is possible to specify TSPAN, Y0 and OPTIONS in the ODE file (see
%   ODEFILE).  If TSPAN or Y0 is empty, then ODE15S calls the ODE file
%   [TSPAN,Y0,OPTIONS] = F([],[],'init') to obtain any values not supplied
%   in the ODE15S argument list.  Empty arguments at the end of the call
%   list may be omitted, e.g. ODE15S('F').   
%
%   The Jacobian matrix dF/dy is critical to reliability and efficiency.
%   Use ODESET to set JConstant 'on' if dF/dy is constant.  Set Vectorized
%   'on' if the ODE file is coded so that F(T,[Y1 Y2 ...]) returns
%   [F(T,Y1) F(T,Y2) ...].  Set JPattern 'on' if dF/dy is a sparse matrix
%   and the ODE file is coded so that F([],[],'jpattern') returns a sparsity
%   pattern matrix of 1's and 0's showing the nonzeros of dF/dy.  Set
%   Jacobian 'on' if the ODE file is coded so that F(T,Y,'jacobian') returns
%   dF/dy.   
%
%   ODE15S can solve problems M(t,y)*y' = F(t,y) with a mass matrix M that
%   is nonsingular.  Use ODESET to set Mass to 'M', 'M(t)', or 'M(t,y)' if
%   the ODE file is coded so that F(T,Y,'mass') returns a constant,
%   time-dependent, or time- and state-dependent mass matrix, respectively.
%   The default value of Mass is 'none'.  
%
%   If M is singular, M(t)*y' = F(t,y) is a differential-algebraic equation
%   (DAE).  DAEs have solutions only when y0 is consistent, i.e. there is a
%   vector yp0 such that M(t0)*yp0 = f(t0,y0).  ODE15S can solve DAEs of
%   index 1 provided that M is not state dependent and y0 is sufficiently
%   close to being consistent.  You can use ODESET to set MassSingular to
%   'yes', 'no', or 'maybe'.  The default of 'maybe' causes ODE15S to test
%   whether the problem is a DAE. If it is, ODE15S treats y0 as a guess,
%   attempts to compute consistent initial conditions that are close to y0,
%   and goes on to solve the problem. When solving DAEs, it is very
%   advantageous to formulate the problem so that M is diagonal (a
%   semi-explicit DAE).  
%
%   [T,Y,TE,YE,IE] = ODE15S('F',TSPAN,Y0,OPTIONS) with the Events property
%   in OPTIONS set to 'on', solves as above while also locating zero
%   crossings of an event function defined in the ODE file.  The ODE file
%   must be coded so that F(T,Y,'events') returns appropriate information.
%   See ODEFILE for details.  Output TE is a column vector of times at which
%   events occur, rows of YE are the corresponding solutions, and indices in
%   vector IE specify which event occurred.
%   
%   See also ODEFILE.

%   ODE15S is a quasi-constant step size implementation in terms of backward
%   differences of the Klopfenstein-Shampine family of Numerical
%   Differentiation Formulas of orders 1-5. The natural "free" interpolants
%   are used. Local extrapolation is not done. By default, Jacobians are
%   generated numerically.  

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997, and in
%   Solving Index-1 DAEs in MATLAB and Simulink, L. F. Shampine,
%   M. W. Reichelt, and J. A. Kierzenka, SIAM Review, 41-3, 1999. 

%   Mark W. Reichelt, Lawrence F. Shampine, and Jacek Kierzenka, 12-18-97
%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 1.76 $  $Date: 2000/06/02 00:11:24 $

true = 1;

false = ~true;

stats = struct('nsteps',0,'nfailed',0,'nfevals',0,...
               'npds',0,'ndecomps',0,'nsolves',0); 

if nargin < 4
  options = [];
  if nargin < 3
    y0 = [];
    if nargin < 2
      tspan = [];
      if nargin < 1
        error('Not enough input arguments.  See ODE15s.');
      end  
    end
  end
end

FcnHandlesUsed = isa(ode,'function_handle');

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, ...
          options, atol, rtol, threshold, normcontrol, normy, hmax, htry]  ...
         = odearguments(FcnHandlesUsed,'ode15s', ode, tspan, y0, options, varargin);
stats.nfevals = stats.nfevals + 1;
one2neq = (1:neq)';

% Handle the output
if nargout > 0
  outfun = odeget(options,'OutputFcn',[],'fast');
else
  outfun = odeget(options,'OutputFcn',@odeplot,'fast');
end
if isempty(outfun)
  haveoutfun = false;
else
  haveoutfun = true;
  outputs = odeget(options,'OutputSel',one2neq,'fast');
  if isa(outfun,'function_handle')
    outputArgs = [{''},varargin];
    outputArgs1 = varargin;
  else       % With v5 syntax do not pass additional input arguments to outfun  
    outputArgs = {};      
    outputArgs1 = {};
  end  
end
refine = odeget(options,'Refine',1,'fast');
if ntspan > 2
  outflag = 1;                          % output only at tspan points
elseif refine <= 1
  outflag = 2;                          % computed points, no refinement
else
  outflag = 3;                          % computed points, with refinement
  S = (1:refine-1)' / refine;
end
printstats = strcmp(odeget(options,'Stats','off','fast'),'on');

% Handle the event function 
[haveeventfun,eventfun,eventargs,valt,teout,yeout,ieout] = ...
    odeevents(FcnHandlesUsed,ode,t0,y0,options,varargin);

% Handle the mass matrix
[Mtype, Mfun, Margs, Mt, Mvs, Mfac] =  odemass(FcnHandlesUsed,ode,t0,y0,...
                                         options,varargin);

% Handle the Jacobian
[Jconstant,Janalytic,Jac,Jargs,Js,jthresh,vectorized] = ...
    odejacobian(FcnHandlesUsed,ode,t0,y0,options,varargin);
% if not set via 'options', initialize constant Jacobian here
if Jconstant 
  if isempty(Jac) % use numjac
    [Jac,fac,g,nF] = numjac(ode,t0,y0,f0,jthresh,[],...
                          vectorized,Js,[],args{:});
    stats.nfevals = stats.nfevals + nF;
    stats.npds = stats.npds + 1;
  elseif ~isa(Jac,'numeric')  % not been set via 'options'  
    Jac = feval(Jac,t0,y0,Jargs{:}); % replace by its value
    stats.npds = stats.npds + 1;
  end
end
    
t = t0;
y = y0;

yp0_OK = false;
DAE = false;
RowScale = [];
if Mtype > 0
  nz = nnz(Mt);
  if nz == 0
    error('The mass matrix must have some non-zero entries.')
  end
   
  Msingular = odeget(options,'MassSingular','maybe','fast');
  switch Msingular
    case 'no',     DAE = false;
    case 'yes',    DAE = true;
    case 'maybe',  DAE = (eps*nz*condest(Mt) > 1);       
  end
   
  if DAE
    yp0 = odeget(options,'InitialSlope',[],'fast');
    if isempty(yp0)
      yp0_OK = false;
      yp0 = zeros(neq,1);  
    else
      yp0 = yp0(:);
      if length(yp0) ~= neq
        error('y0 and yp0 have different lengths.');
      end
      % Test if (y0,yp0) are consistent enough to accept.
      yp0_OK = (norm(Mt*yp0 - f0) <= 1e-3*rtol*max(norm(Mt*yp0),norm(f0)));
    end   
    if ~yp0_OK           % Must compute ICs, so classify them.
      if Mtype >= 3  % state dependent
        ICtype = 3;
      else  % M, M(t)
        % Test for a diagonal mass matrix.
        [r,c] = find(Mt);
        if isequal(r,c)   % diagonal
          ICtype = 1;
        elseif ~issparse(Mt) % not diagonal but full
          ICtype = 2;
        else  % sparse, not diagonal
          ICtype = 3;
        end
      end      
    end
  end
end
Mcurrent = true;
Mtnew = Mt;

maxk = odeget(options,'MaxOrder',5,'fast');
bdf = strcmp(odeget(options,'BDF','off','fast'),'on');

% Initialize method parameters.
G = [1; 3/2; 11/6; 25/12; 137/60];
if bdf
  alpha = [0; 0; 0; 0; 0];
else
  alpha = [-37/200; -1/9; -0.0823; -0.0415; 0];
end
invGa = 1 ./ (G .* (1 - alpha));
erconst = alpha .* G + (1 ./ (2:6)');
difU = [ -1, -2, -3, -4,  -5;           % difU is its own inverse!
          0,  1,  3,  6,  10;
          0,  0, -1, -4, -10;
          0,  0,  0,  1,   5;
          0,  0,  0,  0,  -1 ];
maxK = 1:maxk;
[kJ,kI] = meshgrid(maxK,maxK);
difU = difU(maxK,maxK);
maxit = 4;

warnstat = warning;

% Get the initial slope yp. For DAEs the default is to compute
% consistent initial conditions.
if DAE & ~yp0_OK
  if ICtype < 3
    [y,yp,f0,dfdy,nFE,nPD,fac,g] = daeic12(ode,args,t,ICtype,Mt,y,yp0,f0,rtol,...
         Jconstant,Janalytic,Jac,Jargs,jthresh,vectorized,Js); 
  else    
    [y,yp,f0,dfdy,nFE,nPD,fac,g] = daeic3(ode,args,tspan,htry,Mtype,Mt,Mfun,...
         Margs,Mvs,y,yp0,f0,rtol,Jconstant,Janalytic,Jac,Jargs,jthresh,vectorized,Js);   
  end  
  stats.nfevals = stats.nfevals + nFE;
  stats.npds = stats.npds + nPD;
  if Mtype >= 3
    Mt = feval(Mfun,t,y,Margs{:});
    Mtnew = Mt;
    Mcurrent = true;
  end
else
  if Mtype == 0 
    yp = f0;
  elseif DAE & yp0_OK
    yp = yp0;
  else
    [L,U] = lu(Mt);
    yp = U \ (L \ f0);
    stats.ndecomps = stats.ndecomps + 1;              
    stats.nsolves = stats.nsolves + 1;                
  end
    
  if Jconstant
    dfdy = Jac;
  elseif Janalytic
    dfdy = feval(Jac,t,y,Jargs{:});     
    stats.npds = stats.npds + 1;                            
  else
    [dfdy,fac,g,nF] = numjac(ode,t,y,f0,jthresh,[],...
                            vectorized,Js,[],args{:});
    stats.nfevals = stats.nfevals + nF;    
    stats.npds = stats.npds + 1;                            
  end     
end
Jcurrent = true;

% hmin is a small number such that t + hmin is clearly different from t in
% the working precision, but with this definition, it is 0 if t = 0.
hmin = 16*eps*abs(t);

if isempty(htry)
  % Compute an initial step size h using yp = y'(t).
  if normcontrol
    wt = max(normy,threshold);
    rh = 1.25 * (norm(yp) / wt) / sqrt(rtol);  % 1.25 = 1 / 0.8
  else
    wt = max(abs(y),threshold);
    rh = 1.25 * norm(yp ./ wt,inf) / sqrt(rtol);
  end
  absh = min(hmax, abs(tspan(next) - t));
  if absh * rh > 1
    absh = 1 / rh;
  end
  absh = max(absh, hmin);
  
  if ~DAE
    % The error of BDF1 is 0.5*h^2*y''(t), so we can determine the optimal h.
    h = tdir * absh;
    tdel = (t + tdir*min(sqrt(eps)*max(abs(t),abs(t+h)),absh)) - t;
    f1 = feval(ode,t+tdel,y,args{:});
    stats.nfevals = stats.nfevals + 1;                
    dfdt = (f1 - f0) ./ tdel;
    if normcontrol
      if Mtype > 0
        rh = 1.25 * sqrt(0.5 * (norm(U \ (L \ (dfdt + dfdy*yp))) / wt) / rtol);
      else
        rh = 1.25 * sqrt(0.5 * (norm(dfdt + dfdy*yp) / wt) / rtol);
      end
    else
      if Mtype > 0
        rh = 1.25*sqrt(0.5*norm((U \ (L \ (dfdt+dfdy*yp))) ./ wt,inf) / rtol);
      else
        rh = 1.25 * sqrt(0.5 * norm((dfdt + dfdy*yp) ./ wt,inf) / rtol);
      end
    end
    absh = min(hmax, abs(tspan(next) - t));
    if absh * rh > 1
      absh = 1 / rh;
    end
    absh = max(absh, hmin);
  end
else
  absh = min(hmax, max(hmin, htry));
end
h = tdir * absh;

% Initialize.
k = 1;                                  % start at order 1 with BDF1
K = 1;                                  % K = 1:k
klast = k;
abshlast = absh;

dif = zeros(neq,maxk+2);
dif(:,1) = h * yp;

hinvGak = h * invGa(k);
nconhk = 0;                             % steps taken with current h and k

Miter = Mt - hinvGak * dfdy;

% Account for strongly state-dependent mass matrix.
if Mtype == 4
  psi = dif(:,K) * (G(K) * invGa(k));
  [dMpsidy,Mfac] = numjac(@odemxv,t,y,Mt*psi,jthresh,[],0,Mvs,[],Mfun,Margs,psi);
  Miter = Miter + dMpsidy;
end

% Use explicit scaling of the equations when solving DAEs.
if DAE
  RowScale = 1 ./ max(abs(Miter),[],2);
  Miter = sparse(one2neq,one2neq,RowScale) * Miter;
end
[L,U] = lu(Miter);
stats.ndecomps = stats.ndecomps + 1;                
havrate = false;

% Initialize the output function.
if haveoutfun
  feval(outfun,[t tfinal],y(outputs),'init',outputArgs1{:});
end

% Allocate memory if we're generating output.
if nargout > 0
  if ntspan > 2                         % output only at tspan points
    tout = zeros(ntspan,1);
    yout = zeros(ntspan,neq);
  else                                  % alloc in chunks
    chunk = min(max(100,50*refine),floor((2^13)/neq));
    tout = zeros(chunk,1);
    yout = zeros(chunk,neq);
  end
  nout = 1;
  tout(nout) = t;
  yout(nout,:) = y.';
end

% THE MAIN LOOP

done = false;
while ~done
  
  hmin = 16*eps*abs(t);
  absh = min(hmax, max(hmin, absh));
  h = tdir * absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  if (absh ~= abshlast) | (k ~= klast)
    difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
    dif(:,K) = dif(:,K) * difRU(K,K);

    hinvGak = h * invGa(k);
    nconhk = 0;
    Miter = Mt - hinvGak * dfdy;
    if Mtype == 4
      Miter = Miter + dMpsidy;
    end    
    if DAE
      RowScale = 1 ./ max(abs(Miter),[],2);
      Miter = sparse(one2neq,one2neq,RowScale) * Miter;
    end
    [L,U] = lu(Miter);
    stats.ndecomps = stats.ndecomps + 1;            
    havrate = false;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true                            % Evaluate the formula.
      
    gotynew = false;                    % is ynew evaluated yet?
    while ~gotynew

      % Compute the constant terms in the equation for ynew.
      psi = dif(:,K) * (G(K) * invGa(k));

      % Predict a solution at t+h.
      tnew = t + h;
      if done
        tnew = tfinal;   % Hit end point exactly.
      end
      pred = y + sum(dif(:,K),2);
      ynew = pred;
      
      % The difference, difkp1, between pred and the final accepted 
      % ynew is equal to the backward difference of ynew of order
      % k+1. Initialize to zero for the iteration to compute ynew.
      difkp1 = zeros(neq,1); 
      if normcontrol
        normynew = norm(ynew);
        invwt = 1 / max(max(normy,normynew),threshold);
        minnrm = 100*eps*(normynew * invwt);
      else
        invwt = 1 ./ max(max(abs(y),abs(ynew)),threshold);
        minnrm = 100*eps*norm(ynew .* invwt,inf);
      end

      % Mtnew is required in the RHS function evaluation.
      if Mtype == 2  % state-independent
        if FcnHandlesUsed
          Mtnew = feval(Mfun,tnew,Margs{:}); % mass(t,p1,p2...)
        else                                     
          Mtnew = feval(Mfun,tnew,ynew,Margs{:}); % mass(t,y,'mass',p1,p2...)
        end
      end
      
      % Iterate with simplified Newton method.
      tooslow = false;
      for iter = 1:maxit
        if Mtype >= 3 
          Mtnew = feval(Mfun,tnew,ynew,Margs{:}); % state-dependent
        end
        rhs = hinvGak*feval(ode,tnew,ynew,args{:}) -  Mtnew*(psi+difkp1);
        if DAE                          % Account for row scaling.
          rhs = RowScale .* rhs;
        end
        warning('off');
        del = U \ (L \ rhs);
        warning(warnstat);
        if normcontrol
          newnrm = norm(del) * invwt;
        else
          newnrm = norm(del .* invwt,inf);
        end
        difkp1 = difkp1 + del;
        ynew = pred + difkp1;
        
        if newnrm <= minnrm
          gotynew = true;
          break;
        elseif iter == 1
          if havrate
            errit = newnrm * rate / (1 - rate);
            if errit <= 0.05*rtol       % More stringent when using old rate.
              gotynew = true;
              break;
            end
          else
            rate = 0;
          end
        elseif newnrm > 0.9*oldnrm
          tooslow = true;
          break;
        else
          rate = max(0.9*rate, newnrm / oldnrm);
          havrate = true;                 
          errit = newnrm * rate / (1 - rate);
          if errit <= 0.5*rtol             
            gotynew = true;
            break;
          elseif iter == maxit            
            tooslow = true;
            break;
          elseif 0.5*rtol < errit*rate^(maxit-iter)
            tooslow = true;
            break;
          end
        end
        
        oldnrm = newnrm;
      end                               % end of Newton loop
      stats.nfevals = stats.nfevals + iter;         
      stats.nsolves = stats.nsolves + iter;         
      
      if tooslow
        stats.nfailed = stats.nfailed + 1;          
        % Speed up the iteration by forming new linearization or reducing h.
        if ~Jcurrent
          if Janalytic
            dfdy = feval(Jac,t,y,Jargs{:});
          else
            f0 = feval(ode,t,y,args{:});
            [dfdy,fac,g,nF] = ...
                numjac(ode,t,y,f0,jthresh,fac,vectorized,Js,g,args{:});
            stats.nfevals = stats.nfevals + nF + 1; 
          end             
          stats.npds = stats.npds + 1;            
          Jcurrent = true;
          if ~Mcurrent
            Mt = feval(Mfun,t,y,Margs{:});
            Mcurrent = true;
            if Mtype == 4
              [dMpsidy,Mfac] = numjac(@odemxv,t,y,Mt*psi,jthresh,...
                                      Mfac,0,Mvs,[],Mfun,Margs,psi);  
            end
          end                       
        elseif absh <= hmin
          msg = sprintf(['Failure at t=%e.  Unable to meet integration ' ...
                         'tolerances without reducing the step size below ' ...
                         'the smallest value allowed (%e) at time t.\n'], ...
                        t,hmin);
          warning(msg);
          if haveoutfun
            feval(outfun,[],[],'done',outputArgs1{:});
          end
          if printstats                 % print cost statistics
            fprintf('%g successful steps\n', stats.nsteps);
            fprintf('%g failed attempts\n', stats.nfailed);
            fprintf('%g function evaluations\n', stats.nfevals);
            fprintf('%g partial derivatives\n', stats.npds);
            fprintf('%g LU decompositions\n', stats.ndecomps);
            fprintf('%g solutions of linear systems\n', stats.nsolves);
          end
          if nargout > 0
            tout = tout(1:nout);
            yout = yout(1:nout,:);
            if haveeventfun
              varargout{1} = teout;
              varargout{2} = yeout;
              varargout{3} = ieout;
              if ~FcnHandlesUsed
                varargout{4} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                                stats.npds; stats.ndecomps; stats.nsolves];
              end
            else
              if ~FcnHandlesUsed
                varargout{1} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                                stats.npds; stats.ndecomps; stats.nsolves];   
              end
            end
          end
          return;
        else
          abshlast = absh;
          absh = max(0.3 * absh, hmin);
          h = tdir * absh;
          done = false;

          difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
          dif(:,K) = dif(:,K) * difRU(K,K);
          
          hinvGak = h * invGa(k);
          nconhk = 0;
        end
        Miter = Mt - hinvGak * dfdy;
        if Mtype == 4
          Miter = Miter + dMpsidy;
        end
        if DAE
          RowScale = 1 ./ max(abs(Miter),[],2);
          Miter = sparse(one2neq,one2neq,RowScale) * Miter;
        end
        [L,U] = lu(Miter);
        stats.ndecomps = stats.ndecomps + 1;        
        havrate = false;
      end   
    end     % end of while loop for getting ynew
    
    % difkp1 is now the backward difference of ynew of order k+1.
    if normcontrol
      err = (norm(difkp1) * invwt) * erconst(k);
    else
      err = norm(difkp1 .* invwt,inf) * erconst(k);
    end
    
    if err > rtol                       % Failed step
      stats.nfailed = stats.nfailed + 1;            
      if absh <= hmin
        msg = sprintf(['Failure at t=%e.  Unable to meet integration ' ...
                       'tolerances without reducing the step size below ' ...
                       'the smallest value allowed (%e) at time t.\n'],t,hmin);
        warning(msg);
        if haveoutfun
          feval(outfun,[],[],'done',outputArgs1{:});
        end
        if printstats                   % print cost statistics
          fprintf('%g successful steps\n', stats.nsteps);
          fprintf('%g failed attempts\n', stats.nfailed);
          fprintf('%g function evaluations\n', stats.nfevals);
          fprintf('%g partial derivatives\n', stats.npds);
          fprintf('%g LU decompositions\n', stats.ndecomps);
          fprintf('%g solutions of linear systems\n', stats.nsolves);
        end
        if nargout > 0
          tout = tout(1:nout);
          yout = yout(1:nout,:);
          if haveeventfun
            varargout{1} = teout;
            varargout{2} = yeout;
            varargout{3} = ieout;
            if ~FcnHandlesUsed
              varargout{4} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                              stats.npds; stats.ndecomps; stats.nsolves];
            end  
          else
            if ~FcnHandlesUsed
              varargout{1} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                              stats.npds; stats.ndecomps; stats.nsolves];
            end
          end
        end
        return;
      end
      
      abshlast = absh;
      if nofailed
        nofailed = false;
        hopt = absh * max(0.1, 0.833*(rtol/err)^(1/(k+1))); % 1/1.2
        if k > 1
          if normcontrol
            errkm1 = (norm(dif(:,k) + difkp1) * invwt) * erconst(k-1);
          else
            errkm1 = norm((dif(:,k) + difkp1) .* invwt,inf) * erconst(k-1);
          end
          hkm1 = absh * max(0.1, 0.769*(rtol/errkm1)^(1/k)); % 1/1.3
          if hkm1 > hopt
            hopt = min(absh,hkm1);      % don't allow step size increase
            k = k - 1;
            K = 1:k;
          end
        end
        absh = max(hmin, hopt);
      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      if absh < abshlast
        done = false;
      end
      
      difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
      dif(:,K) = dif(:,K) * difRU(K,K);
      
      hinvGak = h * invGa(k);
      nconhk = 0;
      Miter = Mt - hinvGak * dfdy;
      if Mtype == 4
        Miter = Miter + dMpsidy;
      end      
      if DAE
        RowScale = 1 ./ max(abs(Miter),[],2);
        Miter = sparse(one2neq,one2neq,RowScale) * Miter;
      end
      [L,U] = lu(Miter);
      stats.ndecomps = stats.ndecomps + 1;          
      havrate = false;
      
    else                                % Successful step
      break;
      
    end
  end % while true
  stats.nsteps = stats.nsteps + 1;                  
  
  dif(:,k+2) = difkp1 - dif(:,k+1);
  dif(:,k+1) = difkp1;
  for j = k:-1:1
    dif(:,j) = dif(:,j) + dif(:,j+1);
  end
  
  tstep = tnew;
  %ystep = ynew;
  
  % REZ
   kt        = kt+1;
   tst(kt)   = tstep;
  Dtst       = tstep-tst(kt-1);
   %yst(kt,:) = ystep';
   %stflag = 1;
   %dYst(:,kt) = PEboxdif(tstep,ystep);
   %stflag = 0;
   
   if(tstep > kPr*DtPr)
       tstep
       kPr = kPr+1;
   end;       
       
   
  if haveeventfun
    [te,ye,ie,valt,stop] = odezero(@ntrp15s,eventfun,eventargs,valt,...
                                   t,y,tnew,ynew,t0,h,dif,k);

    nte = length(te);
    if nte > 0
      if nargout > 2
        teout = [teout; te];
        yeout = [yeout; ye.'];
        ieout = [ieout; ie];
      end
      if stop                           % stop on a terminal event
        tnew = te(nte);
        ynew = ye(:,nte);
        done = true;
      end
    end
  end
  
  if nargout > 0
    oldnout = nout;
    if outflag == 2                     % computed points, no refinement
      nout = nout + 1;
      if nout > length(tout)
        tout = [tout; zeros(chunk,1)];
        yout = [yout; zeros(chunk,neq)];
      end
      tout(nout) = tnew;
      yout(nout,:) = ynew.';
    elseif outflag == 3                 % computed points, with refinement
      nout = nout + refine;
      if nout > length(tout)
        tout = [tout; zeros(chunk,1)];  % requires chunk >= refine
        yout = [yout; zeros(chunk,neq)];
      end
      i = oldnout+1:nout-1;
      tout(i) = t + (tnew-t)*S;
      yout(i,:) = ntrp15s(tout(i),[],[],tstep,ystep,h,dif,k).';
      tout(nout) = tnew;
      yout(nout,:) = ynew.';
    elseif outflag == 1                 % output only at tspan points
      while next <= ntspan
        if tdir * (tnew - tspan(next)) < 0
          if haveeventfun & done
            nout = nout + 1;
            tout(nout) = tnew;
            yout(nout,:) = ynew.';
          end
          break;
        elseif tnew == tspan(next)
          nout = nout + 1;
          tout(nout) = tnew;
          yout(nout,:) = ynew.';
          next = next + 1;
          break;
        end
        nout = nout + 1;                % tout and yout are already allocated
        tout(nout) = tspan(next);
        yout(nout,:) = ntrp15s(tspan(next),[],[],tstep,ystep,h,dif,k).';
        next = next + 1;
      end
    end
    
    if haveoutfun
      i = oldnout+1:nout;
      if ~isempty(i) & (feval(outfun,tout(i),yout(i,outputs).',outputArgs{:}) == 1)
        tout = tout(1:nout);
        yout = yout(1:nout,:);
        if haveeventfun
          varargout{1} = teout;
          varargout{2} = yeout;
          varargout{3} = ieout;
          if ~FcnHandlesUsed
            varargout{4} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                            stats.npds; stats.ndecomps; stats.nsolves];
          end
        else
          if ~FcnHandlesUsed
            varargout{1} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                            stats.npds; stats.ndecomps; stats.nsolves];
          end  
        end
        return;
      end
    end
    
  elseif haveoutfun
    if outflag == 2
      if feval(outfun,tnew,ynew(outputs),outputArgs{:}) == 1
        return;
      end
    elseif outflag == 3                 % computed points, with refinement
      tinterp = t + (tnew-t)*S;
      yinterp = ntrp15s(tinterp,[],[],tstep,ystep,h,dif,k);
      if feval(outfun,[tinterp; tnew],[yinterp(outputs,:), ynew(outputs)],outputArgs{:}) == 1
        return;
      end
    elseif outflag == 1                 % output only at tspan points
      ninterp = 0;
      while next <= ntspan 
        if tdir * (tnew - tspan(next)) < 0
          if haveeventfun & done
            ninterp = ninterp + 1;
            tinterp(ninterp,1) = tnew;
            yinterp(:,ninterp) = ynew;
          end
          break;
        elseif tnew == tspan(next)
          ninterp = ninterp + 1;
          tinterp(ninterp,1) = tnew;
          yinterp(:,ninterp) = ynew;
          next = next + 1;
          break;
        end
        ninterp = ninterp + 1;
        tinterp(ninterp,1) = tspan(next);
        yinterp(:,ninterp) = ntrp15s(tspan(next),[],[],tstep,ystep,h,dif,k);
        next = next + 1;
      end
      if ninterp > 0
        if feval(outfun,tinterp(1:ninterp),yinterp(outputs,1:ninterp),outputArgs{:}) == 1
          return;
        end
      end
    end
  end
  
  klast = k;
  abshlast = absh;
  nconhk = min(nconhk+1,maxk+2);
  if nconhk >= k + 2
    temp = 1.2*(err/rtol)^(1/(k+1));
    if temp > 0.1
      hopt = absh / temp;
    else
      hopt = 10*absh;
    end
    kopt = k;
    if k > 1
      if normcontrol
        errkm1 = (norm(dif(:,k)) * invwt) * erconst(k-1);
      else
        errkm1 = norm(dif(:,k) .* invwt,inf) * erconst(k-1);
      end
      temp = 1.3*(errkm1/rtol)^(1/k);
      if temp > 0.1
        hkm1 = absh / temp;
      else
        hkm1 = 10*absh;
      end
      if hkm1 > hopt 
        hopt = hkm1;
        kopt = k - 1;
      end
    end
    if k < maxk
      if normcontrol
        errkp1 = (norm(dif(:,k+2)) * invwt) * erconst(k+1);
      else
        errkp1 = norm(dif(:,k+2) .* invwt,inf) * erconst(k+1);
      end
      temp = 1.4*(errkp1/rtol)^(1/(k+2));
      if temp > 0.1
        hkp1 = absh / temp;
      else
        hkp1 = 10*absh;
      end
      if hkp1 > hopt 
        hopt = hkp1;
        kopt = k + 1;
      end
    end
    if hopt > absh
      absh = hopt;
      if k ~= kopt
        k = kopt;
        K = 1:k;
      end
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  if normcontrol
    normy = normynew;
  end
  Jcurrent = Jconstant;
  switch Mtype
  case {0,1}
    Mcurrent = true;                    % Constant mass matrix I or M.
  case 2
    % M(t) has already been evaluated at tnew in Mtnew.
    Mt = Mtnew;
    Mcurrent = true;
  case {3,4}  % state dependent
    % M(t,y) has not yet been evaluated at the accepted ynew.
    Mcurrent = false;
  end
  
end % while ~done

if haveoutfun
  feval(outfun,[],[],'done',outputArgs1{:});
end

if printstats                           % print cost statistics
  fprintf('%g successful steps\n', stats.nsteps);
  fprintf('%g failed attempts\n', stats.nfailed);
  fprintf('%g function evaluations\n', stats.nfevals);
  fprintf('%g partial derivatives\n', stats.npds);
  fprintf('%g LU decompositions\n', stats.ndecomps);
  fprintf('%g solutions of linear systems\n', stats.nsolves);
end

if nargout > 0
  tout = tout(1:nout);
  yout = yout(1:nout,:);
  if haveeventfun
    varargout{1} = teout;
    varargout{2} = yeout;
    varargout{3} = ieout;
    if ~FcnHandlesUsed
      varargout{4} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                      stats.npds; stats.ndecomps; stats.nsolves];
    end  
  else
    if ~FcnHandlesUsed
      varargout{1} = [stats.nsteps; stats.nfailed; stats.nfevals; ...
                      stats.npds; stats.ndecomps; stats.nsolves];
    end  
  end
end

% --------------------------------------------------------------------------

function yinterp = ntrp15s(tinterp,t,y,tnew,ynew,h,dif,k)
%NTRP15S Interpolation helper function for ODE15S.
%   YINTERP = NTRP15S(TINTERP,T,Y,TNEW,YNEW,H,DIF,K) uses data computed in
%   ODE15S to approximate the solution at time TINTERP.

s = ((tinterp - tnew) / h)';        % may be a row vector

if k == 1
  yinterp = ynew(:,ones(length(tinterp),1)) + dif(:,1) * s;
else                    % cumprod collapses vectors
  K = (1:k)';
  kI = K(:,ones(length(tinterp),1));
  yinterp = ynew(:,ones(length(tinterp),1)) + ...
      dif(:,K) * cumprod((s(ones(k,1),:)+kI-1)./kI);
end
