function s = funstring(fun)
% Yield a string representing fun.

%   Mike Karr, Jacek Kierzenka, 11-19-99
%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 1.7 $  $Date: 2000/06/02 00:11:21 $

if isa(fun, 'function_handle')
  s = upper(func2str(fun));
elseif isstr(fun)
  s = upper(fun);
elseif isa(fun, 'inline')
  s = formula(fun);
else
  s = 'unknown';
end
