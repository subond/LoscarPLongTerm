%% Newton-Raphson
% Calculates new d13c required for GEOCARB to produce Fbg equal to that
% predicted by LOSCAR
% myNR(x0, h, t, fbg0)
% x0 is the d13c guess
% h is the accuracy, 0.001 usually enough
% t (integer) is the time step in million years 1 = present
% fbg0 is the LOSCAR predicted fbg in 10^12 moles C

% x0 is first d13c guess (can be anything)
% fbg is the Fbg calculated by GEOCARB using a 
% new d13c in each iteration until 
% the condition is satisfied
% so GEOCARB module (function) should be returnig Fbg value
% based on d13c supplied
% we could supply annonymus function e.g. myNR(f, x0, h)
% or do it like below where the function (gcmodule) is coded already

function nr = myNR(x0, h, t, fbg0) 

    lastX= x0;
    nextX = lastX + 10*h;
    fbg = fbg0 + 10*h;
    
    while (abs(fbg0 - fbg)>h)
        fbg = gcfun12(nextX, t);
        lastX = nextX;
        nextX = lastX - (fbg-fbg0);% / myDer(@gcfun11, lastX, h);
    end
    
    %     while (abs(lastX - nextX)>h)
%         newY = gcfun11(nextX);
%         lastX = nextX;
%         nextX = lastX - newY / myDer(@gcfun11, lastX, h);
%     end
    
    nr = nextX;
end
% function yp = myDer(f,x,h)
%     yp=(f(x+h) - f(x-h)) / (2.0*h);
% 
% end
% 
% function fgb = gcmodulemock(dbc)
% % global myT
% % myT
% fgb = 2*dbc*dbc-5*dbc+1;
% end




