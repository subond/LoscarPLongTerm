function [An16,Ed15,Pe09,F761,F926,Ba11,Se10] = fpHAuth(X)

k = [001 005 006 015 028 031 073];
N = [004 001 009 013 003 042 008]-1; i=1;
An16.Age = X(k(i):k(i)+N(i),1); 
An16.pH  = X(k(i):k(i)+N(i),2);
An16.CO2 = X(k(i):k(i)+N(i),3);
An16.T   = X(k(i):k(i)+N(i),4); i=2;
Ed15.Age = X(k(i):k(i)+N(i),1); 
Ed15.pH  = X(k(i):k(i)+N(i),2);
Ed15.CO2 = X(k(i):k(i)+N(i),3);
Ed15.T   = X(k(i):k(i)+N(i),4); i=3;
Pe09.Age = X(k(i):k(i)+N(i),1); 
Pe09.pH  = X(k(i):k(i)+N(i),2);
Pe09.CO2 = X(k(i):k(i)+N(i),3);
Pe09.T   = X(k(i):k(i)+N(i),4); i=4;
F761.Age = X(k(i):k(i)+N(i),1); 
F761.pH  = X(k(i):k(i)+N(i),2);
F761.CO2 = X(k(i):k(i)+N(i),3);
F761.T   = X(k(i):k(i)+N(i),4); i=5;
F926.Age = X(k(i):k(i)+N(i),1); 
F926.pH  = X(k(i):k(i)+N(i),2);
F926.CO2 = X(k(i):k(i)+N(i),3);
F926.T   = X(k(i):k(i)+N(i),4); i=6;
Ba11.Age = X(k(i):k(i)+N(i),1); 
Ba11.pH  = X(k(i):k(i)+N(i),2);
Ba11.CO2 = X(k(i):k(i)+N(i),3);
Ba11.T   = X(k(i):k(i)+N(i),4); i=7;
Se10.Age = X(k(i):k(i)+N(i),1); 
Se10.pH  = X(k(i):k(i)+N(i),2);
Se10.CO2 = X(k(i):k(i)+N(i),3);
Se10.T   = X(k(i):k(i)+N(i),4); i=8;

%Ba11.Age = Ba11.Age*1e-3;

end

