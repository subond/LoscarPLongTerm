function [tv,ccd] = interpCCDs(timev,data)
%INTERPCCDS Put all ccd values on the common time interval

% Take out duplicate time points   
[x, index] = unique(timev);
tv = [0:0.5:65];
ccd = interp1(x,data(index),tv);
end

