% Campbell et al 2018 CCD data
% extracting min and max for each simulation to create a range
% of CCD value

% Pacific data
function [x2,inBetween] = campbell18CCD()
campCCDdata_p=csvread('dat/Cenozoicd13c/CenozoicCCD/campbell18pac.csv');
%time
campT_p=campCCDdata_p(:,1);

% Indian data
campCCDdata_i=csvread('dat/Cenozoicd13c/CenozoicCCD/campbell18ind.csv');
%time
campT_i=campCCDdata_i(:,1);

ccdRange=zeros(length(campT_p),3);
ccdRange_i=zeros(length(campT_p),3);
for i=1:length(campT_p)
    ccdMin=min(campCCDdata_p(i,[4,6,8,10]));
    ccdMax=max(campCCDdata_p(i,[4,6,8,10]));
    ccdMin_i=min(campCCDdata_i(i,[4,6,8,10]));
    ccdMax_i=max(campCCDdata_i(i,[4,6,8,10]));
ccdRange(i,:)=[campT_p(i), ccdMin, ccdMax];
ccdRange_i(i,:)=[campT_i(i), ccdMin_i, ccdMax_i];
end

% figure
% plot(campT_p,ccdRange(:,2))
% hold on
% plot(campT_p,ccdRange(:,3),'r')
% hold off
% 
% figure
% plot(campT_i,ccdRange_i(:,2))
% hold on
% plot(campT_i,ccdRange_i(:,3),'r')
% hold off


x = [campT_p,campT_p];        % repeat x values
yy = [ccdRange(:,2),ccdRange(:,3)];   % vector of upper & lower boundaries
% figure
% fill(x,yy,'b')    % fill area defined by x & yy in blue

x = campT_p';
curve1 = ccdRange(:,2)';
curve2 = ccdRange(:,3)';

% figure
% plot(x, curve1, 'r', 'LineWidth', 2);
% hold on;
% plot(x, curve2, 'b', 'LineWidth', 2);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

% figure
% fill(x2, inBetween, [0 1 0 0.5]);
