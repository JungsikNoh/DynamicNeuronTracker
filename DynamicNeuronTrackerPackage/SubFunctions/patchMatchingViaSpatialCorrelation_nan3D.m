function [dX, dY, dZ, condMeanCor, totFramesSatisfied, frIndSatisfied, maxCors] = ...
    patchMatchingViaSpatialCorrelation_nan3D(avgPatch, fieldCylinder, upbdDeformX, ...
    upbdDeformZ, corrThreshold)
% patchMatchingViaSpatialCorrelation_nan3D implements local patch-matching
% for a given local volume and a wider volume across time frames while
% computing maximum rolling spatial correlations and exact jittering
% amounts in X, Y and Z-directions.

 
%% corrfilter3D computes spatial corr using a filter

totFr = size(fieldCylinder, 4);

patchFilter = (avgPatch - mean(avgPatch(:))) ./ std(avgPatch(:));
maxCors = nan(totFr, 1);
rs = nan(totFr, 1);
cs = nan(totFr, 1);
ls = nan(totFr, 1);

for k = 1:totFr
    cf0 = corrfilter3D(patchFilter, fieldCylinder(:,:,:,k));
    [mc, mi] = max(cf0(:));
    [r, c, l] = ind2sub(size(cf0), mi);
    
    maxCors(k) = mc;
    rs(k) = r;
    cs(k) = c;
    ls(k) = l;
end

%% identify jitters via maximal spatial correlations

r0 = (size(cf0, 1) + 1) / 2;
c0 = (size(cf0, 2) + 1) / 2;
l0 = (size(cf0, 3) + 1) / 2;

trRow = rs;
trCol = cs;
trLay = ls;

% corr threshold 
ind1 = (maxCors < corrThreshold);
trRow(ind1) = r0;
trCol(ind1) = c0;
trLay(ind1) = l0;

% new centers need to be within a 3D-eclipse defined by upbdDeformX and Z
distSq = (trRow - r0).^2./upbdDeformX^2 + (trCol - c0).^2 ./ upbdDeformX^2 ...
            + (trLay - l0).^2 ./ upbdDeformZ^2;
ind2 = (distSq > 1); 

trRow(ind2) = r0;
trCol(ind2) = c0;
trLay(ind2) = l0; 

ind3 = ~ind1 & ~ind2;
condMeanCor = mean(maxCors(ind3));  % conditional mean of corr when it is above the threshold
totFramesSatisfied = sum(ind3);
frIndSatisfied = ind3;

%%

dY = trRow - r0;
dX = trCol - c0;
dZ = trLay - l0;

dY(~frIndSatisfied) = nan;
dX(~frIndSatisfied) = nan;
dZ(~frIndSatisfied) = nan;

end
