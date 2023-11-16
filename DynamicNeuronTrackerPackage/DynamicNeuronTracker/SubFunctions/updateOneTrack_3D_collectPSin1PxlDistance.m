function [trackPost, frs] = updateOneTrack_3D_collectPSin1PxlDistance(xvec, yvec, zvec, ...
        PSourMat, imgArray, upbdDeformX, upbdDeformZ, bandWidthX, bandWidthZ, corrThreshold)
% updateOneTrack_3D_collectPSin1PxlDistance updates one neuron-firing
% trajectory by using a local mean image of firing events of the neuron.
%
% J Noh, 2017/11/22
% updated to use nan
% updated not to adjust to ps detected, 2018/11/20

xvec = xvec(:);
yvec = yvec(:);
zvec = zvec(:);

trackPost = struct('x', []);

trackPost.segmentStartFrame = [];
trackPost.segmentEndFrame = [];
trackPost.f = [];
trackPost.x = [];
trackPost.y = [];    
trackPost.z = [];  
trackPost.condMeanCor = [];
trackPost.totFramesSatisfied = [];    
trackPost.frIndSatisfied = [];

%% gather neighboring (<= 1 pxl) point sources and obtain avg img of the PS

xcenter = round(nanmean(xvec));
ycenter = round(nanmean(yvec));
zcenter = round(nanmean(zvec));

ptInd = nan(numel(xvec), 1);
xvec2 = xvec;
yvec2 = yvec;
zvec2 = zvec;

for fr = 1:numel(xvec)
    if ~isnan(yvec(fr)) 
        psArea = PSourMat(yvec(fr)-1:yvec(fr)+1, xvec(fr)-1:xvec(fr)+1,zvec(fr)-1:zvec(fr)+1, fr);
        if (sum(psArea(:)) == 1)            % select clear firing events

            ptInd(fr) = 1;
        else
            ptInd(fr) = sum(psArea(:));
        end

    end
end

indIntersect = (ptInd == 1);
xs = xvec2(indIntersect);
ys = yvec2(indIntersect);
zs = zvec2(indIntersect);
frs = find(indIntersect);

[avgPatch, ~] = localPatchSampling3D(xs,ys,zs,frs, bandWidthX, bandWidthZ, imgArray);

%%  patchMatchingViaSpatialCorrelation_nan3D     
 
xwin = xcenter-upbdDeformX-bandWidthX: xcenter+upbdDeformX+bandWidthX;
ywin = ycenter-upbdDeformX-bandWidthX: ycenter+upbdDeformX+bandWidthX;
zwin = zcenter-upbdDeformZ-bandWidthZ: zcenter+upbdDeformZ+bandWidthZ;

fieldCylinder = imgArray(ywin, xwin, zwin, :);

%
[dX, dY,dZ, condMeanCor, totFramesSatisfied, frIndSatisfied, maxCors] = ...
 patchMatchingViaSpatialCorrelation_nan3D(avgPatch, fieldCylinder, upbdDeformX, ...
                                          upbdDeformZ, corrThreshold);

xupdated = xcenter + dX;
yupdated = ycenter + dY;
zupdated = zcenter + dZ;

% output
trackPost.x = xupdated;
trackPost.y = yupdated;
trackPost.z = zupdated;
trackPost.condMeanCor = condMeanCor;
trackPost.totFramesSatisfied = totFramesSatisfied;
trackPost.frIndSatisfied = frIndSatisfied;
trackPost.maxCors = maxCors;
 
end
