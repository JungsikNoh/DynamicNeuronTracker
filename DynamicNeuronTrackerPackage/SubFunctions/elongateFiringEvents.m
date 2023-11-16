function [tracksIterX1, tracksIterY1, tracksIterZ1] = ...
    elongateFiringEvents(tracksX1_Br4, tracksY1_Br4, tracksZ1_Br4, ...
                        PSourMat2, imgArray, par2, corrThreshold, outFolder)
% elongateFiringEvents Elongates single firing events over the whole time
% frames. Then it selects tracks that are unique and have at least 
% the minimum number of bright point sources.


%% parse

upbdDeformX = par2.upbdDeformX ;
upbdDeformZ = par2.upbdDeformZ ;
bandWidthX = par2.bandWidthX ;
bandWidthZ = par2.bandWidthZ ;
minFiringFramesOfNeurons = par2.minFiringFramesOfNeurons;
imgMarginSizeX = par2.imgMarginSizeX;
imgMarginSizeZ = par2.imgMarginSizeZ;

%% updateOneTrack_3D_collectPSin1PxlDistance

totTr = size(tracksX1_Br4, 1);
disp(['== Total number of firing events: ', num2str(totTr)])
numPS = sum(PSourMat2(:));
disp(['== num Point Sources: ', num2str(numPS)])
%
tracksX1_out = nan(size(tracksX1_Br4));
tracksY1_out = nan(size(tracksX1_Br4));
tracksZ1_out = nan(size(tracksX1_Br4));

% assessment of the update
iterStatus0 = struct('condMeanCor', [], 'totFramesSatisfied', [], ...
        'iterMerged', 0, 'iterConverged', 0, 'L2distPerfr', nan, 'hamdist', nan);
for k = 1:totTr
    iterStatus0(k).condMeanCor = nan; 
    iterStatus0(k).totFramesSatisfied = nan; 
    iterStatus0(k).iterMerged = 0;
    iterStatus0(k).iterConverged = 0;
    iterStatus0(k).L2distPerfr = NaN;
    iterStatus0(k).hamdist = NaN;
end

% boundary coordinates
min0 = upbdDeformX + bandWidthX + 1 + imgMarginSizeX;
minZ0 = upbdDeformZ + bandWidthZ + 1 + imgMarginSizeZ;

xmax0 = size(imgArray, 2) - upbdDeformX - bandWidthX - imgMarginSizeX;
ymax0 = size(imgArray, 1) - upbdDeformX - bandWidthX - imgMarginSizeX;
zmax0 = size(imgArray, 3) - upbdDeformZ - bandWidthZ - imgMarginSizeZ;

parfor itr = 1:totTr
    
    fprintf('%05d ', itr);
    if (mod(itr, 50) == 0); fprintf('\n'); end
    
    triX = tracksX1_Br4(itr, :)';
    triY = tracksY1_Br4(itr, :)';
    triZ = tracksZ1_Br4(itr, :)';    
    
    if (iterStatus0(itr).iterMerged == 0) && (iterStatus0(itr).iterConverged == 0)

        xcenter = round(nanmean(triX));
        ycenter = round(nanmean(triY));
        zcenter = round(nanmean(triZ));

        if (min(xcenter, ycenter) >= min0) && (zcenter >= minZ0) ...
                && (xcenter<=xmax0) && (ycenter<=ymax0) && (zcenter<=zmax0)         

            [trackPost, ~] = updateOneTrack_3D_collectPSin1PxlDistance(...
                triX, triY, triZ, PSourMat2, imgArray, ...
                upbdDeformX, upbdDeformZ, bandWidthX, bandWidthZ, corrThreshold);
            
            dvec = sqrt( (triX - trackPost.x).^2 + (triY-trackPost.y).^2 + (triZ-trackPost.z).^2 );
            trackPost.L2distPerfr = nanmean(dvec);

            nanhamXYZ = nanhamdistXYZ_ignoreAllnan([triX', triY', triZ'], ...
                [trackPost.x', trackPost.y', trackPost.z']);
            trackPost.hamdist = nanhamXYZ;
        
            % record convergence 
            trackPost.merged = 0;           % will be computed below
            trackPost.converged = double(trackPost.hamdist == 0);
            
            % out
            tracksX1_out(itr, :) = trackPost.x';
            tracksY1_out(itr, :) = trackPost.y';
            tracksZ1_out(itr, :) = trackPost.z';

            iterStatus0(itr).condMeanCor = [iterStatus0(itr).condMeanCor, trackPost.condMeanCor];
            iterStatus0(itr).totFramesSatisfied = [iterStatus0(itr).totFramesSatisfied, trackPost.totFramesSatisfied];
            iterStatus0(itr).L2distPerfr = [iterStatus0(itr).L2distPerfr, trackPost.L2distPerfr];
            iterStatus0(itr).hamdist = [iterStatus0(itr).hamdist, trackPost.hamdist];
            
        end
    end
end

%% iteration summary
 
L = size(tracksX1_out, 1);
condMeanCorVec1 = nan(L,1);
totFramesSatisfiedVec1 = nan(L,1);
hamdist = nan(L,1);

for k = 1:L
    condMeanCorVec1(k) = iterStatus0(k).condMeanCor(end);
    totFramesSatisfiedVec1(k) = iterStatus0(k).totFramesSatisfied(end);
    hamdist(k) = iterStatus0(k).hamdist(end);
end

%% save

if ~isfolder(outFolder); mkdir(outFolder); end
save(fullfile(outFolder, 'elongateFiringEvents_afterPatchMatching.mat'), ...
            'tracksX1_out', 'tracksY1_out', 'tracksZ1_out', 'iterStatus0')

%% select tracks that are unique and have at least the min number of bright point sources

[~, ~, ~, ~, indVec] = ...
    selectTracks_Unique_WithMinNumBrightPS(tracksX1_out, tracksY1_out, tracksZ1_out, ...
                                            PSourMat2, imgArray, minFiringFramesOfNeurons);

%% summarize condMeanCorVec, totFramesSatisfiedVec

condMeanCorVec = condMeanCorVec1(indVec);
totFramesSatisfiedVec = totFramesSatisfiedVec1(indVec); 

disp('== Summary of mean of correlations above the threshold of selected tracks')
disp(summaryStatistics(condMeanCorVec))

disp('== Summary of total frames with correlation above the threshold of selected tracks')
disp(summaryStatistics(totFramesSatisfiedVec))

%% tracksIterX1

tracksIterX1 = tracksX1_out(indVec, :);
tracksIterY1 = tracksY1_out(indVec, :);
tracksIterZ1 = tracksZ1_out(indVec, :);

tracksIterXYZ1 = [tracksIterX1, tracksIterY1, tracksIterZ1];

%%

save(fullfile(outFolder, 'tracks_before_Iteration.mat'), ...
    'tracksIterXYZ1', 'corrThreshold', 'upbdDeformX', 'upbdDeformZ', ...
    'bandWidthX', 'bandWidthZ')
 
end
