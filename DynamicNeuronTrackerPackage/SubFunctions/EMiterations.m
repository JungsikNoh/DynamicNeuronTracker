function [tracksIterXmerged0s, tracksIterYmerged0s, tracksIterZmerged0s, iterStatus4s] = ...
    EMiterations(tracksIterX1, tracksIterY1, tracksIterZ1, ...
                imgArray, PSourMat2, MD, par2, corrThreshold1, outFolder1)
% EMiterations implements local patch-matching iterations to compute
% dynamic ROIs. The EM-like iterations are done by tracksXYZ_EMiteration().
% The rest parts are to post-process the converged trajactories.
%
% J Noh.


%% parse

upbdDeformX = par2.upbdDeformX ;
upbdDeformZ = par2.upbdDeformZ ;
bandWidthX = par2.bandWidthX ;
bandWidthZ = par2.bandWidthZ ;
minFiringFramesOfNeurons = par2.minFiringFramesOfNeurons;
imgMarginSizeX = par2.imgMarginSizeX;
imgMarginSizeZ = par2.imgMarginSizeZ;

%% tracksXYZ_EMiteration 

totTr = size(tracksIterX1, 1);
disp(['== total num of Tracks: ', num2str(totTr)])
numPS = sum(PSourMat2(:));
disp(['== numPS: ', num2str(numPS)])
disp(['== corrThreshold: ', num2str(corrThreshold1)])

% tracksXYZ_EMiteration
[tracksIterXYZCell, iterStatus] = tracksXYZ_EMiteration(...
    tracksIterX1, tracksIterY1, tracksIterZ1, ...
    imgArray, PSourMat2, upbdDeformX, upbdDeformZ, bandWidthX, bandWidthZ, ...
    corrThreshold1, imgMarginSizeX, imgMarginSizeZ);

%% save iteration result

tracksIterXYZCell_last3 = tracksIterXYZCell(end-2:end, :);

if ~isfolder(outFolder1); mkdir(outFolder1); end
save(fullfile(outFolder1, 'tracksXYZ_Iteration_result.mat'), 'tracksIterXYZCell_last3', ...
    'corrThreshold1', 'upbdDeformX', 'upbdDeformZ', 'bandWidthX', ...
    'bandWidthZ', 'iterStatus',  ...
    'numPS')

%% iter summary

iterConverged = [iterStatus.iterConverged];
iterMerged = [iterStatus.iterMerged];

ind = max(iterConverged,  iterMerged);

disp('trIter merged or converged')
tabulate(ind)

%
disp('== ID of oscillatory tracks')
nonconvInd = find(ind == 0);
disp(nonconvInd)

disp('== num of converged tracks after EM')
disp(tracksIterXYZCell{end, 5} )

%% output after iteration, final_TrCurr2

numTracks = [tracksIterXYZCell{:, 5}];
ftmp = figure('Visible', par2.figFlag);
plot(numTracks); ylim([0, max(numTracks)]);
xlabel('Iteration'); ylabel('Num of Tracks')

tracksIterX = tracksIterXYZCell{end, 1};
tracksIterY = tracksIterXYZCell{end, 2};
tracksIterZ = tracksIterXYZCell{end, 3};

final_TrCurr = tracksIterXYZCell{end, 4};
final_TrCurr2 = setdiff(final_TrCurr, nonconvInd);

% tracksIterX
tracksIterX = tracksIterX(final_TrCurr2, :);
tracksIterY = tracksIterY(final_TrCurr2, :);
tracksIterZ = tracksIterZ(final_TrCurr2, :);

%% saveas

saveas(ftmp, fullfile(outFolder1, 'numTracks_IterationR1.png'), 'png')

%% iterStatus2

iterStatus2 = iterStatus(final_TrCurr2);    % converged

%% compute correlation along tracks

numTr2 = numel(iterStatus2);

tracksIterXi = tracksIterX;
tracksIterYi = tracksIterY;
tracksIterZi = tracksIterZ;
allFrs = 1:size(tracksIterX, 2);

for k = 1:numTr2    
    fprintf('%05d ', k);
    if (mod(k, 50) == 0); fprintf('\n'); end
    
    x = tracksIterX(k, :);
    y = tracksIterY(k, :);
    z = tracksIterZ(k, :);
    
    xc = imputeByEndPoints(x);
    yc = imputeByEndPoints(y);
    zc = imputeByEndPoints(z);
    tracksIterXi(k, :) = xc;
    tracksIterYi(k, :) = yc;
    tracksIterZi(k, :) = zc;    
    
    % gathering neighboring (<= 1 pxl) PSourMat + localPatchSampling3D()    
    
    ptInd = nan(numel(x), 1);
    xvec2 = x;
    yvec2 = y;
    zvec2 = z;
    
    for fr = 1:numel(x)
        if ~isnan(y(fr)) 
            psArea = PSourMat2(y(fr)-1:y(fr)+1, x(fr)-1:x(fr)+1,z(fr)-1:z(fr)+1, fr);
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
    % clear firing images
    [avgPatch, ~] = localPatchSampling3D(xs,ys,zs,frs, bandWidthX, bandWidthZ, imgArray);
    
    % corr along the updated track
    [~, Patches] = localPatchSampling3D(xc,yc,zc,allFrs, bandWidthX, bandWidthZ, imgArray);
    patchFilter = (avgPatch - mean(avgPatch(:))) ./ std(avgPatch(:));
    corAlongTrack = nan(numel(x), 1);
    for fr = 1:numel(x)
        corAlongTrack(fr) = corrfilter3D(patchFilter, Patches(:,:,:,fr));
    end
    
    % output
    iterStatus2(k).corAlongTrack = corAlongTrack;
    frIndSatisfied = find(~isnan(x));
    condMeanCorHigh = mean(corAlongTrack(frIndSatisfied));  % double check with condMeanCor
    condMeanCorLow = mean(corAlongTrack(isnan(x))); 
    diffMeanCor = condMeanCorHigh - condMeanCorLow;
    
    iterStatus2(k).condMeanCorHigh = condMeanCorHigh;
    iterStatus2(k).condMeanCorLow = condMeanCorLow;
    iterStatus2(k).diffMeanCor = diffMeanCor;
end
 
%% diffMeanCorVec

diffMeanCorVec = [iterStatus2.diffMeanCor];
fDiffMeanCorVec = figure('Visible', par2.figFlag); 
histogram(diffMeanCorVec)
title('condMeanCor High - condMeanCor Low')

%%

saveas(fDiffMeanCorVec, fullfile(outFolder1, 'fDiffMeanCorVec_R1.fig'), 'fig')
saveas(fDiffMeanCorVec, fullfile(outFolder1, 'fDiffMeanCorVec_R1.png'), 'png')

%% totFrSatisfiedVec

totFrSatisfiedVec = nan(numel(iterStatus2), 1);
for i = 1:numel(iterStatus2)
    tmp = iterStatus2(i).totFramesSatisfied(end);
    totFrSatisfiedVec(i) = tmp;
end
fhistTotFrSat = figure('Visible', par2.figFlag); 
histogram(totFrSatisfiedVec, 'BinWidth', 1)
title('Num of total frames satisfying correlation criterion per track')

%%  saveas

saveas(fhistTotFrSat, fullfile(outFolder1, 'fhistTotFrSat_R1.fig'), 'fig')
saveas(fhistTotFrSat, fullfile(outFolder1, 'fhistTotFrSat_R1.png'), 'png')

%% select tracks that are unique and have at least the min number of bright point sources

[~, ~, ~, ~, ind0] = ...
    selectTracks_Unique_WithMinNumBrightPS(tracksIterX, tracksIterY, tracksIterZ, ...
                                            PSourMat2, imgArray, minFiringFramesOfNeurons);
          
%% converged, unique, long enough tracks

iterStatus3 = iterStatus2(ind0);                        

%% scatter again with long ones

tracksIterXl = tracksIterX(ind0, :);
tracksIterYl = tracksIterY(ind0, :);
tracksIterZl = tracksIterZ(ind0, :);

size(tracksIterXl)

xcc = nanmean(tracksIterXl, 2);
ycc = nanmean(tracksIterYl, 2);
zcc = nanmean(tracksIterZl, 2);

%% non-overlapping tracks

tracksIterXYZl = [tracksIterXl, tracksIterYl, tracksIterZl];
%  
DD = pdist(tracksIterXYZl, @nanhamdistXYZ_ignoreAllnan);
summaryStatistics(DD)
sum(DD==0)
%
i = find(DD<1);
rr = nan(size(i));
cc = nan(size(i));
Lt = size(tracksIterXl, 1);

tmp = Lt-1:-1:1;
tmp1= cumsum(tmp);

for j = 1:numel(i)
    rr(j) = find( i(j) <= tmp1, 1);
    if (rr(j) == 1)
        cc(j) = rr(j) + i(j);
    else
        cc(j) = rr(j) + i(j) - tmp1(rr(j)-1);
    end
end

trWithOverlappingId = unique([rr, cc]);

trNonOverlapId = setdiff(1:Lt, trWithOverlappingId);

% scatter plots removed

%% Merge overlapping tracks to a track with the best diffMeanCor score

% numel(trWithOverlappingId)

% merging: overlaggping tracks
[iterStatus3.mergedTo] = deal(0);
    for i = 1:numel(rr)
        if (iterStatus3(rr(i)).diffMeanCor >= iterStatus3(cc(i)).diffMeanCor)
            iterStatus3(cc(i)).mergedTo = rr(i);
        else
            iterStatus3(rr(i)).mergedTo = cc(i);
        end
    end
    
% merged  (tracksIterXlong -> tracksIterXmerged)

final_TrCurr30 = find([iterStatus3.mergedTo] == 0);

tracksIterXmerged0 = tracksIterXl(final_TrCurr30, :);
tracksIterYmerged0 = tracksIterYl(final_TrCurr30, :);
tracksIterZmerged0 = tracksIterZl(final_TrCurr30, :);

tracksIterXYZlmerged0 = [tracksIterXmerged0, tracksIterYmerged0, tracksIterZmerged0];

DD = pdist(tracksIterXYZlmerged0, @nanhamdistXYZ_ignoreAllnan);
summaryStatistics(DD)
sum(DD==0)

% Output of Iteration: tracksIterXmerged0

%% plot tracksIterXmerged0

xccm = nanmean(tracksIterXmerged0, 2);
yccm = nanmean(tracksIterYmerged0, 2);
zccm = nanmean(tracksIterZmerged0, 2);

fig3 = figure('Visible', par2.figFlag); 

fs = scatter3(xccm, yccm, zccm, 15);
ax = gca;
ax.YDir = 'reverse'; ax.ZDir = 'reverse'; 
xlabel('x'); ylabel('y'); zlabel('z')
xlim([1, MD.imSize_(2)]);
ylim([1, MD.imSize_(1)]);
zlim([1, MD.zSize_]);

% trace which merged to which
hold on

mergedto = [iterStatus3.mergedTo]';
m1 = mergedto(mergedto~=0);
mergedfrom = nan(numel(mergedto), 1);
mergedfrom(m1) = 1;

mergetrace = [mergedto, mergedfrom];

indmergedfrom = (mergedto == 0) & (mergedfrom == 1);
xmfrom = nanmean(tracksIterXl(indmergedfrom, :), 2);
ymfrom = nanmean(tracksIterYl(indmergedfrom, :), 2);
zmfrom = nanmean(tracksIterZl(indmergedfrom, :), 2);

indnonmerged = (mergedto == 0) & isnan(mergedfrom);
xmnomerge = nanmean(tracksIterXl(indnonmerged, :), 2);
ymnomerge = nanmean(tracksIterYl(indnonmerged, :), 2);
zmnomerge = nanmean(tracksIterZl(indnonmerged, :), 2);

%
fs = scatter3(xmfrom, ymfrom, zmfrom, [], 'x', 'b');
fs = scatter3(xmnomerge, ymnomerge, zmnomerge, [], '+', 'r');

title(['Merged (b) #', num2str(sum(indmergedfrom)), ', No-overlap (r) #', ...
    num2str(sum(indnonmerged)), ', total #', num2str(numel(xccm))])

saveas(fs, fullfile(outFolder1, '3dscatter_R1.fig'), 'fig')

%% scatterXYmerged0_R1

fig2 = figure('Visible', par2.figFlag);
fs = scatter(xmfrom, ymfrom, [], 'x', 'b');
hold on
fs = scatter(xmnomerge, ymnomerge, [], '+', 'r');
xlabel('x'); ylabel('y');
axis ij
xlim([1, MD.imSize_(2)]); ylim([1, MD.imSize_(1)]); 

title(['Merged (b) #', num2str(sum(indmergedfrom)), ', No-overlap (r) #', ...
    num2str(sum(indnonmerged)), ', total #', num2str(numel(xccm))])

%
saveas(fig2, fullfile(outFolder1, 'scatterXYmerged0_R1.fig'), 'fig')
saveas(fig2, fullfile(outFolder1, 'scatterXYmerged0_R1.png'), 'png')

%% status after merging

ind = find([iterStatus3.mergedTo] == 0);
iterStatus4 = iterStatus3(ind);

%% interpolation for merged/ processed movie

size(tracksIterXmerged0)

tracksIterXmerged0c = tracksIterXmerged0;
tracksIterYmerged0c = tracksIterYmerged0;
tracksIterZmerged0c = tracksIterZmerged0;

for k = 1:size(tracksIterXmerged0, 1)
    x = tracksIterXmerged0(k, :);
    y = tracksIterYmerged0(k, :);
    z = tracksIterZmerged0(k, :);
    
    tracksIterXmerged0c(k, :) = imputeByEndPoints(x);
    tracksIterYmerged0c(k, :) = imputeByEndPoints(y);
    tracksIterZmerged0c(k, :) = imputeByEndPoints(z);
end   

%% order index
 
mXc = mean(tracksIterXmerged0c, 2);
mYc = mean(tracksIterYmerged0c, 2);
mZc = mean(tracksIterZmerged0c, 2);

mx0 = min(mXc);
my0 = min(mYc);
mz0 = min(mZc);

dvec = sqrt( (mXc-mx0).^2 + (mYc - my0).^2 +(mZc-mz0).^2 );

[sdvec, is] = sort(dvec);

mXc2 = mXc(is);
mYc2 = mYc(is);
mZc2 = mZc(is);

% imputedByEnds, sorted
tracksIterXmerged0cs = tracksIterXmerged0c(is, :);
tracksIterYmerged0cs = tracksIterYmerged0c(is, :);
tracksIterZmerged0cs = tracksIterZmerged0c(is, :);

tracksIterXmerged0s = tracksIterXmerged0(is, :);
tracksIterYmerged0s = tracksIterYmerged0(is, :);
tracksIterZmerged0s = tracksIterZmerged0(is, :);

iterStatus4s = iterStatus4(is);

%% save output

save(fullfile(outFolder1, 'tracks_Output_EMiter_postProcessed.mat'), 'tracksIterXYZlmerged0', ...
    'tracksIterXmerged0cs', 'tracksIterYmerged0cs', 'tracksIterZmerged0cs', ...
    'tracksIterXmerged0s', 'tracksIterYmerged0s', 'tracksIterZmerged0s', ...
    'iterStatus', 'iterStatus3', 'iterStatus4', 'iterStatus4s')

%%

if par2.makeMov_trackedNeurons
    savePath = fullfile(outFolder1, 'MIP_EM_tracksIterXYZmerged0s');
    ca3D_plotTracksOnMIP(MD, imgArray, tracksIterXmerged0s, tracksIterYmerged0s, ...
                  tracksIterZmerged0s, savePath, 'figFlag', par2.figFlag)
end

%%   

if par2.makeMov_trackedNeurons
    savePath = fullfile(outFolder1, 'MIP_EM_tracksIterXYZmerged0imputed_');
    ca3D_plotTracksOnMIP(MD, imgArray, tracksIterXmerged0cs, tracksIterYmerged0cs, ...
                  tracksIterZmerged0cs, savePath, 'figFlag', par2.figFlag)
end

end
