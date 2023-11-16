function TrackJitteringFlickering(movieData, imgArray, paramStruct)
% TrackJitteringFlickering identifies groups of the detected firing events 
% that belong to the identical neurons, through local patch-matching
% iterations.
%
% paramStruct contains the following parameters.
%
% patchSizeX, patchSizeZ (in pixel)
%   - Edge size or length in X/Y- and Z-direction, respectively, 
%   of a local volume around each detected point source, which is utilized
%   to track the jittering and flickering neurons via correlation-based
%   patch-matching.
%
% bandWidthX, bandWidthZ (in pixel)
%   - The radius of the local volume or 3D-patch around each detected point source, 
%   that is, 2*bandWidthX+1 = patchSizeX.
%
% upbdDeformX, upbdDeformZ (in pixel)
%   - Upper bound of local deformation in X/Y- and Z-direction,
%   respectively. 
%   - Over the whole time frames, we assume that neuron locations can 
%   deviate from hypothetical central locations by upbdDeformX and 
%   upbdDeformZ (pixels) in X/Y- and Z-direction, respectively. 
%
% corrThreshold (0 < r <1)
%   - Two or multi-dim'l vector specifying multiple thresholds of spatial
%   correlations utilized for patch-matching.
%   - As an ensemble approach, tracking neurons are implemented using each
%   correlation threshold, and then the results are integrated to account
%   for heterogeneity of local image dynamics of individual neurons. 
%
% distBtwTracks_toBeMerged (in pixel)
%   - Distance between tracks to be merged. 
%   - If two or multiple neuron trajectories are too close, then their 
%   segmented ROIs would overlap. Even though the trajectories are correct,
%   overlapping ROIs is what we want to avoid. Thus, the algorithm computes
%   the minimun of euclidean (L2) distance between neuron locations over the
%   frames. If the minimum distance <= distBtwTracks_toBeMerged for two
%   neurons, then we select only one trajectory showing a better
%   measure of tracking. 
%   
% minFramesOfFiringEvents (in # of frames)
%   - Minimal number of frames for bright point sources linked in
%   consecutive time frames to be counted as a single firing event.
%
% minFiringFramesOfNeurons (in # of frames)
%   - For a neuron trajectory over the whole time frames to be valid, it 
%   should contain bright point sources at the specified minimun number of 
%   frames.
%
% imgMarginSizeX, imgMarginSizeZ (in pixel)
%   - Specify a marginal area in X/Y- and Z-direction, which is
%   pre-excluded from the analysis to reduce computational time. 
%   - Detected point sources in the specified marginal area will not be
%   tracked.
% 
% makeMov_trackedNeurons (true or false)
%   - Flag whether to generate MIP videos displaying obtained neuron
%   trajectories on MIP images.
%
% closeFigs (true or false)
%   - Whether to close output plots while running because the pipeline
%   generates many figures. Output figures are always saved in the output
%   directory in .fig and .png formats. It simply executes 'close all;'
%   after saving the figures. 
%
% Output are saved under MD.outputDirectory_/TrackJitteringFlickering.
%
% Jungsik Noh, UTSW, 2023/08
 

%% Input/Set up outputDir

p.outputDir = fullfile(movieData.outputDirectory_, 'TrackJitteringFlickering'); 

MD = movieData;
par2 = paramStruct;
par2.figFlag = 'on';

Part1OutputDir = fullfile(MD.outputDirectory_, 'DetectBrightPointSources');

if ~isfolder(p.outputDir); mkdir(p.outputDir); end
plotDir = fullfile(p.outputDir, 'detailedPlots');
if ~isfolder(plotDir); mkdir(plotDir); end

%% load output from DetectBrightPointSources()

load(fullfile(Part1OutputDir, 'PSourMat2.mat'))
load(fullfile(Part1OutputDir, 'tracks_Bright3.mat'))
load(fullfile(Part1OutputDir, 'tracksX1_Br3.mat'))
load(fullfile(Part1OutputDir, 'tracksY1_Br3.mat'))
load(fullfile(Part1OutputDir, 'tracksZ1_Br3.mat'))

%% filter out firing events longer than or equal to minFramesOfFiringEvents (=3)
% to reduce computation time

lt0 = nan(numel(tracks_Bright3, 1));
for i=1:numel(tracks_Bright3)
    lt0(i) = numel(~isnan(tracks_Bright3(i).x));
end

fig1 = figure('Visible', par2.figFlag);
histogram(lt0, 'BinWidth', 1)
xlabel('Frame'); ylabel('Counts')
title('Length of unfiltered bright tracks')

% save plots
saveas(fig1, fullfile(plotDir, 'lengthOfUnfilteredBrightTracks.fig'), 'fig')
saveas(fig1, fullfile(plotDir, 'lengthOfUnfilteredBrightTracks.png'), 'png')

% par2.minFramesOfFiringEvents
indL = (lt0 >= par2.minFramesOfFiringEvents);

tracks_Bright4 = tracks_Bright3(indL);
tracksX1_Br4 = tracksX1_Br3(indL, :);
tracksY1_Br4 = tracksY1_Br3(indL, :);
tracksZ1_Br4 = tracksZ1_Br3(indL, :);

disp('== Percent of single firing events of which lengths >= 3 frames')
disp(sum(indL)/numel(indL))

%% check uniqueness of the firing events (again)

tracksXYZ1_Br4 = [tracksX1_Br4, tracksY1_Br4, tracksZ1_Br4];
DD = pdist(tracksXYZ1_Br4, @nanhamdistXYZ_ignoreAllnan);

disp('== Number of pairs of non-unique firing events (it is required to be 0.)')
disp(sum(DD==0))

%% describe mean intensity of (unique) firing events

[tracks_Bright4.maxInt] = deal(0);
maxIntVec = nan(numel(tracks_Bright4), 1);
meanIntVec = nan(numel(tracks_Bright4), 1);
for k=1:numel(tracks_Bright4)
    x = round(tracks_Bright4(k).x);
    y = round(tracks_Bright4(k).y);
    z = round(tracks_Bright4(k).z);
    f = tracks_Bright4(k).f;
    tmp = nan(1, numel(x));
    for i=1:numel(x)
        if ~isnan(x(i))
            tmp(i) = imgArray(y(i),x(i),z(i),f(i));
        end
    end
    tracks_Bright4(k).Int = tmp;
    tracks_Bright4(k).maxInt = max(tmp(:), [], 'omitnan');
    maxIntVec(k) = max(tmp(:));
    tracks_Bright4(k).meanInt = nanmean(tmp(:));
    meanIntVec(k) = nanmean(tmp(:));
end

fig2 = figure('Visible', par2.figFlag);
histogram(meanIntVec)
xlabel('Mean Intensity')
title('Mean Intensity of Single Firing Events')

% save plots
saveas(fig2, fullfile(plotDir, 'meanIntOfSingleFiringEvents.fig'), 'fig')
saveas(fig2, fullfile(plotDir, 'meanIntOfSingleFiringEvents.png'), 'png')

if par2.closeFigs; close all; end

%% save
save(fullfile(p.outputDir, 'tracks_Bright4.mat'), 'tracks_Bright4')
save(fullfile(p.outputDir, 'tracksX1_Br4.mat'), 'tracksX1_Br4')
save(fullfile(p.outputDir, 'tracksY1_Br4.mat'), 'tracksY1_Br4')
save(fullfile(p.outputDir, 'tracksZ1_Br4.mat'), 'tracksZ1_Br4')

save(fullfile(p.outputDir, 'par2.mat'), 'par2')

%% Elongation and EM iterations for tracking
% Single firing events are first elongated over the whole time frames.
% For each given correlation threshold, the elongated initial tracks of
% neurons are fed into Expectation-Maximization iterations, where firing
% events of the same neuron are supposed to be linked together. 

numCorrs = numel(par2.corrThreshold);

% tracks X.merged.sorted.pooled.
tracksIterXmerged0sP = cell(1,numCorrs);
tracksIterYmerged0sP = cell(1,numCorrs);
tracksIterZmerged0sP = cell(1,numCorrs);
iterStatus4sP = cell(1,numCorrs);

for K = 1:numCorrs
    
    corrThreshold1 = par2.corrThreshold(K);
    outFolder1 = fullfile(p.outputDir, ['corrThreshold', num2str(corrThreshold1)]);
    
    % elongateFiringEvents
    disp('==')
    disp(['== corrThreshold: ', num2str(corrThreshold1)])
    disp('==')
    
    [tracksIterX1, tracksIterY1, tracksIterZ1] = ...
        elongateFiringEvents(tracksX1_Br4, tracksY1_Br4, tracksZ1_Br4, ...
                            PSourMat2, imgArray, par2, ...
                            corrThreshold1, outFolder1);
    
    [tracksIterXmerged0s, tracksIterYmerged0s, tracksIterZmerged0s, iterStatus4s] = ...
        EMiterations(tracksIterX1, tracksIterY1, tracksIterZ1, ...
                        imgArray, PSourMat2, MD, par2, ...
                        corrThreshold1, outFolder1);
    
    % EM tracks postprocessed, sorted, pooled
    tracksIterXmerged0sP{K} = tracksIterXmerged0s;
    tracksIterYmerged0sP{K} = tracksIterYmerged0s;
    tracksIterZmerged0sP{K} = tracksIterZmerged0s;
    iterStatus4sP{K} = iterStatus4s;
    
end

%% save

save(fullfile(p.outputDir, 'tracks_PooledOutput_EMiter_postProcessed.mat'), ...
    'tracksIterXmerged0sP', 'tracksIterYmerged0sP', 'tracksIterZmerged0sP', 'iterStatus4sP');

%% Combine multiple EM output with different corrThresholds

% merged, sorted, combined
tracksXmsc = tracksIterXmerged0sP{1};
tracksYmsc = tracksIterYmerged0sP{1};
tracksZmsc = tracksIterZmerged0sP{1};
iterStatus4sc = iterStatus4sP{1};

for K = 2:numCorrs
    tracksXmsc = [tracksXmsc; tracksIterXmerged0sP{K}];
    tracksYmsc = [tracksYmsc; tracksIterYmerged0sP{K}];
    tracksZmsc = [tracksZmsc; tracksIterZmerged0sP{K}];
    iterStatus4sc = [iterStatus4sc, iterStatus4sP{K}];
end

tracksXYZmsc = [tracksXmsc, tracksYmsc, tracksZmsc];
tracksXYZmsc_tmp = tracksXYZmsc;
tracksXYZmsc_tmp(isnan(tracksXYZmsc_tmp)) = 0;
[C, ia, ic] = unique(tracksXYZmsc_tmp, 'rows', 'stable');

%%  Make tracks unique

tracksX1 = tracksXmsc(ia, :);
tracksY1 = tracksYmsc(ia, :);
tracksZ1 = tracksZmsc(ia, :);
iterStatus4scu = iterStatus4sc(ia);

%% iter summary

iterConverged = [iterStatus4scu.iterConverged];
iterMerged = [iterStatus4scu.iterMerged];
ind = max(iterConverged,  iterMerged);

disp('trIter merged or converged')
tabulate(ind)

disp('== ID of oscillatory tracks')
nonconvInd = find(ind == 0);
disp(nonconvInd)

%% diffMeanCorVec

diffMeanCorVec = [iterStatus4scu.diffMeanCor];
fDiffMeanCorVec = figure('Visible', par2.figFlag); 
histogram(diffMeanCorVec)
title('condMeanCor High - condMeanCor Low')

%% save plots

saveas(fDiffMeanCorVec, fullfile(plotDir, 'fDiffMeanCorVec_R1Pooled.fig'), 'fig')
saveas(fDiffMeanCorVec, fullfile(plotDir, 'fDiffMeanCorVec_R1Pooled.png'), 'png')

%% select tracks that are unique and have at least the min number of bright point sources

[~, ~, ~, ~, ind0] = ...
    selectTracks_Unique_WithMinNumBrightPS(tracksX1, tracksY1, tracksZ1, ...
                                            PSourMat2, imgArray, par2.minFiringFramesOfNeurons);

%% tracks that are pooled, unique, and having minimal number of bright PSs 

iterStatus5 = iterStatus4scu(ind0);                        % converged, long
tracksIterXl = tracksX1(ind0, :);
tracksIterYl = tracksY1(ind0, :);
tracksIterZl = tracksZ1(ind0, :);

%% non-overlapping tracks
 
tracksIterXYZl = [tracksIterXl, tracksIterYl, tracksIterZl]; 
DD = pdist(tracksIterXYZl, @nanhamdistXYZ_ignoreAllnan);
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

%% Merge overlapping tracks to a track with the best diffMeanCor score

% merging: overlaggping tracks
[iterStatus5.mergedTo] = deal(0);
    for i = 1:numel(rr)
        if (iterStatus5(rr(i)).diffMeanCor >= iterStatus5(cc(i)).diffMeanCor)        
            iterStatus5(cc(i)).mergedTo = rr(i);
        else
            iterStatus5(rr(i)).mergedTo = cc(i);
        end
    end
    
% merged  (tracksIterXlong -> tracksIterXmerged)
final_TrCurr30 = find([iterStatus5.mergedTo] == 0);
tracksXcmerged0 = tracksIterXl(final_TrCurr30, :);
tracksYcmerged0 = tracksIterYl(final_TrCurr30, :);
tracksZcmerged0 = tracksIterZl(final_TrCurr30, :);
tracksXYZcmerged0 = [tracksXcmerged0, tracksYcmerged0, tracksZcmerged0];

DD = pdist(tracksXYZcmerged0, @nanhamdistXYZ_ignoreAllnan);

%% plot tracksXcmerged0

xccm = nanmean(tracksXcmerged0, 2);
yccm = nanmean(tracksYcmerged0, 2);
zccm = nanmean(tracksZcmerged0, 2);

fig3 = figure('Visible', par2.figFlag); 
scatter3(xccm, yccm, zccm, 15);
ax = gca;
ax.YDir = 'reverse'; ax.ZDir = 'reverse'; 
xlabel('x'); ylabel('y'); zlabel('z')
xlim([1, MD.imSize_(2)]);
ylim([1, MD.imSize_(1)]);
zlim([1, MD.zSize_]);

% trace which merged to which
hold on
mergedto = [iterStatus5.mergedTo]';
m1 = mergedto(mergedto~=0);
mergedfrom = nan(numel(mergedto), 1);
mergedfrom(m1) = 1;

indmergedfrom = (mergedto == 0) & (mergedfrom == 1);
xmfrom = nanmean(tracksIterXl(indmergedfrom, :), 2);
ymfrom = nanmean(tracksIterYl(indmergedfrom, :), 2);
zmfrom = nanmean(tracksIterZl(indmergedfrom, :), 2);

indnonmerged = (mergedto == 0) & isnan(mergedfrom);
xmnomerge = nanmean(tracksIterXl(indnonmerged, :), 2);
ymnomerge = nanmean(tracksIterYl(indnonmerged, :), 2);
zmnomerge = nanmean(tracksIterZl(indnonmerged, :), 2);

scatter3(xmfrom, ymfrom, zmfrom, [], 'x', 'b');
scatter3(xmnomerge, ymnomerge, zmnomerge, [], '+', 'r');

title(['Merged (b) #', num2str(sum(indmergedfrom)), ', No-overlap (r) #', ...
    num2str(sum(indnonmerged)), ', total #', num2str(numel(xccm))])

%% saveas

saveas(fig3, fullfile(plotDir, '3dscatter_R1Pooled.fig'), 'fig')

%% scatterXYmerged0_R1Pooled

fig2 = figure('Visible', par2.figFlag);
scatter(xmfrom, ymfrom, [], 'x', 'b');
hold on
scatter(xmnomerge, ymnomerge, [], '+', 'r');
xlabel('x'); ylabel('y');
axis ij
xlim([1, MD.imSize_(2)]); ylim([1, MD.imSize_(1)]); 
title(['Merged (b) #', num2str(sum(indmergedfrom)), ', No-overlap (r) #', ...
    num2str(sum(indnonmerged)), ', total #', num2str(numel(xccm))])

%% saveas
 
saveas(fig2, fullfile(plotDir, 'scatterXYmerged0_R1Pooled.fig'), 'fig')
saveas(fig2, fullfile(plotDir, 'scatterXYmerged0_R1Pooled.png'), 'png')

%% status after merging

ind = find([iterStatus5.mergedTo] == 0);
iterStatus6 = iterStatus5(ind);

%% interpolation for merged/ processed movie

tracksXcmerged0c = tracksXcmerged0;
tracksYcmerged0c = tracksYcmerged0;
tracksZcmerged0c = tracksZcmerged0;

for k = 1:size(tracksXcmerged0, 1)
    x = tracksXcmerged0(k, :);
    y = tracksYcmerged0(k, :);
    z = tracksZcmerged0(k, :);    
    tracksXcmerged0c(k, :) = imputeByEndPoints(x);
    tracksYcmerged0c(k, :) = imputeByEndPoints(y);
    tracksZcmerged0c(k, :) = imputeByEndPoints(z);
end   

%% order index
 
mXc = mean(tracksXcmerged0c, 2);
mYc = mean(tracksXcmerged0c, 2);
mZc = mean(tracksXcmerged0c, 2);
mx0 = min(mXc);
my0 = min(mYc);
mz0 = min(mZc);

dvec = sqrt( (mXc-mx0).^2 + (mYc - my0).^2 +(mZc-mz0).^2 );
[~, is] = sort(dvec);
mXc2 = mXc(is);
mYc2 = mYc(is);
mZc2 = mZc(is);
 
tracksXcmerged0cs = tracksXcmerged0c(is, :);
tracksYcmerged0cs = tracksYcmerged0c(is, :);
tracksZcmerged0cs = tracksZcmerged0c(is, :);
tracksXcmerged0s = tracksXcmerged0(is, :);
tracksYcmerged0s = tracksYcmerged0(is, :);
tracksZcmerged0s = tracksZcmerged0(is, :);
iterStatus6s = iterStatus6(is);

%% save after pooling and merging overlapping tracks

save(fullfile(p.outputDir, 'tracks_PooledEMoutput_MergedOverlaps.mat'),  ...
    'tracksXcmerged0cs', 'tracksYcmerged0cs', 'tracksZcmerged0cs', ...
    'tracksXcmerged0s', 'tracksYcmerged0s', 'tracksZcmerged0s', ...
    'iterStatus4scu', 'iterStatus5', 'iterStatus6', 'iterStatus6s')
     
%% Merge too close tracks, based on minL2 distance btn imputed tracks.

distBtwTracks_toBeMerged = par2.distBtwTracks_toBeMerged;

%% identify too close tracks

tracksXYZc = [tracksXcmerged0cs, tracksYcmerged0cs, tracksZcmerged0cs]; 

% Using minL2distXY
DD = pdist(tracksXYZc, @minL2distXYZ);
L = size(tracksXcmerged0cs, 1);
L*(L-1)/2;
i = find(DD <= distBtwTracks_toBeMerged);
rr = nan(size(i));
cc = nan(size(i));

tmp = L-1:-1:1;
tmp1= cumsum(tmp);

for j = 1:numel(i)
    rr(j) = find( i(j) <= tmp1, 1);
    if (rr(j) == 1)
        cc(j) = rr(j) + i(j);
    else
        cc(j) = rr(j) + i(j) - tmp1(rr(j)-1);
    end
end

%% Merge to a track with the best diffMeanCor score

for i = 1:numel(rr)
    if (iterStatus6s(rr(i)).diffMeanCor >= iterStatus6s(cc(i)).diffMeanCor)
        iterStatus6s(cc(i)).mergedTo = rr(i);
    else
        iterStatus6s(rr(i)).mergedTo = cc(i);
    end
end    

%% merged (tracksXcmerged0s -> tracksXcmerged0sL2)

final_TrCurr3 = find([iterStatus6s.mergedTo] == 0);
tracksXcmerged0sL2 = tracksXcmerged0s(final_TrCurr3, :);
tracksYcmerged0sL2 = tracksYcmerged0s(final_TrCurr3, :);
tracksZcmerged0sL2 = tracksZcmerged0s(final_TrCurr3, :);

%% plot tracksXcmerged0sL2

xccm = nanmean(tracksXcmerged0sL2, 2);
yccm = nanmean(tracksYcmerged0sL2, 2);
zccm = nanmean(tracksZcmerged0sL2, 2);

fig3 = figure('Visible', par2.figFlag); 
scatter3(xccm, yccm, zccm, 15);
ax = gca;
ax.YDir = 'reverse'; ax.ZDir = 'reverse'; 
xlabel('x'); ylabel('y'); zlabel('z')
xlim([1, MD.imSize_(2)]);
ylim([1, MD.imSize_(1)]);
zlim([1, MD.zSize_]);

% trace which merged to which
hold on
mergedto = [iterStatus6s.mergedTo]';
m1 = mergedto(mergedto~=0);
mergedfrom = nan(numel(mergedto), 1);
mergedfrom(m1) = 1;

indmergedfrom = (mergedto == 0) & (mergedfrom == 1);
xmfrom = nanmean(tracksXcmerged0s(indmergedfrom, :), 2);
ymfrom = nanmean(tracksYcmerged0s(indmergedfrom, :), 2);
zmfrom = nanmean(tracksZcmerged0s(indmergedfrom, :), 2);

indnonmerged = (mergedto == 0) & isnan(mergedfrom);
xmnomerge = nanmean(tracksXcmerged0s(indnonmerged, :), 2);
ymnomerge = nanmean(tracksYcmerged0s(indnonmerged, :), 2);
zmnomerge = nanmean(tracksZcmerged0s(indnonmerged, :), 2);

scatter3(xmfrom, ymfrom, zmfrom, [], 'x', 'b');
scatter3(xmnomerge, ymnomerge, zmnomerge, [], '+', 'r');

title(['Merged (b) #', num2str(sum(indmergedfrom)), ', No-overlap (r) #', ...
    num2str(sum(indnonmerged)), ', total #', num2str(numel(xccm))])

%%

saveas(fig3, fullfile(plotDir, '3dscatter_R1Pooled_L2Merged.fig'), 'fig')

%%

fig2 = figure('Visible', par2.figFlag);
scatter(xmfrom, ymfrom, [], 'x', 'b');
hold on
scatter(xmnomerge, ymnomerge, [], '+', 'r');
axis ij
xlim([1, MD.imSize_(2)]); ylim([1, MD.imSize_(1)]); 
xlabel('x'); ylabel('y');

title(['Merged (b) #', num2str(sum(indmergedfrom)), ', No-overlap (r) #', ...
    num2str(sum(indnonmerged)), ', total #', num2str(numel(xccm))])

%%

saveas(fig2, fullfile(plotDir, 'scatterXYmerged0_R1Pooled_L2Merged.fig'), 'fig')
saveas(fig2, fullfile(plotDir, 'scatterXYmerged0_R1Pooled_L2Merged.png'), 'png')

%% status after minL2 merging

tracksXcmerged0csL2 = tracksXcmerged0cs(final_TrCurr3, :);
tracksYcmerged0csL2 = tracksYcmerged0cs(final_TrCurr3, :);
tracksZcmerged0csL2 = tracksZcmerged0cs(final_TrCurr3, :);

iterStatus7 = iterStatus6s(final_TrCurr3);

save(fullfile(p.outputDir, 'tracks_PooledEMoutput_MergedOverlaps_L2Merged.mat'), ...
    'tracksXcmerged0sL2', 'tracksYcmerged0sL2', 'tracksZcmerged0sL2', ...
    'tracksXcmerged0csL2', 'tracksYcmerged0csL2', 'tracksZcmerged0csL2', ...
    'final_TrCurr3', 'iterStatus6s', 'iterStatus7')

if par2.closeFigs; close all; end

%% video output  

if par2.makeMov_trackedNeurons
    savePath = fullfile(p.outputDir, 'MIP_neuronTracks_whenFiring_tracksXcmerged0sL2');
    ca3D_plotTracksOnMIP(MD, imgArray, tracksXcmerged0sL2, tracksYcmerged0sL2, ...
                  tracksZcmerged0sL2, savePath, 'figFlag', par2.figFlag)
end
          
%%  

if par2.makeMov_trackedNeurons
    savePath = fullfile(p.outputDir, 'MIP_neuronTracks_imputedByEnds_tracksXcmerged0csL2');
    ca3D_plotTracksOnMIP(MD, imgArray, tracksXcmerged0csL2, tracksYcmerged0csL2, ...
                  tracksZcmerged0csL2, savePath, 'figFlag', par2.figFlag)
end

%%

disp('== done! ==')
disp('== End of TrackJitteringFlickering (Part 2/4) ==')

end
