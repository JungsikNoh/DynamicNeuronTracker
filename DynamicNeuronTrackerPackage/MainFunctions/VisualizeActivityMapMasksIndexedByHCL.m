function VisualizeActivityMapMasksIndexedByHCL(movieData, imgArray, paramStruct)
% VisualizeActivityMapMasksIndexedByHCL visualizes segmented dynamic ROIs.
% Weighted (if chosen) normalized calcium activity time courses are
% visualized in a heatmap after hierarchical clustering (HCL). ROIs are
% visualized on MIP images when they are active.
%
% paramStruct contains the following parameters.
%
% upSamplingFactor_forMaskContourImg (>= 1)
%   - Specify an image-resizing factor to visualize the masks of tracked
%   neurons, particularly when neuron masks are too small for contour
%   visualization. 
%   - The mask contours when the neurons are firing are visualized on top 
%   of maximum intensity projection images.
%
% makeMov_ROIsAtHighActivities (true or false)
%   - Flag whether to generate a MIP video displaying ROIs. ROIs are shown
%   only when activities are determined to be high by K-means clustering.
%
% makeMov_ROIsFiringAnnotated (true or false)
%   - Flag whether to generate a MIP video displaying ROIs. ROIs are shown
%   only when ROIs are determined to be firing via spatial correlation.
%
% makeMov_singleROIs_AllFrames_firingAnnotated (true or false)
%   - Flag whether to generate all of the single ROI local MIP videos. They
%   display all time frames and ROIs are shown when they are determined to 
%   be firing. 
%
% makeMov_singleROIs_Snapshots_atHighActivities (true or false)
%   - Flag whether to generate all of the single ROI local MIP videos. They
%   display a subset of time frames when the target ROI's activities are
%   determined to be high by K-means clustering.
%
% closeFigs (true or false)
%   - Whether to close output plots while running because the pipeline
%   generates many figures. Output figures are always saved in the output
%   directory in .fig and .png formats. It simply executes 'close all;'
%   after saving the figures. 
%
% Output are saved under MD.outputDirectory_/VisualizeActivityMapMasksIndexedByHCL_weighted
% if par4.activityWeightedByImgCorr = true. If not, the output directory is
% MD.outputDirectory_/VisualizeActivityMapMasksIndexedByHCL_unweighted.
%
% Jungsik Noh, UTSW, 2023/08


%% Input

par4 = paramStruct; 
par4.figFlag = 'on';

MD = movieData;
if par4.activityWeightedByImgCorr
    p.outputDir = fullfile(MD.outputDirectory_, 'VisualizeActivityMapMasksIndexedByHCL_weighted'); 
else
    p.outputDir = fullfile(MD.outputDirectory_, 'VisualizeActivityMapMasksIndexedByHCL_unweighted'); 
end

Part2OutputDir = fullfile(MD.outputDirectory_, 'TrackJitteringFlickering');
Part3OutputDir = fullfile(MD.outputDirectory_, 'SegmentDynamicROIs');
if ~isfolder(p.outputDir); mkdir(p.outputDir); end

%% Load output from Part 2 & 3

% load ROIcell, tracksXcmerged0csL2, psSig, normSig
load(fullfile(Part3OutputDir, 'extractedSignals_indexedByLocations.mat'))
load(fullfile(Part3OutputDir, 'roicenterXYZ.mat'))
load(fullfile(Part3OutputDir, 'ptMaskVol_ROI.mat'))

% load iterStatus_final (iterStatus9)
load(fullfile(Part2OutputDir, 'tracks_PooledEMoutput_MergedOverlaps_L2Merged.mat'))
                                            % include iterStatus7.mat or
                                            % iterStatus9.mat (when subvols
                                            % combined)
load(fullfile(Part2OutputDir, 'par2.mat'))  % load upbdDeformX, etc, for single ROI movies

% change the object name from iterStatus7 or iterStatus9 to
% iterStatus_final
tmp = who('iterStatus*');
if any(strcmp(tmp, 'iterStatus7'))
    iterStatus_final = iterStatus7;
elseif any(strcmp(tmp, 'iterStatus9'))
    iterStatus_final = iterStatus9;
end    

%% determine neural activity matrix

trmax = size(normSig, 1);
frmax = size(normSig, 2);
spatCorrSig = nan(trmax, frmax);
for i = 1:trmax
    tmp = iterStatus_final(i).corAlongTrack;
    spatCorrSig(i, :) =  tmp;
end   

if ~par4.activityWeightedByImgCorr
    actSig = normSig;
    actSig_log = normSig_log;
else
    weightMat = max(spatCorrSig, 0.01);
    actSig = normSig .* weightMat;
    actSig_log = normSig_log .* weightMat;
end

%% when weights are used, the img corr TS is described

if par4.activityWeightedByImgCorr
    
    % Avg activity
    m2 = mean(spatCorrSig, 1);

    f01 = figure('Visible', par4.figFlag);
    plot(m2)
    refline([0,0]) 
    title(['Avgeraged spatial correlation per frame'])
    xlabel('Frame')
    ylabel('Activity')

    % Hierarchical CLustering (HCL) to compute similarities btwn standardized normalized signals
    trmax = size(spatCorrSig, 1);
    inputMat = zscore(spatCorrSig')';

    D = pdist(inputMat);
    tree = linkage(inputMat, 'centroid');
    leafOrder = optimalleaforder(tree, D);

    f02 = figure('Visible', par4.figFlag);
    [H,T,outpermCorr] = dendrogram(tree, trmax, 'Reorder', leafOrder, 'Orientation', 'left');
    xlabel('Distance')
    ylabel('Neuron ID')

    % Activity map indexed by HCL
    spatCorrSigHCL = spatCorrSig(outpermCorr, :);
    q = quantile(spatCorrSig(:), [0.005, 0.995]);

    f03 = figure('Visible', par4.figFlag); 
    imagesc(spatCorrSigHCL, [q(1), q(2)])
    colormap(jet), colorbar 
    title('Spatial correlation TS used as weights')
    xlabel('Frame')
    ylabel('Neuron Index ordered by correlations') 

    % saveas
    saveas(f01, fullfile(p.outputDir, 'spatCorr_avgActivityTS.fig'), 'fig')
    saveas(f02, fullfile(p.outputDir, 'spatCorr_dendrogramHCL.fig'), 'fig')
    saveas(f03, fullfile(p.outputDir, 'spatCorr_HCL-indexedActivityMap.fig'), 'fig') 

    saveas(f01, fullfile(p.outputDir, 'spatCorr_avgActivityTS.png'), 'png')
    saveas(f02, fullfile(p.outputDir, 'spatCorr_dendrogramHCL.png'), 'png')
    saveas(f03, fullfile(p.outputDir, 'spatCorr_HCL-indexedActivityMap.png'), 'png') 

end

if par4.closeFigs; close all; end

%% Avg activity

m2 = mean(actSig, 1);

f1 = figure('Visible', par4.figFlag);
plot(m2)
refline([0,0]) 
title(['Avgeraged activity per frame'])
xlabel('Frame')
ylabel('Activity')

%% Hierarchical CLustering (HCL) to compute similarities btwn standardized normalized signals
 
trmax = size(actSig, 1);
inputMat = zscore(actSig')';

D = pdist(inputMat);
tree = linkage(inputMat, 'centroid');
leafOrder = optimalleaforder(tree, D);

f2 = figure('Visible', par4.figFlag);
[H,T,outperm] = dendrogram(tree, trmax, 'Reorder', leafOrder, 'Orientation', 'left');
xlabel('Distance')
ylabel('Neuron ID')

%% Activity map indexed by HCL

actSigHCL = actSig(outperm, :);
q = quantile(actSig(:), [0.005, 0.995]);

f3 = figure('Visible', par4.figFlag); 
imagesc(actSigHCL, [q(1), q(2)])
cmap = redwhiteblue(q(1), q(2), 256);
colormap(cmap), colorbar 
if par4.activityWeightedByImgCorr
    title('Neural Activity (\DeltaF/F * weight)')
else
    title('Neural Activity (\DeltaF/F)')
end
xlabel('Frame')
ylabel('Neuron Index') 

%% Activity (log) map indexed by HCL

actSigHCL_log = actSig_log(outperm, :);
q = quantile(actSig_log(:), [0.005, 0.995]);

f3log = figure('Visible', par4.figFlag); 
imagesc(actSigHCL_log, [q(1), q(2)])
cmap = redwhiteblue(q(1), q(2), 256);
colormap(cmap), colorbar 
if par4.activityWeightedByImgCorr
    title('Neural Activity (\DeltaF * weight)')
else
    title('Neural Activity (\DeltaF)')
end
xlabel('Frame')
ylabel('Neuron Index') 

%% when using weights, unweighted map is also plotted for comparison at the same scale

if par4.activityWeightedByImgCorr
    
    normSigHCL = normSig(outperm, :);
    q = quantile(actSig(:), [0.005, 0.995]);
    
    f32 = figure('Visible', par4.figFlag); 
    imagesc(normSigHCL, [q(1), q(2)])
    cmap = redwhiteblue(q(1), q(2), 256);
    colormap(cmap), colorbar 
    title('Neural Activity (\DeltaF/F, unweight)') 
    xlabel('Frame')
    ylabel('Neuron Index') 
    
    pause(0.5)
    saveas(f32, fullfile(p.outputDir, 'HCL-indexedActivityMap-unweighted.fig'), 'fig')     
    saveas(f32, fullfile(p.outputDir, 'HCL-indexedActivityMap-unweighted.png'), 'png')     
end

%% Raw activity map indexed by HCL

psSigHCL = psSig(outperm, :);
q = quantile(psSig(:), [0.005, 0.995]);

f4 = figure('Visible', par4.figFlag); 
imagesc(psSigHCL, [q(1), q(2)])
colormap(jet), colorbar 
title('Neural Activity (raw)')
xlabel('Frame')
ylabel('Neuron Index') 

%% HCL-indexed ROIs on MIP

meanVol = mean(imgArray, 4);
roicenterXmergedcsPerm = roicenterXmergedcs(outperm, :);
roicenterYmergedcsPerm = roicenterYmergedcs(outperm, :);
roicenterZmergedcsPerm = roicenterZmergedcs(outperm, :);

mX = mean(roicenterXmergedcsPerm, 2);
mY = mean(roicenterYmergedcsPerm, 2);
mZ = mean(roicenterZmergedcsPerm, 2);

figROI = figure('Visible', par4.figFlag); 
qq = quantile(reshape(meanVol, 1,[]), [0.5, 0.9999]); 
ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;    
[~,~,~,cthree] = j_computeMIPs(meanVol, ZXRatio, qq(1), qq(2));

imagesc(cthree)
colormap gray % default
colorbar
axis off

trmax = size(roicenterXmergedcsPerm, 1);
hold on
maxXBorder=MD.getDimensions('X');
maxYBorder=MD.getDimensions('Y');
maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);

strSize0 = 4;
x = mX(:);
y = mY(:);
z = mZ(:);
z = z .* ZXRatio - ZXRatio/2;

hold on
scatter(x, y, 10, 'MarkerEdgeColor', 'r')
scatter(x, z+maxYBorder+strSize0, 10,  'MarkerEdgeColor', 'r')
scatter(z+maxXBorder+strSize0, y, 10,  'MarkerEdgeColor', 'r')

for k = 1:trmax    
    x = mX(k)+1;
    y = mY(k)+1;
    z = mZ(k)+1;    
    z = z .* ZXRatio - ZXRatio/2;

    text(x, y, num2str(k), 'FontSize', 8, 'Color', 'r');
    text(x, z+maxYBorder+strSize0, num2str(k), 'FontSize', 8, 'Color', 'r')
    text(z+maxXBorder+strSize0, y, num2str(k), 'FontSize', 8, 'Color', 'r')
end

%% saveas

saveas(f1, fullfile(p.outputDir, 'avgActivityTS.fig'), 'fig')
saveas(f2, fullfile(p.outputDir, 'dendrogramHCL.fig'), 'fig')
saveas(f3, fullfile(p.outputDir, 'HCL-indexedActivityMap.fig'), 'fig') 
saveas(f4, fullfile(p.outputDir, 'HCL-indexedRawActivityMap.fig'), 'fig') 
saveas(f3log, fullfile(p.outputDir, 'HCL-indexedActivityMap-log.fig'), 'fig') 

saveas(f1, fullfile(p.outputDir, 'avgActivityTS.png'), 'png')
saveas(f2, fullfile(p.outputDir, 'dendrogramHCL.png'), 'png')
saveas(f3, fullfile(p.outputDir, 'HCL-indexedActivityMap.png'), 'png') 
saveas(f4, fullfile(p.outputDir, 'HCL-indexedRawActivityMap.png'), 'png') 
saveas(f3log, fullfile(p.outputDir, 'HCL-indexedActivityMap-log.png'), 'png') 

saveas(figROI, fullfile(p.outputDir, ['HCL-indexedROIsOnMIP', '.fig']), 'fig')
saveas(figROI, fullfile(p.outputDir, ['HCL-indexedROIsOnMIP', '.png']), 'png')

if par4.closeFigs; close all; end

%% threshold activity by kmeans per-cell (roicenterXmergedcsPerm_kmeans)

indHigh = ones(size(actSigHCL));

for i = 1:size(actSigHCL, 1)
    f = actSigHCL(i, :);
    [idx, C] = kmeans(f(:), 2, 'Replicates', 3);
    % order C(1) =< C(2)
    if C(1) > C(2)
        idx = 3 - idx;
    end
    indHigh(i, idx == 1) = NaN;
end

roicenterXmergedcsPerm_kmeans = roicenterXmergedcsPerm;
roicenterYmergedcsPerm_kmeans = roicenterYmergedcsPerm;
roicenterZmergedcsPerm_kmeans = roicenterZmergedcsPerm;

roicenterXmergedcsPerm_kmeans(isnan(indHigh)) = NaN;
roicenterYmergedcsPerm_kmeans(isnan(indHigh)) = NaN;
roicenterZmergedcsPerm_kmeans(isnan(indHigh)) = NaN;

%% save

roicenterXmergedsPerm = roicenterXmergeds(outperm, :);
roicenterYmergedsPerm = roicenterYmergeds(outperm, :);
roicenterZmergedsPerm = roicenterZmergeds(outperm, :);

save(fullfile(p.outputDir, 'roicenterXYZ_HCLPerm.mat'), ...
    'roicenterXmergedcs', 'roicenterXmergedcsPerm', 'roicenterXmergeds', 'roicenterXmergedsPerm', ...
    'roicenterYmergedcs', 'roicenterYmergedcsPerm', 'roicenterYmergeds', 'roicenterYmergedsPerm', ...
    'roicenterZmergedcs', 'roicenterZmergedcsPerm', 'roicenterZmergeds', 'roicenterZmergedsPerm', ...
    'roicenterXmergedcsPerm_kmeans', 'roicenterYmergedcsPerm_kmeans', 'roicenterZmergedcsPerm_kmeans', ...
    'outperm')

save(fullfile(p.outputDir, 'par4.mat'), 'par4')

save(fullfile(p.outputDir, 'actSig_psSig_HCL-indexed.mat'), ...
    'actSigHCL', 'psSigHCL', 'outperm', 'actSigHCL_log') 

writematrix(actSigHCL, fullfile(p.outputDir, 'actSig_HCLindexed.csv'))
writematrix(actSigHCL_log, fullfile(p.outputDir, 'actSig_log_HCLindexed.csv'))

%% when weights are used, the img corr TS is described (2)

if par4.activityWeightedByImgCorr
    
    % Activity map indexed by HCL
    spatCorrSigHCL2 = spatCorrSig(outperm, :);
    q = quantile(spatCorrSig(:), [0.005, 0.995]);

    f05 = figure('Visible', par4.figFlag); 
    imagesc(spatCorrSigHCL2, [q(1), q(2)])
    colormap(jet), colorbar 
    title('Spatial correlation TS used as weights (2)')
    xlabel('Frame')
    ylabel('Neuron Index ordered by correlations') 

    % saveas
    saveas(f05, fullfile(p.outputDir, 'spatCorr_HCL-indexedActivityMap2.fig'), 'fig') 

    saveas(f05, fullfile(p.outputDir, 'spatCorr_HCL-indexedActivityMap2.png'), 'png') 
    
    save(fullfile(p.outputDir, 'spatCorrSigHCL2.mat'), 'spatCorrSig', 'spatCorrSigHCL2', 'outperm')

end

%% allTSplots

allTSout = fullfile(p.outputDir, 'allTSplots');
if ~isfolder(allTSout); mkdir(allTSout); end

%%

nCells = size(actSigHCL, 1);
frmax = size(actSigHCL, 2);
frMargin = round(0.01 * frmax);
frStart = 1 - frMargin;
frEnd = frmax + frMargin;

for j = 1:5:nCells

    f5 = figure('Position', [10 10 1000 700], 'Visible', par4.figFlag);

    imax = min(j+4, nCells);
    for i = j:imax
        f = actSigHCL(i, :);
        rpos = mod(i, 5);
        if rpos == 0; rpos = 5; end

        subplot(5, 1, rpos)
        plot(f, 'k')
        refline([0, 0])
        ylabel(['ROI ', num2str(i)])
        xlim([frStart, frEnd])
        
        ax = gca;
        if (ax.YLim(1) > -0.02)
            tmp = ax.YLim(2); ax.YLim = [-0.02, tmp];
        end

        fprintf(1, '%g ', i); 
        if (mod(i,50) == 0); fprintf('\n'); end    
    end

    xlabel('Frame')
    
    saveas(f5, fullfile(allTSout, ['allTS_', num2str(j), '.png']), 'png')
    pause(0.2)
    close(f5)
end

%% ordering ptMaskVol_ROI

ptMaskVol_ROI_HCLindex = ptMaskVol_ROI; 

for i = 1:trmax
    ind = find(ptMaskVol_ROI == i);
    permID = find(outperm == i);
    ptMaskVol_ROI_HCLindex(ind) = permID;
end

%%
save(fullfile(p.outputDir, 'ptMaskVol_ROI_HCLindex.mat'), 'ptMaskVol_ROI_HCLindex', '-v7.3')

%% Visualize Contours of masks indexed by HCL 
% masks are shown when firing

if par4.makeMov_ROIsFiringAnnotated    % using roicenterXmergedsPerm

    upSamplingFactor = par4.upSamplingFactor_forMaskContourImg;
    ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;

    savePath = fullfile(p.outputDir, 'maskContoursOnMIP_firingAnnotated');

    ca3D_maskContoursOnMIP_whenFiring(imgArray, ZXRatio, ptMaskVol_ROI_HCLindex, savePath, ...
        upSamplingFactor, roicenterXmergedsPerm, roicenterYmergedsPerm, roicenterZmergedsPerm, ...
        'figFlag', par4.figFlag)
end

%% Visualize Contours of masks indexed by HCL 
% masks at high activities

if par4.makeMov_ROIsAtHighActivities    % using roicenterXmergedcsPerm_kmeans

    upSamplingFactor = par4.upSamplingFactor_forMaskContourImg;
    ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;

    savePath = fullfile(p.outputDir, 'maskContoursOnMIP_atHighActivities');

    ca3D_maskContoursOnMIP_whenFiring(imgArray, ZXRatio, ptMaskVol_ROI_HCLindex, savePath, ...
        upSamplingFactor, roicenterXmergedcsPerm_kmeans, roicenterYmergedcsPerm_kmeans, ...
        roicenterZmergedcsPerm_kmeans, ...
        'figFlag', par4.figFlag)
end

%% Make all movies of ROIs with firing annotated

if par4.makeMov_singleROIs_AllFrames_firingAnnotated   % using roicenterXmergedsPerm

    allMoviesPerROIDir = fullfile(p.outputDir, 'allFramesOfSingleROIs_firingAnnotated');
    if ~isfolder(allMoviesPerROIDir); mkdir(allMoviesPerROIDir); end

    % load roicenterXmergedsPerm, etc
    load(fullfile(Part2OutputDir, 'par2.mat'))

    % Load imgArray
    if ~exist('imgArray')
        imgArray = uint16(zeros(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_));

        disp('======')
        disp('Loading frames:')

        parfor fr=1:MD.nFrames_
            currImage = uint16(MD.channels_(1).loadStack(fr));
            imgArray(:,:,:,fr) = currImage;
            fprintf(1, '%g ', fr); 
            if (mod(fr,50) == 0); fprintf('\n'); end    
        end
        fprintf(1, '\n')
        
        % handle 0 intensities: replace 0 intensity with minimum
        m0 = min(imgArray(imgArray > 0));
        imgArray(imgArray == 0) =  m0; 
    end

    %
    upSamplingFactor = 2;     %par4.upSamplingFactor_forMaskContourImg;
    ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;

    upbdDeformX = par2.upbdDeformX;
    upbdDeformZ = par2.upbdDeformZ;
    bandWidthX = par2.bandWidthX;
    bandWidthZ = par2.bandWidthZ;

    mX = round(mean(roicenterXmergedsPerm, 2, 'omitnan'));
    mY = round(mean(roicenterYmergedsPerm, 2, 'omitnan'));
    mZ = round(mean(roicenterZmergedsPerm, 2, 'omitnan'));
    xmax0 = MD.imSize_(2);
    ymax0 = MD.imSize_(1);
    zmax0 = MD.zSize_;

    subFrs = 1:MD.nFrames_;     % all frames
    for k = 1:numel(mX)
        disp(['== Movie for ROI ', num2str(k)])

        x0 = max(mX(k)-upbdDeformX-bandWidthX, 1);
        x1 = min(mX(k)+upbdDeformX+bandWidthX, xmax0);
        y0 = max(mY(k)-upbdDeformX-bandWidthX, 1);
        y1 = min(mY(k)+upbdDeformX+bandWidthX, ymax0);
        z0 = max(mZ(k)-upbdDeformZ-bandWidthZ, 1);
        z1 = min(mZ(k)+upbdDeformZ+bandWidthZ, zmax0);

        fieldCylinder = imgArray(y0:y1, x0:x1, z0:z1, :);
        ptMaskCylinder = ptMaskVol_ROI_HCLindex(y0:y1, x0:x1, z0:z1, :);
        % effective neurons only
        effNrids = unique(ptMaskCylinder(ptMaskCylinder ~= 0));
        effNrids = [k; setdiff(effNrids, k)];
        roiXtmp = roicenterXmergedsPerm(effNrids, subFrs);
        roiYtmp = roicenterYmergedsPerm(effNrids, subFrs);
        roiZtmp = roicenterZmergedsPerm(effNrids, subFrs);

        savePath = fullfile(allMoviesPerROIDir, ...
            ['hclROI_', num2str(k), '_', ...
            '_centerXYZ_', num2str(mX(k)), '_', ...
            num2str(mY(k)), '_', num2str(mZ(k))]);

        ca3D_maskContoursOnMIP_local_wSubFrames(fieldCylinder, ZXRatio, ptMaskCylinder, savePath, ...
            upSamplingFactor, roiXtmp, ...
            roiYtmp, ...
            roiZtmp, ...
            subFrs, effNrids, 'figFlag', 'on', 'textFlag', false)

        pause(0.2)
        close(gcf)   
    end
end

%% Make snapshotsOfSingleROIs_at_HighActivities_classifiedByKmeans

if par4.makeMov_singleROIs_Snapshots_atHighActivities
    
    snapshotsPerROIDir = fullfile(p.outputDir, 'snapshotsOfSingleROIs_atHighActivities');
    if ~isfolder(snapshotsPerROIDir); mkdir(snapshotsPerROIDir); end
    
    % Load imgArray
    if ~exist('imgArray')
        imgArray = uint16(zeros(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_));

        disp('======')
        disp('Loading frames:')

        parfor fr=1:MD.nFrames_
            currImage = uint16(MD.channels_(1).loadStack(fr));
            imgArray(:,:,:,fr) = currImage;
            fprintf(1, '%g ', fr); 
            if (mod(fr,50) == 0); fprintf('\n'); end    
        end
        fprintf(1, '\n')

        % handle 0 intensities: replace 0 intensity with minimum
        m0 = min(imgArray(imgArray > 0));
        imgArray(imgArray == 0) =  m0; 
    end
    
    %
    upSamplingFactor = 2;     %par4.upSamplingFactor_forMaskContourImg;
    ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;

    upbdDeformX = par2.upbdDeformX;
    upbdDeformZ = par2.upbdDeformZ;
    bandWidthX = par2.bandWidthX;
    bandWidthZ = par2.bandWidthZ;

    mX = round(mean(roicenterXmergedsPerm, 2, 'omitnan'));
    mY = round(mean(roicenterYmergedsPerm, 2, 'omitnan'));
    mZ = round(mean(roicenterZmergedsPerm, 2, 'omitnan'));
    xmax0 = MD.imSize_(2);
    ymax0 = MD.imSize_(1);
    zmax0 = MD.zSize_;

    for k = 1:numel(mX)
        disp(['== Movie for ROI ', num2str(k)])
        
        % select sub-frames when the activity is above the threshold
        act = actSigHCL(k, :);
        subFrs = find(act > par4.thresholdOfNormalizedActivity); 
        if isempty(subFrs)
            continue
        end

        x0 = max(mX(k)-upbdDeformX-bandWidthX, 1);
        x1 = min(mX(k)+upbdDeformX+bandWidthX, xmax0);
        y0 = max(mY(k)-upbdDeformX-bandWidthX, 1);
        y1 = min(mY(k)+upbdDeformX+bandWidthX, ymax0);
        z0 = max(mZ(k)-upbdDeformZ-bandWidthZ, 1);
        z1 = min(mZ(k)+upbdDeformZ+bandWidthZ, zmax0);

        fieldCylinder = imgArray(y0:y1, x0:x1, z0:z1, subFrs);
        ptMaskCylinder = ptMaskVol_ROI_HCLindex(y0:y1, x0:x1, z0:z1, subFrs);
        % effective neurons only
        effNrids = unique(ptMaskCylinder(ptMaskCylinder ~= 0));
        effNrids = [k; setdiff(effNrids, k)];
        roiXtmp = roicenterXmergedcsPerm_kmeans(effNrids, subFrs);
        roiYtmp = roicenterYmergedcsPerm_kmeans(effNrids, subFrs);
        roiZtmp = roicenterZmergedcsPerm_kmeans(effNrids, subFrs);
        
        savePath = fullfile(snapshotsPerROIDir, ...
            ['hclROI_', num2str(k), '_', 'centerXYZ_', num2str(mX(k)), '_', ...
            num2str(mY(k)), '_', num2str(mZ(k))]);

        ca3D_maskContoursOnMIP_local_wSubFrames(fieldCylinder, ZXRatio, ptMaskCylinder, savePath, ...
            upSamplingFactor, roiXtmp, ...
            roiYtmp, ...
            roiZtmp, ...
            subFrs, effNrids, 'figFlag', 'on', 'textFlag', false)

        pause(0.2)
        close(gcf)   
    end 
end

if par4.closeFigs; close all; end

%% 

disp('== done! ==')
disp('== End of VisualizeActivityMapMasksIndexedByHCL (Part 4/4) ==')

end
