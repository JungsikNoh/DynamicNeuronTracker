function SegmentDynamicROIs(movieData, imgArray, paramStruct)
% SegmentDynamicROIs fits a 3D Gaussian function to each local mean image
% of the tracked neurons and defines its central mass as a mask. This mask
% is combined with the neuron trajectory to generate a dynamic ROI.
%
%
% paramStruct contains the following parameters.
%
% PSFsigma_forROI (Point Spread Function Sigma, in pixel)
%   - 2 dim'l vectors, (sigma_XY, sigma_Z), 
%   which specify the standard deviation parameters in X/Y- and Z-direction 
%   of 3D Gaussian functions. 
%   - This sigma vector is just an initial value in fitting 3D Gaussian
%   function for neuron segmentation purpose. 
%   - Set to be a more typical value among the specified PSFsigma values in
%   Part 1.
% 
% levelOf3DGaussianDist_toSegment (0< X <1)
%   - To get a mask for a tracked neuron, images when the neuron is firing
%   (detected as bright point sources) are first averaged. After fitting a
%   3D Gaussian function to the averaged neuron firing image, a central X
%   portion of the Gaussian function is defined to be the mask of the
%   neuron (X = levelOf3DGaussianDist_toSegment).
%   - Larger X leads to a larger mask.
% 
% levelOf3DGaussianDist_toSegment_small (0< X <1)
%   - The algorithm automatically checks and reports overlaps between the 
%   neuron masks segmented using levelOf3DGaussianDist_toSegment. 
%   - If some pairs of masks overlap with each other in at least 10 frames,
%   then those neurons are automatically re-segmented using a smaller threshold
%   (= levelOf3DGaussianDist_toSegment_small). 
%   - After this adjustment, overlaps between the masks are reported again.
%   If it is severe, users can adjust these parameters.
%
% frameLengthForMovingMedian (in # of frames)
%   - Extracted Ca2+ activity time courses of tracked neurons (mean intensities within
%   masks) are normalized, by following (F - F_base)/F_base.
%   - This algorithm utilizes moving medians as the base activity (F_base).
%   - 'frameLengthForMovingMedian' specifies a time period length to
%   compute the moving medians. 
% 
% allROIsOutput (true or false)
%   - Whether or not to generate/save visualization output of maximum 
%   intensity projection of averaged images of tracked neurons and their
%   segmented masks.
%   - Setting 'true' may cost some computational time for plotting.
%   
% makeMov_MIPofROIs (true or false)
%   - Flag whether to generate three MIP videos visualizaing detailed
%   outcomes of the dynamic ROIs.
%
% closeFigs (true or false)
%   - Whether to close output plots while running because the pipeline
%   generates many figures. Output figures are always saved in the output
%   directory in .fig and .png formats. It simply executes 'close all;'
%   after saving the figures. 
%
% Output are saved under MD.outputDirectory_/SegmentDynamicROIs.
%
% Jungsik Noh, UTSW, 2023/08


%% Input/Set up outputDir

p.outputDir = fullfile(movieData.outputDirectory_, 'SegmentDynamicROIs'); 

MD = movieData;
par3 = paramStruct;
par3.figFlag = 'on';

bandWidthX = par3.bandWidthX;
bandWidthZ = par3.bandWidthZ;

Part1OutputDir = fullfile(MD.outputDirectory_, 'DetectBrightPointSources');
Part2OutputDir = fullfile(MD.outputDirectory_, 'TrackJitteringFlickering');

if ~isfolder(p.outputDir); mkdir(p.outputDir); end

plotDir = fullfile(p.outputDir, 'detailedPlots');
if ~isfolder(plotDir); mkdir(plotDir); end

if par3.allROIsOutput
    allROIDir = fullfile(p.outputDir, 'allROIsOutput');
    if ~isfolder(allROIDir)
        mkdir(allROIDir); 
    else
        delete(fullfile(allROIDir, '*'));
    end
end

%% Load output from Part 1 & 2

load(fullfile(Part1OutputDir, 'PSourMat2.mat'))
load(fullfile(Part2OutputDir, 'tracks_PooledEMoutput_MergedOverlaps_L2Merged.mat')) 

%% Segment a mask for each track 

ROIcell = cell(1);
trmax = size(tracksXcmerged0sL2, 1);
convInds = zeros(trmax, 1);

for itr = 1:trmax 
    xvec = tracksXcmerged0sL2(itr, :);
    yvec = tracksYcmerged0sL2(itr, :);
    zvec = tracksZcmerged0sL2(itr, :);

    fprintf(1, '%g ', itr); 
    if (mod(itr, 50) == 0 || (itr == trmax)); fprintf('\n'); end 
    
    [currMask, convId, figtmp] = segmentMaskOfPtSource_3D_anisoGaussian3steps(xvec, yvec, zvec, ...
        par3.PSFsigma_forROI, par3.bandWidthX, par3.bandWidthZ, imgArray, PSourMat2, ...
        'ConfRegionLevel', par3.levelOf3DGaussianDist_toSegment, ...
        'figVisible', 'off', 'figOut', par3.allROIsOutput);
    
    convInds(itr) = convId;
    ROIcell{itr} = currMask;
    
    if par3.allROIsOutput && ~isempty(figtmp)
        figfname = ['ROI_', num2str(itr), '.png'];
        saveas(figtmp, fullfile(allROIDir, figfname), 'png');
    end    
end

%% Exclude rare cases of convId==NaN (cannot collect a clear avg image)

excldInd = isnan(convInds);
ROIcell = ROIcell(~excldInd);
tracksXcmerged0csL2 = tracksXcmerged0csL2(~excldInd, :);
tracksYcmerged0csL2 = tracksYcmerged0csL2(~excldInd, :);
tracksZcmerged0csL2 = tracksZcmerged0csL2(~excldInd, :);
tracksXcmerged0sL2 = tracksXcmerged0sL2(~excldInd, :);
tracksYcmerged0sL2 = tracksYcmerged0sL2(~excldInd, :);
tracksZcmerged0sL2 = tracksZcmerged0sL2(~excldInd, :);

%% Check overlaps of dynamic ROIs and signal extraction

[ROIPairsOverlappingMoreThan10Frames, ROIsOverlappingMoreThan10Frames] = ...
    checkOverlapsOfDynamicROIs(MD, ROIcell, tracksXcmerged0csL2, tracksYcmerged0csL2, ...
                                tracksZcmerged0csL2, par3.bandWidthX, par3.bandWidthZ);

%%  Adjust ROI confRegionLevel for smaller ones

if ~isempty(ROIsOverlappingMoreThan10Frames)

    disp('== Dynamic ROI IDs overlapping with another in >=10 time frames ==')
    disp(ROIsOverlappingMoreThan10Frames')
    disp('== Masks of these neurons are adjusted to smaller ones to reduce the overlaps.')

    for k = 1:numel(ROIsOverlappingMoreThan10Frames)
        itr = ROIsOverlappingMoreThan10Frames(k);
        xvec = tracksXcmerged0csL2(itr, :);
        yvec = tracksYcmerged0csL2(itr, :);
        zvec = tracksZcmerged0csL2(itr, :);

        fprintf(1, '%g ', itr); 
        if (mod(k, 50) == 0 || (k == numel(ROIsOverlappingMoreThan10Frames))); fprintf('\n'); end 

        % levelOf3DGaussianDist_toSegment_small
        [currMask, convId, figtmp] = segmentMaskOfPtSource_3D_anisoGaussian3steps(xvec, yvec, zvec, ...
            par3.PSFsigma_forROI, par3.bandWidthX, par3.bandWidthZ, imgArray, PSourMat2, ...
            'ConfRegionLevel', par3.levelOf3DGaussianDist_toSegment_small, ...
            'figVisible', 'off', 'figOut', true);

        ROIcell{itr} = currMask;

        if par3.allROIsOutput
            figfname = ['ROI_', num2str(itr), '.png'];
            saveas(figtmp, fullfile(allROIDir, figfname), 'png');
        end
    end    
end

%% Check ROI overlaps after adjustment

[ROIPairsOverlappingMoreThan10Frames, ROIsOverlappingMoreThan10Frames] = ...
    checkOverlapsOfDynamicROIs(MD, ROIcell, tracksXcmerged0csL2, tracksYcmerged0csL2, ...
                                tracksZcmerged0csL2, par3.bandWidthX, par3.bandWidthZ);

if ~isempty(ROIsOverlappingMoreThan10Frames)
    disp('== Dynamic ROI IDs, still overlapping with another in >=10 time frames ==')
    disp(ROIsOverlappingMoreThan10Frames')
    disp('== You may try a smaller value for ''levelOf3DGaussianDist_toSegment_small''. ==')
    disp('== Algorithm proceeds to the next step. ==')
end

%% Signal extraction!

psSig = nan(size(tracksXcmerged0csL2));
frmax = size(tracksXcmerged0csL2, 2);
trmax = size(tracksXcmerged0csL2, 1);

for k = 1:trmax
    xvec = tracksXcmerged0csL2(k, :);
    yvec = tracksYcmerged0csL2(k, :);
    zvec = tracksZcmerged0csL2(k, :);
    frs = 1:frmax;

    [~, Patches] = localPatchSampling3D(xvec, yvec, zvec, frs, bandWidthX, ...
        bandWidthZ, imgArray);

    mask0 = double(ROIcell{k});
    mask0(mask0 == 0) = nan;
    maskCube = repmat(mask0, 1, 1, 1, frmax);
    tmp = maskCube .* Patches;
    tmp2 = reshape(tmp, [], frmax);
    psSig(k, :) = nanmean(tmp2, 1);    
end

%% Minor adjustment of tracks using centers of ROIs

roicenterXmergedcs = nan(size(tracksXcmerged0csL2));
roicenterYmergedcs = nan(size(tracksXcmerged0csL2));
roicenterZmergedcs = nan(size(tracksXcmerged0csL2));

[ny, nx, nz] = size(ROIcell{1});
ys = (-ny+1)/2:(ny-1)/2;
xs = (-nx+1)/2:(nx-1)/2;
zs = (-nz+1)/2:(nz-1)/2;
coorys = repmat(repmat(ys(:), 1, nx), 1, 1, nz);
coorxs = repmat(repmat(xs, ny, 1), 1, 1, nz);
coorzs0 = repmat(repmat(zs(:), 1, ny), 1, 1, nx);
coorzs = shiftdim(coorzs0, 1);

disp('== center location adjustment in x y z')
for itr = 1:trmax    
    currMask = ROIcell{itr};
    ROIsize = sum(currMask(:));
    massCenter1 = round(sum(coorys(:) .* currMask(:)) / ROIsize);
    massCenter2 = round(sum(coorxs(:) .* currMask(:)) / ROIsize);
    massCenter3 = round(sum(coorzs(:) .* currMask(:)) / ROIsize);

    mcenter = [massCenter2, massCenter1, massCenter3];
    disp([itr, mcenter])
    
    roicenterXmergedcs(itr, :) = tracksXcmerged0csL2(itr, :) + massCenter2;
    roicenterYmergedcs(itr, :) = tracksYcmerged0csL2(itr, :) + massCenter1;
    roicenterZmergedcs(itr, :) = tracksZcmerged0csL2(itr, :) + massCenter3;
end

%% output

roicenterXmergeds = roicenterXmergedcs;
roicenterYmergeds = roicenterYmergedcs;
roicenterZmergeds = roicenterZmergedcs;
roicenterXmergeds(isnan(tracksXcmerged0sL2)) = NaN;
roicenterYmergeds(isnan(tracksYcmerged0sL2)) = NaN;
roicenterZmergeds(isnan(tracksZcmerged0sL2)) = NaN;

%% Indexed ROIs

fig4 = figure('Visible', par3.figFlag);
meanVol = mean(imgArray, 4);
qq = quantile(reshape(meanVol, 1,[]), [0.5, 0.9999]);
ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;
[~,~,~,cthree] = j_computeMIPs(meanVol,ZXRatio, qq(1), qq(2));

imagesc(cthree)
colormap  gray % default
colorbar
axis off

trmax = size(roicenterXmergedcs, 1); 
maxXBorder=MD.getDimensions('X');
maxYBorder=MD.getDimensions('Y');
maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);

strSize0 = 4;
X = round(mean(roicenterXmergedcs, 2));
Y = round(mean(roicenterYmergedcs, 2));
Z = round(mean(roicenterZmergedcs, 2));
Z = Z .* ZXRatio - ZXRatio/2;

hold on
scatter(X, Y, 10,  'MarkerEdgeColor', 'r')
scatter(X, Z+maxYBorder+strSize0, 10,  'MarkerEdgeColor', 'r')
scatter(Z+maxXBorder+strSize0, Y, 10,  'MarkerEdgeColor', 'r')

for k = 1:numel(ROIcell)       
    x = X(k)+1; y = Y(k)+1; z = Z(k)+1;
    text(x, y, num2str(k), 'FontSize', 8, 'Color', 'r');    
    text(x, z+maxYBorder+strSize0, num2str(k), 'FontSize', 8, 'Color', 'r')
    text(z+maxXBorder+strSize0, y, num2str(k), 'FontSize', 8, 'Color', 'r')
end

%% save output

levelOf3DGaussianDist_toSegment = par3.levelOf3DGaussianDist_toSegment;
levelOf3DGaussianDist_toSegment_small = par3.levelOf3DGaussianDist_toSegment_small;

save(fullfile(p.outputDir, 'ROIcell_3DGaussian.mat'), 'ROIcell', ...
    'levelOf3DGaussianDist_toSegment', 'levelOf3DGaussianDist_toSegment_small', ...
    'ROIPairsOverlappingMoreThan10Frames', 'ROIsOverlappingMoreThan10Frames')

save(fullfile(p.outputDir, 'roicenterXYZ.mat'), 'roicenterXmergedcs', ...
    'roicenterYmergedcs', 'roicenterZmergedcs', 'roicenterXmergeds', ...
    'roicenterYmergeds', 'roicenterZmergeds')   

save(fullfile(p.outputDir, 'par3.mat'), 'par3')

%% a metric of the amount of jittering via roicenterXYZmergeds

difx = diff(roicenterXmergeds, 1, 2);
dify = diff(roicenterYmergeds, 1, 2);
difz = diff(roicenterZmergeds, 1, 2);
difx(isnan(difx)) = 0;
dify(isnan(dify)) = 0;
difz(isnan(difz)) = 0;
totJitters = [NaN, sum(abs(difx), 1) + sum(abs(dify), 1) + sum(abs(difz), 1)];

fig11 = figure('Visible', par3.figFlag); 
plot(totJitters)
xlabel('Frame'); ylabel('Amount of Jitters')

saveas(fig11, fullfile(plotDir, ['totalJitters', '.fig']), 'fig')
saveas(fig11, fullfile(plotDir, ['totalJitters', '.png']), 'png')


%% MIP of ROI centers

if par3.makeMov_MIPofROIs

    savePath = fullfile(p.outputDir, 'MIP_ROIcenter_imputedByEnds_roicenterXmergedcs');
    ca3D_plotTracksOnMIP(MD, imgArray, roicenterXmergedcs, roicenterYmergedcs, ...
                  roicenterZmergedcs, savePath, 'figFlag', par3.figFlag)

    savePath = fullfile(p.outputDir, 'MIP_ROIcenter_whenFiring_roicenterXmergeds');
    ca3D_plotTracksOnMIP(MD, imgArray, roicenterXmergeds, roicenterYmergeds, ...
                  roicenterZmergeds, savePath, 'figFlag', par3.figFlag)

end

%% heatmap of signals

numSig = size(psSig, 1);
vals = psSig(:);
q = quantile(vals, [0.001, 0.995]);

fig5 = figure('Visible', par3.figFlag); 
imagesc(psSig, [q(1), q(2)])
colorbar; colormap(jet)
title(['Raw signals of point sources: ', num2str(numSig) ' neurons'])
xlabel('Frame')

%% heatmap of log-transformed-signals

psSig_log = log(psSig + 1);

numSig = size(psSig_log, 1);
vals = psSig_log(:);
q = quantile(vals, [0.001, 0.995]);

fig52 = figure('Visible', par3.figFlag); 
imagesc(psSig_log, [q(1), q(2)])
colorbar; colormap(jet)
title(['Raw signals (log) of point sources: ', num2str(numSig) ' neurons'])
xlabel('Frame')

%% F-F0/F0 (moving median normalization)

tmp = nan(size(psSig));
normSig = nan(size(psSig));

for k = 1:size(psSig, 1)
    x = psSig(k, :);
    xmed = movmedian(x, par3.frameLengthForMovingMedian);
    tmp(k,:) = xmed;
    normx = (x - xmed) ./ xmed;
    normSig(k, :) = normx;
end

vals = normSig(:);
q = quantile(vals, [0.001, 0.995]);

fig7 = figure('Visible', par3.figFlag); 
imagesc(normSig, [q(1), q(2)])
colorbar; colormap(jet)
title('Moving-median Normalized Signals (F-F0/F0)')
xlabel('Frame')

%% F-F0/F0 (moving median normalization) for psSig_log

tmp = nan(size(psSig_log));
normSig_log = nan(size(psSig_log));

for k = 1:size(psSig_log, 1)
    x = psSig_log(k, :);
    xmed = movmedian(x, par3.frameLengthForMovingMedian);
    tmp(k,:) = xmed;
    normx = (x - xmed) ;
    normSig_log(k, :) = normx;
end

vals = normSig_log(:);
q = quantile(vals, [0.001, 0.995]);

fig72 = figure('Visible', par3.figFlag); 
imagesc(normSig_log, [q(1), q(2)])
colorbar; colormap(jet)
title('Moving-median Normalized Signals (log) (F-F0)')
xlabel('Frame')

%% Make mask vol movie

ptMaskVol = zeros(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_);  % row, col
ptMaskVol = uint16(ptMaskVol);

frmax = size(tracksXcmerged0csL2, 2);
trmax = size(tracksXcmerged0csL2, 1);
for k = 1:trmax
    xvec = tracksXcmerged0csL2(k, :);
    yvec = tracksYcmerged0csL2(k, :);
    zvec = tracksZcmerged0csL2(k, :);    
    mask0 = double(ROIcell{k});
    
    for fr = 1:frmax
        xwin = xvec(fr)-bandWidthX:xvec(fr)+bandWidthX;
        ywin = yvec(fr)-bandWidthX:yvec(fr)+bandWidthX;        
        zwin = zvec(fr)-bandWidthZ:zvec(fr)+bandWidthZ; 
        maskVoltmp = ptMaskVol(ywin, xwin, zwin, fr);       
        maskVoltmp(mask0 == 1) = mod(k,100) + 1;
        ptMaskVol(ywin, xwin, zwin, fr) = maskVoltmp;
    end
end

%% mask volume

ptMaskVol2 = ptMaskVol + 10;
ptMaskVol2(ptMaskVol == 0) = 0;
ptMaskVol3 = mat2gray(ptMaskVol2);

clear ptMaskVol
clear ptMaskVol2

%% vol to mip

if par3.makeMov_MIPofROIs
    savePath = fullfile(p.outputDir, 'MIP_dynamicROIMasks');
    ca3D_Vol2MIP(MD, ptMaskVol3, savePath)
end

%% tif output for dynamic ROI Masks

savePath = fullfile(p.outputDir, 'tif_volume_per_frame_ROIMasks');
if ~isfolder(savePath); mkdir(savePath); end

ptMaskVol4 = uint8((2^8-1)*ptMaskVol3);

for fr=1:MD.nFrames_
    fprintf(1, '%g ', fr)
    if (mod(fr,50) == 0); fprintf('\n'); end
    fname = ['Maskimg_', sprintf('%04d', fr), '.tif'];
    tmp = ptMaskVol4(:,:,:,fr);
    for k = 1:size(tmp,3)
        imwrite(tmp(:,:,k), fullfile(savePath, fname), 'WriteMode', 'append');
    end
end

clear ptMaskVol3

%% Save mask vol in a .mat format

ptMaskVol_ROI = zeros(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_);  % row, col
ptMaskVol_ROI = uint16(ptMaskVol_ROI);

frmax = size(tracksXcmerged0csL2, 2);
trmax = size(tracksXcmerged0csL2, 1);
for k = 1:trmax
    xvec = tracksXcmerged0csL2(k, :);
    yvec = tracksYcmerged0csL2(k, :);
    zvec = tracksZcmerged0csL2(k, :);
    
    mask0 = double(ROIcell{k});
    
    for fr = 1:frmax
        xwin = xvec(fr)-bandWidthX:xvec(fr)+bandWidthX;
        ywin = yvec(fr)-bandWidthX:yvec(fr)+bandWidthX;        
        zwin = zvec(fr)-bandWidthZ:zvec(fr)+bandWidthZ; 
        maskVoltmp = ptMaskVol_ROI(ywin, xwin, zwin, fr);       
        maskVoltmp(mask0 == 1) = k;
        ptMaskVol_ROI(ywin, xwin, zwin, fr) = maskVoltmp;
    end
end

save(fullfile(p.outputDir, 'ptMaskVol_ROI.mat'), 'ptMaskVol_ROI', '-v7.3')

%% save figure

saveas(fig4, fullfile(plotDir, ['indexedROIsOnMIP', '.fig']), 'fig')
saveas(fig4, fullfile(plotDir, ['indexedROIsOnMIP', '.png']), 'png')

saveas(fig5, fullfile(plotDir, ['RawSignalMap', '.fig']), 'fig')
saveas(fig5, fullfile(plotDir, ['RawSignalMap', '.png']), 'png')

saveas(fig7, fullfile(plotDir, ['normalizedSignalMap', '.fig']), 'fig')
saveas(fig7, fullfile(plotDir, ['normalizedSignalMap', '.png']), 'png')

saveas(fig52, fullfile(plotDir, ['RawSignalMap-log', '.fig']), 'fig')
saveas(fig52, fullfile(plotDir, ['RawSignalMap-log', '.png']), 'png')

saveas(fig72, fullfile(plotDir, ['normalizedSignalMap-log', '.fig']), 'fig')
saveas(fig72, fullfile(plotDir, ['normalizedSignalMap-log', '.png']), 'png')

if par3.closeFigs; close all; end

%% save output

save(fullfile(p.outputDir, 'extractedSignals_indexedByLocations.mat'), ...
    'psSig', 'normSig', 'ROIcell', 'tracksXcmerged0csL2', 'tracksYcmerged0csL2', ...
    'tracksZcmerged0csL2', 'psSig_log', 'normSig_log')

%% ROI volume distribution

ptMaskVol_ROI1 = double(ptMaskVol_ROI(:,:,:,1));
numPixels = nan(trmax, 1);
roiVolume = nan(trmax, 1);

for k = 1:trmax
    numPixels(k) = sum(ptMaskVol_ROI1(:) == k);
end

roiVolume = numPixels .* (MD.pixelSize_/1000)^2 .* (MD.pixelSizeZ_/1000);  % vol in um^3

ROIvolumeTable = table(numPixels, roiVolume);

%% save

save(fullfile(p.outputDir, 'ROIvolumeTable.mat'), 'ROIvolumeTable')

fig9 = figure('Visible', par3.figFlag); 
histogram(numPixels); title('Number of voxels of ROIs')

fig10 = figure('Visible', par3.figFlag); 
histogram(roiVolume); title('Volumes of ROIs in um^3')

saveas(fig9, fullfile(p.outputDir, 'numPixels.png'), 'png')
saveas(fig10, fullfile(p.outputDir, 'roiVolume.png'), 'png')

saveas(fig9, fullfile(p.outputDir, 'numPixels.fig'), 'fig')
saveas(fig10, fullfile(p.outputDir, 'roiVolume.fig'), 'fig')

%% example time series showing movMedian normalization

if (trmax <= 3)
    exId = 1:trmax;
else
    tmp = randsample(trmax, 3);
    exId = sort(tmp);
end

fig8 = cell(1, 3);
for i = 1:numel(exId)
    k = exId(i);
    x = psSig(k, :);
    xmed = movmedian(x, par3.frameLengthForMovingMedian);
    normx = normSig(k, :);

    fig8{i} = figure('Visible', par3.figFlag);
    subplot(2, 1, 1)
    plot(x)
    hold on
    plot(xmed)
    title(['ROI: ', num2str(k)])
    legend('single trajectory', 'moving median')

    subplot(2, 1, 2)
    plot(normx)
    xlabel('Frames')
    title('Normalized Sig = (F-F0)/F0, F0 = movmedian')
    refline([0 0])
end

%% example TS plots

for i = 1:3
    if ~isempty(fig8{i})
        saveas(fig8{i}, fullfile(plotDir, ['exampleTimeSeries_normalization_', num2str(i), '.fig']), 'fig')
        saveas(fig8{i}, fullfile(plotDir, ['exampleTimeSeries_normalization_', num2str(i), '.png']), 'png')
    end
end

if par3.closeFigs; close all; end

%%

disp('== done! ==')
disp('== End of SegmentDynamicROIs (Part 3/4) ==')

end
