function DetermineResponsiveness(movieData, paramStruct)
% DetermineResponsiveness implements statistical analyses to determine 
% the combinatorial responsiveness of each segmented neuron to given
% stimuli for one movie. It assumes that multiple stimuli are applied 
% randomly and repeatedly. Stimulation information needs to be provided by 
% a .csv file.
%
% paramStruct contains the following parameters.
%
% stimulationFile
%   - A file path and name for the stimulation information (.csv file). The file
%   is typically located at the same folder with the movieData object (.mat)
%   file.
%
% stimulusInterval (in # of frames)
%   - A time interval (in frames) between consecutive two different stimuli. 
%
% negCtrlStimulusName
%   - 'stimulusLabel' for a negative control condition. To determine the
%   responsiveness, the statistical test requires a negative control.
% 
% activityWeightedByImgCorr (true or false)
%   - Whether to analyze 'weighted' or 'unweighted' calcium activities. 
%   - This option depends on which type was specified in Part 4
%   (VisualizeActivityMapMasksIndexedByHCL).
%
% make_AllTSPlots (true or false)
%   - Flag whether to generate time series plots for all ROIs.
%
% make_BoxplotsOfTtests (true or false)
%   - Flag whether to generate boxplots of t-tests for all ROIs.
%   
% make_AllResponseCurves_perStimulus (true or false)
%   - Flag whether to generate plots of across-trial response curves and their mean
%   curves per-stimulus for all ROIs.
%
% makeROIMovieWithStimulusLabel_atHighActivities (true or false)
%   - Flag whether to generate a MIP video displaying ROIs and stimulation 
%   status. ROIs are shown only when activities are determined to be high 
%   by K-means clustering.
%
% makeROIMovieWithStimulusLabel_firingAnnotated (true or false)
%   - The same as the above except that ROIs are shown only when ROIs are 
%   determined to be firing via spatial correlation.
%
% makeAllFramesOfSingleROIs_firingAnnotated (true or false)
%   - Flag whether to generate all of the single ROI local MIP videos that
%   also show stimulation status. The MIP videos display all time frames 
%   and ROIs are shown when they are determined to be firing. 
%
% makeSnapshotsOfSingleROIs_atHighActivities (true or false)
%   - Flag whether to generate all of the single ROI local MIP videos that
%   also show stimulation status. The MIP videos display a subset of time 
%   frames when the target ROI's activities are determined to be high by 
%   K-means clustering.
%
% FalseDiscoveryRate_threshold
%   - The level of significance for adjusted P-values.
%
% thresholdForMinimalNumberOfNeuronsForEachStimulusCombination (X: integer)
%   - After t-tests and FDR controls, we select significant stimulus
%   subsets that activated at least a certain number of neurons (= X). See
%   the paper.
%
% numNaNRows_btwnNeuronClusters
%   - When mean activities per stimulus are visualized along with
%   identified stimulus subsets, numNaNRows_btwnNeuronClusters specifies
%   empty spaces between the stimulus subsets for visualization purposes.
%
% closeFigs (true or false)
%   - Whether to close output plots while running because the pipeline
%   generates many figures. Output figures are always saved in the output
%   directory in .fig and .png formats. It simply executes 'close all;'
%   after saving the figures. 
%
% Jungsik Noh, UTSW, 2023/08

%% Set up params

MD = movieData;
par5 = paramStruct;
par5.figFlag = 'on';

if par5.activityWeightedByImgCorr
    p.outputDir = fullfile(MD.outputDirectory_, 'DetermineResponsiveness_weighted'); 
    p.outputDir2 = fullfile(MD.outputDirectory_, 'DetermineResponsiveness_weighted', 'detailedAnalysis');
else
    p.outputDir = fullfile(MD.outputDirectory_, 'DetermineResponsiveness_unweighted'); 
    p.outputDir2 = fullfile(MD.outputDirectory_, 'DetermineResponsiveness_unweighted', 'detailedAnalysis');
end

if par5.activityWeightedByImgCorr
    Part4OutputDir = fullfile(MD.outputDirectory_, 'VisualizeActivityMapMasksIndexedByHCL_weighted');
else
    Part4OutputDir = fullfile(MD.outputDirectory_, 'VisualizeActivityMapMasksIndexedByHCL_unweighted');
end

Part2OutputDir = fullfile(MD.outputDirectory_, 'TrackJitteringFlickering');

if ~isfolder(p.outputDir); mkdir(p.outputDir); end
if ~isfolder(p.outputDir2); mkdir(p.outputDir2); end
save(fullfile(p.outputDir, 'par5.mat'), 'par5')

%% Load input data

% load psSig, actSigHVL
load(fullfile(Part4OutputDir, 'actSig_psSig_HCL-indexed.mat'))

%% 1. stimuliTable

stimuliTable0 = readtable(par5.stimulationFile);
analFrameMax = size(actSigHCL, 2);

% startFrame <= analFrameMax
ind = (stimuliTable0{:, 2} > analFrameMax);
stimuliTable = stimuliTable0(~ind, :);
stimuliNames = unique(stimuliTable{:, 1});        % cell object
numStimuli = numel(stimuliNames);
stimuliInfo = struct();

for i = 1:numStimuli
    stimuliInfo(i).stimuliNames = stimuliNames{i};

    ind = strcmp(stimuliTable{:, 1}, stimuliNames{i});
    stimuliInfo(i).startFrames = reshape(stimuliTable{ind, 2}, 1, []);
    stimuliInfo(i).endFrames = reshape(stimuliTable{ind, 3}, 1, []);
end

% stimuliMatrix
frmax = size(actSigHCL, 2);
stimuliMatrix = zeros(numStimuli, frmax);

for i = 1:size(stimuliTable, 1)
    stimName = stimuliTable{i, 1};
    stFr = stimuliTable{i, 2};
    endFr = stimuliTable{i, 3};
    id0 = find(strcmp(stimuliNames, stimName{:}));
    stimuliMatrix(id0, stFr:endFr) = 1;
end

%% 2. Make a ROI movie with stimulus labels

if (par5.makeROIMovieWithStimulusLabel_atHighActivities || par5.makeROIMovieWithStimulusLabel_firingAnnotated)
    
    load(fullfile(Part4OutputDir, 'par4.mat'))
    load(fullfile(Part4OutputDir, 'ptMaskVol_ROI_HCLindex.mat'))
    load(fullfile(Part4OutputDir, 'roicenterXYZ_HCLPerm.mat'))
    
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
    
    % Visualize Contours of masks indexed by HCL 
    upSamplingFactor = par4.upSamplingFactor_forMaskContourImg;
    ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;

    if par5.makeROIMovieWithStimulusLabel_atHighActivities
        savePath = fullfile(p.outputDir, 'maskContoursOnMIP_WithStimulusLabel_atHighActivities');
        ca3D_maskContoursOnMIP_whenFiring_wStimulus(imgArray, ZXRatio, ptMaskVol_ROI_HCLindex, savePath, ...
            upSamplingFactor, roicenterXmergedcsPerm_kmeans, roicenterYmergedcsPerm_kmeans, ...
            roicenterZmergedcsPerm_kmeans, ...
            stimuliMatrix, stimuliNames, 'figFlag', par5.figFlag)
    end
    
    if par5.makeROIMovieWithStimulusLabel_firingAnnotated
        savePath = fullfile(p.outputDir, 'maskContoursOnMIP_WithStimulusLabel_firingAnnotated');
        ca3D_maskContoursOnMIP_whenFiring_wStimulus(imgArray, ZXRatio, ptMaskVol_ROI_HCLindex, savePath, ...
            upSamplingFactor, roicenterXmergedsPerm, roicenterYmergedsPerm, ...
            roicenterZmergedsPerm, ...
            stimuliMatrix, stimuliNames, 'figFlag', par5.figFlag)
    end
end

%% 4. Avgeraged Activity

m2 = mean(actSigHCL, 1);

f1 = figure('Position', [100 100 1140 413], 'Visible', par5.figFlag);
plot(m2, 'Color', 'k', 'DisplayName', 'Avg Activity')
refline([0,0])
hold on
ax=gca;

% plot startFrame
colmap = jet(numStimuli);

for s = 1:numStimuli
    if ~isempty(stimuliInfo(s).startFrames)
        for i = 1:numel(stimuliInfo(s).startFrames)
            h(s) = plot([stimuliInfo(s).startFrames(i), stimuliInfo(s).startFrames(i)], ...
                ax.YLim, 'Color', colmap(s, :), 'DisplayName', stimuliInfo(s).stimuliNames);
        end
    end    
end

lg = legend(h); 
lg.Location = 'eastoutside'; 

title(['Avgeraged Activity'])
xlabel('Frame')
ylabel('Activity')

%% save

saveas(f1, fullfile(p.outputDir, 'avgActivy_stimuli.fig'), 'fig')
saveas(f1, fullfile(p.outputDir, 'avgActivy_stimuli.png'), 'png')

save(fullfile(p.outputDir, 'stimuliInfo.mat'), 'stimuliInfo')
save(fullfile(p.outputDir, 'stimuliTable.mat'), 'stimuliTable')

%% 5. activity_stimuli_map

q = quantile(actSigHCL(:), [0.005, 0.995]);

f2 = figure('Position', [100 100 1140 613], 'Visible', par5.figFlag);

subplot('Position', [0.1 0.45, 0.8, 0.45])
imagesc(actSigHCL, [q(1), q(2)])
ylabel('Neuron Index')
cmap = redwhiteblue(q(1), q(2), 256);
colorbar; colormap(cmap)
title('Neural Activity (\DeltaF/F)')

subplot('Position', [0.1, 0.1, 0.8, 0.25])
imtmp = imagesc(stimuliMatrix);
bwmap = gray(2);
colorbar; colormap(gca, bwmap)
title('Stimulation')
xlabel('Frames')

ax = gca;
ax.YTick = 1:numel(stimuliNames);
ax.YTickLabel = stimuliNames;

%% save

saveas(f2, fullfile(p.outputDir, 'activity_stimuli_map.fig'), 'fig')
saveas(f2, fullfile(p.outputDir, 'activity_stimuli_map.png'), 'png')

%% 6. allTSplots

nCells = size(actSigHCL, 1);
frmax = size(actSigHCL, 2);
frMargin = round(0.01 * frmax);
frStart = 1 - frMargin;
frEnd = frmax + frMargin;
    
if par5.make_AllTSPlots

    allTSout = fullfile(p.outputDir, 'allTSplots');
    if ~isfolder(allTSout); mkdir(allTSout); end

    % All TS plots
    for j = 1:5:nCells

        f3 = figure('Position', [10 10 1000 700], 'Visible', par5.figFlag);

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

        saveas(f3, fullfile(allTSout, ['allTS_', num2str(j), '.png']), 'png')
        pause(0.2)
        close(f3)
    end

end

%% 7. meanActivityOfNeuronsPerStimulus
%% colect activity blocks per stimulus

activityPerStimulusInfo = stimuliInfo;
addedBlock = nan(nCells, par5.stimulusInterval);
actSigHCL1 = [actSigHCL, addedBlock];

for s = 1:numStimuli
    stFrs = activityPerStimulusInfo(s).startFrames;
    tmp = nan(nCells, par5.stimulusInterval, numel(stFrs));
    
    for i = 1:numel(stFrs)
        tmp1 = actSigHCL1(:, stFrs(i):(stFrs(i) + par5.stimulusInterval - 1) );
        tmp(:, :, i) = tmp1;
    end    
    activityPerStimulusInfo(s).blockactSigHCL = tmp;
end

% activity TS per stimulus 
stimulusMeanMat = nan(nCells, par5.stimulusInterval);
[activityPerStimulusInfo.stimulusMeanMat] = deal(stimulusMeanMat);

for s = 1:numStimuli
    tmp = activityPerStimulusInfo(s).blockactSigHCL;
    avgMat = mean(tmp, 3, 'omitnan');
    activityPerStimulusInfo(s).stimulusMeanMat = avgMat;
end

%% maximizingFrame_PerNeuronStimulus (fr0) 

maximizingFrame_PerNeuronStimulus = nan(nCells, 1);
[activityPerStimulusInfo.maximizingFrame_PerNeuronStimulus] = deal(maximizingFrame_PerNeuronStimulus);

for s = 1:numStimuli
    stimulusMeanMat = activityPerStimulusInfo(s).stimulusMeanMat;   % nCells-by-stimulusInterval
    
    for i = 1:nCells
        [~, fr] = max(stimulusMeanMat(i, :));
        activityPerStimulusInfo(s).maximizingFrame_PerNeuronStimulus(i) = fr;
    end    
end

% activity at the maximizing frame
activityAtMaximizingFrame = [];
[activityPerStimulusInfo.activityAtMaximizingFrame] = deal(activityAtMaximizingFrame);

for s = 1:numStimuli
    activityBlock = activityPerStimulusInfo(s).blockactSigHCL;     % nCells*Interval*repeats
    fr0 = activityPerStimulusInfo(s).maximizingFrame_PerNeuronStimulus;
    tmp = nan( nCells, size(activityBlock, 3) );
    
    for i = 1:nCells
        tmp(i, :) = reshape(activityBlock(i, fr0(i), :), 1, []);
    end
    
    activityPerStimulusInfo(s).activityAtMaximizingFrame = tmp;    
end

%% mean activity over the 1st half of par5.stimulusInterval 

halfIntervalLength = round(par5.stimulusInterval / 2);

% mean activity over the 1st half of stimulusInterval
for s = 1:numStimuli
    activityBlock = activityPerStimulusInfo(s).blockactSigHCL;     % nCells*Interval*repeats
    actBlockHalfInterval = activityBlock(:, 1:halfIntervalLength, :);
    tmp = squeeze(mean(actBlockHalfInterval, 2));
    activityPerStimulusInfo(s).meanActivityOverHalfInterval = tmp;    
end

%% save

save(fullfile(p.outputDir, 'activityPerStimulusInfo.mat'), 'activityPerStimulusInfo')

%% meanActivityOfNeuronsPerStimulus

meanActivityOfNeuronsPerStimulus = nan(nCells, numStimuli);

for s = 1:numStimuli
    for i = 1:nCells
        tmp = activityPerStimulusInfo(s).meanActivityOverHalfInterval(i, :);
        meanActivityOfNeuronsPerStimulus(i, s) = mean(tmp, 'omitnan');
    end
end

f5 = figure('Visible', par5.figFlag);
q = quantile(meanActivityOfNeuronsPerStimulus(:), [0.005, 0.995]);
imagesc(meanActivityOfNeuronsPerStimulus, [q(1), q(2)])
ylabel('Neuron Index')
cmap = redwhiteblue(q(1), q(2), 256);
colorbar; colormap(cmap)
title('Neural Activity (\DeltaF/F)')

ax = gca;
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;

%% save

save(fullfile(p.outputDir, 'meanActivityOfNeuronsPerStimulus.mat'), ...
    'meanActivityOfNeuronsPerStimulus')
saveas(f5, fullfile(p.outputDir, 'meanActivityOfNeuronsPerStimulus.png'), 'png')
saveas(f5, fullfile(p.outputDir, 'meanActivityOfNeuronsPerStimulus.fig'), 'fig')

if par5.closeFigs
    close all;
end

%% myHCL for meanActivityOfNeuronsPerStimulus

inputMat = meanActivityOfNeuronsPerStimulus';

linkageMethod = 'complete';

[~, fig2, fig3, outperm, hclmap, fig4, fig5, outpermCol, hclmap2] = ...
    myHCL(inputMat, linkageMethod, ...
    'rowNames', stimuliNames, 'xlabel', 'Neuron Index', 'ylabel', "", ...
    'title', ['Neural Activity - ', linkageMethod], 'figFlag', par5.figFlag);

%% save

save(fullfile(p.outputDir2, 'myHCL_meanActivityOfNeuronsPerStimulus.mat'), ...
                'inputMat', 'outperm', 'hclmap', 'outpermCol', 'hclmap2')
            
saveas(fig2, fullfile(p.outputDir2, 'myHCL_dendro1_meanActivityOfNeuronsPerStimulus.png'), 'png')
saveas(fig3, fullfile(p.outputDir2, 'myHCL_rowSorted_meanActivityOfNeuronsPerStimulus.png'), 'png')
saveas(fig4, fullfile(p.outputDir2, 'myHCL_dendroCol_meanActivityOfNeuronsPerStimulus.png'), 'png')
saveas(fig5, fullfile(p.outputDir2, 'myHCL_rowColumnSorted_meanActivityOfNeuronsPerStimulus.png'), 'png')

saveas(fig2, fullfile(p.outputDir2, 'myHCL_dendro1_meanActivityOfNeuronsPerStimulus.fig'), 'fig')
saveas(fig3, fullfile(p.outputDir2, 'myHCL_rowSorted_meanActivityOfNeuronsPerStimulus.fig'), 'fig')
saveas(fig4, fullfile(p.outputDir2, 'myHCL_dendroCol_meanActivityOfNeuronsPerStimulus.fig'), 'fig')
saveas(fig5, fullfile(p.outputDir2, 'myHCL_rowColumnSorted_meanActivityOfNeuronsPerStimulus.fig'), 'fig')

if par5.closeFigs; close all; end

%% 8. PCA visualization

inputMat = zscore(meanActivityOfNeuronsPerStimulus');
numComp0 = min(min(size(inputMat, 1), size(inputMat, 2)), 6);

if (numComp0 >= 3)
    
    % PCA
    [coeff,score,latent,tsquared,explained,mu] = pca(inputMat, 'NumComponents', numComp0);

    f6 = figure('Visible', par5.figFlag);
    scatter(score(:, 1), score(:, 2), 'filled')
    hold on, text(score(:, 1)+0.5, score(:, 2) , stimuliNames, 'FontSize', 8)

    ax = gca;
    xlim(ax.XLim * 1.2)
    ylim(ax.YLim * 1.2)
    refline([0,0])
    line([0, 0], ax.YLim)
    xlabel('PC1')
    ylabel('PC2')

    f62 = figure('Visible', par5.figFlag);
    plot(explained)
    xlabel('Number of Principal Components')
    ylabel('Variance explained by each component')

    f63 = figure('Visible', par5.figFlag);
    scatter(score(:, 1), score(:, 3), 'filled')
    hold on, text(score(:, 1)+1, score(:, 3) , stimuliNames)
    ax = gca;
    xlim(ax.XLim * 1.2)
    ylim(ax.YLim * 1.2)
    refline([0,0])
    line([0, 0], ax.YLim)
    xlabel('PC1')
    ylabel('PC3')

    f64 = figure('Visible', par5.figFlag);
    scatter(score(:, 2), score(:, 3), 'filled')
    hold on, text(score(:, 2)+1, score(:, 3) , stimuliNames)
    ax = gca;
    xlim(ax.XLim * 1.2)
    ylim(ax.YLim * 1.2)
    refline([0,0])
    line([0, 0], ax.YLim)
    xlabel('PC2')
    ylabel('PC3')

    %% save

    save(fullfile(p.outputDir2, 'PCA_score.mat'), 'score', 'coeff', 'explained')
    saveas(f6, fullfile(p.outputDir2, 'PCA12.png'), 'png')
    saveas(f6, fullfile(p.outputDir2, 'PCA12.fig'), 'fig')

    saveas(f62, fullfile(p.outputDir2, 'PCA_varExplained.png'), 'png')
    saveas(f62, fullfile(p.outputDir2, 'PCA_varExplained.fig'), 'fig')
    saveas(f63, fullfile(p.outputDir2, 'PCA13.png'), 'png')
    saveas(f63, fullfile(p.outputDir2, 'PCA13.fig'), 'fig')
    saveas(f64, fullfile(p.outputDir2, 'PCA23.png'), 'png')
    saveas(f64, fullfile(p.outputDir2, 'PCA23.fig'), 'fig')

    if par5.closeFigs; close all; end
    
end

%% 9. ttest with negative control
%% ttest2 with unequal variance
% at first, activityAtMaximizingFrame were tested. Changed to meanActivityOverHalfInterval

[PvaluesOfNeuronsPerStimulus, BHFDROfNeuronsPerStimulus] = ...
    testResponses_ttest2FDR_BH(activityPerStimulusInfo, par5, p, nCells, numStimuli, stimuliNames);


%% Boxplots

if par5.make_BoxplotsOfTtests

    BoxplotTestingDir = fullfile(p.outputDir, 'Boxplots_t-Testing');
    if ~isfolder(BoxplotTestingDir); mkdir(BoxplotTestingDir); end

    % Boxplot of activityAttheFrame

    tmp = nan(1, numStimuli);
    for s = 1:numStimuli
        tmp(s) = size(activityPerStimulusInfo(s).meanActivityOverHalfInterval, 2);
    end
    maxRepeats = max(tmp);

    for j = 1:3:nCells

        f10 = figure('Position', [10 10 1500 600], 'Visible', par5.figFlag);

        imax = min(j+2, nCells);
        for i = j:imax  
            matOut = nan(maxRepeats, numStimuli);

            for s = 1:numStimuli            
                tmp = activityPerStimulusInfo(s).meanActivityOverHalfInterval(i, :);
                matOut(1:numel(tmp), s) = tmp(:);
            end

            rpos = mod(i, 3);
            if rpos == 0; rpos = 3; end
            subplot(3, 1, rpos) 

            condNames = stimuliNames;
            f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames);
            hold on
            mattmp = 0.05*randn(size(matOut));
            mat2 = mattmp + [1:size(matOut, 2)];
            s = scatter(mat2(:), matOut(:));
            s.LineWidth = 0.6;
            s.MarkerEdgeColor = 'w';
            s.MarkerFaceColor = 'b';
            meanVec = nanmean(matOut, 1);
            s2 = scatter(1:size(matOut,2), meanVec, 72, 'r+');

            % mark significance
            hold on
            ax = gca;
            for s = 1:numStimuli
               if (PvaluesOfNeuronsPerStimulus(i, s) < 0.05)
                   tt = text(s - 0.2, ax.YLim(2) * 0.9, '*', 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'k'); 
               end

               if (BHFDROfNeuronsPerStimulus(i, s) < par5.FalseDiscoveryRate_threshold)
                   tt = text(s + 0.1, ax.YLim(2) * 0.9, '*', 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'r'); 
               end           
            end

            ylabel('\DeltaF/F')
            refline([0, 0])
            title(['ROI ', num2str(i)])

        end

        saveas(f10, fullfile(BoxplotTestingDir, ['BoxplotTesting_', num2str(j), '.png']), 'png')
        pause(0.2)
        close(f10)
    end

end

%% 10. Significance annotated - activity TS per stimulus 

if par5.make_AllResponseCurves_perStimulus

    responseCurves_perStimulus_withPvalueAndFDR = fullfile(p.outputDir, 'responseCurves_perStimulus_withPvalueAndFDR');
    if ~isfolder(responseCurves_perStimulus_withPvalueAndFDR); mkdir(responseCurves_perStimulus_withPvalueAndFDR); end

    % One figure with plots for 3 ROIs
    frameMiddle = round(par5.stimulusInterval / 2);

    for j = 1:3:nCells

        f11 = figure('Position', [10 10 1500 600], 'Visible', par5.figFlag);

        imax = min(j+2, nCells);
        for i = j:imax  

            y = actSigHCL(i, :);
            ymax = max(y)*1.2; ymin = min(min(y), -0.02);

            for s = 1:numStimuli            
                rpos = mod(i, 3);
                if rpos == 0; rpos = 3; end

                subplot(3, numStimuli, (rpos-1)*numStimuli + s)            
                tmp = activityPerStimulusInfo(s).blockactSigHCL(i, :, :);
                fs = squeeze(tmp);   
                plot(fs, 'k')
                refline([0, 0])

                avgTS = activityPerStimulusInfo(s).stimulusMeanMat(i, :);
                hold on
                plot(avgTS, 'r', 'LineWidth', 2)

                xlabel(activityPerStimulusInfo(s).stimuliNames)
                ylim([ymin, ymax])
                title(['ROI ', num2str(i)])

                % P-value annotation
                ax = gca;
                if (PvaluesOfNeuronsPerStimulus(i, s) < 0.05)
                   tt = text(frameMiddle - 2, ax.YLim(2) * 0.9, '*', 'FontSize', 15, ...
                       'FontWeight', 'bold', 'Color', 'k'); 
                end
                if (BHFDROfNeuronsPerStimulus(i, s) < par5.FalseDiscoveryRate_threshold)
                   tt = text(frameMiddle + 2, ax.YLim(2) * 0.9, '*', 'FontSize', 15, ...
                       'FontWeight', 'bold', 'Color', 'r'); 
                end 
            end
        end

        saveas(f11, fullfile(responseCurves_perStimulus_withPvalueAndFDR, ['meanPerStimulus_', num2str(j), '.png']), 'png')
        saveas(f11, fullfile(responseCurves_perStimulus_withPvalueAndFDR, ['meanPerStimulus_', num2str(j), '.fig']), 'fig')
        pause(0.2)
        close(f11)
    end

end

%% 11. Determine responsiveness matrix by FDR < threshold

responseMat = (BHFDROfNeuronsPerStimulus < par5.FalseDiscoveryRate_threshold);

% numStimuli-bit number representation
responseCodes = nan(nCells, 1);
for i = 1:nCells
    resVec = responseMat(i, :);
    numOnes = sum(resVec);
    bases = (numStimuli + 1) .^((numOnes - 1):-1:0);
    stId = find(resVec);
    code0 = sum(stId .* bases);
    responseCodes(i) = code0;
end

resLabels = cell(nCells, 2);
for i = 1:nCells
    resVec = responseMat(i, :);
    resNames = stimuliNames(resVec);
    tmp = strjoin(resNames, " + ");
    if isempty(tmp); tmp = 'None'; end
    resLabels{i, 1} = tmp;
    resLabels{i, 2} = responseCodes(i);
end

%% bar graphs

resOrNot = (responseCodes > 0);
% tabulate resOrNot;
tbl = cell(2, 3);
tbl(:, 1) = {0; 1};
freq1 = sum(resOrNot);
pct1 = freq1 / numel(resOrNot) * 100;
tbl(:, 2) = {numel(resOrNot) - freq1; freq1};
tbl(:, 3) = {100 - pct1; pct1};

table_CountsResponsiveOrNot = table(tbl(:,1), cell2mat(tbl(:,2)), ...
    cell2mat(tbl(:,3)), 'VariableNames', {'Value', 'Count', 'Percent'} );

f12 = figure('Visible', par5.figFlag);
X = categorical({'Non-responsive', 'Responsive'});
X = reordercats(X, {'Non-responsive', 'Responsive'});
bar(X, cell2mat(tbl(:, 2)))
title('Neuron Response to Stimulus')
ylabel('Count')
ax = gca;
ylim(ax.YLim .* 1.1);
hold on
text(1, tbl{1, 2} * 1.1, num2str(tbl{1, 2}), 'FontSize', 12);
text(2, tbl{2, 2} * 1.1, num2str(tbl{2, 2}), 'FontSize', 12);

f13 = figure('Visible', par5.figFlag);
X = categorical({'Non-responsive', 'Responsive'});
X = reordercats(X, {'Non-responsive', 'Responsive'});
bar(X, cell2mat(tbl(:, 3)))
title('Neuron Response to Stimulus')
ylabel('Percent')
ax = gca;
ylim(ax.YLim .* 1.1);
hold on
text(1, tbl{1, 3} * 1.1, num2str(round(tbl{1, 3})), 'FontSize', 12);
text(2, tbl{2, 3} * 1.1, num2str(round(tbl{2, 3})), 'FontSize', 12);

%% counts per-stimulus

resPerStim = sum(responseMat, 1);
table_neuronCountsPerStimulus = table(stimuliNames, resPerStim(:), 'VariableNames', {'Stimulus', 'Count'});

f16 = figure('Visible', par5.figFlag);
X = categorical(stimuliNames);
bar(X, resPerStim)
title('Neuron Response to Stimulus')
ylabel('Count')
xlabel('Stimulus')

%% Stimulus combinations 

indSpecific = find(responseCodes > 0);
responseCodes1 = responseCodes(indSpecific);
resLabel1 = resLabels(indSpecific, 1);
tbl = tabulate(categorical(responseCodes1));
tblValue = tbl(:, 1);

val = cellfun(@(x) str2double(x), tblValue);
tmp = arrayfun(@(x) find(responseCodes1 == x, 1), val);
tblX = resLabel1(tmp);
nn = size(tbl, 1);

% table output
if isempty(tbl)
    tbl = cell.empty(0, 3);
end

table_neuronCounts_responseToMultipleStimuli = table(tbl(:,1), tblX, ...
    cell2mat(tbl(:,2)), ...
    cell2mat(tbl(:,3)), ...
    'VariableNames', {'Code of Stimuli combination', 'MultipleStimuli', 'Count', 'Percent'} );

f14 = figure('Position', [100 100 1070 413], 'Visible', par5.figFlag);
X = categorical(tblX(nn:-1:1));
X = reordercats(X, tblX(nn:-1:1));
b = barh(X, cell2mat(tbl(nn:-1:1, 2)));
title('Detected Stimulus Subsets')
xlabel('Count')
ax = gca;
ax.YAxis.FontSize = 6;

f15 = figure('Position', [100 100 1070 413], 'Visible', par5.figFlag);
X = categorical(tblX(nn:-1:1));
X = reordercats(X, tblX(nn:-1:1));
b = barh(X, cell2mat(tbl(nn:-1:1, 3)));
title('Detected Stimulus Subsets')
xlabel('Percent')
ax = gca;
ax.YAxis.FontSize = 6;

%% save
writetable(table_neuronCounts_responseToMultipleStimuli, ...
    fullfile(p.outputDir, 'table_neuronCounts_responseToMultipleStimuli.csv'))

%% Stimulus combinations plot 2 - dot representation of all detected stimulus subsets

thr = par5.thresholdForMinimalNumberOfNeuronsForEachStimulusCombination;

% combination code -> response vector
stimCombMat = zeros(nn, numStimuli);
for k = 1:nn
    stimComb = val(k);
    btmp = mod(stimComb, numStimuli + 1);
    if (btmp ~= 0); stimCombMat(k, btmp) = 1; end
    while stimComb > btmp
        stimComb = (stimComb - btmp) / (numStimuli + 1);
        btmp = mod(stimComb, numStimuli + 1);
        if (btmp ~= 0); stimCombMat(k, btmp) = 1; end
    end    
end

[krow, kcol] = find(stimCombMat);

f14dot = figure('Position', [100 100 1100 500], 'Visible', par5.figFlag);

subplot('Position', [0.4 0.2, 0.55, 0.7])
b = barh(cell2mat(tbl(nn:-1:1, 2)));
title('Detected Stimulus Subsets')
xlabel('Neuron Count')
yticklabels([])
ax = gca;
plot2ylim = ax.YLim;
plot2FontSize = ax.FontSize;
h = line([thr, thr], plot2ylim); h.LineStyle = '--'; h.Color = 'r';
h.LineWidth = 2;

subplot('Position', [0.1 0.2, 0.28, 0.7])
scatter(kcol, krow, 18, 'k', 'filled')
xlim([0, numStimuli+1])
ylim(plot2ylim)                   
ylabel('Stimulus Subsets')
ax = gca; 
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45; 
grid on;
ax.Box = 'on';
ax.FontSize = plot2FontSize;
axis ij;

%% Stimulus combinations plot 2 (in pct) - dot representation of all detected stimulus subsets

f15dot = figure('Position', [100 100 1100 500], 'Visible', par5.figFlag);

subplot('Position', [0.4 0.2, 0.55, 0.7])
b = barh(cell2mat(tbl(nn:-1:1, 3)));
title('Detected Stimulus Subsets')
xlabel('Neuron Percentage')
yticklabels([])
ax = gca;
plot2ylim = ax.YLim;
plot2FontSize = ax.FontSize;

subplot('Position', [0.1 0.2, 0.28, 0.7])
scatter(kcol, krow,18, 'k', 'filled')
xlim([0, numStimuli+1])
ylim(plot2ylim)                 
ylabel('Stimulus Subsets')
ax = gca;
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;
grid on;
ax.Box = 'on';
ax.FontSize = plot2FontSize;
axis ij;

%% Histogram of activated neuron coutns for each detected stimulus subset

cnt = table_neuronCounts_responseToMultipleStimuli.Count;
cnttbl = tabulate(cnt);

f150 = figure('Position', [100 100 1070 413], 'Visible', par5.figFlag);
if ~isempty(cnttbl); bar(cnttbl(:, 2)); end
xlabel('No. of neurons triggered by one stimulus subset')
ylabel('No. of stimulus subsets')
title('Distribution of activated neuron counts for each stimuls subset')

%% Stimulus subsets that activated enough number of neurons
% After cutting off low counts in table_neuronCounts_responseToMultipleStimuli

thr = par5.thresholdForMinimalNumberOfNeuronsForEachStimulusCombination; 

indHighCounts = find(cell2mat(tbl(:, 2)) >= thr);
tbl2 = tbl(indHighCounts, :);
pct2 = cell2mat(tbl2(:,2))./ sum(cell2mat(tbl2(:,2))) .* 100;

tblValue2 = tbl2(:, 1);
val2 = cellfun(@(x) str2double(x), tblValue2);
tmp2 = arrayfun(@(x) find(responseCodes1 == x, 1), val2);
tblX2 = resLabel1(tmp2);

nn2 = size(tbl2, 1);

% table output
table2_neuronCounts_responseToMultipleStimuli = table(tbl2(:,1), tblX2, ...
    cell2mat(tbl2(:,2)), ...
    cell2mat(tbl2(:,3)), ...
    'VariableNames', {'Code of Stimuli combination', 'MultipleStimuli', 'Count', 'Percent'} );

f14thr = figure('Position', [100 100 1070 413], 'Visible', par5.figFlag);
X = categorical(tblX2(nn2:-1:1));
X = reordercats(X, tblX2(nn2:-1:1));
b = barh(X, cell2mat(tbl2(nn2:-1:1, 2)));
title(['Stimulus subsets that triggered enough # of neurons (>= ', num2str(ceil(thr)), ')'])
xlabel('Count')
ax = gca;
ax.YAxis.FontSize = 10;

f14thrPct = figure('Position', [100 100 1070 413], 'Visible', par5.figFlag);
X = categorical(tblX2(nn2:-1:1));
X = reordercats(X, tblX2(nn2:-1:1));
b = barh(X, pct2(nn2:-1:1));
title(['Stimulus subsets that triggered enough # of neurons'])
xlabel('Percent')
ax = gca;
ax.YAxis.FontSize = 10;

%% Stimulus subsets that activated enough number of neurons - dot representation
% After cutting off low counts in table_neuronCounts_responseToMultipleStimuli

% combination code -> response vector
stimCombMatThr = zeros(nn2, numStimuli);
for k = 1:nn2
    stimComb = val2(k);
    btmp = mod(stimComb, numStimuli + 1);
    if (btmp ~= 0); stimCombMatThr(k, btmp) = 1; end
    while stimComb > btmp
        stimComb = (stimComb - btmp) / (numStimuli + 1);
        btmp = mod(stimComb, numStimuli + 1);
        if (btmp ~= 0); stimCombMatThr(k, btmp) = 1; end
    end    
end

[krow, kcol] = find(stimCombMatThr);

f14thrdot = figure('Position', [100 100 1100 500], 'Visible', par5.figFlag);

subplot('Position', [0.4 0.2, 0.55, 0.7])
b = barh(cell2mat(tbl2(nn2:-1:1, 2)));
title(['Stimulus subsets that triggered enough # of neurons (>= ', num2str(ceil(thr)), ')'])
xlabel('Neuron Count')
yticklabels([])
ax = gca;
plot2ylim = ax.YLim;
plot2FontSize = ax.FontSize;
h = line([thr, thr], plot2ylim); h.LineStyle = '--'; h.Color = 'r';
h.LineWidth = 2;

subplot('Position', [0.1 0.2, 0.28, 0.7])
scatter(kcol, krow, 18, 'k', 'filled')
xlim([0, numStimuli+1])
ylim(plot2ylim)                  
ylabel('Stimulus Subsets')
ax = gca;
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;
grid on;
ax.Box = 'on';
ax.FontSize = plot2FontSize;
axis ij;

%% percentage

f14thrPctdot = figure('Position', [100 100 1100 500], 'Visible', par5.figFlag);

subplot('Position', [0.4 0.2, 0.55, 0.7])
b = barh(pct2(nn2:-1:1));
title(['Stimulus subsets that triggered enough # of neurons'])
xlabel('Neuron Percentage')
yticklabels([])
ax = gca;
plot2ylim = ax.YLim;
plot2FontSize = ax.FontSize;

subplot('Position', [0.1 0.2, 0.28, 0.7])
scatter(kcol, krow, 18, 'k', 'filled')
xlim([0, numStimuli+1])
ylim(plot2ylim)                   
ylabel('Stimulus Subsets')
ax = gca;
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;
grid on;
ax.Box = 'on';
ax.FontSize = plot2FontSize;
axis ij;

%% save table2
writetable(table2_neuronCounts_responseToMultipleStimuli, ...
    fullfile(p.outputDir, 'table2_neuronCounts_responseToMultipleStimuli.csv'))

%% save

saveas(f12, fullfile(p.outputDir2,  'barCnt_CountsResponsiveOrNot.png'), 'png')
saveas(f13, fullfile(p.outputDir2,  'barPercent_CountsResponsiveOrNot.png'), 'png')
saveas(f14, fullfile(p.outputDir,  'barCnt_allDetectedStimulusSubsets.png'), 'png')
saveas(f15, fullfile(p.outputDir,  'barPercent_allDetectedStimulusSubsets.png'), 'png')
saveas(f150, fullfile(p.outputDir,  'Activated_Neuron_Counts_per_category.png'), 'png')
saveas(f16, fullfile(p.outputDir2,  'barCnt_CountsPerStimulus.png'), 'png')
saveas(f14thr, fullfile(p.outputDir,  'barCnt_identifiedStimulusSubsets_AfterThresholding.png'), 'png')
saveas(f14thrPct, fullfile(p.outputDir,  'barPercent_identifiedStimulusSubsets_AfterThresholding.png'), 'png')
saveas(f14dot, fullfile(p.outputDir,  'barCnt_dotsPlot_allDetectedStimulusSubsets.png'), 'png')
saveas(f15dot, fullfile(p.outputDir,  'barPercent_dotPlot_allDetectedStimulusSubsets.png'), 'png')
saveas(f14thrdot, fullfile(p.outputDir,  'barCnt_dotPlot_identifiedStimulusSubsets_AfterThresholding.png'), 'png')
saveas(f14thrPctdot, fullfile(p.outputDir,  'barPercent_dotPlot_identifiedStimulusSubsets_AfterThresholding.png'), 'png')
    
saveas(f12, fullfile(p.outputDir2,  'barCnt_CountsResponsiveOrNot.fig'), 'fig')
saveas(f13, fullfile(p.outputDir2,  'barPercent_CountsResponsiveOrNot.fig'), 'fig')
saveas(f14, fullfile(p.outputDir,  'barCnt_allDetectedStimulusSubsets.fig'), 'fig')
saveas(f15, fullfile(p.outputDir,  'barPercent_allDetectedStimulusSubsets.fig'), 'fig')
saveas(f150, fullfile(p.outputDir,  'Activated_Neuron_Counts_per_category.fig'), 'fig')
saveas(f16, fullfile(p.outputDir2,  'barCnt_CountsPerStimulus.fig'), 'fig')
saveas(f14thr, fullfile(p.outputDir,  'barCnt_identifiedStimulusSubsets_AfterThresholding.fig'), 'fig')
saveas(f14thrPct, fullfile(p.outputDir,  'barPercent_identifiedStimulusSubsets_AfterThresholding.fig'), 'fig')
saveas(f14dot, fullfile(p.outputDir,  'barCnt_dotsPlot_allDetectedStimulusSubsets.fig'), 'fig')
saveas(f15dot, fullfile(p.outputDir,  'barPercent_dotPlot_allDetectedStimulusSubsets.fig'), 'fig')
saveas(f14thrdot, fullfile(p.outputDir,  'barCnt_dotPlot_identifiedStimulusSubsets_AfterThresholding.fig'), 'fig')
saveas(f14thrPctdot, fullfile(p.outputDir,  'barPercent_dotPlot_identifiedStimulusSubsets_AfterThresholding.fig'), 'fig')

if par5.closeFigs; close all; end

%% 12. Heatmaps with neuron-clustered by responsiveness
%% FDR_positiveNeuronMap_clustered

[~, ord1] = sort(responseCodes);

% BHFDROfNeuronsPerStimulus - neuron clustered
posMap = double(BHFDROfNeuronsPerStimulus < par5.FalseDiscoveryRate_threshold);
posMap(isnan(BHFDROfNeuronsPerStimulus)) = NaN;
posMap1 = posMap(ord1, :);

f17 = figure('Visible', par5.figFlag);
figtmp = imagesc(posMap1);
colormap(gray)
ylabel('Neuron Index') 
title('Test positives - one-sided-t-test')

figtmp.AlphaData = 1 - isnan(posMap1) .* 0.5;
ax = gca;
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;

%% save

saveas(f17, fullfile(p.outputDir2, 'FDR_positiveNeuronMap_clustered.png'), 'png')
saveas(f17, fullfile(p.outputDir2, 'FDR_positiveNeuronMap_clustered.fig'), 'fig')

%% activityMap_Annotated

[rc, ord1] = sort(responseCodes);
tmp = rle(rc);

runs = tmp(2:2:numel(tmp));
rcunq = unique(rc);

actSigHCL_tested = [];
cumsumRuns = [0, cumsum(runs)];

for k = 1:numel(rcunq)
    indBlock = (cumsumRuns(k) + 1):cumsumRuns(k+1);
    blockOrd = ord1(indBlock);
    tmp = actSigHCL(blockOrd, :);
    actSigHCL_tested = [actSigHCL_tested; tmp; ...
        nan(par5.numNaNRows_btwnNeuronClusters, size(actSigHCL, 2))];
end

runs2 = runs + par5.numNaNRows_btwnNeuronClusters;
cumsumRuns2 = [0, cumsum(runs2)];
locForLabels = cumsumRuns2(1:end-1) + 1;

% labels
multiStimuliLabel = resLabels(ord1, 1);
uniqueMultiStimuliLabel = multiStimuliLabel(cumsumRuns(1:end-1) + 1);

f18 = figure('Position', [100 100 900 700], 'Visible', par5.figFlag);
q = quantile(actSigHCL_tested(:), [0.005, 0.995]);
figtmp = imagesc(actSigHCL_tested, [q(1), q(2)]);
colorbar
cmap = redwhiteblue(q(1), q(2), 256);
colormap(cmap)

figtmp.AlphaData = 1 - isnan(actSigHCL_tested);
ax = gca;
ax.YTick = locForLabels;
ax.YTickLabel = uniqueMultiStimuliLabel;
ax.YAxis.FontSize = 6;

xlabel('Frame') 
title('Neural Activity (\DeltaF/F)')

%% save

saveas(f18, fullfile(p.outputDir2, 'activityMap_Annotated.png'), 'png')
saveas(f18, fullfile(p.outputDir2, 'activityMap_Annotated.fig'), 'fig')

%% meanActivityMap_Annotated

[rc, ord1] = sort(responseCodes);
tmp = rle(rc);

runs = tmp(2:2:numel(tmp));
rcunq = unique(rc);

meanActivityOfNeuronsPerStimulus_tested = [];
cumsumRuns = [0, cumsum(runs)];

for k = 1:numel(rcunq)
    indBlock = (cumsumRuns(k) + 1):cumsumRuns(k+1);
    blockOrd = ord1(indBlock);
    tmp = meanActivityOfNeuronsPerStimulus(blockOrd, :);
    meanActivityOfNeuronsPerStimulus_tested = [meanActivityOfNeuronsPerStimulus_tested; ...
        tmp; nan(par5.numNaNRows_btwnNeuronClusters, size(meanActivityOfNeuronsPerStimulus, 2))];
end

runs2 = runs + par5.numNaNRows_btwnNeuronClusters;
cumsumRuns2 = [0, cumsum(runs2)];
locForLabels = cumsumRuns2(1:end-1) + 1;

% labels
multiStimuliLabel = resLabels(ord1, 1);
uniqueMultiStimuliLabel = multiStimuliLabel(cumsumRuns(1:end-1) + 1);

f19 = figure('Position', [100 100 900 700], 'Visible', par5.figFlag);
q = quantile(meanActivityOfNeuronsPerStimulus_tested(:), [0.005, 0.995]);
figtmp = imagesc(meanActivityOfNeuronsPerStimulus_tested, [q(1), q(2)]);
colorbar
cmap = redwhiteblue(q(1), q(2), 256);
colormap(cmap)

figtmp.AlphaData = 1 - isnan(meanActivityOfNeuronsPerStimulus_tested);
ax = gca;
ax.YTick = locForLabels;
ax.YTickLabel = uniqueMultiStimuliLabel;
ax.YAxis.FontSize = 6;
 
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;
xlabel('Stimulus') 

title('Mean Activity Per-Stimulus (\DeltaF/F)')

%% save

saveas(f19, fullfile(p.outputDir2, 'meanActivityMap_Annotated.png'), 'png')
saveas(f19, fullfile(p.outputDir2, 'meanActivityMap_Annotated.fig'), 'fig')

%% activityMap_Annotated 2 - categories with low counts are replaced with 'inconclusive'

[rc, ord1] = sort(responseCodes);
[~, ord2] = sort(ord1);
tmp = rle(rc);

vals = tmp(1:2:numel(tmp));
runs = tmp(2:2:numel(tmp));
ind = find(runs < thr);
vals2 = vals; 
vals2(ind) = -1;

tmp2 = tmp;
tmp2(1:2:numel(tmp2)) = vals2;
rc2 = irle(tmp2);
responseCodesReplaced = rc2(ord2)';

resLabelsReplaced = resLabels;
ind1 = find(responseCodesReplaced == -1);
for i = 1:numel(ind1)
    disp(i)
    resLabelsReplaced{ind1(i), 1} = 'Inconclusive profiles';
end

%% plot resorted activities  

[rc, ord1] = sort(responseCodesReplaced);
tmp = rle(rc);

runs = tmp(2:2:numel(tmp));
rcunq = unique(rc);

actSigHCL_tested = [];
cumsumRuns = [0, cumsum(runs)];

for k = 1:numel(rcunq)
    indBlock = (cumsumRuns(k) + 1):cumsumRuns(k+1);
    blockOrd = ord1(indBlock);
    tmp = actSigHCL(blockOrd, :);
    actSigHCL_tested = [actSigHCL_tested; tmp; ...
        nan(par5.numNaNRows_btwnNeuronClusters, size(actSigHCL, 2))];
end

runs2 = runs + par5.numNaNRows_btwnNeuronClusters;
cumsumRuns2 = [0, cumsum(runs2)];
locForLabels = cumsumRuns2(1:end-1) + 1;

% labels
multiStimuliLabel = resLabelsReplaced(ord1, 1);
uniqueMultiStimuliLabel = multiStimuliLabel(cumsumRuns(1:end-1) + 1);

f182 = figure('Position', [100 100 900 700], 'Visible', par5.figFlag);
q = quantile(actSigHCL_tested(:), [0.005, 0.995]);
figtmp = imagesc(actSigHCL_tested, [q(1), q(2)]);
colorbar
cmap = redwhiteblue(q(1), q(2), 256);
colormap(cmap)

figtmp.AlphaData = 1 - isnan(actSigHCL_tested);
ax = gca;
ax.YTick = locForLabels;
ax.YTickLabel = uniqueMultiStimuliLabel;
ax.YAxis.FontSize = 8;

xlabel('Frame') 
title('Neural Activity (\DeltaF/F)')

%% save

saveas(f182, fullfile(p.outputDir, 'activityMap_Annotated_withInconclusive.png'), 'png')
saveas(f182, fullfile(p.outputDir, 'activityMap_Annotated_withInconclusive.fig'), 'fig')

%% meanActivityMap_Annotated 2 with 'inconclusive'

[rc, ord1] = sort(responseCodesReplaced);
tmp = rle(rc);

runs = tmp(2:2:numel(tmp));
rcunq = unique(rc);

meanActivityOfNeuronsPerStimulus_tested = [];
cumsumRuns = [0, cumsum(runs)];

for k = 1:numel(rcunq)
    indBlock = (cumsumRuns(k) + 1):cumsumRuns(k+1);
    blockOrd = ord1(indBlock);
    tmp = meanActivityOfNeuronsPerStimulus(blockOrd, :);
    meanActivityOfNeuronsPerStimulus_tested = [meanActivityOfNeuronsPerStimulus_tested; ...
        tmp; nan(par5.numNaNRows_btwnNeuronClusters, size(meanActivityOfNeuronsPerStimulus, 2))];
end

runs2 = runs + par5.numNaNRows_btwnNeuronClusters;
cumsumRuns2 = [0, cumsum(runs2)];
locForLabels = cumsumRuns2(1:end-1) + 1;

% labels
multiStimuliLabel = resLabelsReplaced(ord1, 1);
uniqueMultiStimuliLabel = multiStimuliLabel(cumsumRuns(1:end-1) + 1);

f192 = figure('Position', [100 100 900 700], 'Visible', par5.figFlag);
q = quantile(meanActivityOfNeuronsPerStimulus_tested(:), [0.005, 0.995]);
figtmp = imagesc(meanActivityOfNeuronsPerStimulus_tested, [q(1), q(2)]);
colorbar
cmap = redwhiteblue(q(1), q(2), 256);
colormap(cmap)

figtmp.AlphaData = 1 - isnan(meanActivityOfNeuronsPerStimulus_tested);
ax = gca;
ax.YTick = locForLabels;
ax.YTickLabel = uniqueMultiStimuliLabel;
ax.YAxis.FontSize = 8;
 
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;
xlabel('Stimulus') 

title('Mean Activity Per-Stimulus (\DeltaF/F)')

%% save

saveas(f192, fullfile(p.outputDir, 'meanActivityMap_Annotated_withInconclusive.png'), 'png')
saveas(f192, fullfile(p.outputDir, 'meanActivityMap_Annotated_withInconclusive.fig'), 'fig')

if par5.closeFigs; close all; end

%% save table data

significantStimuli_perHCLindexedNeurons =  resLabels;
significantStimuli_perHCLindexedNeurons_withInconclusive =  resLabelsReplaced;

save(fullfile(p.outputDir,  'tables_neuronCounts_upon_stimulus.mat'), ...
    'significantStimuli_perHCLindexedNeurons', 'responseCodes', 'table_CountsResponsiveOrNot', ...
    'table_neuronCounts_responseToMultipleStimuli', ...
    'table_neuronCountsPerStimulus', ...
    'table2_neuronCounts_responseToMultipleStimuli', ...
    'significantStimuli_perHCLindexedNeurons_withInconclusive')

%% 13. Make all movies Annotated per ROI

if par5.makeAllFramesOfSingleROIs_firingAnnotated

    allMoviesPerROIDir = fullfile(p.outputDir, 'allFramesOfSingleROIs_firingAnnotated');
    if ~isfolder(allMoviesPerROIDir); mkdir(allMoviesPerROIDir); end

    % load roicenterXmergedsPerm, etc
    load(fullfile(Part2OutputDir, 'par2.mat'))
    load(fullfile(Part4OutputDir, 'roicenterXYZ_HCLPerm.mat'))
    load(fullfile(Part4OutputDir, 'ptMaskVol_ROI_HCLindex.mat'))
    load(fullfile(Part4OutputDir, 'par4.mat'))

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
    upSamplingFactor = par4.upSamplingFactor_forMaskContourImg;
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
            ['hclROI_', num2str(k), '_', significantStimuli_perHCLindexedNeurons{k,1}, ...
            '_centerXYZ_', num2str(mX(k)), '_', ...
            num2str(mY(k)), '_', num2str(mZ(k))]);

        ca3D_maskContoursOnMIP_local_wStimulusSubFrames(fieldCylinder, ZXRatio, ptMaskCylinder, savePath, ...
            upSamplingFactor, roiXtmp, ...
            roiYtmp, ...
            roiZtmp, ...
            stimuliMatrix(:, subFrs), stimuliNames, ...
            subFrs, effNrids, 'figFlag', 'on', 'textFlag', false)

        pause(0.2)
        close(gcf)   
    end
end

%% 14. Make snapshotsOfSingleROIs_atLocalMaximumActivity
% a movie with sub-frames when the activity is above the user-specified threshold
% ROI marked with roicenterXmergedsPerm_kmeans (high, low activities)

if par5.makeSnapshotsOfSingleROIs_atHighActivities
    
    snapshotsPerROIDir = fullfile(p.outputDir, 'snapshotsOfSingleROIs_atHighActivities');
    if ~isfolder(snapshotsPerROIDir); mkdir(snapshotsPerROIDir); end
    
    % load roicenterXmergedsPerm, etc
    load(fullfile(Part2OutputDir, 'par2.mat'))
    load(fullfile(Part4OutputDir, 'roicenterXYZ_HCLPerm.mat'))
    load(fullfile(Part4OutputDir, 'ptMaskVol_ROI_HCLindex.mat'))
    load(fullfile(Part4OutputDir, 'par4.mat'))

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
    upSamplingFactor = par4.upSamplingFactor_forMaskContourImg;
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
        subFrs = find(act > par5.thresholdOfNormalizedActivity); 
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
            ['hclROI_', num2str(k), '_', significantStimuli_perHCLindexedNeurons{k,1}, ...
            '_centerXYZ_', num2str(mX(k)), '_', ...
            num2str(mY(k)), '_', num2str(mZ(k))]);

        ca3D_maskContoursOnMIP_local_wStimulusSubFrames(fieldCylinder, ZXRatio, ptMaskCylinder, savePath, ...
            upSamplingFactor, roiXtmp, ...
            roiYtmp, ...
            roiZtmp, ...
            stimuliMatrix(:, subFrs), stimuliNames, ...
            subFrs, effNrids, 'figFlag', 'on', 'textFlag', false)

        pause(0.2)
        close(gcf)   
    end 
end

%%

disp('== done! ==')
disp('== End of DetermineResponsiveness ==')

end
