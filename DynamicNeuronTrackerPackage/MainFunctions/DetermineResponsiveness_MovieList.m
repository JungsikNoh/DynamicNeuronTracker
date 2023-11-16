function DetermineResponsiveness_MovieList(ML, paramStruct)
% DetermineResponsiveness_MovieList implements statistical analyses to determine 
% the combinatorial responsiveness of each segmented neuron to given
% stimuli based on multiple moviedata. 
% It assumes that each movie data has been already processed by Parts 1-5.
%
% paramStruct contains the following parameters.
%
% activityWeightedByImgCorr (true or false)
%   - Whether to analyze 'weighted' or 'unweighted' calcium activities. 
%   - This option depends on which type was specified in Part 4
%   (VisualizeActivityMapMasksIndexedByHCL).
%    
% stimulusInterval (in # of frames)
%   - A time interval (in frames) between consecutive two different stimuli. 
%
% negCtrlStimulusName
%   - 'stimulusLabel' for a negative control condition. To determine the
%   responsiveness, the statistical test requires a negative control.
%
% make_AllResponseCurves_perStimulus_MovieList (true or false)
%   - Flag whether to generate plots of across-trial response curves and their mean
%   curves per-stimulus for all ROIs and all MD.
%
% make_responsesOfNeurons_activatedByDeterminedStimulusSubsets (true or false)
%   - Flag whether to generate plots of across-trial response curves and
%   their mean curves for the selected neurons which are activated by the 
%   determined stimulus subsets. 
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

par5 = paramStruct;
par5.figFlag = 'on';

if par5.activityWeightedByImgCorr
    Part5OutName = 'DetermineResponsiveness_weighted';
    Part4OutName = 'VisualizeActivityMapMasksIndexedByHCL_weighted';
    p.outputDir2 = fullfile(ML.outputDirectory_, 'DetermineResponsiveness_weighted_MovieList', 'detailedAnalysis');
else
    Part5OutName = 'DetermineResponsiveness_unweighted';
    Part4OutName = 'VisualizeActivityMapMasksIndexedByHCL_unweighted';
    p.outputDir2 = fullfile(ML.outputDirectory_, 'DetermineResponsiveness_unweighted_MovieList', 'detailedAnalysis');
end

p.outputDir = fullfile(ML.outputDirectory_, [Part5OutName, '_MovieList']);

if ~isfolder(p.outputDir); mkdir(p.outputDir); end
if ~isfolder(p.outputDir2); mkdir(p.outputDir2); end
save(fullfile(p.outputDir, 'par5.mat'), 'par5')

%% load activityPerStimulusInfo.mat for each MD

numMDs = numel(ML.movieDataFile_);

MDs = cell(numMDs, 1);
activityPerStimulusInfoCell = cell(numMDs, 1);
actSigHCLCell = cell(numMDs, 1);
for i = 1:numMDs
    tmp = load(fullfile(ML.movieDataFile_{i}));
    MDs{i} = tmp.MD;
    MDoutPath = fullfile(MDs{i}.outputDirectory_, Part5OutName);
    tmp =  load(fullfile(MDoutPath, 'activityPerStimulusInfo.mat'));
    activityPerStimulusInfoCell{i} = tmp.activityPerStimulusInfo;
    
    MDoutPath2 = fullfile(MDs{i}.outputDirectory_, Part4OutName);
    tmp2 =  load(fullfile(MDoutPath2, 'actSig_psSig_HCL-indexed.mat'));
    actSigHCLCell{i} = tmp2.actSigHCL;
end

%% save

save(fullfile(p.outputDir, 'activityPerStimulusInfoCell.mat'), 'activityPerStimulusInfoCell')

%% set up input variables

numStimuli = numel(activityPerStimulusInfoCell{1});
nCellsML = cell(numMDs, 1);

for k = 1:numMDs
    tmp = activityPerStimulusInfoCell{k}(1).activityAtMaximizingFrame;
    nCellsML{k} = size(tmp, 1);
end

totNumCells = sum(cell2mat(nCellsML));   

% get stimuliNames from the 1st MD output
stimuliNames = cell(numStimuli, 1);
for s = 1:numStimuli
    stimuliNames{s} = activityPerStimulusInfoCell{1}(s).stimuliNames;
end

%% meanActivityOfNeuronsPerStimulus

meanActivityOfNeuronsPerStimulus = nan(totNumCells, numStimuli);

for k = 1:numMDs
    tmp = cell2mat(nCellsML);
    rid = sum(tmp(1: (k-1)));
    for s = 1:numStimuli
        for i = 1:nCellsML{k}
            tmp = activityPerStimulusInfoCell{k}(s).meanActivityOverHalfInterval(i, :);
            meanActivityOfNeuronsPerStimulus(rid + i, s) = mean(tmp, 'omitnan');
        end
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

%% myHCL for meanActivityOfNeuronsPerStimulus

inputMat = meanActivityOfNeuronsPerStimulus';

linkageMethod = 'complete';

[~, fig2, fig3, outpermRow, hclmap, fig4, fig5, outpermCol, hclmap2] = ...
    myHCL(inputMat, linkageMethod, ...
    'rowNames', stimuliNames, 'xlabel', 'Neuron Index', 'ylabel', "", ...
    'title', ['Neural Activity - ', linkageMethod], 'figFlag', par5.figFlag);

%% save

save(fullfile(p.outputDir2, 'myHCL_meanActivityOfNeuronsPerStimulus.mat'), ...
                'inputMat', 'outpermRow', 'hclmap', 'outpermCol', 'hclmap2')
            
saveas(fig2, fullfile(p.outputDir2, 'myHCL_dendro1_meanActivityOfNeuronsPerStimulus.png'), 'png')
saveas(fig3, fullfile(p.outputDir2, 'myHCL_rowSorted_meanActivityOfNeuronsPerStimulus.png'), 'png')
saveas(fig4, fullfile(p.outputDir2, 'myHCL_dendroCol_meanActivityOfNeuronsPerStimulus.png'), 'png')
saveas(fig5, fullfile(p.outputDir2, 'myHCL_rowColumnSorted_meanActivityOfNeuronsPerStimulus.png'), 'png')

saveas(fig2, fullfile(p.outputDir2, 'myHCL_dendro1_meanActivityOfNeuronsPerStimulus.fig'), 'fig')
saveas(fig3, fullfile(p.outputDir2, 'myHCL_rowSorted_meanActivityOfNeuronsPerStimulus.fig'), 'fig')
saveas(fig4, fullfile(p.outputDir2, 'myHCL_dendroCol_meanActivityOfNeuronsPerStimulus.fig'), 'fig')
saveas(fig5, fullfile(p.outputDir2, 'myHCL_rowColumnSorted_meanActivityOfNeuronsPerStimulus.fig'), 'fig')

if par5.closeFigs; close all; end

%% 1. PCA visualization

inputMat = zscore(meanActivityOfNeuronsPerStimulus');
numComp0 = min(min(size(inputMat, 1), size(inputMat, 2)), 6);

if (numComp0 >= 3)
    
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

%% 2. ttest with negative control
%% ttest2 with unequal variance
% This result should be the same as the one of the single movie. Run this
% again just for convenience. 

PvaluesOfNeuronsPerStimulus = nan(totNumCells, numStimuli);
% neg control
s0 = find(strcmp(stimuliNames, par5.negCtrlStimulusName));

for k = 1:numMDs
    tmp = cell2mat(nCellsML);
    rid = sum(tmp(1: (k-1)));
    for s = 1:numStimuli
        if (s ~= s0)
            for i = 1:nCellsML{k}
                tmp = activityPerStimulusInfoCell{k}(s).meanActivityOverHalfInterval(i, :);
                ctrltmp = activityPerStimulusInfoCell{k}(s0).meanActivityOverHalfInterval(i, :);
                tmp = tmp(~isnan(tmp));
                ctrltmp = ctrltmp(~isnan(ctrltmp));
                [~, pval, ~, ~] = ttest2(tmp(:), ctrltmp(:), 'Tail', 'right', 'Vartype', 'unequal');
                
                PvaluesOfNeuronsPerStimulus(rid + i, s) = pval;
            end
        end
    end
end

%% PvaluesOfNeuronsPerStimulus

posMap = double(PvaluesOfNeuronsPerStimulus < 0.05);
posMap(isnan(PvaluesOfNeuronsPerStimulus)) = NaN;

f7 = figure('Visible', par5.figFlag);
figtmp = imagesc(posMap);
colormap(gray)
ylabel('Neuron Index') 
title('Test positives - one-sided-t-test')

figtmp.AlphaData = 1 - isnan(posMap) .* 0.5;
ax = gca;
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;

%% p-value histogram

f72 = figure('Visible', par5.figFlag);
histogram(PvaluesOfNeuronsPerStimulus(:), 'Normalization', 'pdf')
ax = gca;
ax.FontSize = 14;
ylabel('Density') 
xlabel('P-values') 

%% csv output of P-values

PvaluesOfNeuronsPerStimulus2 = array2table(PvaluesOfNeuronsPerStimulus, ...
    'VariableNames', stimuliNames);
writetable(PvaluesOfNeuronsPerStimulus2, ...
    fullfile(p.outputDir, 'PvaluesOfNeuronsPerStimulus.csv'))

%% FDR (Benjamini & Hochberg)
% Though Pvalues should be the same, the adjusted P-values will be
% (slightly) different from the separate analysis of each MD.

pvals = reshape(PvaluesOfNeuronsPerStimulus, [], 1);

fFDRdiagplot = figure('Visible', par5.figFlag);
fdrAll = mafdr(pvals,'BHFDR', true, 'Showplot',true);

BHFDROfNeuronsPerStimulus = reshape(fdrAll, totNumCells, numStimuli);

posMap = double(BHFDROfNeuronsPerStimulus < par5.FalseDiscoveryRate_threshold);
posMap(isnan(BHFDROfNeuronsPerStimulus)) = NaN;

f8 = figure('Visible', par5.figFlag);
figtmp = imagesc(posMap);
colormap(gray)
ylabel('Neuron Index') 
title(['Test positives - one-sided-t-test - FDR < ', ...
    num2str(par5.FalseDiscoveryRate_threshold)])

figtmp.AlphaData = 1 - isnan(posMap) .* 0.5;
ax = gca;
ax.XTick = 1:numel(stimuliNames);
ax.XTickLabel = stimuliNames;
ax.XTickLabelRotation = 45;

%% save

save(fullfile(p.outputDir, 'oneSidedTwoSampTTests.mat'), ...
    'PvaluesOfNeuronsPerStimulus', 'BHFDROfNeuronsPerStimulus')

saveas(f7, fullfile(p.outputDir, 'Pvalue_positiveNeuronMap.png'), 'png')
saveas(f72, fullfile(p.outputDir, 'Pvalue_histogram.png'), 'png')
saveas(f8, fullfile(p.outputDir, 'FDR_positiveNeuronMap.png'), 'png')
saveas(fFDRdiagplot, fullfile(p.outputDir, 'BHFDR_pqvalue_diagPlot.png'), 'png')

saveas(f7, fullfile(p.outputDir, 'Pvalue_positiveNeuronMap.fig'), 'fig')
saveas(f72, fullfile(p.outputDir, 'Pvalue_histogram.fig'), 'fig')
saveas(f8, fullfile(p.outputDir, 'FDR_positiveNeuronMap.fig'), 'fig')
saveas(fFDRdiagplot, fullfile(p.outputDir, 'BHFDR_pqvalue_diagPlot.fig'), 'fig')

if par5.closeFigs; close all; end

%% 3. Significance annotated - activity TS per stimulus 
% Because FDR were recalculated, across-trail response curve plots are made again with
% different annotation.

if par5.make_AllResponseCurves_perStimulus_MovieList

    meanTS_perStimulus_withPvalueFDR = fullfile(p.outputDir, 'responseCurves_perStimulus_withPvalueAndFDR_MovieList');
    if ~isfolder(meanTS_perStimulus_withPvalueFDR); mkdir(meanTS_perStimulus_withPvalueFDR); end

    % One figure with plots for 3 ROIs
    frameMiddle = round(par5.stimulusInterval / 2);

    for k = 1:numMDs
        tmp = cell2mat(nCellsML);
        rid = sum(tmp(1: (k-1)));

        for j = 1:3:nCellsML{k}

            f11 = figure('Position', [10 10 1500 600], 'Visible', par5.figFlag);

            imax = min(j+2, nCellsML{k});
            for i = j:imax  

                y = actSigHCLCell{k}(i, :);
                ymax = max(y)*1.2; ymin = min(min(y), -0.02);

                for s = 1:numStimuli            
                    rpos = mod(i, 3);
                    if rpos == 0; rpos = 3; end

                    subplot(3, numStimuli, (rpos-1)*numStimuli + s)            
                    tmp = activityPerStimulusInfoCell{k}(s).blockactSigHCL(i, :, :);
                    fs = squeeze(tmp);   
                    plot(fs, 'k')
                    refline([0, 0])

                    avgTS = activityPerStimulusInfoCell{k}(s).stimulusMeanMat(i, :);
                    hold on
                    plot(avgTS, 'r', 'LineWidth', 2)

                    xlabel(activityPerStimulusInfoCell{k}(s).stimuliNames)
                    ylim([ymin, ymax])
                    title(['ROI ', num2str(i)])

                    % P-value annotation
                    ax = gca;
                    if (PvaluesOfNeuronsPerStimulus(rid + i, s) < 0.05)
                       tt = text(frameMiddle - 2, ax.YLim(2) * 0.9, '*', 'FontSize', 15, ...
                           'FontWeight', 'bold', 'Color', 'k'); 
                    end
                    if (BHFDROfNeuronsPerStimulus(rid + i, s) < par5.FalseDiscoveryRate_threshold)
                       tt = text(frameMiddle + 2, ax.YLim(2) * 0.9, '*', 'FontSize', 15, ...
                           'FontWeight', 'bold', 'Color', 'r'); 
                    end 
                end
            end

            % Output filename shows MD index
            saveas(f11, fullfile(meanTS_perStimulus_withPvalueFDR, ['meanPerStimulus_MD_', num2str(k), ...
                    '_cell_', num2str(j), '.png']), 'png')
            saveas(f11, fullfile(meanTS_perStimulus_withPvalueFDR, ['meanPerStimulus_MD_', num2str(k), ...
                    '_cell_', num2str(j), '.fig']), 'fig')
            pause(0.2)
            close(f11)
        end
    end
    
end

%% 4. Determine responsiveness matrix by FDR < threshold

responseMat = (BHFDROfNeuronsPerStimulus < par5.FalseDiscoveryRate_threshold);

% numStimuli-bit number representation
responseCodes = nan(totNumCells, 1);
for i = 1:totNumCells
    resVec = responseMat(i, :);
    numOnes = sum(resVec);
    bases = (numStimuli + 1) .^((numOnes - 1):-1:0);
    stId = find(resVec);
    code0 = sum(stId .* bases);
    responseCodes(i) = code0;
end

resLabels = cell(totNumCells, 2);
for i = 1:totNumCells
    resVec = responseMat(i, :);
    resNames = stimuliNames(resVec);
    tmp = strjoin(resNames, " + ");
    if isempty(tmp); tmp = 'None'; end
    resLabels{i, 1} = tmp;
    resLabels{i, 2} = responseCodes(i);
end

%% bar graphs

resOrNot = (responseCodes > 0);
tbl = tabulate(resOrNot);
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
scatter(kcol, krow, 3, 'k', 'filled')
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
title('Distribution of neuron counts for each stimuls subset')

%% Stimulus subsets that activated enough number of neurons
% After cutting off low counts in table_neuronCounts_responseToMultipleStimuli

thr = par5.thresholdForMinimalNumberOfNeuronsForEachStimulusCombination;
par5.thresholdForResponsivenessInProportion = thr / totNumCells;

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
h = line([thr, thr], plot2ylim); h.LineStyle = '--'; h.Color = 'k';
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

%% Extract per-movie tables of response profiles to draw stacked bar charts
% Stimulus subsets that activated enough number of neurons - dot representation
% After cutting off low counts in table_neuronCounts_responseToMultipleStimuli

responseCodesPerMD = cell(numMDs, 1);
perMDTable_neuronCounts_responseToMultipleStimuli = cell(numMDs, 1);
perMDTable2_neuronCounts_responseToMultipleStimuli = cell(numMDs, 1);

cumsumId = 0;
for i = 1:numMDs
    responseCodesPerMD{i} = responseCodes((cumsumId+1):(cumsumId+nCellsML{i}));
    cumsumId = cumsumId + nCellsML{i};
end

tableNRow = size(table_neuronCounts_responseToMultipleStimuli, 1);
for i = 1:numMDs
    % drop 'percentage' variable and keep resCode, resLabel and count
    perMDTable_neuronCounts_responseToMultipleStimuli{i} = ...
        table_neuronCounts_responseToMultipleStimuli(:, 1:3);   
    for k = 1:tableNRow
        tmp_resCode = double(string(table_neuronCounts_responseToMultipleStimuli.("Code of Stimuli combination"){k}));
        perMDCnt = sum(responseCodesPerMD{i} == tmp_resCode);
        perMDTable_neuronCounts_responseToMultipleStimuli{i}{k, 3} = perMDCnt;
    end    
end

tableNRow2 = size(table2_neuronCounts_responseToMultipleStimuli, 1);
for i = 1:numMDs
    % drop 'percentage' variable and keep resCode, resLabel and count
    perMDTable2_neuronCounts_responseToMultipleStimuli{i} = ...
        table2_neuronCounts_responseToMultipleStimuli(:, 1:3);   
    for k = 1:tableNRow2
        tmp_resCode = double(string(table2_neuronCounts_responseToMultipleStimuli.("Code of Stimuli combination"){k}));
        perMDCnt = sum(responseCodesPerMD{i} == tmp_resCode);
        perMDTable2_neuronCounts_responseToMultipleStimuli{i}{k, 3} = perMDCnt;
    end    
end

% count matrix, tableNRow2 * numMDs
perMDNeuronCounts = nan(tableNRow2, numMDs);
for i = 1:numMDs
    perMDNeuronCounts(:, i) = perMDTable2_neuronCounts_responseToMultipleStimuli{i}{:, 3};
end
stackedBarLgd = strcat('Movie', num2str(transpose(1:numMDs)));

%% Determined stimulus combinations as a statcked bar chart

% combination code -> response vector
stimCombMatThr = zeros(nn2, numStimuli);    % nn2 = tableNRow2
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

f14thrdotStacked = figure('Position', [100 100 1100 500], 'Visible', par5.figFlag);

subplot('Position', [0.4 0.2, 0.55, 0.7])
b = barh( perMDNeuronCounts(nn2:-1:1, :), 'stacked' );
b(2).FaceColor = [0    0.7410   0.4470];
legend(stackedBarLgd, 'Location', 'southeast');

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
saveas(f14thrdotStacked, fullfile(p.outputDir,  'barCnt_dotPlot_Stacked_identifiedStimulusSubsets_AfterThresholding.png'), 'png')

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
saveas(f14thrdotStacked, fullfile(p.outputDir,  'barCnt_dotPlot_Stacked_identifiedStimulusSubsets_AfterThresholding.fig'), 'fig')

if par5.closeFigs; close all; end    
   
%% 12. Heatmaps with neuron-clustered by responsiveness
%% FDR_positiveNeuronMap_clustered

[~, ord1] = sort(responseCodes);

% neuron clustered
posMap = double(BHFDROfNeuronsPerStimulus < par5.FalseDiscoveryRate_threshold);
posMap(isnan(PvaluesOfNeuronsPerStimulus)) = NaN;
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

saveas(f17, fullfile(p.outputDir, 'FDR_positiveNeuronMap_clustered.png'), 'png')
saveas(f17, fullfile(p.outputDir, 'FDR_positiveNeuronMap_clustered.fig'), 'fig')

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

%% meanActivityMap_Annotated 2 - categories with low counts are replaced with 'inconclusive'

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
    'significantStimuli_perHCLindexedNeurons_withInconclusive', ...
    'responseCodesPerMD', ...
    'perMDTable_neuronCounts_responseToMultipleStimuli', ...
    'perMDTable2_neuronCounts_responseToMultipleStimuli')


%% Collect neurons for identified stimulus-subsets after thresholding

if par5.make_responsesOfNeurons_activatedByDeterminedStimulusSubsets

    frameMiddle = round(par5.stimulusInterval / 2);

    % rc = response code
    [rc, ord1] = sort(responseCodesReplaced);
    resLabelsReplaced_sort = resLabelsReplaced(ord1, 1);

    % Remove 'low counts' and 'none' categories
    ind_conclusive = find(rc > 0);
    rc_conclusive = rc(ind_conclusive);
    ord1_conclusive = ord1(ind_conclusive);
    resLabelsReplaced_sort_conclusive = resLabelsReplaced_sort(ind_conclusive);

    tmp = [];
    if ~isempty(rc_conclusive); tmp = rle(rc_conclusive); end
    runs = tmp(2:2:numel(tmp));
    rcunq = unique(rc_conclusive);
    cumsumRuns = [0, cumsum(runs)];

    % labels
    uniqueMultiStimuliLabel = resLabelsReplaced_sort_conclusive(cumsumRuns(1:end-1) + 1);

    % Make subdirs and draw response curves in each identified stimulus-combination
    nCellsmat = cell2mat(nCellsML);
    id_firstOfEachMD = [1; 1 + cumsum(nCellsmat(1:(end-1)))];

    for k = 1:numel(rcunq)
        % Make a dir for each identified stimulus-combination
        meanTS_perStimulus_Dir = fullfile(p.outputDir, ...
            'IdentifiedStimulusSubsets-responseCurvesOfNeurons', uniqueMultiStimuliLabel{k});
        if ~isfolder(meanTS_perStimulus_Dir); mkdir(meanTS_perStimulus_Dir); end

        % Collect neurons triggered by the combination
        indBlock = (cumsumRuns(k) + 1):cumsumRuns(k+1);
        % Row ids for the neurons in combined dataset
        extendedIDs_k = ord1_conclusive(indBlock);               

        MDIds = arrayfun(@(x) sum(x >= id_firstOfEachMD), extendedIDs_k);
        cellIds = arrayfun(@(x) x - id_firstOfEachMD(sum(x >= id_firstOfEachMD)) + 1, extendedIDs_k);

        j = 1;
        while (j <= numel(extendedIDs_k))

            figts = figure('Position', [10 10 1500 600], 'Visible', par5.figFlag);
            m0 = MDIds(j);
            mvec = m0;
            c0 = cellIds(j);
            cvec = c0;
            
            % check if the next two ROIs belong to the same MD or not
            if (j+1 <= numel(extendedIDs_k)) && (m0 == MDIds(j+1) )
                mvec = [mvec, MDIds(j+1) ];
                cvec = [cvec, cellIds(j+1) ];
                if (j+2 <= numel(extendedIDs_k)) && (m0 == MDIds(j+2) )
                    mvec = [mvec, MDIds(j+2) ];
                    cvec = [cvec, cellIds(j+2) ];
                end
            end
                        
            for i = 1:numel(mvec)
                m = mvec(i);
                c = cvec(i);

                % Activity TS of c-th cell of m-th MD
                y = actSigHCLCell{m}(c, :);
                ymax = max(y)*1.2; ymin = min(min(y), -0.02);

                for s = 1:numStimuli            
                    rpos = i;                   % row-position in subplots

                    subplot(3, numStimuli, (rpos-1)*numStimuli + s)            
                    tmp = activityPerStimulusInfoCell{m}(s).blockactSigHCL(c, :, :);
                    fs = squeeze(tmp);   
                    plot(fs, 'k')
                    refline([0, 0])

                    avgTS = activityPerStimulusInfoCell{m}(s).stimulusMeanMat(c, :);
                    hold on
                    plot(avgTS, 'r', 'LineWidth', 2)

                    xlabel(activityPerStimulusInfoCell{m}(s).stimuliNames)
                    ylim([ymin, ymax])
                    title(['ROI ', num2str(c)])

                    % P-value annotation
                    ax = gca;
                    if (PvaluesOfNeuronsPerStimulus(extendedIDs_k(j + i - 1), s) < 0.05)
                       tt = text(frameMiddle - 2, ax.YLim(2) * 0.9, '*', 'FontSize', 15, ...
                           'FontWeight', 'bold', 'Color', 'k'); 
                    end
                    if (BHFDROfNeuronsPerStimulus(extendedIDs_k(j + i - 1), s) < par5.FalseDiscoveryRate_threshold)
                       tt = text(frameMiddle + 2, ax.YLim(2) * 0.9, '*', 'FontSize', 15, ...
                           'FontWeight', 'bold', 'Color', 'r'); 
                    end 
                end
            end

            % Filename shows MD index and cell index
            saveas(figts, fullfile(meanTS_perStimulus_Dir, ['meanPerStimulus_MD_', num2str(m0), ...
                    '_cell_', num2str(c0), '.png']), 'png')
            saveas(figts, fullfile(meanTS_perStimulus_Dir, ['meanPerStimulus_MD_', num2str(m0), ...
                    '_cell_', num2str(c0), '.fig']), 'fig')
            pause(0.2)
            close(figts)
            
            % Move j to the next position
            j  = j + numel(mvec);
        end
    end
    
end

%% 

disp('== done! ==')
disp('== End of DetermineResponsiveness_MovieList ==')

end
