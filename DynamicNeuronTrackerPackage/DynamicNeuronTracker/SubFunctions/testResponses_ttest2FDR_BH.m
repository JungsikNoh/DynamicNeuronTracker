function [PvaluesOfNeuronsPerStimulus, BHFDROfNeuronsPerStimulus] = ...
            testResponses_ttest2FDR_BH(activityPerStimulusInfo, par5, p, ...
                                        nCells, numStimuli, stimuliNames)
% testResponses_ttest2FDR_BH implements one-sided t-tests followed by
% BH-FDR controls to compare mean activities under each stimulus and a
% negative control condition.


%% ttest2 unequal 
% at first, activityAtMaximizingFrame were tested. Changed to meanActivityOverHalfInterval

PvaluesOfNeuronsPerStimulus = nan(nCells, numStimuli);
% neg control
s0 = find(strcmp(stimuliNames, par5.negCtrlStimulusName));

for s = 1:numStimuli
    if (s ~= s0)
        for i = 1:nCells
            tmp = activityPerStimulusInfo(s).meanActivityOverHalfInterval(i, :);
            ctrltmp = activityPerStimulusInfo(s0).meanActivityOverHalfInterval(i, :);
            tmp = tmp(~isnan(tmp));
            ctrltmp = ctrltmp(~isnan(ctrltmp));
            [~, pval, ~, ~] = ttest2(tmp(:), ctrltmp, 'Tail', 'right', 'Vartype', 'unequal');
            PvaluesOfNeuronsPerStimulus(i, s) = pval;
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

%% csv output for BUM FDR

PvaluesOfNeuronsPerStimulus2 = array2table(PvaluesOfNeuronsPerStimulus, ...
    'VariableNames', stimuliNames);
writetable(PvaluesOfNeuronsPerStimulus2, ...
    fullfile(p.outputDir, 'PvaluesOfNeuronsPerStimulus.csv'))

%% FDR all columns  (Benjamini & Hochberg)

pvals = reshape(PvaluesOfNeuronsPerStimulus, [], 1);

fFDRdiagplot = figure('Visible', par5.figFlag);
fdrAll = mafdr(pvals,'BHFDR', true, 'Showplot',true);

BHFDROfNeuronsPerStimulus = reshape(fdrAll, nCells, numStimuli);

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


end