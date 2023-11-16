%% Part5_DetermineResponsiveness_forMovieList
%
% This pipeline is to determine combinatorial responsiveness of segmented
% neurons to given stimuli based on multiple calcium videos.
% It assumes that each movie data has been already processed by Parts 1-5.
% The multiple movie data need to be registered into a movieList object.%
%
% Author: Jungsik Noh, UTSW, Dallas, TX  

%% Load input into workspace: 
%   - ML (a movieList object that contains multiple movieData objects. Each
%   movieData represents one calcium video taken while applying the same
%   set of stimuli.)

%% Specify Parameters for DetermineResponsiveness_MovieList()
%
% Parameter description
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

par5 = struct();
par5.activityWeightedByImgCorr = true;    
par5.stimulusInterval = 20;                     % in frame
par5.negCtrlStimulusName = "Ringers";
par5.make_AllResponseCurves_perStimulus_MovieList = false;
par5.make_responsesOfNeurons_activatedByDeterminedStimulusSubsets = true;
par5.FalseDiscoveryRate_threshold = 0.10;    
par5.thresholdForMinimalNumberOfNeuronsForEachStimulusCombination = 6;
par5.numNaNRows_btwnNeuronClusters = 10;        
par5.closeFigs = true;

%% DetermineResponsiveness_MovieList

DetermineResponsiveness_MovieList(ML, par5)

%% EOF