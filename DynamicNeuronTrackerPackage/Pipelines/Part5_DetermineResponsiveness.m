%% Part5_DetermineResponsiveness
%
% Part5 of DyNT pipeline is a statistical module to determine the combinatorial 
% responsiveness of each segmented neuron to a given set of stimuli.
% It assumes that multiple stimuli are applied randomly and repeatedly.
% Stimulation information needs to be provided by a .csv file in the folder 
% having the movieData object. The .csv file needs to contain three columns 
% with the information for 'stimulusLable', 'startFrame', and 'endFrame'.
% For example, 
%       stimulusLabel	startFrame	endFrame
%        10 uM CDCA      6           10
%        10 uM Q1570     26          30
%        MF 1:300        46          50
%           ...          ...         ...        
%
% Author: Jungsik Noh, UTSW, Dallas, TX   


%% Load input into workspace: 
%   - MD (a movieData object which has been already processed by Parts 1-4) 

%% Specify Parameters for Part 5 (DetermineResponsiveness)
%
% Parameter description
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

par5 = struct();
par5.stimulationFile = ...
    fullfile(MD.outputDirectory_, '2018_02_06_sample1a_BA_VNO_reduced_registered_rigid_stimuliDataset.csv');
par5.stimulusInterval = 20;                     % in frame
par5.negCtrlStimulusName = "Ringers";
par5.activityWeightedByImgCorr = true;
par5.make_AllTSPlots = false;
par5.make_BoxplotsOfTtests = false;
par5.make_AllResponseCurves_perStimulus = false;
par5.makeROIMovieWithStimulusLabel_atHighActivities = true;
par5.makeROIMovieWithStimulusLabel_firingAnnotated = false;
par5.makeAllFramesOfSingleROIs_firingAnnotated = false;
par5.makeSnapshotsOfSingleROIs_atHighActivities = false;
    % only when par5.makeSnapshotsOfSingleROIs_atHighActivities = true
    % subframes chosen when the activities > threshold below
    par5.thresholdOfNormalizedActivity = 0.02;
par5.FalseDiscoveryRate_threshold = 0.10;    
% The following threshold parameter needs to be determined after looking at
% a 'activated neuron count' histogram that will be generated in this
% Part5.
par5.thresholdForMinimalNumberOfNeuronsForEachStimulusCombination = 6;  
par5.numNaNRows_btwnNeuronClusters = 10;
par5.closeFigs = true;

%% DetermineResponsiveness

DetermineResponsiveness(MD, par5)

%% End