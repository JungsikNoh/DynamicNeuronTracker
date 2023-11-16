%% Main pipeline of DynamicNeuronTracker
%
% DynamicNeuronTracker (DyNT) is a computational tool to segment single neurons
% from 3D calcium imaging videos of deforming tissues.
% For the case when tissue deformation is complicated and the images cannot
% be well registered, DyNT tracks and segments jittering single neurons by
% generating dynamic ROIs. 
%
% This main pipeline consists of four steps:
%   (1) DetectBrightPointSources()
%   (2) TrackJitteringFlickering()
%   (3) SegmentDynamicROIs()
%   (4) VisualizeActivityMapMasksIndexedByHCL()
%
% Author: Jungsik Noh, UTSW, Dallas, TX (2023/08)  


%% Load input into workspace: 
%   (1) MD (movieData object) 
%   (2) TrackingParams_init (a Matlab struct with pre-specified parameters for tracking point sources)


%% Specify Parameters for Part 1 (DetectBrightPointSources)
%
% Parameter description
% detectionAlpha       
%   - In point source (firing neurons) detection, cadidate point sources  
%   are tested to be detected by comparing background intensities and maximal 
%   intensities at the point source. 
%   - 'detectionAlpha' defines the significance threshold in such tests. 
%   - A smaller 'detectionAlpha' means more strict detection criterion leading 
%   to a smaller number of point sources. 
%
% PSFsigma (Point Spread Function Sigma, in pixel)
%   - A cell with three or multiple 2 dim'l vectors, (sigma_XY, sigma_Z), 
%   which specify the standard deviation parameters in X/Y- and Z-direction 
%   of 3D Gaussian functions. The Gaussian functions are fitted to candidate 
%   point sources. 
%   - PSFsigma is proportional to the size or radius of point sources. Larger
%   PSFsigma is optimal for detection of bigger point sources.
%   - Typically three PSFsigma vectors are used to detect neurons with 
%   different sizes. It can take any number of different sigma vectors. 
%
% TopX_ThreshodForBrightness (0 < X < 1)
%   - Among detected point sources, top 100*X% of point sources in brightness 
%   or point source intensity are selected as neuron firing events, 
%   which are fed into the next tracking step.
%
% makeMov_MIP (true or false)
%   - Flag whether to generate two Maximum Intensity Projection (MIP) videos
%   that display the raw 3D images and detected point source on the images.
%
% figFlag ('on' or 'off')
%   - Whether to display output plots, which are always saved in the output
%   directory in .fig and .png formats. 

par1 = struct();
par1.detectionAlpha = {0.001; 0.001; 0.001};
par1.PSFsigma = cell(1, 3);
par1.PSFsigma{1} = [3.5; 1.5];
par1.PSFsigma{2} = [2; 1];
par1.PSFsigma{3} = [1; 1];
par1.TopX_ThreshodForBrightness = 0.2;
par1.makeMov_MIP = true;
par1.figFlag = 'off';

%% DetectBrightPointSources

imgArray = DetectBrightPointSources(MD, par1, TrackingParams_init);


%% Specify Parameters for Part 2 (TrackJitteringFlickering)
%
% Parameter description
% patchSizeX, patchSizeZ (in pixel)
%   - Edge size or length in X/Y- and Z-direction, respectively, 
%   of a local volume around each detected point source, which is utilized
%   to track the jittering and flickering neurons via correlation-based
%   patch-matching.
%
% bandWidthX, bandWidthZ (in pixel)
%   - The radius of the local volume or 3D-patch around each detected point source. 
%   It determines the patchSize via 2*bandWidthX+1 = patchSizeX.
%
% upbdDeformX, upbdDeformZ (in pixel)
%   - Upper bound of local deformation or jittering in X/Y- and Z-direction,
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
%   neurons, then we select only one trajectory that shows a better
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
%   trajectories (or firing events) on MIP images.
%
% closeFigs (true or false)
%   - Whether to close output plots while running because the pipeline
%   generates many figures. Output figures are always saved in the output
%   directory in .fig and .png formats. It simply executes 'close all;'
%   after saving the figures. 

ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;
par2 = struct; 
par2.bandWidthX = 5;
par2.bandWidthZ = ceil(par2.bandWidthX / ZXRatio);
par2.patchSizeX = 1 + 2 * par2.bandWidthX;
par2.patchSizeZ = 1 + 2 * par2.bandWidthZ;
par2.upbdDeformX = 3;
par2.upbdDeformZ = 2;
par2.corrThreshold = [0.9 0.8];
par2.distBtwTracks_toBeMerged = 3;   
par2.imgMarginSizeX = 0; 
par2.imgMarginSizeZ = 0;
par2.minFramesOfFiringEvents = 3;
par2.minFiringFramesOfNeurons = 3;               
par2.makeMov_trackedNeurons = false;
par2.closeFigs = true;

%% TrackJitteringFlickering

TrackJitteringFlickering(MD, imgArray, par2)

%% Specify Parameters for Part 3 (SegmentDynamicROIs)
%
% Parameter description
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

par3 = struct; 
par3.PSFsigma_forROI = par1.PSFsigma{2};
par3.levelOf3DGaussianDist_toSegment = 0.4;
par3.levelOf3DGaussianDist_toSegment_small = 0.2;
par3.bandWidthX = par2.bandWidthX;
par3.bandWidthZ = par2.bandWidthZ;
par3.frameLengthForMovingMedian = 81;
par3.allROIsOutput = true;
par3.makeMov_MIPofROIs = false;
par3.closeFigs = true;

%% SegmentDynamicROIs

SegmentDynamicROIs(MD, imgArray, par3)

%% Specify Parameters for Part 4 (VisualizeActivityMapMasksIndexedByHCL)
%
% Parameter description
% activityWeightedByImgCorr (true or false)
%   - Whether to apply a weighting scheme when collecting calcium
%   activities from the ROIs. If true, spatial correlation-based weights
%   are multiplied to normalized calcium activities to reduce spurious
%   peaks in the signals. 
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

par4 = struct;
par4.activityWeightedByImgCorr = true;
par4.upSamplingFactor_forMaskContourImg = 2;
par4.makeMov_ROIsAtHighActivities = true;
par4.makeMov_ROIsFiringAnnotated = false;
par4.makeMov_singleROIs_AllFrames_firingAnnotated = false;
par4.makeMov_singleROIs_Snapshots_atHighActivities = false;
    % only when par4.makeMov_singleROIs_Snapshots_atHighActivities = true
    % subframes chosen when the activities > threshold below
    par4.thresholdOfNormalizedActivity = 0.02;
par4.closeFigs = true;

%% VisualizeActivityMapMasksIndexedByHCL

VisualizeActivityMapMasksIndexedByHCL(MD, imgArray, par4)

%% End of pipeline

disp("== Pipeline_DyNT() is completed.")
