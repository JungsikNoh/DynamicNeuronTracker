function imgArray = DetectBrightPointSources(movieData, paramStruct, ...
                            TrackingParams_init)
% DetectBrightPointSources detects small bright objects, which are firing 
% neurons, and tracking their movements within consecutive time frames.
%
% paramStruct contains the following parameters.
%
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
% 
% Output are saved under MD.outputDirectory_/DetectBrightPointSources.
%
% Jungsik Noh, UTSW, 2023/08


%% Input/Set up outputDir

p.outputDir = fullfile(movieData.outputDirectory_, 'DetectBrightPointSources'); 

MD = movieData;
par1 = paramStruct;

if ~isfolder(p.outputDir); mkdir(p.outputDir); end

%% Load imgArray

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

%%  MIP if necessary

if par1.makeMov_MIP
    ca3D_plotPSonMIP(MD, false(size(imgArray)), imgArray, ...
        'DetectBrightPointSources/MIP_rawImg', 'figFlag', 'on')
end

%% Correlations betn consecutive frames (CC) metric to check overall jittering degree

CCmetric = repmat(NaN, 1, MD.nFrames_);
for k = 2:MD.nFrames_
    x = double(reshape(imgArray(:,:,:,k-1), [], 1));
    y = double(reshape(imgArray(:,:,:,k), [], 1));
    x1 = zscore(x);
    y1 = zscore(y);
    CCmetric(k) = mean(x1.*y1);
end

fig1 = figure('Visible', par1.figFlag);
plot(CCmetric)
xlabel('Frame'); 
title('Correlation between Consecutive time frames (CC metric)')

% save plots
saveas(fig1, fullfile(p.outputDir, 'CCmetric.fig'), 'fig')
saveas(fig1, fullfile(p.outputDir, 'CCmetric.png'), 'png')


%% Step 1: PointSourceDetectionProcess3D, TrackingProcess

K = numel(par1.detectionAlpha);

for i = 1:K
    
    % get PointSourceDetection Process Index
    itmp = MD.getProcessIndex('PointSourceDetectionProcess3D');
    if isempty(itmp)
        PSDproc = PointSourceDetectionProcess3D(MD);
        MD.addProcess(PSDproc);
        itmp = MD.getProcessIndex('PointSourceDetectionProcess3D');
    else
        PSDproc = MD.getProcess(itmp);
    end

    % param
    PSDparams_ = PSDproc.funParams_;
    PSDparams_.OutputDirectory = fullfile(MD.outputDirectory_, 'TrackingPackage/pointsource3D_detect');
    PSDparams_.algorithmType = {'pointSourceAutoSigmaFit'};
    PSDparams_.filterSigma = par1.PSFsigma{i};        % user input
    PSDparams_.alpha = par1.detectionAlpha{i};         % usef input
    PSDparams_.isoCoord = 0;
    
    % Resave the parameters
    parseProcessParams(MD.processes_{itmp}, PSDparams_);
    
    % Run the process 
    MD.processes_{itmp}.run();
    MD.sanityCheck(); 

    % Set up tracking process
    itmp2 = MD.getProcessIndex('TrackingProcess');
    if isempty(itmp2)
        TRproc = TrackingProcess(MD);
        MD.addProcess(TRproc)
        itmp2 = MD.getProcessIndex('TrackingProcess');
        
    end

    % Parse TrackingParams_init
    TRparams_ = TrackingParams_init;
    TRparams_.OutputDirectory = fullfile(MD.outputDirectory_, 'TrackingPackage/tracks');
    TRparams_.DetProcessIndex = itmp;   
    
    % Resave the parameters
    parseProcessParams(MD.processes_{itmp2}, TRparams_);
    
    % Run the process
    MD.processes_{itmp2}.run();
    MD.sanityCheck();

    fdname = fullfile(p.outputDir, ...
        ['TrackingPackage_', num2str(par1.PSFsigma{i}(1)), '_', ...
        num2str(par1.PSFsigma{i}(2)), '_', num2str(par1.detectionAlpha{i})]);
    copyfile(fullfile(MD.outputDirectory_, 'TrackingPackage'), fdname);  

end

%% Step 2: Select top X% bright point sources and tracks

tracks_Bright = cell(K,1);
 
trackDir = cell(K,1);
for k = 1:K
    trackDir{k} = fullfile(p.outputDir, ...
        ['TrackingPackage_', num2str(par1.PSFsigma{k}(1)), '_', ...
        num2str(par1.PSFsigma{k}(2)), '_', num2str(par1.detectionAlpha{k})]);
end

% movieInfo to PSourMat
PSourMat = false(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_);  % row, col

% set up PSourMatIntensity
PSourMatIntensity = zeros(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_);  % row, col
volSiz = [MD.imSize_(1), MD.imSize_(2), MD.zSize_];
    
for detectInd = 1:K
    
    % load movieInfo
    S = load(fullfile(trackDir{detectInd}, 'pointsource3D_detect', 'channel_1.mat'));
    movieInfo = S.movieInfo;
    
    for fr = 1:MD.nFrames_       
        if ~isempty(movieInfo(fr).xCoord)
            x = movieInfo(fr).xCoord(:, 1);
            y = movieInfo(fr).yCoord(:, 1);
            z = movieInfo(fr).zCoord(:, 1);            
            x = round(x); y = round(y);
            z = round(z);

            for i=1:numel(x)
                if (mod(x(i), 1) == 0) && (mod(y(i), 1) == 0) && (mod(z(i), 1) == 0)
                    PSourMat(y(i), x(i), z(i), fr) = 1;
                else
                    disp('== Non-integer x, y, z is found ==')
                    break
                end
            end
        end
    end

    for fr = 1:MD.nFrames_
        lind = find(PSourMat(:,:,:,fr) == 1);
        [y, x, z] = ind2sub(volSiz, lind);

        for i = 1:numel(y)
            ix = (x(i)>1) & (x(i)<MD.imSize_(2) );
            iy = (y(i)>1) & (y(i)<MD.imSize_(1) );
            iz = (z(i)>1) & (z(i)<MD.zSize_ );
            if ix && iy && iz     
                xwin = x(i)-1:x(i)+1;       % simple 3*3*3 cube mean intensity
                ywin = y(i)-1:y(i)+1;
                zwin = z(i)-1:z(i)+1;
                tmp = imgArray(ywin, xwin, zwin, fr);
                m0 = mean(tmp(:));
                PSourMatIntensity(y(i), x(i), z(i), fr) = m0;
            end
        end
    end
end
    
tmp = PSourMatIntensity(PSourMatIntensity > 0);

% Select 
thr_Bright = quantile(tmp, 1 - par1.TopX_ThreshodForBrightness);
disp(['== Threshold (intensity) for bright point sources (top ', ...
    num2str(par1.TopX_ThreshodForBrightness*100), '%)'])
disp(thr_Bright)

PSourMat2 = PSourMat;
PSourMat2(PSourMatIntensity < thr_Bright) = 0;
disp('== Total number of bright point sources:')
disp(sum(PSourMat2(:)))

% tracks only having bright PS
for detectInd = 1:K
    S = load(fullfile(trackDir{detectInd}, 'tracks', 'Channel_1_tracking_result.mat'));
    tracks = TracksHandle(S.tracksFinal);
    
    tracks_Bright{detectInd} = struct('x', []);
    
    k = 1;
       
    for i = 1:numel(tracks)
        tri = tracks(i);
        x = round(tri.x);
        y = round(tri.y);
        z = round(tri.z);
        f = round(tri.f);           
        ind0 = zeros(1, numel(x));
        
        for j = 1:numel(x)
            if ~isnan(y(j))
                ind0(j) = PSourMat2(y(j), x(j), z(j), f(j));
            end
        end

        if (mean(ind0) > 0.9)      % 90% of track PS is bright, then keep the track
            tracks_Bright{detectInd}(k).x = x;
            tracks_Bright{detectInd}(k).y = y;
            tracks_Bright{detectInd}(k).z = z;
            tracks_Bright{detectInd}(k).f = f;
            k = k + 1;
        end
    end
end

disp('== Number of tracks with bright point sources')
for detectInd = 1:K
    disp(['== PSFsigma type = ', num2str(detectInd)])
    disp(numel(tracks_Bright{detectInd}))
end

%% MIP for bright PS
%  to reduce memory burden
clear PSourMatIntensity
clear PSourMat 

if par1.makeMov_MIP
    outname = 'DetectBrightPointSources/MIP_BrightPointSources'; 
    ca3D_plotPSonMIP(MD, PSourMat2, imgArray, outname, 'figFlag', par1.figFlag)
end
    
%% Step 3: Pooling unique tracks_Bright

tracks_Bright2 = [];
if ~isempty(tracks_Bright{1}(1).x)
    tracks_Bright2 = tracks_Bright{1};
end
    
for k = 2:K
    if ~isempty(tracks_Bright{k}(1).x)
        tracks_Bright2 = [tracks_Bright2, tracks_Bright{k}];
    end
end

L = numel(tracks_Bright2);       
frmax = MD.nFrames_;
tracksX1 = nan(L, frmax);
tracksY1 = nan(L, frmax);
tracksZ1 = nan(L, frmax);
 
% here, tracks are double, and tracksXYZ is integer
for k = 1:L
    tri = tracks_Bright2(k);
    x = round(tri.x);
    y = round(tri.y);
    z = round(tri.z);
    f = tri.f;
    
    tracksX1(k,f) = reshape(x, 1, []);
    tracksY1(k,f) = reshape(y, 1, []);
    tracksZ1(k,f) = reshape(z, 1, []);
end 

tracksXYZ1 = [tracksX1, tracksY1, tracksZ1]; 
tracksXYZ1_tmp = tracksXYZ1;
tracksXYZ1_tmp(isnan(tracksXYZ1)) = 0;
[C, ia, ~] = unique(tracksXYZ1_tmp, 'rows', 'stable');
disp('== Number of unique bright tracks')
disp(size(C, 1))
tracks_Bright3 = tracks_Bright2(ia);
tracksX1_Br3 = tracksX1(ia, :);
tracksY1_Br3 = tracksY1(ia, :);
tracksZ1_Br3 = tracksZ1(ia, :);

%%  save tracks_Bright3, PSourMat2, imgArray

save(fullfile(p.outputDir, 'tracks_Bright3.mat'), 'tracks_Bright3')
save(fullfile(p.outputDir, 'tracksX1_Br3.mat'), 'tracksX1_Br3')
save(fullfile(p.outputDir, 'tracksY1_Br3.mat'), 'tracksY1_Br3')
save(fullfile(p.outputDir, 'tracksZ1_Br3.mat'), 'tracksZ1_Br3')

save(fullfile(p.outputDir, 'PSourMat2.mat'), 'PSourMat2', '-v7.3')
save(fullfile(p.outputDir, 'par1.mat'), 'par1') 
save(fullfile(p.outputDir, 'TrackingParams_init.mat'), 'TrackingParams_init') 

%%

disp('== done! ==')
disp('== End of DetectBrightPointSources (Part 1/4) ==')

end
