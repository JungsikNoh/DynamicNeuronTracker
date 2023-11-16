function [imgOut, mask, vI, fvI, bfvI] = segCellLCH_MD(movieDataOrProcess, varargin)% img, varargin)
%segCellLCH_MD segments cell outlines in phase constract images
% simple script to pre-process/segment LCH cells for Deep learning.
% use example: [mask vI fvI bfvI]=segCellLCH(imread('./14-May-2017_atcc_s06_t120_x998_y1586_t130_f6.png'));
%
%  INPUT: image
%          UPDATE....
%             
%  OUTPUT: binary mask
%    (optional outputs include image processing steps for debugging)
%
%   
%
% by Andrew R. Jamieson, Dec 2017
% see also:f segCellLCH.m

%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.addParameter('tagVerify',[], @ischar);
ip.KeepUnmatched = true;
ip.parse(movieDataOrProcess, varargin{:});
p = ip.Results;
paramsIn = p.paramsIn; % extra?

if isa(movieDataOrProcess,'MovieData') && ~isempty(p.tagVerify)
    thisProc = movieDataOrProcess.findProcessTag(p.tagVerify);
    movieData = movieDataOrProcess;
else
    [movieData, thisProc] = getOwnerAndProcess(movieDataOrProcess,'LCHCellSegmentationProcess', true);
end

if ~isempty(p.tagVerify)
    assert(strcmp(p.tagVerify, thisProc.tag_), ['Failed to verify process tag! thisProc:' thisProc.tag_ ' vs. tagVerify:' p.tagVerify])
end

MD = movieData;

pDirCheck = thisProc.funParams_;
pDirCheck.OutputDirectory = [pDirCheck.OutputDirectory filesep thisProc.tag_];
thisProc.setParameters(pDirCheck);

p = parseProcessParams(thisProc, paramsIn); % in case more parameters passed in.

% Sanity Checks
nChan = numel(movieData.channels_);
if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

% Input paths
% Set up the input directories (input images)
inFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex
    if isempty(p.ProcessIndex)
        inFilePaths{1,i} = movieData.getChannelPaths{i};
    else
        inFilePaths{1,i} = movieData.processes_{p.ProcessIndex}.outFilePaths_{1,i}; 
    end
end
thisProc.setInFilePaths(inFilePaths);

% Output paths
dName = 'chan';%String for naming the mask directories for each channel
outFilePaths = cell(2, numel(movieData.channels_));
mkClrDir(p.OutputDirectory);
for iChan = p.ChannelIndex;
    % Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(iChan)];
    outFilePaths{1,iChan} = currDir;
    thisProc.setOutMaskPath(iChan, currDir);
    mkClrDir(outFilePaths{1,iChan});

    outFilePaths{2, iChan} = [p.OutputDirectory filesep 'LCHSegmentMetaData.mat'];
    
end
thisProc.setOutFilePaths(outFilePaths);
% prep for parfor linearization 
pRunSet = {}; cell(MD.nFrames_ * length(p.ChannelIndex));
outDir = {}; 
outDir2 = {}; 

i = 1;
for chIdx = p.ChannelIndex
    imFileName = p.imageFileNamePrefix;
    outDir{chIdx} = thisProc.outFilePaths_{1, chIdx};
%     outDir2{chIdx} = outFilePaths_{2, chIdx};
    for frameIdx = 1:MD.nFrames_
        % frameNum = num2str(frameIdx);
        frameNum = sprintf('%05d', frameIdx);

        pRunSet{i}.chIdx = chIdx; %struct('chIdx', chIdx)
        pRunSet{i}.frameIdx = frameIdx; %struct('chIdx', chIdx)
        imFileNames_pRunSet{chIdx, frameIdx} = ['segLCH_' imFileName '_f' frameNum]; 
        i = i + 1;
    end
end

% for pfi = 1:length(pRunSet)
numDetect_all = {}; %cell(1, MD.nFrames_);
numDetect_allChan = {}; % per channel

% for pfi = 1:length(pRunSet)
parfor pfi = 1:length(pRunSet)

    chIdx = pRunSet{pfi}.chIdx;
    frameIdx = pRunSet{pfi}.frameIdx;

    %% segCellLCH
    % disp('=====');
    disp(['Channel ' , num2str(chIdx), ' Frame: ', num2str(frameIdx)]);    
    im = MD.channels_(chIdx).loadImage(frameIdx);
    [imgSeg, mask, numDetect] = segCellLCH_im(im, 'frameIdx', frameIdx, p);
    numDetect_all{pfi} = numDetect;
    
    % write out the mask
    imwrite(mask, fullfile(outFilePaths{chIdx}, [imFileNames_pRunSet{chIdx, frameIdx} '.tif']));
end


% unmix the numDetect (due to parfor restrictions)
for pfi = 1:length(pRunSet)
    
    chIdx = pRunSet{pfi}.chIdx;
    frameIdx = pRunSet{pfi}.frameIdx;
    numDetect_allChan{frameIdx, chIdx} = numDetect_all{pfi};
    
end

save(outFilePaths{2,1}, 'numDetect_allChan');
   

function [imgOut, mask, numDetect, vI, fvI, bfvI] = segCellLCH_im(img, varargin)
    % modified from function segCellLCH.m
    ip = inputParser;
    ip.addRequired('img', @isnumeric);
    ip.addParameter('closureRadius', 0, @isnumeric); % from use of phasecontrastSeg.m
    ip.addParameter('algorithm', 'basic', @ischar); 
    ip.addParameter('preview', false, @islogical);
    ip.addParameter('leverParams', [], @isstruct);
    ip.addParameter('frameIdx', [], @isnumeric);
    ip.KeepUnmatched = true;
    ip.parse(img, varargin{:});
    p = ip.Results;
    
    frameIdx = p.frameIdx;
    I = mat2gray(img);
    numDetect = 0;
    
    switch p.algorithm
        case 'basic' % ARJ's homebrew hack method
            % Core image processing steps
            %%%%%%%%%%%%%%%%%%%%%%%
        %         I = mat2gray(img);
            gI = imfilter(I, fspecial('gaussian', 5, .25));
            vI = stdfilt(gI);
            fvI = imfilter(vI, fspecial('gaussian', 5, 1));
            bfvI = imbinarize(fvI, .02);
            bfvI = imdilate(bfvI, strel('disk', 5));
            maskAll = imclose(bfvI, strel('disk', 5));
            % maskAll = imerode(maskAll, strel('disk', 4));
            % maskAll = imerode(maskAll, strel('disk', 2));
            maskAll = bwfill(maskAll, 'holes');
            %%%%%%%%%%%%%%%%%%%%%%%

            % See how many objects
            CC = bwconncomp(maskAll);

            % check if in center
            % if yes, keep just this object
            % if multiple, select one that overlaps with center
            if CC.NumObjects > 1
                centerPts = round(size(maskAll)./2);
                mask = bwselect(maskAll,centerPts,centerPts);

                % check if mask present at center.
                if isempty(find(mask,1))
                    % if not, just take the larget
                    mask = bwareafilt(maskAll,1);
                end
            else
                mask = maskAll;
            end
        case 'phasecontrast' % some old method I found on the repo...
            
            mask = phasecontrastSeg(I, p.closureRadius);
        
        case 'lever' % based on 
            % https://git-bioimage.coe.drexel.edu/opensource/leverjs/tree/master
            % note, mask will take the cumulative logical union of
            % all detected cell masks.
            % this will later be optional. (for example to take the largest)

            [mask, Cells] = lever_FrameSegment_texture(I, frameIdx, p);
            numDetect = length(Cells);
        
        otherwise
            error('No segmentation algorithm selected!')
    end

    imgOut = mask.*I;  
    
    if p.preview && strcmp(p.algorithm,'basic')  
        % preview results
        figure; imshow(imgOut);
        figure;
        subplot(3,2,1);imshow(bfvI,[]); title('after thre & dilation');
        subplot(3,2,2);imshow(vI,[]);title('variance filter');
        subplot(3,2,3);imshow(fvI,[]); title('smoothing');
        subplot(3,2,4);imshow(imgOut,[]);title('final');
        subplot(3,2,5);imshow(maskAll,[]);title('all object mask');
        subplot(3,2,6);imshow(I,[]); title('original');
    end