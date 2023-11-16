function generate_8SubVolumes(movieData, varargin)
% generate_8SubVolumes Generates 8 subvolumes of the original volume for
% the given movieData. The subvolumes overlap at the boundary. 
%
% Usage:    generate_8SubVolumes(MD, 'ratioSubVolumeEdge', 0.65)


%% input

ip = inputParser;
ip.addParameter('ratioSubVolumeEdge', 0.65); 
parse(ip, varargin{:})
p = ip.Results;

MD = movieData;
outDir = MD.outputDirectory_;

%% load mov

imgArray = nan(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_);

parfor fr=1:MD.nFrames_
    currImage = MD.channels_(1).loadStack(fr);
    imgArray(:,:,:,fr) = currImage;
    fprintf(1, '%g ', fr); 
    if (mod(fr,50) == 0); fprintf('\n'); end
end
imgArray = uint16(imgArray);

%%

param_generate_8SubVolumes = struct();
subArray = cell(8,1); 

param_generate_8SubVolumes.ratioSubVolumeEdge = p.ratioSubVolumeEdge;
x0 = MD.imSize_(2);
y0 = MD.imSize_(1);
z0 = MD.zSize_;
x1 = round(MD.imSize_(2) * param_generate_8SubVolumes.ratioSubVolumeEdge);
y1 = round(MD.imSize_(1) * param_generate_8SubVolumes.ratioSubVolumeEdge);
z1 = round(MD.zSize_ * param_generate_8SubVolumes.ratioSubVolumeEdge);

% set x, y, z partitioning
x2 = MD.imSize_(2) - x1 + 1;
y2 = MD.imSize_(1) - y1 + 1;
z2 = MD.zSize_ - z1 + 1;

param_generate_8SubVolumes.x0 = x0; param_generate_8SubVolumes.y0 = y0; param_generate_8SubVolumes.z0 = z0;
param_generate_8SubVolumes.x1 = x1; param_generate_8SubVolumes.y1 = y1; param_generate_8SubVolumes.z1 = z1;

%        --------> x
%       -- 1  2 
%      - -  
%     - 3- 4
%    -   - 5  6
%   -    -
%  -    7  8
% -      -
% y      z 

subArray{1} = imgArray(1:y1, 1:x1, 1:z1, :);
subArray{2} = imgArray(1:y1, x2:x0, 1:z1, :);
subArray{3} = imgArray(y2:y0, 1:x1, 1:z1, :);
subArray{4} = imgArray(y2:y0, x2:x0, 1:z1, :);
subArray{5} = imgArray(1:y1, 1:x1, z2:z0, :);
subArray{6} = imgArray(1:y1, x2:x0, z2:z0, :);
subArray{7} = imgArray(y2:y0, 1:x1, z2:z0, :);
subArray{8} = imgArray(y2:y0, x2:x0, z2:z0, :);

%
save(fullfile(outDir, 'param_generate_8SubVolumes.mat'), 'param_generate_8SubVolumes')

%% write volumes

for i = 1:8

    dname = strcat('subvol', sprintf('%01d', i));

    outputDir = fullfile(outDir, dname, 'rawImages');
    if ~isdir(outputDir); mkdir(outputDir); end

    for fr = 1:MD.nFrames_
        outfname = fullfile(outputDir, ['img', '_', sprintf('%04d', fr), '.tif']);
        imwrite(subArray{i}(:,:,1,fr), outfname);
        for slc = 2:z1        
            imwrite(subArray{i}(:,:,slc,fr), outfname, 'WriteMode', 'append');
        end
        fprintf(1, '%g ', fr);
        if (mod(fr,50) == 0); fprintf('\n'); end
    end

    % make sub-movieData objects
    mdDir = fullfile(pwd, dname);
    channel = Channel(outputDir);

    mdnew = MovieData(channel, mdDir);
    % Set the path where to store the MovieData object.
    mdnew.setPath(mdDir);
    mdnew.setFilename('movieData.mat');
    % Run sanityCheck on MovieData.
    % Check image size and number of frames are consistent.
    % Save the movie if successfull
    mdnew.sanityCheck; 

    % Set some additional movie properties
    mdnew.pixelSize_ = MD.pixelSize_;      % in nm after binning
    mdnew.pixelSizeZ_ = MD.pixelSizeZ_;
    mdnew.timeInterval_= MD.timeInterval_;      % in sec 

    % Save the movieData
    mdnew.sanityCheck()
    mdnew.save;

end

end
