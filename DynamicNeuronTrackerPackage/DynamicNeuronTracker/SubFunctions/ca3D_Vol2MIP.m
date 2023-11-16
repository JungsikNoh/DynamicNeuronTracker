function ca3D_Vol2MIP(MD, imgArray, savePath)
% ca3D_Vol2MIP


%%

ip = inputParser;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.parse(MD);

% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

if(MD.zSize_==1)
    warning('This seems to be a 2D movie, No MIP produced.');
    return;
end

%savePath=[MD.outputDirectory_ filesep outFolderName];
if ~isdir(savePath); mkdir(savePath); end


XYFilesPattern=[savePath filesep 'XY' filesep 'XY_frame_nb%04d.png'];
YZFilesPattern=[savePath filesep 'ZY' filesep 'ZY_frame_nb%04d.png'];
XZFilesPattern=[savePath filesep 'ZX' filesep 'ZX_frame_nb%04d.png'];
ThreeFilesPattern=[savePath filesep 'three' filesep 'Three_frame_nb%04d.png'];

maxXBorder=MD.getDimensions('X');
maxYBorder=MD.getDimensions('Y');
maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
minXBorder=1;
minYBorder=1;
minZBorder=1;
frameNb=MD.nFrames_;

%%

ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;
minIntensityNorm = min(imgArray(:));
maxIntensityNorm = max(imgArray(:));

threeCell = cell(1, MD.nFrames_);
parfor frameIdx=1:MD.nFrames_
    vol = imgArray(:,:,:, frameIdx);
    [cmaxXY,cmaxZY,cmaxZX,cthree]= j_computeMIPs(vol,ZXRatio,minIntensityNorm,maxIntensityNorm);
    threeCell{frameIdx} = cthree;
end

% color range
tmp = double(cell2mat(threeCell));
a1 = min(tmp(:));
a2 = max(tmp(:));

volSiz = [MD.imSize_(1), MD.imSize_(2), MD.zSize_];
        
%
ftmp = figure('Visible', 'off');
for frameIdx = 1:MD.nFrames_
    
    fprintf(1, '%g ', frameIdx)
    if (mod(frameIdx, 30) == 0); fprintf('\n'); end
    
    imwrite(threeCell{frameIdx}, fullfile(savePath, ['three_', sprintf('%04d', frameIdx), '.tif']), 'tif')
    hold off
end

end

