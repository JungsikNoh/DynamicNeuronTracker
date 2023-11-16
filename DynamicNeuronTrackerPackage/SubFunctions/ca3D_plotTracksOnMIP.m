function ca3D_plotTracksOnMIP(MD, imgArray, tracksIterXmergedcs, ...
    tracksIterYmergedcs, tracksIterZmergedcs, savePath, varargin)
% ca3D_plotTracksOnMIP Plot tracked neurons on MIP.


%% Input

ip = inputParser;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('figFlag', 'on');

ip.parse(MD, varargin{:});
p = ip.Results;

% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

if(MD.zSize_==1)
    warning('This seems to be a 2D movie, No MIP produced.');
    return;
end

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
strSize0 = 4;
parfor frameIdx=1:MD.nFrames_
    vol = imgArray(:,:,:, frameIdx);
    [cmaxXY,cmaxZY,cmaxZX,cthree] = ...
        j_computeMIPs(vol,ZXRatio,minIntensityNorm,maxIntensityNorm,'stripeSize',strSize0);
    threeCell{frameIdx} = cthree;
end

% color range
tmp = double(cell2mat(threeCell));
a2 = quantile(tmp(tmp~=0), 0.999);
a1 = quantile(tmp(tmp~=0), 0.03);

L = size(tracksIterXmergedcs, 1);
colmap = jet(L);  
rng(2021)
colmap2 = colmap(randsample(L,L), :);
%
fwidth = size(threeCell{1}, 2);
fheight = size(threeCell{1}, 1);
pause(1)
ftmp = figure('Visible', p.figFlag, 'Position', [251, 251, fwidth, fheight]);
pause(1)
for fr = 1:MD.nFrames_
    
    fprintf(1, '%g ', fr)
    if (mod(fr, 30) == 0); fprintf('\n'); end
    
    imshow(threeCell{fr}, [a1, a2]);
    
    x = tracksIterXmergedcs(:, fr);
    y = tracksIterYmergedcs(:, fr);
    z = tracksIterZmergedcs(:, fr);

    z = z .* ZXRatio - ZXRatio/2;

    hold on
    scatter(x, y, [], colmap2)
    scatter(x, z+maxYBorder+strSize0, [], colmap2)
    scatter(z+maxXBorder+strSize0, y, [], colmap2)

    tmpframe = getframe(ftmp);
    tmpim = frame2im(tmpframe);

    imwrite(tmpim, fullfile(savePath, ['MIPtracksIter_', sprintf('%04d', fr), '.tiff']), 'tiff')

    hold off
    
end

disp('== done ==')

end

