function ca3D_maskContoursOnMIP_local_wSubFrames(imgArray, ZXRatio, ptMaskVol_ROI, savePath, ...
    upSamplingFactor, xx, yy, zz, subFrs, effNrids, varargin)
% ca3D_maskContoursOnMIP_local_wSubFrames generates a MIP video of a local
% volume around a given ROI, while annotating active neighborhood ROIs,
% allowing up-sampling, and allowing sub-frame sampling.


%% input

ip = inputParser;
ip.addParameter('figFlag', 'on');
ip.addParameter('textFlag', true);
ip.parse(varargin{:});
p = ip.Results;

if (numel(size(imgArray))<4)
    warning('This seems to be a 2D movie, No MIP produced.');
    return;
end

if ~isfolder(savePath); mkdir(savePath); end

maxXBorder= size(imgArray, 2);
maxYBorder= size(imgArray, 1);
maxZBorder= ceil(size(imgArray, 3) *ZXRatio);
minXBorder=1;
minYBorder=1;
minZBorder=1;
frameNb= size(imgArray,4);

%%

minIntensityNorm = min(imgArray(:));
maxIntensityNorm = max(imgArray(:));

threeCell = cell(1, frameNb);
strSize0 = 4;
parfor frameIdx=1:frameNb
    vol = imgArray(:,:,:, frameIdx);
    [~,~,~,cthree] = ...
        j_computeMIPs(vol, ZXRatio, minIntensityNorm, maxIntensityNorm, 'stripeSize', strSize0);
    threeCell{frameIdx} = imresize(cthree, upSamplingFactor);
end

% color range
tmp = double(cell2mat(threeCell));
a2 = quantile(tmp(tmp~=0), 0.997); 
a1 = quantile(tmp(tmp~=0), 0.03);

d1 = size(imgArray,1);
d2 = size(imgArray,2);
d3 = size(imgArray,3);
volSiz = [d1 d2 d3]; 

%
fwidth = size(threeCell{1}, 2);
fheight = size(threeCell{1}, 1);

yy = yy .* upSamplingFactor - upSamplingFactor/2;
xx = xx .* upSamplingFactor - upSamplingFactor/2;
zz = zz .* upSamplingFactor - upSamplingFactor/2;
zz = zz .* ZXRatio;

pause(.1)
ftmp = figure('Visible', p.figFlag, 'Position', [51, 51, fwidth, fheight]);
pause(.1)

K = size(xx,1);
colmap0 = autumn(K);
colmap = colmap0(end:-1:1, :);

for frameIdx = 1:frameNb
    
    fprintf(1, '%g ', frameIdx)
    if (mod(frameIdx, 30) == 0); fprintf('\n'); end
    imshow(threeCell{frameIdx}, [a1, a2])
    hold on
        
    roiCell = cell(1, K);    

    for nrid = 1:K
        if ~isnan(xx(nrid, frameIdx))
            maski = (ptMaskVol_ROI(:,:,:,frameIdx) == effNrids(nrid));
            if (sum(maski(:)) == 0)
                continue
            end
            
            [~,~,~,cthree2] = j_computeMIPs(maski,ZXRatio,0,1,'stripeSize',strSize0);
            roiCell{nrid} = imresize(cthree2, upSamplingFactor);
            C = contourc(double(roiCell{nrid}), [.5 .5]);
            ca3D_plotTransparentContourLines(C, 0.4, colmap(nrid,:));
            
            if p.textFlag 
                text(xx(nrid, frameIdx)+5, yy(nrid, frameIdx), num2str(effNrids(nrid)), 'Color', 'y')
                text(xx(nrid, frameIdx)+5, zz(nrid, frameIdx)+(maxYBorder+strSize0)*upSamplingFactor, num2str(effNrids(nrid)), 'Color', 'y')
                text(zz(nrid, frameIdx)+5+(maxXBorder+strSize0)*upSamplingFactor, yy(nrid, frameIdx), num2str(effNrids(nrid)), 'Color', 'y')
            end
        end
    end
    
    tmpframe = getframe(gca);
    tmpim = frame2im(tmpframe);
    oriFr = subFrs(frameIdx);
    imwrite(tmpim, fullfile(savePath, ['three_', sprintf('%04d', oriFr), '.png']), 'png')
    hold off
end

disp('== done ==')

end

