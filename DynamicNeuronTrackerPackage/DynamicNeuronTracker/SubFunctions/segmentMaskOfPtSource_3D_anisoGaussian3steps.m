function [currMask, convId, fig] = segmentMaskOfPtSource_3D_anisoGaussian3steps(xvec, yvec, ...
    zvec, PSFsig, bandWidthX, bandWidthZ, imgArray, PSourMat2, varargin)
% segmentMaskOfPtSource_3D_anisoGaussian3steps fits a 3D-Gaussian function
% with anisotropic sigma values to segment a central area of local volume
% image of a point source. It reports various fitting cases via convId.
%
% Updates:
% J Noh. psArea on track => +-1 area due to using bright PSs
% instead of corPSourMat.


%% Input

ip = inputParser;
ip.addParameter('ConfRegionLevel', 0.6);
ip.addParameter('figVisible', 'off');
ip.addParameter('figOut', false);

ip.parse(varargin{:});
p = ip.Results;

fig = [];

%% get avgPatch image

ptInd = nan(numel(xvec), 1);
xvec2 = xvec;
yvec2 = yvec;
zvec2 = zvec;

for fr = 1:numel(xvec)
    if ~isnan(yvec(fr)) 
        psArea = PSourMat2(yvec(fr)-1:yvec(fr)+1, xvec(fr)-1:xvec(fr)+1,zvec(fr)-1:zvec(fr)+1, fr);
        ptInd(fr) = sum(psArea(:));
    end
end

indIntersect = (ptInd == 1);
xs = xvec2(indIntersect);
ys = yvec2(indIntersect);
zs = zvec2(indIntersect);
frs = find(indIntersect);

[avgPatch, ~] = localPatchSampling3D(xs,ys,zs,frs, bandWidthX, bandWidthZ, imgArray);

currImage = avgPatch;
if (sum(indIntersect) == 0)
    convId = NaN;    
    currMask = [];
    return; % error('No ps')
end
   
%% mild center adjustment based on imgaussfilt3

filtCurrImg = imgaussfilt3(currImage); 

[~, i] = max(filtCurrImg(:));
[my, mx, mz] = ind2sub(size(currImage), i);

if (abs(bandWidthX+1-mx) > ceil(PSFsig(1))) ...
        || (abs(bandWidthX+1-my) > ceil(PSFsig(1))) ...
        || (abs(bandWidthZ+1-mz) > ceil(PSFsig(2)))
    xin = bandWidthX;
    yin = bandWidthX;
    zin = bandWidthZ;
else
    xin = mx - 1;     % fitGaussian3D xyz start with [0 0 0]
    yin = my - 1;
    zin = mz - 1;
end

%% fitGaussian3D, center=0,0,0 fixed

convId = 0;

currImage2 = currImage;
%[prmVect0, ~, ~, res0, ~] = fitAnisoGaussian2D(currImage, ...
%    [0 0 max(currImage(:)) 1.5 1.5 pi/6 min(currImage(:))], 'arstc');

amp_init = currImage(xin+1, yin+1, zin+1) - min(currImage(:));
[prmVect0, ~, ~, res0, ~] = fitGaussian3D(currImage, [xin yin ...
    zin amp_init PSFsig(1) PSFsig(2)  min(currImage(:))], 'asrc');

if (prmVect0(5) > 2*PSFsig(1)) || (prmVect0(6) > 2*PSFsig(2))
    disp('=== one of sigma0 exceeds upperBound. ===')
    sxinit = min(prmVect0(5), 2*PSFsig(1));
    szinit = min(prmVect0(6), 2*PSFsig(2));
    [prmVect0, ~, ~, res0, ~] = fitGaussian3D(currImage, ...
    [xin yin zin amp_init sxinit szinit min(currImage(:))], 'ac');
end

% outliers mainly from another ptSource
outlierInd0 = (abs(res0.data) > 3*res0.std);
if any(outlierInd0(:))
    %prmVectAdj = [0,0,0, prmVect0(4:7)];
    prmVectAdj = [xin - bandWidthX, yin - bandWidthX, zin - bandWidthZ, prmVect0(4:7)];
    F0 = ftnGaussian3D(prmVectAdj, -bandWidthX:bandWidthX, -bandWidthX:bandWidthX, ...
    -bandWidthZ:bandWidthZ);
    currImage2(outlierInd0) = F0(outlierInd0);
    convId = 10;
end

% step2
amp_init = currImage2(xin+1, yin+1, zin+1) - min(currImage2(:));
[prmVect, ~, ~, res1, ~] = fitGaussian3D(currImage2, ...
    [xin yin zin amp_init PSFsig(1) PSFsig(2)  min(currImage2(:))], ...
    'asrc');

outlierInd1 = (abs(res1.data) > 3*res1.std);
if any(outlierInd1(:))
    prmVectAdj = [xin - bandWidthX, yin - bandWidthX, zin - bandWidthZ, prmVect(4:7)];
    F1 = ftnGaussian3D(prmVectAdj, -bandWidthX:bandWidthX, -bandWidthX:bandWidthX, ...
    -bandWidthZ:bandWidthZ);

    currImage2(outlierInd1) = F1(outlierInd1);
    convId = convId + 1;

    % step3
    [prmVect, ~, ~, ~, ~] = fitGaussian3D(currImage2, ...
        [xin yin zin amp_init PSFsig(1) PSFsig(2)  min(currImage2(:))], ...
        'asrc');
end

% no outlier in all steps: convId = 0
%convId

%%
if (prmVect(5) > 2*PSFsig(1)) || (prmVect(6) > 2*PSFsig(2))
    convId = convId + 10;
    disp('=== one of sigmas exceeds upperBound. ===')
    sxinit = min(prmVect(5), 2*PSFsig(1));
    szinit = min(prmVect(6), 2*PSFsig(2));
    [prmVect, ~, ~, res2, ~] = fitGaussian3D(currImage2, ...
    [xin yin zin amp_init sxinit szinit min(currImage2(:))], 'ac');
    
    outlierInd2 = (abs(res2.data) > 3*res2.std);
    if any(outlierInd2(:))
        prmVectAdj = [xin - bandWidthX, yin - bandWidthX, zin - bandWidthZ, prmVect(4:7)];
        F2 = ftnGaussian3D(prmVectAdj, -bandWidthX:bandWidthX, -bandWidthX:bandWidthX, ...
        -bandWidthZ:bandWidthZ);

        currImage2(outlierInd2) = F2(outlierInd2);
        convId = convId + 1;
    end

    % step3
    amp_init = currImage2(xin+1, yin+1, zin+1) - min(currImage2(:));
    [prmVect, ~, ~, ~, ~] = fitGaussian3D(currImage2, ...
        [xin yin zin amp_init PSFsig(1) PSFsig(2)  min(currImage2(:))], ...
        'asrc');

    if (prmVect(5) > 2*PSFsig(1)) || (prmVect(6) > 2*PSFsig(2))
        convId = convId + 1;
        sxinit2 = min(prmVect(5), 2*PSFsig(1));
        szinit2 = min(prmVect(6), 2*PSFsig(2));
        [prmVect, ~, ~, ~, ~] = fitGaussian3D(currImage2, ...
        [xin yin zin amp_init sxinit2 szinit2 min(currImage2(:))], 'ac');        
    end
end

%% chisquare(2, 1-alpha)

prmVectAdj = [xin - bandWidthX, yin - bandWidthX, zin - bandWidthZ, prmVect(4:7)];
F3 = ftnGaussian3D(prmVectAdj, -bandWidthX:bandWidthX, -bandWidthX:bandWidthX, ...
-bandWidthZ:bandWidthZ);

% p.ConfRegionAlpha = 0.2 (default)
dimDF = 3;
amp = prmVect(4);
F3_alpha = amp * exp((-1)/2 * gaminv(p.ConfRegionLevel, dimDF/2, 2)) + prmVect(7);

% mask
currMask = (F3 > F3_alpha);

%% figure

if p.figOut

    vol = avgPatch;
    vol(~currMask) = min(avgPatch(:));
    outMIPmat = simpleVol2MIP(vol);

    %
    outMIPmat0 = simpleVol2MIP(avgPatch);
    outMIPmat0(isnan(outMIPmat0)) = min(avgPatch(:));

    C2 = imfuse(outMIPmat0, outMIPmat, 'montage');
    fig = figure('Visible', p.figVisible);
    imtmp = imagesc(C2);
    
    colormap(gray)
    title1 = ['Num frs with PS: ', num2str(sum(indIntersect)), ...
        ', sigx:', num2str(round(prmVect(5),1)), ...
        ', sigz:', num2str(round(prmVect(6),1))];
    title2 = ['ConvID: ', num2str(convId), ',', ' center loc in xyz: ', num2str(xin+1),' ', num2str(yin+1), ' ', num2str(zin+1)];
    title({title1; title2})
    imtmp.AlphaData = 1 - isnan(C2);

    colormap(gray); colorbar

end

end

