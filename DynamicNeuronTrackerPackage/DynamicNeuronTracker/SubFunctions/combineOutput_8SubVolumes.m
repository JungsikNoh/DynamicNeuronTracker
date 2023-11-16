function combineOutput_8SubVolumes(MD, param_generate_8SubVolumes, activityWeightedByImgCorr)
% combineOutput_8SubVolumes integrates sub-volume outputs into the full-volume 
% output, while assessing potentially redundant ROIs within overlapping areas 
% and merging the redundant ROIs into unique ones without requiring an 
% additional parameter
%
% Jungsik Noh, UTSW, 2023/08


%% Set up directories

SVolDirName = {'subvol1', 'subvol2','subvol3','subvol4','subvol5','subvol6', ...
    'subvol7', 'subvol8'};
nV = numel(SVolDirName);

outDirName1 = 'DetectBrightPointSources';
outDirName2 = 'TrackJitteringFlickering';
outDirName3 = 'SegmentDynamicROIs';

if activityWeightedByImgCorr
    outDirName4 = 'VisualizeActivityMapMasksIndexedByHCL_weighted'; 
else
    outDirName4 = 'VisualizeActivityMapMasksIndexedByHCL_unweighted'; 
end

outSubDir1 = fullfile(MD.outputDirectory_, SVolDirName{1}, outDirName1);
outSubDir2 = fullfile(MD.outputDirectory_, SVolDirName{1}, outDirName2);
outSubDir3 = fullfile(MD.outputDirectory_, SVolDirName{1}, outDirName3);
outSubDir4 = fullfile(MD.outputDirectory_, SVolDirName{1}, outDirName4);

% load par1, par2
load(fullfile(outSubDir1, 'par1.mat'));
load(fullfile(outSubDir2, 'par2.mat'));
load(fullfile(outSubDir3, 'par3.mat'));
load(fullfile(outSubDir4, 'par4.mat'));

Part1OutputDir = fullfile(MD.outputDirectory_, outDirName1);
Part2OutputDir = fullfile(MD.outputDirectory_, outDirName2);

if ~isfolder(Part1OutputDir); mkdir(Part1OutputDir); end
if ~isfolder(Part2OutputDir); mkdir(Part2OutputDir); end

%% Volume offsets

x0 = param_generate_8SubVolumes.x0; y0 = param_generate_8SubVolumes.y0; 
z0 = param_generate_8SubVolumes.z0;  % MD.imSize
x1 = param_generate_8SubVolumes.x1; y1 = param_generate_8SubVolumes.y1; 
z1 = param_generate_8SubVolumes.z1; 

x2 = x0 - x1 + 1;
y2 = y0 - y1 + 1;
z2 = z0 - z1 + 1;

ofs = zeros(8, 3);  % offsets in row, col, z
ofs(1,:) = [0, 0, 0];
ofs(2,:) = [0, x2-1, 0];
ofs(3,:) = [y2-1, 0, 0];
ofs(4,:) = [y2-1, x2-1, 0];
ofs(5,:) = [0, 0, z2-1];
ofs(6,:) = [0, x2-1, z2-1];
ofs(7,:) = [y2-1, 0, z2-1];
ofs(8,:) = [y2-1, x2-1, z2-1];

%% Combine part 1 results for 8 sub volumes

% combine bright point source matrix
PSourMat2 = false(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_);  % row, col

for i = 1:nV
    disp(SVolDirName{i})
    subdir = fullfile(MD.outputDirectory_, SVolDirName{i}, outDirName1);
    S = load(fullfile(subdir, 'PSourMat2.mat'));
    subPSourMat = S.PSourMat2;
    
    yind = (1:y1) + ofs(i,1); 
    xind = (1:x1) + ofs(i,2);
    zind = (1:z1) + ofs(i,3);    
    SVol = PSourMat2(yind, xind, zind, :);
    SVol2 = max(SVol, subPSourMat);
    PSourMat2(yind, xind, zind, :) = SVol2;
end

%%  save tracks_Bright3, PSourMat2, imgArray
 
save(fullfile(Part1OutputDir, 'PSourMat2.mat'), 'PSourMat2', '-v7.3')
save(fullfile(Part1OutputDir, 'par1.mat'), 'par1')  

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
    
    % MIP for bright PS
    outname = 'DetectBrightPointSources/MIP_BrightPointSources'; 
    ca3D_plotPSonMIP(MD, PSourMat2, imgArray, outname, 'figFlag', par1.figFlag)
end

%%

disp('== done! ==')
disp('== End of DetectBrightPointSources (Part 1/4) ==')

%% Combine part 2 results. Set up min/max boundaries in xyz.

upbdDeformX = par2.upbdDeformX;
upbdDeformZ = par2.upbdDeformZ;

xup = round((MD.imSize_(2) + 1) / 2) + 2*upbdDeformX; 
xlow = floor((MD.imSize_(2) + 1) / 2) - 2*upbdDeformX; 
yup = round((MD.imSize_(1) + 1) / 2) + 2*upbdDeformX;
ylow = floor((MD.imSize_(1) + 1) / 2) - 2*upbdDeformX;
zup = round((MD.zSize_ + 1) / 2) + 2*upbdDeformZ; 
zlow = floor((MD.zSize_ + 1) / 2) - 2*upbdDeformZ; 

bndmin = zeros(8, 3); 
bndmin(1,:) = [1, 1, 1];
bndmin(2,:) = [1, xlow, 1];
bndmin(3,:) = [ylow, 1, 1];
bndmin(4,:) = [ylow, xlow, 1];
bndmin(5,:) = [1, 1, zlow];
bndmin(6,:) = [1, xlow, zlow];
bndmin(7,:) = [ylow, 1, zlow];
bndmin(8,:) = [ylow, xlow, zlow];

bndmax = zeros(8, 3); 
bndmax(1,:) = [yup, xup, zup];
bndmax(2,:) = [yup,  x0, zup];
bndmax(3,:) = [y0, xup, zup];
bndmax(4,:) = [y0, x0, zup];
bndmax(5,:) = [yup, xup, z0];
bndmax(6,:) = [yup, x0, z0];
bndmax(7,:) = [y0, xup, z0];
bndmax(8,:) = [y0, x0, z0];

%% tracks from subvolumes

disp(SVolDirName{1})
subdir = fullfile(MD.outputDirectory_, SVolDirName{1}, outDirName2);
S = load(fullfile(subdir, 'tracks_PooledEMoutput_MergedOverlaps_L2Merged.mat'));

tracksXcmerged0csL2 = S.tracksXcmerged0csL2;
tracksXcmerged0sL2 = S.tracksXcmerged0sL2;
tracksYcmerged0csL2 = S.tracksYcmerged0csL2;
tracksYcmerged0sL2 = S.tracksYcmerged0sL2;
tracksZcmerged0csL2 = S.tracksZcmerged0csL2;
tracksZcmerged0sL2 = S.tracksZcmerged0sL2;

iterStatus7 = S.iterStatus7;

% boundary check
Xind1 = any((tracksXcmerged0csL2 < bndmin(1, 2)), 2);
Xind2 = any((tracksXcmerged0csL2 > bndmax(1, 2)), 2);
Yind1 = any((tracksYcmerged0csL2 < bndmin(1, 1)), 2);
Yind2 = any((tracksYcmerged0csL2 > bndmax(1, 1)), 2);
Zind1 = any((tracksZcmerged0csL2 < bndmin(1, 3)), 2);
Zind2 = any((tracksZcmerged0csL2 > bndmax(1, 3)), 2);
idout = Xind1 | Xind2 | Yind1 | Yind2 | Zind1 | Zind2;
    
tracksXcmerged0csL2 = tracksXcmerged0csL2(~idout, :);
tracksXcmerged0sL2 = tracksXcmerged0sL2(~idout, :);
tracksYcmerged0csL2 = tracksYcmerged0csL2(~idout, :);
tracksYcmerged0sL2 = tracksYcmerged0sL2(~idout, :);
tracksZcmerged0csL2 = tracksZcmerged0csL2(~idout, :);
tracksZcmerged0sL2 = tracksZcmerged0sL2(~idout, :);    

iterStatus7C = iterStatus7(~idout);

% load subvol output
for i = 2:nV
    disp(SVolDirName{i})
    subdir = fullfile(MD.outputDirectory_, SVolDirName{i}, outDirName2);
    S = load(fullfile(subdir, 'tracks_PooledEMoutput_MergedOverlaps_L2Merged.mat'));
    Xcs = S.tracksXcmerged0csL2;
    Xs = S.tracksXcmerged0sL2;
    Ycs = S.tracksYcmerged0csL2;
    Ys = S.tracksYcmerged0sL2;
    Zcs = S.tracksZcmerged0csL2;
    Zs = S.tracksZcmerged0sL2;
    
    iterStatusTmp = S.iterStatus7;
    
    Xcs = Xcs + ofs(i,2); 
    Xs = Xs + ofs(i,2); 
    Ycs = Ycs + ofs(i,1);
    Ys = Ys + ofs(i,1);
    Zcs = Zcs + ofs(i,3);    
    Zs = Zs + ofs(i,3);    
    
    % boundary check
    Xind1 = any((Xcs < bndmin(i, 2)), 2);
    Xind2 = any((Xcs > bndmax(i, 2)), 2);
    Yind1 = any((Ycs < bndmin(i, 1)), 2);
    Yind2 = any((Ycs > bndmax(i, 1)), 2);
    Zind1 = any((Zcs < bndmin(i, 3)), 2);
    Zind2 = any((Zcs > bndmax(i, 3)), 2);    
    idout = Xind1 | Xind2 | Yind1 | Yind2 | Zind1 | Zind2;
    
    Xcs = Xcs(~idout, :);
    Xs = Xs(~idout, :);
    Ycs = Ycs(~idout, :);
    Ys = Ys(~idout, :);
    Zcs = Zcs(~idout, :);
    Zs = Zs(~idout, :);    
    iterStatusTmp = iterStatusTmp(~idout);
    
    tracksXcmerged0csL2 = [tracksXcmerged0csL2; Xcs];
    tracksXcmerged0sL2 = [tracksXcmerged0sL2; Xs];
    tracksYcmerged0csL2 = [tracksYcmerged0csL2; Ycs];
    tracksYcmerged0sL2 = [tracksYcmerged0sL2; Ys];
    tracksZcmerged0csL2 = [tracksZcmerged0csL2; Zcs];
    tracksZcmerged0sL2 = [tracksZcmerged0sL2; Zs];
    iterStatus7C = [iterStatus7C, iterStatusTmp];
end

%% make tracks unique

tracksXYZ1 = [tracksXcmerged0csL2, tracksYcmerged0csL2, tracksZcmerged0csL2];
 
tracksXYZ1_tmp = tracksXYZ1;
tracksXYZ1_tmp(isnan(tracksXYZ1)) = 0;
[C, ia, ~] = unique(tracksXYZ1_tmp, 'rows', 'stable');

tracksXcmerged0csL2 = tracksXcmerged0csL2(ia,:);
tracksYcmerged0csL2 = tracksYcmerged0csL2(ia,:);
tracksZcmerged0csL2 = tracksZcmerged0csL2(ia,:);
size(tracksXcmerged0csL2) 

iterStatus8 = iterStatus7C(ia);
tracksXcmerged0sL2 = tracksXcmerged0sL2(ia,:);
tracksYcmerged0sL2 = tracksYcmerged0sL2(ia,:);
tracksZcmerged0sL2 = tracksZcmerged0sL2(ia,:);
size(tracksXcmerged0sL2) 

%% Merge too close tracks, based on minL2 distance btn imputed tracks.

distBtwTracks_toBeMerged = par2.distBtwTracks_toBeMerged;

[final_TrCurr3, iterStatus8] = ...
    minL2Merging_diffMeanCor(distBtwTracks_toBeMerged, tracksXcmerged0csL2, tracksYcmerged0csL2, ...
    tracksZcmerged0csL2, iterStatus8);

tracksXcmerged0csL2 = tracksXcmerged0csL2(final_TrCurr3, :);
tracksYcmerged0csL2 = tracksYcmerged0csL2(final_TrCurr3, :);
tracksZcmerged0csL2 = tracksZcmerged0csL2(final_TrCurr3, :);

tracksXcmerged0sL2 = tracksXcmerged0sL2(final_TrCurr3, :);
tracksYcmerged0sL2 = tracksYcmerged0sL2(final_TrCurr3, :);
tracksZcmerged0sL2 = tracksZcmerged0sL2(final_TrCurr3, :);

iterStatus9 = iterStatus8(final_TrCurr3);

%% plot tracksXcmerged0sL2

xccm = nanmean(tracksXcmerged0sL2, 2);
yccm = nanmean(tracksYcmerged0sL2, 2);
zccm = nanmean(tracksZcmerged0sL2, 2);

fig3 = figure('Visible', par2.figFlag); 
scatter3(xccm, yccm, zccm, 15);
ax = gca;
ax.YDir = 'reverse'; ax.ZDir = 'reverse'; 
xlabel('x'); ylabel('y'); zlabel('z')
xlim([1, MD.imSize_(2)]);
ylim([1, MD.imSize_(1)]);
zlim([1, MD.zSize_]);

title(['Total #: ', num2str(numel(xccm))])

%%

saveas(fig3, fullfile(Part2OutputDir, '3dscatter_boundary_L2Merged.fig'), 'fig')

%% 2D scatterplot

fig2 = figure('Visible', par2.figFlag); 
scatter(xccm, yccm, 15);
axis ij
xlim([1, MD.imSize_(2)]); ylim([1, MD.imSize_(1)]); 
xlabel('x'); ylabel('y');

title(['Total #:', num2str(numel(xccm))])

%%
saveas(fig2, fullfile(Part2OutputDir, 'scatterXYmerged0_boundary_L2Merged.fig'), 'fig')
saveas(fig2, fullfile(Part2OutputDir, 'scatterXYmerged0_boundary_L2Merged.png'), 'png')

%% save

save(fullfile(Part2OutputDir, 'tracks_PooledEMoutput_MergedOverlaps_L2Merged.mat'), ...
    'tracksXcmerged0sL2', 'tracksYcmerged0sL2', 'tracksZcmerged0sL2', ...
    'tracksXcmerged0csL2', 'tracksYcmerged0csL2', 'tracksZcmerged0csL2', ...
    'final_TrCurr3', 'iterStatus7C', 'iterStatus8', 'iterStatus9')

save(fullfile(Part2OutputDir, 'par2.mat'), 'par2')  

%% video output  

if par2.makeMov_trackedNeurons
    savePath = fullfile(Part2OutputDir, 'MIP_neuronTracks_whenFiring_tracksXcmerged0sL2');
    ca3D_plotTracksOnMIP(MD, imgArray, tracksXcmerged0sL2, tracksYcmerged0sL2, ...
                  tracksZcmerged0sL2, savePath, 'figFlag', par2.figFlag)
end          

%%  

if par2.makeMov_trackedNeurons
    savePath = fullfile(Part2OutputDir, 'MIP_neuronTracks_imputedByEnds_tracksXcmerged0csL2');
    ca3D_plotTracksOnMIP(MD, imgArray, tracksXcmerged0csL2, tracksYcmerged0csL2, ...
              tracksZcmerged0csL2, savePath, 'figFlag', par2.figFlag)
end  

%%

disp('== done! ==')
disp('== End of TrackJitteringFlickering (Part 2/4) ==')

%% run part 3

SegmentDynamicROIs(MD, imgArray, par3)

%% run part 4

VisualizeActivityMapMasksIndexedByHCL(MD, imgArray, par4)

%%

disp('== done! ==')
disp('== End of combineOutput_8SubVolumes ==')

end
