function [ROIPairsOverlappingMoreThan10Frames, ROIsOverlappingMoreThan10Frames] = ...
    checkOverlapsOfDynamicROIs(MD, ROIcell, tracksXcmerged0csL2, tracksYcmerged0csL2, ...
                                tracksZcmerged0csL2, bandWidthX, bandWidthZ)
% checkOverlapsOfDynamicROIs


%% 

ptMaskMark = zeros(MD.imSize_(1), MD.imSize_(2), MD.zSize_, MD.nFrames_);  % row, col
overlapStat = [];

frmax = size(tracksXcmerged0csL2, 2);
trmax = size(tracksXcmerged0csL2, 1);
for k = 1:trmax
    xvec = tracksXcmerged0csL2(k, :);
    yvec = tracksYcmerged0csL2(k, :);
    zvec = tracksZcmerged0csL2(k, :);
    frs = 1:frmax;

    [~, markPatches] = localPatchSampling3D(xvec, yvec, zvec, frs, bandWidthX, ...
        bandWidthZ, ptMaskMark);

    mask0 = double(ROIcell{k});
    mask0(mask0 == 0) = nan;
    maskCube = repmat(mask0, 1, 1, 1, frmax); 
    
    % check overlapping
    tmp = maskCube .* markPatches;
    for fr = 1:frmax
        marks = tmp(:,:,:,fr);
        marks1 = marks(~isnan(marks));
        marks2 = marks1(marks1 ~= 0);
        if ~isempty(marks2)
            marks2u = unique(marks2);
            for i=1:numel(marks2u)
                overlapStat0 = [k, marks2u(i), fr];
                overlapStat = [overlapStat; overlapStat0];
            disp(['== ', num2str(k), ' ROI overlaps with ', num2str(marks2u(i)), ...
            ' at frame of ', num2str(fr), ' ==']) 
            end
        end
        
        % mark
        mask2 = k * double(ROIcell{k});
        xwin = xvec(fr)-bandWidthX:xvec(fr)+bandWidthX;
        ywin = yvec(fr)-bandWidthX:yvec(fr)+bandWidthX;        
        zwin = zvec(fr)-bandWidthZ:zvec(fr)+bandWidthZ;       
        maskMarktmp = ptMaskMark(ywin, xwin, zwin, fr);
        ptMaskMark(ywin, xwin, zwin, fr) = max(maskMarktmp, mask2);        
    end
end
%% 

fprintf('\n')
ROIPairsOverlappingMoreThan10Frames = [];
ROIsOverlappingMoreThan10Frames = [];
if isempty(overlapStat)
    disp('== No dynamic ROI overlaps with each other. ==')
else
    %numel(unique(overlapStat(:,1:2)))   % # of neurons
    [unq_overlapStat, ia, ic] = unique(overlapStat(:, 1:2), 'stable', 'rows');
    %tabulate(ic)                        % # of pairs of neurons
    tbl = tabulate(ic);

    % find neuron pairs overlapping in >=10 frames
    tbl10 = tbl((tbl(:,2) >=10), 1:2);

    disp('== Pairs of dynamic ROIs overlapping in >10= frames ==')
    ROIPairsOverlappingMoreThan10Frames = unq_overlapStat(tbl10(:,1), :);
    disp(ROIPairsOverlappingMoreThan10Frames)
    ROIsOverlappingMoreThan10Frames = unique(ROIPairsOverlappingMoreThan10Frames(:));
end
  
end