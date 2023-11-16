function [tracksIterXYZCell, iterStatus] = tracksXYZ_EMiteration(...
            tracksIterX1, tracksIterY1, tracksIterZ1, ...
            imgArray, PSourMat2, upbdDeformX, upbdDeformZ, bandWidthX, bandWidthZ, ...
            corrThreshold1, imgMarginSizeX, imgMarginSizeZ)
% tracksXYZ_EMiteration implement a set of EM iterations of local
% patch-matching while updating local mean images of firing neurons and 
% their trajectories in order to find the best connections between firing 
% events from the identical neurons.


%% initialize iteration 

iterStatus = struct('condMeanCor', [], 'totFramesSatisfied', [], ...
        'iterMerged', 0, 'iterConverged', 0, 'L2distPerfr', nan, 'hamdist', nan);
for k = 1:size(tracksIterX1, 1)
    iterStatus(k).condMeanCor = nan; %tracksIter(k).condMeanCor;
    iterStatus(k).totFramesSatisfied = nan; %tracksIter(k).totFramesSatisfied;
    iterStatus(k).iterMerged = 0;
    iterStatus(k).iterConverged = 0;
    iterStatus(k).L2distPerfr = NaN;
    iterStatus(k).hamdist = NaN;
end


tracksIterXYZCell = cell(3,5);
tracksIterXYZCell{1, 1} = tracksIterX1;
tracksIterXYZCell{1, 2} = tracksIterY1;
tracksIterXYZCell{1, 3} = tracksIterZ1;
tracksIterXYZCell{1, 4} = 1:size(tracksIterX1, 1);
tracksIterXYZCell{1, 5} = size(tracksIterX1, 1);

% only to handle when iterations end at 2 in an extreme case.
tracksIterXYZCell{3, 1} = tracksIterX1;
tracksIterXYZCell{3, 2} = tracksIterY1;
tracksIterXYZCell{3, 3} = tracksIterZ1;
tracksIterXYZCell{3, 4} = 1:size(tracksIterX1, 1);
tracksIterXYZCell{3, 5} = size(tracksIterX1, 1);

%% Run one EM iteration for all tracks and merge non-unique updated tracks 
%% for efficient iterations

for trIter = 2:20
    %% update tracks
    
    tracksIterX_pre = tracksIterXYZCell{trIter-1, 1};
    tracksIterY_pre = tracksIterXYZCell{trIter-1, 2};
    tracksIterZ_pre = tracksIterXYZCell{trIter-1, 3};

    tracksIterX = tracksIterX_pre;
    tracksIterY = tracksIterY_pre;
    tracksIterZ = tracksIterZ_pre;
    
    min0 = upbdDeformX + bandWidthX + 1 + imgMarginSizeX;
    minZ0 = upbdDeformZ + bandWidthZ + 1 + imgMarginSizeZ;

    xmax0 = size(imgArray, 2) - upbdDeformX - bandWidthX - imgMarginSizeX;
    ymax0 = size(imgArray, 1) - upbdDeformX - bandWidthX - imgMarginSizeX;
    zmax0 = size(imgArray, 3) - upbdDeformZ - bandWidthZ - imgMarginSizeZ;

    totTr = size(tracksIterX_pre, 1); 
    
    disp('== EM iteration ID:')
    disp(trIter)
    %ll = 0;
    
    parfor itr = 1:totTr
        fprintf('%05d ', itr);
        if (mod(itr, 50) == 0); fprintf('\n'); end 
        
        triX = tracksIterX_pre(itr, :)';
        triY = tracksIterY_pre(itr, :)';
        triZ = tracksIterZ_pre(itr, :)';

        if (iterStatus(itr).iterMerged == 0) && (iterStatus(itr).iterConverged == 0)

            xcenter = round(nanmean(triX));
            ycenter = round(nanmean(triY));
            zcenter = round(nanmean(triZ));

            if (min(xcenter, ycenter) >= min0) && (zcenter >= minZ0) ...
                    && (xcenter<=xmax0) && (ycenter<=ymax0) && (zcenter<=zmax0)         
    
                % updateOneTrack_3D_collectPSin1PxlDistance
                [trackPost, ~] = updateOneTrack_3D_collectPSin1PxlDistance(...
                    triX, triY, triZ, PSourMat2, imgArray, ...
                    upbdDeformX, upbdDeformZ, bandWidthX, bandWidthZ, corrThreshold1);


                dvec = sqrt( (triX - trackPost.x).^2 + (triY-trackPost.y).^2 + (triZ-trackPost.z).^2 );
                trackPost.L2distPerfr = nanmean(dvec);

                nanhamXYZ = nanhamdistXYZ_ignoreAllnan([triX', triY', triZ'], ...
                    [trackPost.x', trackPost.y', trackPost.z']);
                trackPost.hamdist = nanhamXYZ;

                % record convergence 
                trackPost.merged = 0;           % will be computed below
                trackPost.converged = double(trackPost.hamdist == 0);

                % out
                tracksIterX(itr, :) = trackPost.x';
                tracksIterY(itr, :) = trackPost.y';
                tracksIterZ(itr, :) = trackPost.z';

                iterStatus(itr).condMeanCor = [iterStatus(itr).condMeanCor, trackPost.condMeanCor];
                iterStatus(itr).totFramesSatisfied = [iterStatus(itr).totFramesSatisfied, trackPost.totFramesSatisfied];
                iterStatus(itr).L2distPerfr = [iterStatus(itr).L2distPerfr, trackPost.L2distPerfr];
                iterStatus(itr).hamdist = [iterStatus(itr).hamdist, trackPost.hamdist];
                if (trackPost.converged == 1)
                    iterStatus(itr).iterConverged = trIter;
                end
            end
        end
    end

    %% summarize condMeanCorVec, totFramesSatisfiedVec

    L = size(tracksIterX, 1);
    condMeanCorVec = nan(L,1);
    totFramesSatisfiedVec = nan(L,1);

    for k = 1:L
        condMeanCorVec(k) = iterStatus(k).condMeanCor(end);
        totFramesSatisfiedVec(k) = iterStatus(k).totFramesSatisfied(end);
    end
    
    disp('==')
    disp('== Summary of per-track averages of correlations above the threshold')
    disp(summaryStatistics(condMeanCorVec))

    disp('== Summary of per-track total frames with correlation above the threshold')
    disp(summaryStatistics(totFramesSatisfiedVec))

    %% identify the uniqueness using hamming dist 

    TrCurr_pre = tracksIterXYZCell{trIter-1, 4};
    TrOut_pre = setdiff(1:size(tracksIterX, 1), TrCurr_pre);
    tracksIterXYZ = [tracksIterX, tracksIterY, tracksIterZ]; 
 
    tracksIterXYZ(TrOut_pre, :) = nan;

    DD = pdist(tracksIterXYZ, @nanhamdistXYZ_ignoreAllnan); 

    %
    i = find(DD==0);
    rr = nan(size(i));
    cc = nan(size(i));

    tmp = L-1:-1:1;
    tmp1= cumsum(tmp);

    for j = 1:numel(i)
        rr(j) = find( i(j) <= tmp1, 1);
        if (rr(j) == 1)
            cc(j) = rr(j) + i(j);
        else
            cc(j) = rr(j) + i(j) - tmp1(rr(j)-1);
        end
    end
    
    %rr
    %cc

    %% merge non-unique tracks & output

    if ~isempty(cc)
        for k = unique(cc)
            iterStatus(k).iterMerged = trIter;
        end
    end

    TrCurr = setdiff(TrCurr_pre, cc);
    numTrCurr = numel(TrCurr);

    tracksIterXYZCell{trIter, 4} = TrCurr;
    tracksIterXYZCell{trIter, 5} = numTrCurr;

    tracksIterXYZCell{trIter, 1} = tracksIterX;
    tracksIterXYZCell{trIter, 2} = tracksIterY;
    tracksIterXYZCell{trIter, 3} = tracksIterZ;

    %% report

    disp('== iteration: ')
    disp(trIter)
    disp('== num of unique tracks')
    disp(numTrCurr)

    %
    iterConverged = [iterStatus.iterConverged];
    iterMerged = [iterStatus.iterMerged]; 
    ind = max(iterConverged,  iterMerged);
    disp('== num of not merged nor converged')
    disp(sum(ind==0))

    disp('== hamming distance average in this iteration')
    iterHamdist = nan(numel(iterStatus), 1);
    for k = 1:numel(iterStatus)
        if (numel(iterStatus(k).hamdist) ==  trIter)
            iterHamdist(k) = iterStatus(k).hamdist(end);
        end
    end
    disp(nanmean(iterHamdist))

    disp('== whether tracks remain the same?')
    iterHamdist1 = iterHamdist(~isnan(iterHamdist));
    if all(iterHamdist1 == 0)
        disp('=the same=')
        break
    else
        disp('==different==')
    end

    % Break the EM iteration in the following conditions:
    % (1) non-merged nor converged tracks are less than 20% of initial # of
    % tracks; (2) iteration steps >= 20; (3) recent 16 iterations didn't
    % make tracks converged. 
    if (sum(ind==0) <= 0.2*tracksIterXYZCell{1,5})
        if (trIter >= 20) && all(diff([tracksIterXYZCell{trIter-16:trIter, 5}]) == 0)
            disp('= Non-converged tracks are less than 20%, and numTrCurr remain the same. =')
            break
        end
    end 

end

end
