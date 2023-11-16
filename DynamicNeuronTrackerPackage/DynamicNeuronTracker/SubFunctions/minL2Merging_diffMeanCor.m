function [final_TrCurr3, iterStatus] = minL2Merging_diffMeanCor(...
                distBtwTracks_toBeMerged, tracksXcmerged0cs, ...
                tracksYcmerged0cs, tracksZcmerged0cs, iterStatus)
% minL2Merging_diffMeanCor implements merging of trajectories to make sure
% that all trajectories are unique and ROIs are not overlapping to much.
% When two tracks are two close or overlapping, one is merged to the other
% that is better in terms of maximum rolling spatial correlations. 


%% identify too close tracks

tracksXYZc = [tracksXcmerged0cs, tracksYcmerged0cs, tracksZcmerged0cs]; 

% minL2distXY
DD = pdist(tracksXYZc, @minL2distXYZ); 

L = size(tracksXcmerged0cs, 1);
L*(L-1)/2;

%
i = find(DD <= distBtwTracks_toBeMerged);
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

%% Merge to a track with the best diffMeanCor score

for i = 1:numel(rr)
    if (iterStatus(rr(i)).diffMeanCor >= iterStatus(cc(i)).diffMeanCor)
        iterStatus(cc(i)).mergedTo = rr(i);
    else
        iterStatus(rr(i)).mergedTo = cc(i);
    end
end    

%% merged (tracksXcmerged0s -> tracksXcmerged0sL2)

final_TrCurr3 = find([iterStatus.mergedTo] == 0);

end
