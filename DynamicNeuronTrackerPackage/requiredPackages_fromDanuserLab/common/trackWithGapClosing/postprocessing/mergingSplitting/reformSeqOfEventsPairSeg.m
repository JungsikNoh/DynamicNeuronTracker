function seqOfEvents = reformSeqOfEventsPairSeg( segmentMerge,segmentSplit,mergeIndex,splitIndex,possibleArtifactsTime,eventTime,firstFrame,seqOfEvents,tracksCoordAmpCG)
% reformSeqOfEventsPairSeg reform seq of events considering intensity and
% moviment to pairwise segments. This fuction is to be used in events of
% merging and spliting.
%
% SYNOPSIS
%
% seqOfEvents = removesimultaneousMergeSplit( seqOfEventsIn )
%
% INPUT
%        INPUT
%               segmentMerge     : the segment that is merging
%               segmentSplit     : the segment that is splitting
%               seqOfEvents      : seqOfEvents coming from the specific track
%        number
%               tracksCoordAmpCG : tracks coordinates coming from the specific track
%        number
%
% OUTPUT
%          seqOfEvents: reform seqOfEvents, considering the pairwise
%          conditions
%
%
% Luciana de Oliveira, Febuary 2018


%%
% pre-alocate the ratios for the calculation of
% costMatrix
ratioIntensity=zeros(length(segmentMerge),length(segmentSplit));
ratioMeanDisp=zeros(length(segmentMerge),length(segmentSplit));

%%
% calculate number of merges and splits
numberMerges=1:length(segmentMerge);
numberSplits=1:length(segmentSplit);


% calculate all the possible combinations between merges and splits

[mergeIndices,splitIndices] = meshgrid(numberMerges,numberSplits);

% Concatenate indices
concatenateIndices=cat(2,mergeIndices',splitIndices');
% final combination of all merges and splits
% indices
indicesMS=reshape(concatenateIndices,[],2);

%%
for indexMS=1: size(indicesMS,1)
    % calculate mean intensity and displacements for all
    % merges and splits
    
    %intensity
    
    meanIntMerge=nanmean(tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),4:8:8*(eventTime-firstFrame+1)-1),2);% take into account that track does not necessarily start at frame 1, in all the times
    
    
    % take the x and y values
    mergeXvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),1:8:8*(eventTime-firstFrame+1)-1);
    mergeYvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),2:8:8*(eventTime-firstFrame+1)-1);
    
    %calculate displacement
    
    dispMergeX=mergeXvals(:,2:end)-mergeXvals(:,1:end-1);
    dispMergeY=mergeYvals(:,2:end)-mergeYvals(:,1:end-1);
    
    
    % calculate mean displacement
    
    meanDispMerge=nanmean(sqrt(dispMergeX.^2+dispMergeY.^2),2);
    
    
    % splits
    
    %intensity
    
    % take into account that track does not necessarily start at frame 1, in all the times
    meanIntSplit=nanmean(tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),8*(eventTime-firstFrame+1)+1:8:end),2);
    
    % take the x and y values
    
    splitXvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),8*(eventTime-firstFrame+1)+1:8:end);
    splitYvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),8*(eventTime-firstFrame+1)+2:8:end);
    
    %calculate displacement
    
    dispSplitX=splitXvals(:,2:end)-splitXvals(:,1:end-1);
    dispSplitY=splitYvals(:,2:end)-splitYvals(:,1:end-1);
    
    % calculate mean displacement
    
    meanDispSplit=nanmean(sqrt(dispSplitX.^2+dispSplitY.^2),2);
    
    % calculate ratio
    % check which mean have the largest value
    
    %intensity
    
    intLarger=max(meanIntMerge,meanIntSplit);
    intSmaller=min(meanIntMerge,meanIntSplit);
    
    %displacement
    
    dispLarger=max(meanDispMerge,meanDispSplit);
    dispSmaller=min(meanDispMerge,meanDispSplit);
    
    %calculate ratio
    ratioIntensity(indicesMS(indexMS,1),indicesMS(indexMS,2))=intLarger/intSmaller;
    ratioMeanDisp(indicesMS(indexMS,1),indicesMS(indexMS,2))=dispLarger/dispSmaller;
    
end


% for some cases where the split occour and it is
% followed by a end we have ratios with value NaN, and
% NaN is not a valid input for the LAP function. To
% avoid this problem, in the cases where we have
% multiple merges and split, we will take the maximum
% value of cost matrix and apply for the NaN value



if  all(all(isnan (ratioMeanDisp))) % if any ratio displacement is NaN
    
    %if all the displacements are NaN we will take only
    %the value of intensity for the comparison
    
    ratioMeanDisp(:)=1;
    
elseif any(any(isnan(ratioMeanDisp)))
    
    indiceNaNDisp= isnan(ratioMeanDisp);
    newValRatioDisp=max(max(ratioMeanDisp));
    ratioMeanDisp(indiceNaNDisp)=newValRatioDisp;
    
end

% calculate the costMatrix

costMatrix=ratioIntensity.*ratioMeanDisp;


% after all the combinations are calculated and saved in costMatrix, we can
% analyse which is the most likely to be a pair of merge and split
[mergePositions, splitPositions] = lap(costMatrix, [],[],1,[]);

%%
%here I need to combine all the possibilities from merge and splits. If for
%example I have multiples merges but only one split, from the output of lap
%I need to have the identification of this split is going with this merge.
%maybe the best solution could be combine only the number of merges with
%the splits, given by the length of merges and splits

if length(numberMerges)<length(numberSplits)
    goodPairs=find(mergePositions<=max(numberMerges)); % for some cases one can have more than only one merge or split, for that in selecting the good
    %pairs need to take the maximum of the number of
    %merges or splits.
else
    goodPairs=find(splitPositions<=max(numberSplits));
end
pairsFinal=[mergePositions(goodPairs),splitPositions(goodPairs)];

% go over all the pairs and reform seqOfEvents and tracks

for indexPair=1:size(pairsFinal,1)
    segMergeTemp=segmentMerge(pairsFinal(indexPair,1));% it gives the segments that is merging based in the list of segments that merge
    segSplitTemp=segmentSplit(pairsFinal(indexPair,2));% it gives the segments that is splitting based in the list of segments that split
    
    % define which are the segments to be removed
    removeIndex=[possibleArtifactsTime(mergeIndex(pairsFinal(indexPair,1))),possibleArtifactsTime(splitIndex(pairsFinal(indexPair,2)))];   % it  finds the position in the list of the artifacts to be removed
    
    %remove the merge and the split
    seqOfEvents(removeIndex,:)=[];
    
    % replace the segment number to have only one continuos
    % segment
    seqOfEvents(seqOfEvents(:,3)==segSplitTemp,3)=segMergeTemp;
    seqOfEvents(seqOfEvents(:,4)==segSplitTemp,4)=segMergeTemp;
    
    % for all segments with index larger than removed
    % segment, reduce their index by 1
    
    seqOfEventsIndex= seqOfEvents(:,3:4);
    segmentNumbersOld=unique(seqOfEventsIndex);
    segmentNumbersOld(isnan(segmentNumbersOld))=[];
    for i=1:length(segmentNumbersOld)
        if segmentNumbersOld(i)~=i
            seqOfEventsIndex(seqOfEventsIndex==segmentNumbersOld(i))=i;
        end
    end
    
    %update seqOfEvents
    seqOfEvents(:,3:4)=seqOfEventsIndex;
end

