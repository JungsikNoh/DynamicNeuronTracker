function tracksReform = removeSimultaneousMergeSplit( tracksIn )
% function  tracksReform = removeSimultaneousMergeSplit( tracksIn )
% This function takes care of the artifacts coming from simultaneous merges
% and splits occourring in the same frame with the same segments. It
% reforms the compTracks.
%
% SYNOPSIS
%
% seqOfEvents = removesimultaneousMergeSplit( seqOfEventsIn )
%
% INPUT
%        INPUT
%               tracksIn       : Output of trackCloseGapsKalman.
% OUTPUT
%          tracksReform     : tracks after removing artifacts coming from
%          the simultaneous merge and split events.
%
%
% Luciana de Oliveira, January 2018
% Modification March 2018
%
% This function is divided in 4 parts:
%
% 1) Check if there are more than one event happening in the same frame;
% 2) Check if these events are happining with the same segments;
% 3) Check if the events are a combination of merges and splits.
% 4) Reform the seqOfEvents and compTracks. In this last step it will
% decouple the segments that are not part of the same compTrack anymore and
% correct starting times for the segments and segment numbers in the
% seqOfEvents.

%%
% load compTracks
compTracks=tracksIn;


%initiate the reform compTrack
tracksReform=tracksIn;

% go over all the tracks

for iTrack = 1: length(compTracks);
  
    
    %load sequence of events
    
    seqOfEvents= compTracks(iTrack).seqOfEvents;
    
    %For experimental data, the first iteration (frame) in seqOfEvents
    %can be some value other than 1.  This is occurs becuase an object
    %can appear or disappear at any time while in simulation, all of
    %the objects are present at frame 1. This affects how aggregState
    %is accessed since the frame # is supposed to correspond to a
    %column in aggregState. So, need to shift frame numbers in
    %seqOfEvents, if necessary, when accessing aggregState. The very
    %first frame # is the shift amount.
    firstFrame= seqOfEvents(1,1);
    
    
    %load tracksFeatIndxCG and tracksCoordAmpCG
    
    tracksFeatIndxCG= tracksReform(iTrack).tracksFeatIndxCG;
    tracksCoordAmpCG= tracksReform(iTrack).tracksCoordAmpCG;
    
    %calculate length seqOfEvents
    
    lengthSeqOfEvents= size(seqOfEvents,1);
    
    
    % if the seqOfEvents have more than one segment
    
    
    if lengthSeqOfEvents > 4
        
        %find indices for which a start of a new track is not due a birth and
        % an end is not due a death
        
        eventRows = ~isnan(seqOfEvents(:,4));
        
        %% Check if the events are happening in the same frame
        
        % find which are the times that appear more than two times in the same
        % seqOfEvents, because these are candidates to have simultaneous merge and
        % split events
        
        % identify the frames where the events are happening
        timeEvents=seqOfEvents(eventRows,1);
        uniqueTimeEvents=unique(timeEvents);
        timeRepetition=uniqueTimeEvents(1<histc(timeEvents,unique(timeEvents)));
        
        % if there are multiple events happening in the same frame
        while timeRepetition
            % check for each time that repeats if the events are simultaneous
            
            eventTime=timeRepetition(1);
            %% check if the events are happening with the same segment
            
            
            % identify the possible segments with artifacts
            
            possibleArtifactsTime=find(seqOfEvents(:,1)==eventTime);
            
            %find if they are interacting with the same segment
            segmentInteraction=seqOfEvents(possibleArtifactsTime,4);
            
            %determine which ones are the same and how much times it appears
            uniqueSegment=unique(segmentInteraction);
            
            % remove NaNs
            uniqueSegment(isnan(uniqueSegment))=[];
            
            %count the repetitions
            segmentRepetition=histc(segmentInteraction,uniqueSegment);
            
            
            for indexSegment=1:length(uniqueSegment)
                
                % if there are multiple events happening with the same segment
                
                if segmentRepetition(indexSegment)>1
                    
                    
                    % check if it is a merge or a split
                    
                    mergeIndex=find(seqOfEvents(possibleArtifactsTime,2)==2 & seqOfEvents(possibleArtifactsTime,4)==uniqueSegment(indexSegment) );
                    
                    splitIndex=find(seqOfEvents(possibleArtifactsTime,4)==uniqueSegment(indexSegment) & seqOfEvents(possibleArtifactsTime,2)==1);
                    
                    % pair merge and split, because there are cases where the number of
                    % merges and splits are not equal, so the solution is pairwise all the
                    % merges and splits and leave the "extra" events how they are.
                    
                    
                    
                    %% reform seqOfevents- replace segment numbers
                    
                    % identify segment merge and split
                    
                    segmentMerge=seqOfEvents(possibleArtifactsTime(mergeIndex),3);
                    segmentSplit=seqOfEvents(possibleArtifactsTime(splitIndex),3);
                    %                     segmentsContinuous=unique(seqOfEvents(possibleArtifactsTime,4));
                    
                    % test if there are multiple merges and split events
                    % happening simultaneously
                    if length(mergeIndex)==1 && length(splitIndex)==1 % if there is only one merge and one split
                        
                        %find rows to be removed
                        rowMerge=possibleArtifactsTime(mergeIndex);
                        rowSplit=possibleArtifactsTime(splitIndex);
                        rowMergeSplit=[rowMerge,rowSplit];
                        
                        % the continuous segment is maintened as it is and
                        % pairwise the other two segments
                        
                        % paiwise tracks
                        [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,segmentMerge,segmentSplit);
                        
                        
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         %find rows to be removed
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         rowMerge=possibleArtifactsTime(1);
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         rowSplit=possibleArtifactsTime(2);
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         % replace the segment number to have only one continuos
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         % segment
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         seqOfEvents(seqOfEvents(:,3)==segmentSplit,3)=segmentMerge;
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         seqOfEvents(seqOfEvents(:,4)==segmentSplit,4)=segmentMerge;
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         %remove the merge and the split
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         seqOfEvents([rowMerge,rowSplit],:)=[];
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         % for all segments with index larger than removed
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         % segment, reduce their index by 1
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         seqOfEventsIndex= seqOfEvents(:,3:4);
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         segmentNumbersOld=unique(seqOfEventsIndex);
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         segmentNumbersOld(isnan(segmentNumbersOld))=[];
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         for i=1:length(segmentNumbersOld)
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             if segmentNumbersOld(i)~=i
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 seqOfEventsIndex(seqOfEventsIndex==segmentNumbersOld(i))=i;
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             end
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         end
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         %update seqOfEvents
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         seqOfEvents(:,3:4)=seqOfEventsIndex;
                        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                        
                        seqOfEvents = reformSeqOfEventsFromMerge2split(seqOfEvents,segmentMerge,segmentSplit,rowMergeSplit);
                        
                        %% save tracks in the final format
                        tracksReform(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
                        tracksReform(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
                        tracksReform(iTrack).seqOfEvents = seqOfEvents;
                        
                    else % if there are multiple simultaneous merges and splits
                        
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % pre-alocate the ratios for the calculation of
                        % costMatrix
                        ratioIntensity=zeros(length(segmentMerge),length(segmentSplit));
                        ratioMeanDisp=zeros(length(segmentMerge),length(segmentSplit));
                        
                        
                        % calculate number of merges and splits
                        numberMerges=length(segmentMerge);
                        numberSplits=length(segmentSplit);
                        
                        
                        % calculate all the possible combinations between merges and splits
                        
                        [mergeIndices,splitIndices] = meshgrid(1:numberMerges,1:numberSplits);
                        
                        %put indices together
                        indicesMS=[mergeIndices,splitIndices];
                        
                        for indexMS=1: size(indicesMS,1)
                            % calculate mean intensity and displacements for all
                            % merges and splits
                            
                            %intensity
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % modification LRO 2018/03/16
                            % included the condition to calculate the
                            % intensity only in the interval between the last event before the
                            % merge and the first after the split.
                            
                            timeBefAftEvent  = timeBeforeAfterEvent(seqOfEvents,segmentMerge(indicesMS(indexMS,1)),segmentSplit(indicesMS(indexMS,2)),eventTime);
                            
                            
                            
                            % check if there is an event before the merge for both merge and split
                            
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % take only the events before the merge
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             possibleBefMerge=eventRows(:,1)<eventTime;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             possibleAftMerge=eventRows(:,1)>eventTime;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % take only the seqOfEvents part that has no NaNs and is before and
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % after the merge 2 split event
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             seqOfEventsBefore = seqOfEvents(notNanEvents(possibleBefMerge),:);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             seqOfEventsAfter  = seqOfEvents(notNanEvents(possibleAftMerge),:);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % find if there is events before and after
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % merge for both merge and continuous segments
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             [eventBefMerge,~]=find(seqOfEventsBefore(:,3:4)==segmentMerge);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             [eventAfterMerge,~]=find(seqOfEventsAfter(:,3:4)==segmentMerge);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             [eventBefSplit,~]=find(seqOfEventsBefore(:,3:4)==segmentSplit);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             [eventAfterSplit,~]=find(seqOfEventsAfter(:,3:4)==segmentSplit);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             if eventBefMerge
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeBefMerge=seqOfEvents(eventBefMerge(end),1);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             else
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeBefMerge=1;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             end
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             if eventAfterMerge
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeAftMerge=seqOfEvents(eventAfterMerge(1),1);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             else
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeAftMerge=1;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             end
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             if eventBefSplit
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeBefSplit=seqOfEvents(eventBefSplit(end),1);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             else
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeBefSplit=1;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             end
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             if eventAfterSplit
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeAfterSplit=seqOfEvents(eventAfterCont(1),1);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             else
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 timeAfterSplit=1;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             end
                            
                            % calculate mean intensity
                            
                            % take into account that track does not necessarily start at frame 1, in all the times
                            
                            % merge segment
                            
                            meanIntMerge=nanmean(tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),[8*(timeBefAftEvent(1)-firstFrame)+4:8:8*(eventTime-firstFrame+1),8*(eventTime+1)+4:8:8*(timeBefAftEvent(2)-firstFrame+1)]),2);
                            
                            %split segment
                            
                            meanIntSplit=nanmean(tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),[8*(timeBefAftEvent(1)-firstFrame)+4:8:8*(eventTime-firstFrame+1),8*(eventTime+1)+4:8:8*(timeBefAftEvent(2)-firstFrame+1)]),2);
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            % take the x and y values
                            % % % % % %                             mergeXvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),1-firstFrame+1:8:end);
                            % % % % % %                             mergeYvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),2-firstFrame+1:8:end);
                            mergeXvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),[8*(timeBefAftEvent(1)-firstFrame+1):8:8*(eventTime-firstFrame+1),8*(eventTime+1)+1:8:8*(timeBefAftEvent(2)-firstFrame+1)]);
                            mergeYvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),[8*(timeBefAftEvent(1)-firstFrame)+2:8:8*(eventTime-firstFrame+1),8*(eventTime+1)+2:8:8*(timeBefAftEvent(2)-firstFrame+1)]);
                            
                            %calculate displacement
                            
                            dispMergeX=diff(mergeXvals,1,2);
                            dispMergeY=diff(mergeYvals,1,2);
                            
                            
                            % calculate mean displacement
                            
                            meanDispMerge=nanmean(sqrt(dispMergeX.^2+dispMergeY.^2),2);
                            
                            
                            % splits
                            
                            
                            
                            % take the x and y values
                            
                            % % % % % % % % % % % %                             splitXvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),1-firstFrame+1:8:end);
                            % % % % % % % % % % % %                             splitYvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),2-firstFrame+1:8:end);
                            splitXvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),[8*(timeBefAftEvent(1)-firstFrame+1):8:8*(eventTime-firstFrame+1),8*(eventTime+1)+1:8:8*(timeBefAftEvent(2)-firstFrame+1)]);
                            splitYvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)), [8*(timeBefAftEvent(1)-firstFrame)+2:8:8*(eventTime-firstFrame+1),8*(eventTime+1)+2:8:8*(timeBefAftEvent(2)-firstFrame+1)]);
                            %calculate displacement
                            
                            dispSplitX=diff(splitXvals,1,2);
                            dispSplitY=diff(splitYvals,1,2);
                            
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
                        
                        
                        
                        if  all(isnan(ratioMeanDisp(:))) % if any ratio displacement is NaN
                            
                            %if all the displacements are NaN we will take only
                            %the value of intensity for the comparison
                            
                            ratioMeanDisp(:)=1;
                            
                        elseif any(isnan(ratioMeanDisp(:)))
                            
                            indiceNaNDisp= isnan(ratioMeanDisp);
                            newValRatioDisp=max(ratioMeanDisp(:));
                            ratioMeanDisp(indiceNaNDisp)=newValRatioDisp;
                            
                        end
                        
                        % calculate the costMatrix
                        
                        costMatrix=ratioIntensity.*ratioMeanDisp;
                        
                        
                        % after all the combinations are calculated and saved in costMatrix, we can
                        % analyse which is the most likely to be a pair of merge and split
                        [linkMerge2Split, ~] = lap(costMatrix, [],[],1,[]);
                        
                        
                        %here I need to combine all the possibilities from merge and splits. If for
                        %example I have multiples merges but only one split, from the output of lap
                        %I need to have the identification of this split is going with this merge.
                        %maybe the best solution could be combine only the number of merges with
                        %the splits, given by the length of merges and splits
                        
                        goodPairsMerge=linkMerge2Split(1:numberMerges);
                        goodPairs=find(goodPairsMerge<=numberSplits);
                        
                        % put the pairs together
                        
                        pairsFinal=[goodPairs',goodPairsMerge(goodPairs)'];
                        
                        % go over all the pairs and reform seqOfEvents and tracks
                        
                        for indexPair=1:size(pairsFinal,1)
                            segMergeTemp=segmentMerge(pairsFinal(indexPair,1));% it gives the segments that is merging based in the list of segments that merge
                            segSplitTemp=segmentSplit(pairsFinal(indexPair,2));% it gives the segments that is splitting based in the list of segments that split
                            
                            % define which are the segments to be removed
                            rowMergeSplit=[possibleArtifactsTime(mergeIndex(pairsFinal(indexPair,1))),possibleArtifactsTime(splitIndex(pairsFinal(indexPair,2)))];   % it  finds the position in the list of the artifacts to be removed
                            
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             %remove the merge and the split
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             seqOfEvents(removeIndex,:)=[];
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % replace the segment number to have only one continuos
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % segment
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             seqOfEvents(seqOfEvents(:,3)==segSplitTemp,3)=segMergeTemp;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             seqOfEvents(seqOfEvents(:,4)==segSplitTemp,4)=segMergeTemp;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % for all segments with index larger than removed
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             % segment, reduce their index by 1
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             seqOfEventsIndex= seqOfEvents(:,3:4);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             segmentNumbersOld=unique(seqOfEventsIndex);
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             segmentNumbersOld(isnan(segmentNumbersOld))=[];
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             for i=1:length(segmentNumbersOld)
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 if segmentNumbersOld(i)~=i
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                     seqOfEventsIndex(seqOfEventsIndex==segmentNumbersOld(i))=i;
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 end
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             end
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             %update seqOfEvents
                            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             seqOfEvents(:,3:4)=seqOfEventsIndex;
                            
                            seqOfEvents = reformSeqOfEventsFromMerge2split(seqOfEvents,segmentMerge,segmentSplit,rowMergeSplit);
                            
                            %% reformCompTracks
                            
                            % pairwise segments
                            
                            [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,segMergeTemp,segSplitTemp);
                            
                            % save tracks in the final format
                            tracksReform(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
                            tracksReform(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
                            tracksReform(iTrack).seqOfEvents = seqOfEvents;
                        end
                    end
                    
                end
                
            end
            
            %increase index
            
            timeRepetition=timeRepetition(2:end);
        end
        
    end
    
end
%% regroup tracks

tracksReform = reformCompTracksAfterRemoveArtifacts(tracksReform);

end

