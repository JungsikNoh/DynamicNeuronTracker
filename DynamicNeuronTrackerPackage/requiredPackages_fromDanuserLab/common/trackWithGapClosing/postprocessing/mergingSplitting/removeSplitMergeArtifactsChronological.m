function   tracksReform = removeSplitMergeArtifactsChronological( tracksIn,thresholdTime)
%REMOVESPLITMERGEARTIFACTSCHRONOLOGICAL remove potentially artifactual merges and
%                           splits, resulting for instance from detection
%
% SYNOPSIS  seqOfEvents  = removeSplitMergeArtifactsChronological(seqOfEventsIn,replaceSegNum,timeThreshold)
%
% INPUT
%       seqOfEventsIn     :   Matrix with number of rows equal to number
%                              of events happening in a compound track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 = start of track segment, 2 = end of track segment;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN = start is a birth and end is a death,
%                                   number = start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%
%      thresholdTime      : structure with the following fiels
%
%                            thersSt2M :   threshold time start to merge
%                            Default=1.
%                            thresS2M  :   threshold time split to merge
%                            Default=1.
%                            thersS2E  :   threshold time split to end
%                            Default=0.
%
%
% OUTPUT
%          seqOfEvents     : the seqOfEvents after the possible artifacts
%          are identified. The situations that are taken into account are:
%          1) a merge followed by a split;
%          2) a merge of a segment that just appeared
%          3) a split followed by a remerge
%          4) a split followed by a vanishing
%
% Luciana de Oliveira, August 2017
% Modified Luciana de Oliveira, December 2017.

%% Input


%get threshold information

if isfield(thresholdTime,'thersSt2M')
    thersSt2M = thresholdTime.thersSt2M;
else
    thersSt2M = 1;
end

if isfield(thresholdTime,'thresS2M')
    thresS2M = thresholdTime.thresS2M;
else
    thresS2M = 1;
end

if isfield(thresholdTime,'thersS2E')
    thersS2E = thresholdTime.thersS2E;
else
    thersS2E = 0;
end

% load compTracks
compTracks=tracksIn;
%initiate the reform compTrack
tracksReform=tracksIn;
% go over all compTracks
% iTrackNew=length(compTracks);
% add a new field in reformTrack to take account of all the modifications
% done in the compTracks
tracksReform(1).oldTracksInfo=[];

 for iTrack = 1: length(compTracks);
%     for iTrack = 15;
    %initiate the flag for segments that were already checked
    
%     segCheck=[];
    
%     fprintf('\nProcessing track=%s ',int2str(iTrack));
    
    % load seqOfEvents
    
    seqOfEventsIn = compTracks(iTrack).seqOfEvents;
    
    % load clustHistory
    %     clustHistory = clustHistoryAll{iTrack,1};
    
    %calculate length seqOfEvents
    
    lengthSeqOfEvents= size(seqOfEventsIn,1);
    
    % if the seqOfEvents have more than one segment
    
    if lengthSeqOfEvents > 2
        
        %% identify artifacts
        
%         tracksFeatIndxCG=compTracks(iTrack).tracksFeatIndxCG;
%         numberOfSegments=size(tracksFeatIndxCG,1);
        
        %Identify which events come from artifacts
        
        %copy sequence of events
        
        seqOfEvents = seqOfEventsIn;
        
        %find indices for which a start of a new track is not due a birth and
        % an end is not due a death
        
        eventRows = find(~isnan(seqOfEvents(:,4)));
        
        % initiate the vector for the information of segments that were modified
        
        segInden=zeros(length(eventRows),1); % segIden is the identification of the
        % segment that was already analysed in modified in seqOfEvents
        
        %initiate the index that will take account the segments that are included
        %in segIden
        indexSeg=1;
        
        % determine if the event is a merge or a split and depending on that make
        % the tests for artifacts
        
        %go over all the events
        
        for indexEvent=1:length(eventRows)
            
            % to be sure of analysing only the rows that some event is happening
            iEvent=eventRows(indexEvent);
            
            % determine if the event was already analysed
            segFlag=find(segInden==iEvent, 1);
            
            if isempty(segFlag)
                
                % determine the kind of event
                kindOfEvent=seqOfEvents(iEvent,2);
                
                %check for merging
                if kindOfEvent==2
                    % merge part
                    
                    % for the merging we will have two possible artifacts:
                    % (a) if there is a split right after the merge
                    % (b) if the segment is just borm. This condition will only be tested  if
                    % there is no split after the merge.
                    
                    
                    %find the two merging segments and the time of merging
                    segment1 = seqOfEvents(iEvent,3);
                    segment2 = seqOfEvents(iEvent,4);
                    timeMerge = seqOfEvents(iEvent,1);
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % check if a split happened right after the merge
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    %             %check whether these two segments split from each other again
                    %             iSplit = find(seqOfEvents(:,4)==segment2 & seqOfEvents(:,2)==1);
                    %
                    %             %get times of split
                    %             splitTime = seqOfEvents(iSplit,1);
                    %
                    %             %keep the earliest time of split that comes after the time of
                    %             %merge
                    %             iSplit = iSplit(splitTime>timeMerge);
                    %             iSplit = min(iSplit);
                    %             splitTime = seqOfEvents(iSplit,1);
                    %
                    %             %if there is a consequent split ...
                    %             if ~isnan(iSplit)
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % %%%%%%%%%%%%%%%%%%%%%Modification LRO 2017/12/14%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % % modification 1: comment the part of the function taking care of merge to split times and taking of all of those that
                    % % are bellow the threshold.
                    % %
                    %
                    %  %calculate merge-to-split time
                    %                 timeMerge2SplitTmp = splitTime - timeMerge;
                    %       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %             % couting merge2 split time
                    %             countMerge2SplitTemp(timeMerge2SplitTmp)=countMerge2SplitTemp(timeMerge2SplitTmp)+1;
                    %
                    %
                    %                 %%%%% add a condition to check only until the values of prob that I have
                    %       if timeMerge2SplitTmp <= length(probArtifactMerge2split)
                    %
                    %                 % index for the probe, the lifetime will be in timeStep
                    %                 % unities and I need it in integer quantities
                    %
                    % %                indexLifeTime= timeMerge2SplitTmp*10;
                    % % % % % % % %
                    % % % % % % % %                 %if time between merge and split is smaller or equal the theresholdTime,
                    % % % % % % % %                 %replace this merge for a end of segment
                    % % % % % % % %
                    % % % % % % % %                 if timeMerge2SplitTmp <=thresM2S
                    %
                    % % modification 2: now if it is merge followed by an split, the test to be
                    % % done is if this is an artifact or not. This test is done using the
                    % % probability to be an artifact calculated previously
                    %
                    % % calculate the prob of this observation be an artifact
                    % probMerge2splitArtifacts=rand;
                    %
                    % % Determine which rows will be removed
                    %
                    % if probMerge2splitArtifacts<probArtifactMerge2split(timeMerge2SplitTmp)
                    %
                    %
                    %
                    %
                    %     %Replace the segment that merge with an end and the time
                    %     %where the segment ends by mergeTime-1
                    %
                    %
                    %     seqOfEvents(iEvent,4) = NaN;
                    %     seqOfEvents(iEvent,1) = timeMerge-1;
                    %
                    %     % Replace the segment that start from the split with a
                    %     % born event
                    %
                    %     seqOfEvents(iSplit,4) = NaN;
                    %
                    %     % save these segments in the vector segInden
                    %     segInden(indexSeg)=iEvent;
                    %     segInden(indexSeg+1)=iSplit;
                    %     %increase indexSeg
                    %     indexSeg=indexSeg+2;
                    % end
                    %       end
                    %             end
                    % else there is no split after merge
                    
                    %          comment merge part for now
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                 check if the merge is for a segment that just appeared
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %                 find when each of the segments started
                    iStart1 = seqOfEvents(:,3)==segment1 & ...
                        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1;
                    timeStart1 = seqOfEvents(iStart1,1);
                    iStart2 = seqOfEvents(:,3)==segment2 & ...
                        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1;
                    timeStart2 = seqOfEvents(iStart2,1);
                    
                    %             get the start-to-merge time for both segments
                    timeStart12Merge = timeMerge - timeStart1;
                    timeStart22Merge = timeMerge - timeStart2;
                    
                    %                 if either time is equal or smaller than the threshold, discard the merge
                    
                    if any([timeStart12Merge timeStart22Merge]<=thersSt2M)
                        
                        %                 replace the merge as the end of the segment
                        seqOfEvents(iEvent,4) = NaN;
                        seqOfEvents(iEvent,1) = timeMerge-1;

                        
                        %                     save thes segment in the vector segInden
                        segInden(indexSeg)=iEvent;
                        
                        %                 increase indexSeg
                        indexSeg=indexSeg+1;
                        
                    end
                    
                    
                    
                end % if there is a split
                
                
                %% split part
                
                %else check for spliting
            else
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % check if the split is followed by a merge
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %find the two splitting segments and the time of splitting
                segment1 = seqOfEvents(iEvent,3);
                segment2 = seqOfEvents(iEvent,4);
                timeSplit = seqOfEvents(iEvent,1);
                
                
                %check whether these two segments merge with each other again
                iMerge = find(any(seqOfEvents(:,3:4)==segment1,2) & ...
                    any(seqOfEvents(:,3:4)==segment2,2) & seqOfEvents(:,2)==2);
                
                %if they merge ...
                if ~isempty(iMerge)
                    
                    %calculate split to merge time
                    timeMerge = seqOfEvents(iMerge,1);
                    timeSplit2Merge = timeMerge - timeSplit;
                    
                    %if the split to merge time is lower or equal timeThreshold
                    if timeSplit2Merge <=thresS2M %== 1%timeThreshold
                        
                        
                        %Replace the segment that split as a new birth
                        seqOfEvents(iEvent,4) = NaN;
                        % Replace the merge as the end of the segment
                        seqOfEvents( iMerge,4) = NaN;
                        
                        
                        % save these segments in the vector segInden
                        segInden(indexSeg)=iEvent;
                        segInden(indexSeg+1)=iMerge;
                        %increase indexSeg
                        indexSeg=indexSeg+2;
                    end % if the split is followed by a merge
                    
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % check if the split is followed by the end of the segment
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %find when each of the segments ends
                    iEnd1 = seqOfEvents(:,3)==segment1 & ...
                        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2;
                    timeEnd1 = seqOfEvents(iEnd1,1);
                    iEnd2 = seqOfEvents(:,3)==segment2 & ...
                        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2;
                    timeEnd2 = seqOfEvents(iEnd2,1);
                    
                    %get the split-to-end time for both segments
                    timeSplit2End1 = timeEnd1 - timeSplit;
                    timeSplit2End2 = timeEnd2 - timeSplit;
                    
                    %if either time is equal to 0 frame, discard the split
                    if any([timeSplit2End1 timeSplit2End2]==thersS2E)
                        
                        %Replace the segment that split as a new birth
                        seqOfEvents(iEvent,4) = NaN;
                        
                        % save thes segment in the vector segInden
                        segInden(indexSeg)=iEvent;
                        
                        %increase indexSeg
                        indexSeg=indexSeg+1;
                        
                    end
                end
                
            end
        end
        tracksReform(iTrack).seqOfEvents= seqOfEvents;
    end
    
end

%% regroup segments
tracksReform = reformCompTracksAfterRemoveArtifacts(tracksReform);

end
