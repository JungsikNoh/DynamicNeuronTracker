function tracksReform = reformCompTracksRemoveArtifacts( tracksIn,thresholdTime)
%reformCompTracksRemoveArtifacts will call the function
%removeSplitMergeArtifactsChronological and reform the compTracks based in
%the output of this function
%
% SYNOPSIS  tracksReform = reformCompTracksRemoveArtifacts( tracksIn,clustHistoryAll,probArtifactMerge2split)
%
% INPUT
%               tracksIn       : Output of trackCloseGapsKalman.
%
%         
%          thresholdTime      : structure with the following fiels
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
%          tracksReform     : tracks after removing all artifacts for
%          merging to splitting times. This structure contains the same fields
%          as the original track plus one addicional field:
%          - oldTracksInfo: it is a cell array where the first field is the
%          track number of this new track in the original data and the
%          second field is the number of segment in the orinal data.
%
%
% Luciana de Oliveira, August 2017
% Modified Luciana de Oliveira, December 2017

% This function is divided in two parts:
% 1) identify which segments are part of the same compTrack
% 2) reform tracksFeatIndxCG, tracksCoordAmpCG and seqOfEvents.

% load compTracks
compTracks=tracksIn;
%initiate the reform compTrack
tracksReform=tracksIn;
% go over all compTracks
iTrackNew=length(compTracks);
% add a new field in reformTrack to take account of all the modifications
% done in the compTracks
tracksReform(1).oldTracksInfo=[];

for iTrack = 1: length(compTracks);
    
    %initiate the flag for segments that were already checked
    
    segCheck=[];
    
    fprintf('\nProcessing track=%s ',int2str(iTrack));
    
    % load seqOfEvents
    
    seqOfEventsIn = compTracks(iTrack).seqOfEvents;
    
    % load clustHistory
%     clustHistory = clustHistoryAll{iTrack,1};
    
    %calculate length seqOfEvents
    
    lengthSeqOfEvents= size(seqOfEventsIn,1);
    
    % if the seqOfEvents have more than one segment
    
    if lengthSeqOfEvents > 2
        
        %% identify artifacts
        
        tracksFeatIndxCG=compTracks(iTrack).tracksFeatIndxCG;
        numberOfSegments=size(tracksFeatIndxCG,1);
        
        %Identify which events come from artifacts
        
       
      seqOfEventsArtifacts  = removeSplitMergeArtifactsChronological(seqOfEventsIn,thresholdTime);

   
% % % % % % % % % % % % % % % % %         %% regroup segments
% % % % % % % % % % % % % % % % %         % initiate the cell array that will contain the segments that are
% % % % % % % % % % % % % % % % %         % part of the same group
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         compTrackSegCell=cell(1,numberOfSegments);
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         % check the interactions for each segment in this compTrack
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         for segIndex=1:numberOfSegments
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % check if this segment was already checked
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             if ~ismember(segIndex,segCheck)
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 %  copy seqOfEventsArtifacts, because I will remove the rows
% % % % % % % % % % % % % % % % %                 %  that are already checked
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 seqOfEventsNew=seqOfEventsArtifacts;
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 %initiate a new group of segments
% % % % % % % % % % % % % % % % %                 compTrackSeg=[];
% % % % % % % % % % % % % % % % %               % load compTracks
compTracks=tracksIn;
%initiate the reform compTrack
tracksReform=tracksIn;
% go over all compTracks
iTrackNew=length(compTracks);
% add a new field in reformTrack to take account of all the modifications
% done in the compTracks
tracksReform(1).oldTracksInfo=[];

for iTrack = 1: length(compTracks);
    
    %initiate the flag for segments that were already checked
    
    segCheck=[];
    
    fprintf('\nProcessing track=%s ',int2str(iTrack));
    
    % load seqOfEvents
    
    seqOfEventsIn = compTracks(iTrack).seqOfEvents;
    
    % load clustHistory
%     clustHistory = clustHistoryAll{iTrack,1};
    
    %calculate length seqOfEvents
    
    lengthSeqOfEvents= size(seqOfEventsIn,1);
    
    % if the seqOfEvents have more than one segment
    
    if lengthSeqOfEvents > 2
        
        %% identify artifacts
        
        tracksFeatIndxCG=compTracks(iTrack).tracksFeatIndxCG;
        numberOfSegments=size(tracksFeatIndxCG,1);
        
        %Identify which events come from artifacts
        
       
      seqOfEventsArtifacts  = removeSplitMergeArtifactsChronological(seqOfEventsIn,thresholdTime);

   
% % % % % % % % % % % % % % % % %         %% regroup segments
% % % % % % % % % % % % % % % % %         % initiate the cell array that will contain the segments that are
% % % % % % % % % % % % % % % % %         % part of the same group
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         compTrackSegCell=cell(1,numberOfSegments);
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         % check the interactions for each segment in this compTrack
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         for segIndex=1:numberOfSegments
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % check if this segment was already checked
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             if ~ismember(segIndex,segCheck)
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 %  copy seqOfEventsArtifacts, because I will remove the rows
% % % % % % % % % % % % % % % % %                 %  that are already checked
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 seqOfEventsNew=seqOfEventsArtifacts;
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 %initiate a new group of segments
% % % % % % % % % % % % % % % % %                 compTrackSeg=[];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % determine in which rows this segment appears
% % % % % % % % % % % % % % % % %                 rowsSegment= seqOfEventsNew(:,3)==segIndex|seqOfEventsNew(:,4)==segIndex;
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % determine the number of segments that are interacting
% % % % % % % % % % % % % % % % %                 segmentNumber= unique (seqOfEventsNew(rowsSegment,3:4));
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 %removing NaNs
% % % % % % % % % % % % % % % % %                 segmentNumber(isnan(segmentNumber))=[];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % save these segments to their group
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 compTrackSeg=[compTrackSeg,segmentNumber];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % remove these rows from seqOfEventsArtifacts, I need to do
% % % % % % % % % % % % % % % % %                 % that, otherwise I have an infinite while loop
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 seqOfEventsNew(rowsSegment,:)=[];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 clear rowsSegment
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % check for the segments that are part of this group which
% % % % % % % % % % % % % % % % %                 % are the other segments that they are interacting
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % initiate the flag for the segments interacting with the
% % % % % % % % % % % % % % % % %                 % segment but were already cheched in the previous loop.
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % determine in which rows this segment appears
% % % % % % % % % % % % % % % % %                 rowsSegment= seqOfEventsNew(:,3)==segIndex|seqOfEventsNew(:,4)==segIndex;
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % determine the number of segments that are interacting
% % % % % % % % % % % % % % % % %                 segmentNumber= unique (seqOfEventsNew(rowsSegment,3:4));
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 %removing NaNs
% % % % % % % % % % % % % % % % %                 segmentNumber(isnan(segmentNumber))=[];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % save these segments to their group
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 compTrackSeg=[compTrackSeg,segmentNumber];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % remove these rows from seqOfEventsArtifacts, I need to do
% % % % % % % % % % % % % % % % %                 % that, otherwise I have an infinite while loop
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 seqOfEventsNew(rowsSegment,:)=[];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 clear rowsSegment
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % check for the segments that are part of this group which
% % % % % % % % % % % % % % % % %                 % are the other segments that they are interacting
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 % initiate the flag for the segments interacting with the
% % % % % % % % % % % % % % % % %                 % segment but were already cheched in the previous loop.
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 iSegChecked=[];
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 while~isempty(segmentNumber)
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     %get the first seg
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     iSeg = segmentNumber(1);
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     % this analysis need to be done only for the segments
% % % % % % % % % % % % % % % % %                     % that are not segIndex
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     if iSeg~=segIndex
% % % % % % % % % % % % % % % % %                         if iSeg~=segCheck
% % % % % % % % % % % % % % % % %                             % determine in which rows this segment appears
% % % % % % % % % % % % % % % % %                             rowsSegment= seqOfEventsNew(:,3)==iSeg|seqOfEventsNew(:,4)==iSeg;
% % % % % % % % % % % % % % % % %                             
% % % % % % % % % % % % % % % % %                             % determine the number of segments that are interacting
% % % % % % % % % % % % % % % % %                             segmentNumberNew= unique (seqOfEventsNew(rowsSegment,3:4));
% % % % % % % % % % % % % % % % %                             %removing NaNs
% % % % % % % % % % % % % % % % %                             segmentNumberNew(isnan(segmentNumberNew))=[];
% % % % % % % % % % % % % % % % %                             
% % % % % % % % % % % % % % % % %                             %update segment number list
% % % % % % % % % % % % % % % % %                             segExtra=setdiff(segmentNumberNew,segmentNumber);
% % % % % % % % % % % % % % % % %                             segmentNumber=[segmentNumber;segExtra];
% % % % % % % % % % % % % % % % %                             
% % % % % % % % % % % % % % % % %                             % if this segment is still not part of the compTrack
% % % % % % % % % % % % % % % % %                             if any(~ismember(segmentNumber,compTrackSeg))
% % % % % % % % % % % % % % % % %                                 
% % % % % % % % % % % % % % % % %                                 
% % % % % % % % % % % % % % % % %                                 % save this segments to group this segments in the compTrack
% % % % % % % % % % % % % % % % %                                 compTrackSeg=[compTrackSeg;segmentNumber];
% % % % % % % % % % % % % % % % %                                 
% % % % % % % % % % % % % % % % %                             end
% % % % % % % % % % % % % % % % %                         end
% % % % % % % % % % % % % % % % %                     end
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     %save the segments that where already checked from the
% % % % % % % % % % % % % % % % %                     %list of segments
% % % % % % % % % % % % % % % % %                     segCheck=[segCheck;setdiff(iSeg,segCheck )];
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     % some segments can interact multiple times with the same segment. so
% % % % % % % % % % % % % % % % %                     % update compTracks
% % % % % % % % % % % % % % % % %                     compTrackSegCell{segIndex}=unique(compTrackSeg);
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     %remove this segment from the list
% % % % % % % % % % % % % % % % %                     segmentNumber=segmentNumber(2:end);
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                     
% % % % % % % % % % % % % % % % %                 end
% % % % % % % % % % % % % % % % %             end
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         % remove the empty cells
% % % % % % % % % % % % % % % % %         compTrackSegCell=compTrackSegCell(~cellfun('isempty',compTrackSegCell));
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         %% reform compTrack
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         %load tracksFeatIndxCG and tracksCoordAmpCG
% % % % % % % % % % % % % % % % %         tracksFeatIndxCG= compTracks(iTrack).tracksFeatIndxCG;
% % % % % % % % % % % % % % % % %         tracksCoordAmpCG= compTracks(iTrack).tracksCoordAmpCG;
% % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % %         for indexGroup=1:length(compTrackSegCell)
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             %load the segments that are part of goup 1
% % % % % % % % % % % % % % % % %             segmentGroup=compTrackSegCell{indexGroup};
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % the group that have segment 1 will stay in the original
% % % % % % % % % % % % % % % % %             % position in compTrack, the other will be add at the end of the
% % % % % % % % % % % % % % % % %             % compTrack
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % check if the group have the first segment
% % % % % % % % % % % % % % % % %             if any (segmentGroup==1)
% % % % % % % % % % % % % % % % %                 iTrackCurr= iTrack;
% % % % % % % % % % % % % % % % %             else
% % % % % % % % % % % % % % % % %                 %increase iTrackNew
% % % % % % % % % % % % % % % % %                 iTrackNew=iTrackNew+1;
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %                 iTrackCurr= iTrackNew;
% % % % % % % % % % % % % % % % %                 
% % % % % % % % % % % % % % % % %             end
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             %% reform tracksFeatIndxCG and tracksCoordAmpCG
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % take only the segments that are part of the group segmentGroup
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             tracksReform( iTrackCurr).tracksFeatIndxCG = tracksFeatIndxCG(segmentGroup,:);
% % % % % % % % % % % % % % % % %             tracksReform( iTrackCurr).tracksCoordAmpCG =tracksCoordAmpCG(segmentGroup,:);
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             %% reform seqOfevents
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % determine which are the rows that remain for this group
% % % % % % % % % % % % % % % % %             rowsGood= ismember(seqOfEventsArtifacts(:,3),segmentGroup)|ismember(seqOfEventsArtifacts(:,4),segmentGroup);
% % % % % % % % % % % % % % % % %             tracksReform(iTrackCurr).seqOfEvents=seqOfEventsArtifacts(rowsGood,:);
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % replace number of segments
% % % % % % % % % % % % % % % % %             if sum(rowsGood)==2
% % % % % % % % % % % % % % % % %                 tracksReform( iTrackCurr).seqOfEvents(:,3)=1;
% % % % % % % % % % % % % % % % %             else
% % % % % % % % % % % % % % % % %                 seqOfEventsIndex= tracksReform(iTrackCurr).seqOfEvents(:,3:4);
% % % % % % % % % % % % % % % % %                 for i=1:length(segmentGroup)
% % % % % % % % % % % % % % % % %                     seqOfEventsIndex(seqOfEventsIndex==segmentGroup(i))=i;
% % % % % % % % % % % % % % % % %                 end
% % % % % % % % % % % % % % % % %                 %update seqOfEvents
% % % % % % % % % % % % % % % % %                 tracksReform( iTrackCurr).seqOfEvents(:,3:4)=seqOfEventsIndex;
% % % % % % % % % % % % % % % % %             end
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             %% reform tracksFeatIndxCG and tracksCoordAmpCG to give only information for the frames that they appear
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % with the reform some tracks can start now at some point
% % % % % % % % % % % % % % % % %             % different of frame 1, but we need the information of tracks
% % % % % % % % % % % % % % % % %             % only after they appear and until they disappear
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             %calculate if there is difference between the time that the
% % % % % % % % % % % % % % % % %             %original and the new track initiates or finishes
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % calculate the initiation times
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             timeStartOri=compTracks(iTrack).seqOfEvents(1,1);
% % % % % % % % % % % % % % % % %             timeStartNew=tracksReform( iTrackCurr).seqOfEvents(1,1);
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % calculate the end times
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             timeEndOri=compTracks(iTrack).seqOfEvents(end,1);
% % % % % % % % % % % % % % % % %             timeEndNew=tracksReform( iTrackCurr).seqOfEvents(end,1);
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % calculate the difference between these times
% % % % % % % % % % % % % % % % %             timeDiffStart= timeStartNew-timeStartOri;
% % % % % % % % % % % % % % % % %             timeDiffEnd= timeEndOri-timeEndNew;
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % reform tracksFeatIndxCG and tracksCoordAmpCG
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             tracksReform( iTrackCurr).tracksFeatIndxCG =tracksReform( iTrackCurr).tracksFeatIndxCG (:,timeDiffStart+1:end-timeDiffEnd);
% % % % % % % % % % % % % % % % %             tracksReform( iTrackCurr).tracksCoordAmpCG =tracksReform( iTrackCurr).tracksCoordAmpCG(:,8*timeDiffStart+1:end-8*timeDiffEnd);
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %             % add a new structure field to the compTracks, that will take
% % % % % % % % % % % % % % % % %             % account of all the alterations from the original compTrack
% % % % % % % % % % % % % % % % %             tracksReform(iTrackCurr).oldTracksInfo={iTrack,segmentGroup};
% % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % %     else
% % % % % % % % % % % % % % % % %         % It was only one track, just copy this info in the field oldTracksInfo
% % % % % % % % % % % % % % % % %         tracksReform(iTrack).oldTracksInfo={iTrack,1};
% % % % % % % % % % % % % % % % %     end
    end
end
    end
end

 %% regroup segments
 tracksReform = reformCompTracksAfterRemoveArtifacts(tracksReform);

