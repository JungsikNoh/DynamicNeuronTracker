function tracksReformat = reformatSplitToMergeEvents( tracksIn )
% This function verify if for an split followed by a merge, the longest
% segment will remain as a unique segment with a short segment interacting
% with it. It reformat the compTracks.

% SYNOPSIS
%
% tracksReform = reformatSplitToMergeEvents( tracksIn )
%
% INPUT
%        INPUT
%               tracksIn       : Output of trackCloseGapsKalman.
% OUTPUT
%          tracksReformat     : tracks after reformating to have a long
%          segment and a short one interacting with it.
%
%
% Luciana de Oliveira, February 2018

%%

%initiate the reform compTrack
tracksReformat=tracksIn;

% go over all the tracks

 for iTrack = 1: length(tracksReformat);

    %initiate variable for reform track, if it is zero, means that there is
    %no need to rearange tracks.
    reformTrack=0;
    
    %load sequence of events
    
    seqOfEvents= tracksReformat(iTrack).seqOfEvents;
    
    %correction for indexing
    seqOfEventsCorr = seqOfEvents(:,1)-seqOfEvents(1,1)+1;
    
    %find indices for which a start of a new track is not due a birth and
    % an end is not due a death
    
    eventRows = find(~isnan(seqOfEvents(:,4)));
    
    % initiate the vector for the information of segments that were modified
    
    segInden=zeros(length(eventRows),1); % segIden is the identification of the
    % segment that was already analysed in modified in seqOfEvents      
    
    % determine if the event is a merge or a split and depending on that make
    % the tests for artifacts
    
    %go over all the events
    
    for indexEvent=1:length(eventRows)
        reformTrack=0;
        % to be sure of analysing only the rows that some event is happening
        iEvent=eventRows(indexEvent);
        
        % determine if the event was already analysed
        segFlag=find(segInden==iEvent, 1);
        
        if isempty(segFlag)
            
            % determine the kind of event
            
            kindOfEvent=seqOfEvents(iEvent,2);
            
            if kindOfEvent==1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % check if the split is followed by a merge
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %find the two splitting segments and the time of splitting
                segment1 = seqOfEvents(iEvent,3);
                segment2 = seqOfEvents(iEvent,4);
                               
                
                %check whether these two segments merge with each other again
                iMerge = find(any(seqOfEvents(:,3:4)==segment1,2) & ...
                    any(seqOfEvents(:,3:4)==segment2,2) & seqOfEvents(:,2)==2);
                
                % determine time merge
              
                timeMerge=seqOfEventsCorr(iMerge,1);
                
                %if they merge ...
                if ~isempty(iMerge)
                    
                 % check if the older segment has been conserved, that is, if the newest segment is that spliting and merging
                 
                    segmentMerge1=seqOfEvents(iMerge,3);
                    segmentMerge2=seqOfEvents(iMerge,4);
                    
                    % if not, change the segment order
                    
                  if segmentMerge1<segmentMerge2
                      reformTrack=1;
                      %replace name of segment 
                     segment1=segmentMerge1;
                     segment2=segmentMerge2;
                      
                     seqOfEvents(iMerge,3)=segment2;
                     seqOfEvents(iMerge,4)=segment1;
                     
                     % replace the segment number in the other interactions
                     seqOfEventsTemp=seqOfEvents(iMerge+1:end,3:4);
%                      segmentReplace= seqOfEventsTemp==segmentMerge2;
                     
                     seqOfEventsTemp(seqOfEventsTemp==segment2)=segment1;    
                     
                     % replace in seqOfEvents
                     seqOfEvents(iMerge+1:end,3:4)=seqOfEventsTemp;
                     
                     
                  end
                 
                  if reformTrack
                  %% reformCompTracks
                
                  
                %load tracksFeatIndxCG and tracksCoordAmpCG
                tracksFeatIndxCGTmp=tracksReformat(iTrack).tracksFeatIndxCG;
                tracksCoordAmpCGTmp= tracksReformat(iTrack).tracksCoordAmpCG;
                
                               
                % take the rows for the new and old segments
                
                tracksFeatIndxCGSeg1=tracksFeatIndxCGTmp(segment1,:);
                tracksFeatIndxCGSeg2=tracksFeatIndxCGTmp(segment2,:);
                
                tracksCoordAmpCGSeg1=tracksCoordAmpCGTmp(segment1,:);
                tracksCoordAmpCGSeg2=tracksCoordAmpCGTmp(segment2,:);
                
                %from the segment new, take the time after merge and
                %replace this time interval in segment old
                
                % put together the two piece of tracks, that is, take the non nan
                % info from the segment that split and attach it in the segment that
                % merge.
                tracksFeatIndxCGSeg1(timeMerge:end)=tracksFeatIndxCGSeg2(timeMerge:end);
                tracksCoordAmpCGSeg1(8*(timeMerge-1)+1:end)=tracksCoordAmpCGSeg2(8*(timeMerge-1)+1:end);
                
                tracksFeatIndxCGSeg2(timeMerge:end)=0;
                tracksCoordAmpCGSeg2(8*(timeMerge-1)+1:end)=NaN;
                
                % replace these values in the full tracks
                tracksFeatIndxCGTmp(segment1,:)=tracksFeatIndxCGSeg1;
                tracksCoordAmpCGTmp(segment1,:)=tracksCoordAmpCGSeg1;
                
                tracksFeatIndxCGTmp(segment2,:)=tracksFeatIndxCGSeg2;
                tracksCoordAmpCGTmp(segment2,:)=tracksCoordAmpCGSeg2;
                
                % save tracks in the final format
                tracksReformat(iTrack).tracksFeatIndxCG= tracksFeatIndxCGTmp;
                tracksReformat(iTrack).tracksCoordAmpCG= tracksCoordAmpCGTmp;
                tracksReformat(iTrack).seqOfEvents=seqOfEvents;  
                  end
                end
            end
        end
    end
    
end

