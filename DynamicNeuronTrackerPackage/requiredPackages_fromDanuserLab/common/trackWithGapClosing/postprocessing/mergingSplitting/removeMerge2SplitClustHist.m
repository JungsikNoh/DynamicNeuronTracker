function  tracksReform = removeMerge2SplitClustHist(clustHistoryAll,tracksIn,probArtifactMerge2split)
% REMOVEMERGE2SPLITCLUSTHIST remove potentially artifactual merges and
%                           splits, using the information of clustHistory
%                           and the probability to be an artifact
%                           calculated for all the outputs.
%
% SYNOPSIS
% seqOfEvents = removeMerge2SplitClustHist(clustHitory,seqOfEvents,probArtifactMerge2split)
%
%
%        INPUT
%               tracksIn       : Output of trackCloseGapsKalman.
%
%
%       probArtifactMerge2split : probability to be an artifact calculated
%                                 considering time merge to split distribution
%                                 coming from all the observations.
%                                 The probability is calculated using the fit
%                                 exponetial giving the
%                                 distribution of real interactions and the
%                                 possible artifacts.
%
%
%OUTPUT
%          tracksReform     : tracks after removing artifacts coming from
%          the merge to split events.
%
%
% Luciana de Oliveira, December 2017.
% Modified February 2018.
%%

% load compTracks
compTracks=tracksIn;
%initiate the reform compTrack
tracksReform=tracksIn;


%% check where merges to splits are happening

 for iTrack = 1: length(compTracks);
% for iTrack = 2;
    
%     fprintf('\nProcessing track=%s ',int2str(iTrack));
    
    % load seqOfEvents
    
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    
    %For experimental data, the first iteration (frame) in seqOfEvents
    %can be some value other than 1.  This is occurs becuase an object
    %can appear or disappear at any time while in simulation, all of
    %the objects are present at frame 1. This affects how aggregState
    %is accessed since the frame # is supposed to correspond to a
    %column in aggregState. So, need to shift frame numbers in
    %seqOfEvents, if necessary, when accessing aggregState. The very
    %first frame # is the shift amount.
    firstFrame= seqOfEvents(1,1);
    
    
    
    % load tracksFeatIndxCG and tracksCoordAmpCG
    
    tracksFeatIndxCG= tracksReform(iTrack).tracksFeatIndxCG;
    tracksCoordAmpCG= tracksReform(iTrack).tracksCoordAmpCG;
    
    % load clustHistory
    clustHistory = clustHistoryAll{iTrack,1};
    %calculate length seqOfEvents
    
    lengthSeqOfEvents= size(seqOfEvents,1);
    
    % if the seqOfEvents have more than one segment
    
    if lengthSeqOfEvents > 2
        
        % identify in which rows merge to split events are happening
        
        indexMerge2Split=find(clustHistory(:,6)==2 & clustHistory(:,7)==1);
        
        % go over all the merge to split events
        
        for indexTime=1:length(indexMerge2Split)
            
            % take info from clustHistory
            
            
            timeMerge=clustHistory(indexMerge2Split(indexTime),3);
            timeSplit=clustHistory(indexMerge2Split(indexTime),4);
            lifetime= clustHistory(indexMerge2Split(indexTime),5);
            segmentContinuous=clustHistory(indexMerge2Split(indexTime),1);
            
            % calculate random number to test if this observation can be an artifact
            
            randNumberArtifacts=rand;
            
            %% decide if the observation is an artifact, based on the probability
            
            if randNumberArtifacts<probArtifactMerge2split(lifetime)
                
                
                % identify the rows where the merge is happining
                
                iMerge= seqOfEvents(:,1)==timeMerge & seqOfEvents(:,2)==2 & seqOfEvents(:,4)==segmentContinuous;
                
                
                %check where these two segments split
                
                iSplit = seqOfEvents(:,4)==segmentContinuous & seqOfEvents(:,2)==1&seqOfEvents(:,1)==timeSplit;
                
                % merge and split segment number
                
                segmentMerge=seqOfEvents(iMerge,3);
                segmentSplit=seqOfEvents(iSplit,3);
                
                %determine merge and split rows
                
                rowMerge=find(iMerge);
                rowSplit=find(iSplit);
                rowMergeSplit=[rowMerge,rowSplit];
                %% pairwise for merge and split
                
                
                % copy the positions in the merge to split time to the merge segment
                % for the intensity need to calculate the fraction of intensity for each
                % segment. For that need to take the intensity before the merge, but if
                % some event happened before the event of merge, need to take only the
                % part of the segment that comes before the merge event.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                timeBefAftEvent  = timeBeforeAfterEvent(seqOfEvents,segmentMerge,segmentSplit,timeMerge,segmentContinuous,timeSplit,firstFrame);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %calculate the intensities before the the merge
                meanIntSegMerge =  nanmean(tracksCoordAmpCG(segmentMerge,[8*(timeBefAftEvent(1)-firstFrame)+4:8:8*(timeMerge-firstFrame+1),8*(timeSplit-firstFrame)+4:8:8*(timeBefAftEvent(2)-firstFrame+1)]));
                meanIntSegCont  =  nanmean(tracksCoordAmpCG(segmentContinuous,[8*(timeBefAftEvent(3)-firstFrame)+4:8:8*(timeMerge-firstFrame+1),8*(timeBefAftEvent-firstFrame)+4:8:8*(timeBefAftEvent(4)-firstFrame+1)]));
                
                % distribute the intensity
                intMergeProp = meanIntSegMerge/(meanIntSegMerge+meanIntSegCont);
                intSpliProp  = meanIntSegCont/(meanIntSegMerge+meanIntSegCont);
                
                %first replace the whole matrix
                % % % % % % % % % % % % % % % % % % % % % % %  check the time intervals
                tracksFeatIndxCG(segmentMerge,timeMerge-firstFrame+1:timeSplit-firstFrame)         = tracksFeatIndxCG(segmentContinuous,timeMerge-firstFrame+1:timeSplit-firstFrame);
                tracksCoordAmpCG(segmentMerge,8*(timeMerge-firstFrame)+1:8*(timeSplit-firstFrame)) = tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame)+1:8*(timeSplit-firstFrame));
                
                % then replace the specific value of intensity
                
                tracksCoordAmpCG(segmentMerge,8*(timeMerge-firstFrame+1)+4:8:8*(timeSplit-firstFrame)+4)      = intMergeProp.*tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame+1)+4:8:8*(timeSplit-firstFrame)+4);
                tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame+1)+4:8:8*(timeSplit-firstFrame)+4) = intSpliProp.*tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame+1)+4:8:8*(timeSplit-firstFrame)+4);
                
                
                
                % paiwise tracks
                [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,segmentMerge,segmentSplit);
                
                
                
                %% reform seqOfEvents and tracks
                
                
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 % replace the segment number to have only one continuos
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 % segment
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 seqOfEvents(seqOfEvents(:,3)==segmentSplit,3)=segmentMerge;
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 seqOfEvents(seqOfEvents(:,4)==segmentSplit,4)=segmentMerge;
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 clustHistory(clustHistory(:,1)==segmentSplit)=segmentMerge;
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 %remove the merge and the split
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 seqOfEvents([rowMerge,rowSplit],:)=[];
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 % for all segments with index larger than removed
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 % segment, reduce their index by 1
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 seqOfEventsIndex= seqOfEvents(:,3:4);
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 segmentNumbersOld=unique(seqOfEventsIndex);
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 segmentNumbersOld(isnan(segmentNumbersOld))=[];
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 for i=1:length(segmentNumbersOld)
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                     if segmentNumbersOld(i)~=i
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         seqOfEventsIndex(seqOfEventsIndex==segmentNumbersOld(i))=i;
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         clustHistory(clustHistory(:,1)==segmentNumbersOld(i))=i;
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                     end
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 end
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 %update seqOfEvents
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                 seqOfEvents(:,3:4)=seqOfEventsIndex;
                
                [ seqOfEvents, clustHistory] = reformSeqOfEventsFromMerge2split(seqOfEvents,segmentMerge,segmentSplit,rowMergeSplit,clustHistory);
                
                
                %% save tracks in the final format
                tracksReform(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
                tracksReform(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
                tracksReform(iTrack).seqOfEvents = seqOfEvents;
                %                     end
            end
        end
    end
    
end
%% regroup tracks

%         if reformTackFlag
tracksReform= reformCompTracksAfterRemoveArtifacts(tracksReform);
%         end

end