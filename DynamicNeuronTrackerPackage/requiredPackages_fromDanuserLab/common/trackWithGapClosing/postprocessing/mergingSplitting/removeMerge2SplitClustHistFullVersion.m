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
    
    fprintf('\nProcessing track=%s ',int2str(iTrack));
    
    % load seqOfEvents
    
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    
    reformTracksFlag=0;
    
    
    %For experimental data, the first iteration (frame) in seqOfEvents
    %can be some value other than 1.  This is occurs becuase an object
    %can appear or disappear at any time while in simulation, all of
    %the objects are present at frame 1. This affects how aggregState
    %is accessed since the frame # is supposed to correspond to a
    %column in aggregState. So, need to shift frame numbers in
    %seqOfEvents, if necessary, when accessing aggregState. The very
    %first frame # is the shift amount.
    %     firstFrame= seqOfEvents(1,1);
    
    %load tracksFeatIndxCG and tracksCoordAmpCG
    
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
            
            
            timeMerge=round(clustHistory(indexMerge2Split(indexTime),3));
            timeSplit=round(clustHistory(indexMerge2Split(indexTime),4));
            lifetime= round(clustHistory(indexMerge2Split(indexTime),5));
            segmentContinuous= round(clustHistory(indexMerge2Split(indexTime),1));
            
            
            % % % % % % % % % % % % % % % % % %             %%%%%%%%%% ?????????????????????????????
            % % % % % % % % % % % % % % % % % %             if length(segmIndex)>2
            % % % % % % % % % % % % % % % % % %                 segmIndexNew=find( seqOfEvents(:,1)==timeSplit&seqOfEvents(:,2)==1);
            % % % % % % % % % % % % % % % % % %
            % % % % % % % % % % % % % % % % % %                 segmIndex=intersect(segmIndexNew,segmIndex);
            % % % % % % % % % % % % % % % % % %             end
            %             ???????????????????????????????????????????????????
            % how to avoid ambiguits when calculating the segment continuous. maybe I
            % should already calculate the segment merge and split and try to find
            % where these are interacting, or try to find this info in the clust
            % history
            
            
            % % % % % % % % % % % %              if sum(segmIndex)~=0
            
            %             segmentContinuous=seqOfEvents(segmIndex,4); % this is the segment of segment merge with and from segment split from
            
            %%%%%%%%%%%%%%%% solve the problem with the empty segments
            
            
            %             if sum(segmIndex) > 1 % for large seqOfEvents if could have multiple events happening in the same frame,
            %                 %this if statment is to prevent these problems.
            %
            %                 segmIndex= seqOfEvents(:,1)==timeSplit&seqOfEvents(:,2)==1;
            %                 segmentNumberSplit=seqOfEvents(segmIndex,4);
            %
            %                 segmentContinuous=intersect(segmentContinuous,segmentNumberSplit);
            %
            %             end
            % calculate random number to test if this observation can be an artifact
            
            randNumberArtifacts=rand;
            
            %% decide if the observation is an artifact, based on the probability
            
            if randNumberArtifacts<probArtifactMerge2split(lifetime)
                
                reformTracksFlag=1; %flag for the reform tracks and regroup
                
                % identify the rows where the merge is happining
                
                iMerge= seqOfEvents(:,1)==timeMerge & seqOfEvents(:,2)==2 & seqOfEvents(:,4)==segmentContinuous;
                
                
                %check where these two segments split
                
                iSplit = seqOfEvents(:,4)==segmentContinuous & seqOfEvents(:,2)==1&seqOfEvents(:,1)==timeSplit;
                
                % merge and split segment number
                
                segmentMerge=seqOfEvents(iMerge,3);
                segmentSplit=seqOfEvents(iSplit,3);
                
                
                %% pairwise for merge and split
                %                     % test for the intensity and moviment of segment merge,
                %                     % continuos and split
                %
                %                     % pre-alocate the ratios for the calculation of
                %                     % costMatrix
                %                     % for this case will have only 2 possibilities of merges and 2 possibilities of splits
                %                     ratioIntensity=zeros(1,2);
                %                     ratioMeanDisp=zeros(1,2);
                
                % calculate all the possible combinations between merges and splits
                
                %                 [mergeIndices,splitIndices] = meshgrid([1,2],[1,2]);
                
                % Concatenate indices
                %                 concatenateIndices=cat(2,mergeIndices',splitIndices');
                % final combination of all merges and splits
                % indices
                %                 indicesMS=reshape(concatenateIndices,[],2);
                
                
                % copy the positions in the merge to split time to the merge segment
                
                
                tracksFeatIndxCG(segmentMerge,timeMerge:timeSplit) = tracksFeatIndxCG(segmentContinuous,timeMerge:timeSplit-1);
                tracksCoordAmpCG(segmentMerge,timeMerge:timeSplit) = tracksCoordAmpCG(segmentContinuous,timeMerge:timeSplit-1);
                
                % paiwise tracks
                [tracksFeatIndxCGTmp,tracksCoordAmpCGTmp] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,segmentMerge,segmentSplit);
                
                % calculate mean intensity and displacements for
                % segment merge, split and continuous
                
                %                     %% intensity
                %
                %                     % segment continuous before merging
                %
                %                     meanIntSegContinuousBefore=nanmean(tracksCoordAmpCG(segmentContinuous,4:8:8*(timeMerge-firstFrame+1)-1),2);
                %
                %                     % segment merge
                %                     %                   meanIntSegMerge=nanmean(tracksCoordAmpCG(segmentMerge,4:8:8*(timeMerge-firstFrame+1)-1),2);
                %
                %                     % segment continuous after splitting
                %
                %                     meanIntSegContinuousAfter=nanmean(tracksCoordAmpCG(segmentContinuous,8*(timeSplit-firstFrame+1)+4:8:end),2);
                %
                %                     % segment split
                %
                %                     meanIntSegSplit=nanmean(tracksCoordAmpCG(segmentSplit,8*(timeSplit-firstFrame+1)+4:8:end),2);
                %
                %
                %                     %% calculate displacements
                %
                %                     %take the x and y values
                %
                %                     % segment continuous before merging
                %
                %                     segContinuousBeforeXvals=tracksCoordAmpCG(segmentMerge,1:8:8*(timeMerge-firstFrame+1)-1);
                %                     segContinuousBeforeYvals=tracksCoordAmpCG(segmentMerge,2:8:8*(timeMerge-firstFrame+1)-1);
                %
                %                     %calculate displacement
                %
                %                     dispSegContinuousBeforeX=segContinuousBeforeXvals(:,2:end)-segContinuousBeforeXvals(:,1:end-1);
                %                     dispSegContinuousBeforeY=segContinuousBeforeYvals(:,2:end)-segContinuousBeforeYvals(:,1:end-1);
                %
                %
                %                     % calculate mean displacement
                %
                %                     meanDispSegContinuousBefore=nanmean(sqrt(dispSegContinuousBeforeX.^2+dispSegContinuousBeforeY.^2),2);
                %
                %                     % segment merge
                %                     %
                %                     %                     segMergeXvals=tracksCoordAmpCG(segmentMerge,1:8:8*(timeMerge-firstFrame+1)-1);
                %                     %                     segMergeYvals=tracksCoordAmpCG(segmentMerge,2:8:8*(timeMerge-firstFrame+1)-1);
                %                     %
                %                     %calculate displacement
                %
                %                     %                     dispSegMergeX=segMergeXvals(:,2:end)-segMergeXvals(:,1:end-1);
                %                     %                     dispSegMergeY=segMergeYvals(:,2:end)-segMergeYvals(:,1:end-1);
                %
                %
                %                     % calculate mean displacement
                %
                %                     %                     meanDispMerge=nanmean(sqrt(dispSegMergeX.^2+dispSegMergeY.^2),2);
                %
                %                     % segment continuous after splitting
                %
                %                     segContinuousAfterXvals=tracksCoordAmpCG(segmentContinuous,8*(timeSplit-firstFrame+1)+1:8:end);
                %                     segContinuousAfterYvals=tracksCoordAmpCG(segmentContinuous,8*(timeSplit-firstFrame+1)+2:8:end);
                %
                %                     %calculate displacement
                %
                %                     dispSegContinuousAfterX=segContinuousAfterXvals(:,2:end)-segContinuousAfterXvals(:,1:end-1);
                %                     dispSegContinuousAfterY=segContinuousAfterYvals(:,2:end)-segContinuousAfterYvals(:,1:end-1);
                %
                %
                %                     % calculate mean displacement
                %
                %                     meanDispSegContinuousAfter=nanmean(sqrt(dispSegContinuousAfterX.^2+dispSegContinuousAfterY.^2),2);
                %
                %                     % segment split
                %
                %                     segSplitXvals=tracksCoordAmpCG(segmentSplit,8*(timeSplit-firstFrame+1)+1:8:end);
                %                     segSplitYvals=tracksCoordAmpCG(segmentSplit,8*(timeSplit-firstFrame+1)+2:8:end);
                %
                %
                %                     %calculate displacement
                %
                %                     dispSegSplitX=segSplitXvals(:,2:end)-segSplitXvals(:,1:end-1);
                %                     dispSegSplitY=segSplitYvals(:,2:end)-segSplitYvals(:,1:end-1);
                %
                %
                %                     % calculate mean displacement
                %
                %                     meanDispSegSplit=nanmean(sqrt(dispSegSplitX.^2+dispSegSplitY.^2),2);
                %
                %
                %                     %% combine intensity and displacements
                %
                %                     intensityMean=[meanIntSegContinuousAfter, meanIntSegSplit];
                %                     displacementMean=[meanDispSegContinuousAfter,meanDispSegSplit];
                %
                %
                %
                %
                %                     for indexMS=1:2
                %
                %                         %% calculate ratio
                %
                %                         % check which mean have the largest value
                %
                %                         %intensity
                %
                %                         intLarger=max(meanIntSegContinuousBefore,intensityMean(indexMS));
                %                         intSmaller=min(meanIntSegContinuousBefore,intensityMean(indexMS));
                %
                %                         %displacement
                %
                %                         dispLarger=max(meanDispSegContinuousBefore,displacementMean(indexMS));
                %                         dispSmaller=min(meanDispSegContinuousBefore,displacementMean(indexMS));
                %
                %                         %calculate ratio
                %                         ratioIntensity(1,indexMS)=intLarger/intSmaller;
                %                         ratioMeanDisp(1,indexMS)=dispLarger/dispSmaller;
                %
                %                     end
                
                
                % for some cases where the split occour and it is
                % followed by a end we have ratios with value NaN, and
                % NaN is not a valid input for the LAP function. To
                % avoid this problem, in the cases where we have
                % multiple merges and split, we will take the maximum
                % value of cost matrix and apply for the NaN value
                
                
                %
                %                     if any(any(isnan(ratioMeanDisp)))||all(all(isnan (ratioMeanDisp)))
                %
                %
                %                         if  all(all(isnan (ratioMeanDisp))) % if any ratio displacement is NaN
                %
                %                             %if all the displacements are NaN we will take only
                %                             %the value of intensity for the comparison
                %
                %                             ratioMeanDisp(:)=1;
                %
                %                         else
                %
                %
                %                             indiceNaNDisp= isnan(ratioMeanDisp);
                %                             newValRatioDisp=max(max(ratioMeanDisp));
                %                             ratioMeanDisp(indiceNaNDisp)=newValRatioDisp;
                %
                %                         end
                %
                %
                %                     end
                %
                %                     % calculate the costMatrix
                %
                %                     costMatrix=ratioIntensity.*ratioMeanDisp;
                %
                %
                %                     % after all the combinations are calculated and saved in costMatrix, we can
                %                     % analyse which is the most likely to be a pair of segment
                %                     % before and after merge
                %                     [beforeMergePositions, afterSplitPositions] = lap(costMatrix, [],[],1,[]);
                %
                %
                %                     % Take only the possible pairs, as we are comparing only
                %                     % segment continuous, it is known that the only possible
                %                     % combination is 1 with one of the other part segments after
                %                     % splitting, so what is need to find is the position of 1 in
                %                     % the list of mergePositions
                %
                %                     goodPairs=find(beforeMergePositions==1);
                %
                %                     pairsFinal=[beforeMergePositions(goodPairs),afterSplitPositions(goodPairs)];
                %
                % define which are the combinations of pairs in terms of
                % segment number
                
                
                
                %% reform seqOfEvents and tracks
                
                %                     %find rows to be removed
                %                     rowMerge=find(iMerge);
                %                     rowSplit=find(iSplit);
                %
                %                     %remove the merge and the split
                %                     seqOfEvents([rowMerge,rowSplit],:)=[];
                %
                %
                %                     %
                %                     if pairsFinal(2)==1 %% segment continuous maintain the same as it is
                %
                %                         % replace the segment number to have only one continuos
                %                         % segment
                seqOfEvents(seqOfEvents(:,3)==segmentSplit,3)=segmentMerge;
                seqOfEvents(seqOfEvents(:,4)==segmentSplit,4)=segmentMerge;
                
                %                     end
                %
                
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
                
                
                
                %% save tracks in the final format
                tracksReform(iTrack).tracksFeatIndxCG = tracksFeatIndxCGTmp;
                tracksReform(iTrack).tracksCoordAmpCG = tracksCoordAmpCGTmp;
                tracksReform(iTrack).seqOfEvents = seqOfEvents;
                %                     end
            end
        end
    end
    
           
        
        %% regroup tracks
        if reformTracksFlag
        tracksReform = reformCompTracksAfterRemoveArtifacts(tracksReform);
        end
      
end



% % % %                 %% reform seqOfEvents
% % % %
% % % %                 % remove the merge and split events from the seqOfEvents
% % % %                 seqOfEvents([iMerge,iSplit],:)=[];
% % % %
% % % %                 % replace segment split number by segment merge number
% % % %
% % % %                 seqOfEvents(seqOfEvents(:,3)== segmentSplit,3)= segmentMerge;
% % % %                 seqOfEvents(seqOfEvents(:,4)== segmentSplit,4)= segmentMerge;
% % % %
% % % %                 % for all segments with index larger than removed
% % % %                 % segment, reduce their index by 1
% % % %
% % % %                 seqOfEventsIndex= seqOfEvents(:,3:4);
% % % %                 segmentNumbersOld=unique(seqOfEventsIndex);
% % % %                 segmentNumbersOld(isnan(segmentNumbersOld))=[];
% % % %                 for i=1:length(segmentNumbersOld)
% % % %                     if segmentNumbersOld(i)~=i
% % % %                         seqOfEventsIndex(seqOfEventsIndex==segmentNumbersOld(i))=i;
% % % %                     end
% % % %                 end
% % % %                 %update seqOfEvents
% % % %                 seqOfEvents(:,3:4)=seqOfEventsIndex;
% % % %
% % % %                 %% reform tracks
% % % %
% % % %                 % take the rows for the segment that merge and split
% % % %                 tracksFeatIndxCGMerge=tracksFeatIndxCG(segmentMerge,:);
% % % %                 tracksFeatIndxCGSplit=tracksFeatIndxCG(segmentSplit,:);
% % % %
% % % %                 tracksCoordAmpCGMerge=tracksCoordAmpCG(segmentMerge,:);
% % % %                 tracksCoordAmpCGSplit=tracksCoordAmpCG(segmentSplit,:);
% % % %
% % % %                 tracksFeatIndxCGMerge2Split=tracksFeatIndxCG(segmentNumber,timeMerge:timeSplit-1);
% % % %                 tracksCoordAmpCGMerge2Split=tracksCoordAmpCG(segmentNumber,8*(timeMerge)+1:8*timeSplit);%%% check that time
% % % %
% % % %                 % put together the two piece of tracks, that is, take the non nan
% % % %                 % info from the segment that split and attach it in the segment that
% % % %                 % merge.
% % % %
% % % %
% % % %                 tracksFeatIndxCGMerge(timeMerge:end)=tracksFeatIndxCGSplit((timeMerge:end));
% % % %                 tracksCoordAmpCGMerge(8*timeMerge+1:end)=tracksCoordAmpCGSplit(8*timeMerge+1:end);
% % % %
% % % %                 %replace the time merge to split by the intensity and
% % % %                 %positions of the segment those merged to
% % % %                 tracksFeatIndxCGMerge(timeMerge:timeSplit-1)=tracksFeatIndxCGMerge2Split;
% % % %                 tracksCoordAmpCGMerge(8*timeMerge+1:8*timeSplit)=tracksCoordAmpCGMerge2Split;
% % % %
% % % %
% % % %                 % update tracksFeatIndxCG and tracksCoordAmpCG
% % % %                 tracksFeatIndxCG(segmentMerge,:) = tracksFeatIndxCGMerge;
% % % %                 tracksCoordAmpCG(segmentMerge,:) = tracksCoordAmpCGMerge;
% % % %                 % remove split segment
% % % %                 tracksFeatIndxCG(segmentSplit,:)=[];
% % % %                 tracksCoordAmpCG(segmentSplit,:)= [];
% % % %
% % % %                 % update reformTracks
% % % %                 tracksReform(iTrack).tracksFeatIndxCG=tracksFeatIndxCG;
% % % %                 tracksReform(iTrack).tracksCoordAmpCG=tracksCoordAmpCG;
% % % %
% % % %             end
% % % %         end
% % % %     end
% % % % end
% % % % end
% % % %
