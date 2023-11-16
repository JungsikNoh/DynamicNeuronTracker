function [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracksMerge2Split( tracksFeatIndxCG,tracksCoordAmpCG,segmentMerge,segmentSplit,segmentContinuous, timeSplit)
% This function pairwise the segments from simultaneos merges and splits
% artifacs
% SYNOPSIS
%
% tracksReform = pairwiseSegTracks( tracksFeatIndxCG,tracksCoordAmpCG, segmentMerge,segmentSplit)
%
% INPUT
%        tracksFeatIndxCG : tracks intensity coming from the specific track
%        number
%        tracksCoordAmpCG : tracks coordinates coming from the specific track
%        number
%        segmentMerge     : segment that is merging in the simultaneous event
%        segmentSplit     : segment that is spliting in the simultaneous event 
%        segmentContinuous: segment that is continued in the merge to split
%        events
%        timeSplit        : time where the split occours   
% OUTPUT
%          tracksReform     : tracks after pairwesing the merges and the
%          splits, creating longer segments
%
%
% Luciana de Oliveira, January 2018 
   %%             
              
                
                % take the rows for the segment that need to be linked
                tracksFeatIndxCGMerge=tracksFeatIndxCG(segmentMerge,:);
                tracksFeatIndxCGContinuous=tracksFeatIndxCG(segmentContinuous,:);
                tracksFeatIndxCGSplit=tracksFeatIndxCG(segmentSplit,:);
               
                
                tracksCoordAmpCGMerge=tracksCoordAmpCG(segmentMerge,:);
                tracksCoordAmpCGContinuous=tracksCoordAmpCG(segmentContinuous,:);
                tracksCoordAmpCGSplit=tracksCoordAmpCG(segmentSplit,:);
                
                % put together the two pieces of tracks
                tracksCoordAmpCGMerge=[tracksCoordAmpCGMerge(1:timeSplit-1),tracksCoordAmpCGContinuous(timeSplit:end)];
                tracksFeatIndxCGMerge=[tracksFeatIndxCGMerge(1:timeSplit-1),tracksFeatIndxCGContinuous(timeSplit:end)];
                
                tracksCoordAmpCGContinuous=[tracksCoordAmpCGContinuous(1:timeSplit-1),tracksCoordAmpCGSplit(timeSplit:end)];
                tracksFeatIndxCGContinuous=[tracksFeatIndxCGContinuous(1:timeSplit-1),tracksFeatIndxCGSplit(timeSplit:end)];
                
                                
                % replace the new segments and remove the
                % split segment
                
                tracksFeatIndxCG(segmentMerge,:)=tracksFeatIndxCGMerge;
                tracksCoordAmpCG(segmentMerge,:)= tracksCoordAmpCGMerge;
                
                tracksFeatIndxCG(segmentContinuous,:)=tracksCoordAmpCGContinuous;
                tracksCoordAmpCG(segmentContinuous,:)=tracksFeatIndxCGContinuous;
                
                % remove split segment
                tracksFeatIndxCG(segmentSplit,:)=[];
                tracksCoordAmpCG(segmentSplit,:)= [];
                
               
          

end

