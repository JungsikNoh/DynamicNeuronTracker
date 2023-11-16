function [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,segmentMerge,segmentSplit)
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
%        segmentMerge     : segment that is merging in the merge event to
%        the contineous segment
%        segmentSplit     : segment that is spliting in the split event
%        from the contineous segment
% OUTPUT
%          tracksReform     : tracks after pairwesing the merges and the
%          splits, creating longer segments
%
%
% Luciana de Oliveira, January 2018 
   %%             
              
   
                % take the rows for the segment that merge and that split
                tracksFeatIndxCGMerge=tracksFeatIndxCG(segmentMerge,:);
                tracksFeatIndxCGSplit=tracksFeatIndxCG(segmentSplit,:);
                
                tracksCoordAmpCGMerge=tracksCoordAmpCG(segmentMerge,:);
                tracksCoordAmpCGSplit=tracksCoordAmpCG(segmentSplit,:);
                
                % put together the two piece of tracks, that is, take the non nan
                % info from the segment that split and replace it in the segment that
                % merge.
                                              
                tracksFeatIndxCGMerge=max(tracksFeatIndxCGMerge,tracksFeatIndxCGSplit);
                tracksCoordAmpCGMerge=max(tracksCoordAmpCGMerge,tracksCoordAmpCGSplit);
                
                % replace the merge segment for the new combined one and remove the
                % split segment
                tracksFeatIndxCG(segmentMerge,:)=tracksFeatIndxCGMerge;
                tracksCoordAmpCG(segmentMerge,:)= tracksCoordAmpCGMerge;
                
                % remove split segment
                tracksFeatIndxCG(segmentSplit,:)=[];
                tracksCoordAmpCG(segmentSplit,:)= [];
                
               
        
end

