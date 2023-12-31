function [ seqOfEventsReform, clustHistoryReform] = reformSeqOfEventsFromMerge2split(seqOfEvents,segmentMerge,segmentSplit,rowMergeSplit,clustHistory)
%reformSeqOfEventsFromMerge2split reform seqOfEvents and clustHistory after
%removing artifacts.
%
% [ seqOfEvents, clustHistory] = reformSeqOfEventsFromMerge2split(seqOfEvents,segmentMerge,segmentSplit,rowMergeSplit,clustHistory)
%
% INPUT
%       seqOfEvents   :        Matrix with number of rows equal to number
%                              of events happening in a compound track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 = start of track segment, 2 = end of track segment;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN = start is a birth and end is a death,
%                              number = start is due to a split, end
%                              is due to a merge, number is the index
%                              of track segment for the merge/split.
%
%       segmentMerge  :        segment that is merging in the merge event to
%                              the contineous segment
%
%       segmentSplit  :        segment that is spliting in the split event
%                              from the contineous segment
%
%       rowMergeSplit :        row numbers of segment merge and split that
%                              need to be removed from seqOfEvents
%
%       clustHistory  :        clustHistoryAll:  a 1D cell with rows = number of tracks in
%                              compTracks.
%                              Each entry contains a clusterHistory table for a
%                              track in compTracks. In this version cluster
%                              history is recorded for all clusters in the
%                              obervation time.
%                              clusterHistory is a 2D array with each row
%                              corresponding to an association or a dissociation
%                              event. The 9 colums give the following information:
%                              1) Track segment number.
%                              2) Cluster size.
%                              3) Start time (same units as input timeStep etc).
%                              4) End time (same units as input timeStep etc).
%                              5) Lifetime (same units as input timeStep etc).
%                              6) Event that started the cluster
%                              (1 = dissociation, 2 = association).
%                              7) Event that ended the cluster
%                              (1 = dissociation, 2 = association).
%                              8) Resulting cluster size.
%                              9) Association flag - 1 indicates the segment
%                              and its partner are both listed and NaN
%                              indicates only the current segment is listed,
%                              i.e. the partner is not listed.
%
%
% OUTPUT
%       seqOfEventsReform   :        reform SeqOfEvents
%
%       clustHistoryReform  :        reform clustHistory
%
% This funciton has the oprional to reform only seqOfEvents or both seqOfEvents and clustHistory 
%
% Luciana de Oliveira, March 2018
%
%%

% check if clustHistory is an input and replace segmentSplit index by
% segmentMerge index to have only one long segment

if nargin==5 % if clustHistory is an input
    
    clustHistory(clustHistory(:,1)==segmentSplit)=segmentMerge;
end

% same thing for seqOfEvents

seqOfEvents(seqOfEvents(:,3)==segmentSplit,3)=segmentMerge;
seqOfEvents(seqOfEvents(:,4)==segmentSplit,4)=segmentMerge;


%remove event from seqOfEvents

seqOfEvents(rowMergeSplit,:)=[];

% for all segments with index larger than removed
% segment, reduce their index by 1

% take only 3trd and 4th coluns from seqOfEvents
seqOfEventsIndex= seqOfEvents(:,3:4); 

% check the list of segments index from the orinal seqOfEvents, since we
% remove segments, we need to update the order of segment index
segmentNumbersOld=unique(seqOfEventsIndex);
segmentNumbersOld(isnan(segmentNumbersOld))=[];

% check if there are discrepancies between the segment index and the
% segment number it should have
for i=1:length(segmentNumbersOld)
    if segmentNumbersOld(i)~=i
        % replace segment index by the right value
        seqOfEventsIndex(seqOfEventsIndex==segmentNumbersOld(i))=i;
        
        if nargin==5 % if clustHistory is an input, correct also the segment index in that.
            clustHistory(clustHistory(:,1)==segmentNumbersOld(i))=i;
        end
    end
end

%update seqOfEvents
seqOfEvents(:,3:4)=seqOfEventsIndex;

%% output

seqOfEventsReform=seqOfEvents;

if nargin==5 % if clustHistory is an input
    clustHistoryReform=clustHistory;
end

if nargout > 0 && ~exist('clustHistoryReform','var') % if clustHistory is an input/output
    clustHistoryReform = [];
end
end

