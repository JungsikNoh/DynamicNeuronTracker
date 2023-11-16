function [statsGeneral,statsMotion,statsMerging,statsSplitting] = ...
    calcStatsMSAlt(tracksCell,diffAnalysisResCell,diffModeAnalysisResCell,minTrackLen,ignoreUn)
%CALCSTATSMSAlt calculates merge/split statistics for different motion types and modes
%
%SYNOPSIS [statsGeneral,statsMotion,statsMerging,statsSplitting] = ...
%    calcStatsMSAlt(tracksCell,diffAnalysisResCell,diffModeAnalysisResCell,minTrackLen)
%
%INPUT  tracksCell     : Cell array of particle tracks in alternative format, as output by
%                    convTrackFormatDefault2Alt. 
%       diffAnalysisResCell:Cell array of diffusion analysis results, as output by
%                    trackDiffusionAnalysis1 when run on the tracks in
%                    alternative format.
%       diffModeAnalysisResCell:Cell array of diffusion mode analysis results, as output by
%                    trackDiffModeAnalysis when run on the tracks in
%                    alternative format.
%       minTrackLen: Minimum length of a compound track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%
%OUTPUT
%
%
%Khuloud Jaqaman, February 2018

%% input

if nargin < 3 || isempty(tracksCell) || isempty(diffAnalysisResCell) || isempty(diffModeAnalysisResCell)
    disp('calcStatsMS: Missing tracks and/or their motion analysis!');
    return
end

if nargin < 4 || isempty(minTrackLen)
    minTrackLen = 5;
end

%% preamble
%take cells out of cell, but record size
%NOTE: Variable structure is being changed here! 
cellNum = length(tracksCell);
tracks = vertcat(tracksCell{:});
diffAnalysisRes = vertcat(diffAnalysisResCell{:});
diffModeAnalysisRes = vertcat(diffModeAnalysisResCell{:});
%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria

tracks = tracks(indx);
diffAnalysisRes = diffAnalysisRes(indx);
diffModeAnalysisRes = diffModeAnalysisRes(indx);

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1))*cellNum;

%put tracks in matrix format
[tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);

%get number of track segments
numTrackSegments = size(tracksMat,1);

%% features

%get average number of features per frame
numFeatTot = length(find(tracksIndxMat(:)));
aveFeatPerFrame = numFeatTot / numFrames;

%% track segment types

%get track segment classification from diffusion analysis
trackSegmentClass = vertcat(diffAnalysisRes.classification);

%get indices of linear/directed, free, confined, immobile and undetermined
%track segments
indxLin    = find( trackSegmentClass(:,1) == 1 | trackSegmentClass(:,2) == 3 );
indxFree   = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 2 );
indxConf   = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 1 );
indxImm    = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 0 );
indxUndet  = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) );

%store track segment type in an array
%-1 = undetermined, 0 = immobile, 1 = confined, 2 = free,
%3 = linear/directed
trackSegmentType = zeros(numTrackSegments,1);
trackSegmentType(indxUndet) = -1;
trackSegmentType(indxConf)  = 1;
trackSegmentType(indxFree)  = 2;
trackSegmentType(indxLin)   = 3;

if isempty(ignoreUn)
    %calculate number of track segments per diffusion type
numSegmentsType = [length(indxLin) length(indxFree) length(indxConf) ...
    length(indxImm) length(indxUndet)]';
else   
    %calculate number of track segments per diffusion type
numSegmentsType = [length(indxLin) length(indxFree) length(indxConf) ...
    length(indxImm)]';
end

%get track segment mode from diffusion mode analysis
trackSegmentMode = vertcat(diffModeAnalysisRes.diffMode);


trackSegmentMode(isnan(trackSegmentMode)) = -1; %replace NaN, indicating undetermined mode, with -1

%get indices of track segments in different modes
%and calculate number of track segments per diffusion mode
maxMode = max(trackSegmentMode);
indxMode = cell(maxMode+1,1);
numSegmentsMode = zeros(maxMode+1,1);
for iMode = 1  : maxMode
    indxMode{iMode} = find( trackSegmentMode == maxMode-iMode+1 );
    numSegmentsMode(iMode) = length(indxMode{iMode});
end

indxMode{maxMode+1} = find(trackSegmentMode == -1);

numSegmentsMode(maxMode+1) = length(indxMode{maxMode+1});

%calculate fraction of track segments falling in each type/mode
if isempty(ignoreUn)
    fracSegmentsType = numSegmentsType / numTrackSegments;
    fracSegmentsMode = numSegmentsMode / numTrackSegments;
else
    fracSegmentsType = numSegmentsType / sum(numSegmentsType);
    fracSegmentsMode = numSegmentsMode / sum(numSegmentsMode);
end

if isempty(ignoreUn)
    %calculate number of features in each motion type
numFeatType = [length(find(tracksIndxMat(indxLin,:))) ...
    length(find(tracksIndxMat(indxFree,:))) ...
    length(find(tracksIndxMat(indxConf,:))) ...
    length(find(tracksIndxMat(indxImm,:))) ...
    length(find(tracksIndxMat(indxUndet,:)))]';
else
    %calculate number of features in each motion type
numFeatType = [length(find(tracksIndxMat(indxLin,:))) ...
    length(find(tracksIndxMat(indxFree,:))) ...
    length(find(tracksIndxMat(indxConf,:))) ...
    length(find(tracksIndxMat(indxImm,:)))]'; 
end
%calculate number of features in each motion mode
if isempty(ignoreUn)
    numFeatMode = zeros(maxMode+1,1);
    for iMode = 1 : maxMode+1
        numFeatMode(iMode) = length(find(tracksIndxMat(indxMode{iMode},:)));
    end
else
    numFeatMode = zeros(maxMode,1);
    for iMode = 1 : maxMode
        numFeatMode(iMode) = length(find(tracksIndxMat(indxMode{iMode},:)));
    end
end

%get fraction of features in each type/mode - this is the probability
%of a feature undergoing a certain motion type/mode
if isempty(ignoreUn)
    probFeatType = numFeatType / numFeatTot;
    probFeatMode = numFeatMode / numFeatTot;
else
    probFeatType = numFeatType / sum(numFeatType);
    probFeatMode = numFeatMode / sum(numFeatMode);
end

%also calculate average number of features in each type/mode per frame
aveFeatTypePerFrame = numFeatType / numFrames;
aveFeatModePerFrame = numFeatMode / numFrames;

%% merges/splits vs. features

%initialize lists of merges/splits
listOfMerges = [];
listOfSplits = [];

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get the sequence of events of this compound track
    seqOfEvents = tracks(iTrack).seqOfEvents;

    %find where this track has merges and splits
    indxSplit = find( seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)) );
    indxMerge = find( seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)) );
    
    %store the merge/split time and the segment after merging/before
    %splitting (remember tracks are in alternative format)
    infoSplit = unique(seqOfEvents(indxSplit,[1 4]),'rows'); %#ok<FNDSB>
    infoMerge = unique(seqOfEvents(indxMerge,[1 4]),'rows'); %#ok<FNDSB>
    
    %go over the splits
    for iMS = 1 : size(infoSplit,1)
        
        %find the rows of this split
        rowSplit = find( seqOfEvents(:,1)==infoSplit(iMS,1) & seqOfEvents(:,4)==infoSplit(iMS,2) );
        
        %make a row vector of the segment before splitting and segments after
        %splitting
        segmentsMS = [infoSplit(iMS,2) seqOfEvents(rowSplit,3)']; %#ok<FNDSB>
        
        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;

        %add this split to the list of splits
        listOfSplits = [listOfSplits; segmentsMS]; %#ok<AGROW>

    end

    %go over the merges
    for iMS = 1 : size(infoMerge,1)
        
        %find the rows of this merge
        rowMerge = find( seqOfEvents(:,1)==infoMerge(iMS,1) & seqOfEvents(:,4)==infoMerge(iMS,2) );
        
        %make a row vector of the segment after merging and segments before
        %merging
        segmentsMS = [infoMerge(iMS,2) seqOfEvents(rowMerge,3)']; %#ok<FNDSB>
        
        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;

        %add this merge to the list of merges
        listOfMerges = [listOfMerges; segmentsMS]; %#ok<AGROW>

    end

end

%get total number of merges and splits
numMergesTot = size(listOfMerges,1);
numSplitsTot = size(listOfSplits,1);

%calculate average number of merges/splits per frame
aveMergePerFrame = numMergesTot / numFrames;
aveSplitPerFrame = numSplitsTot / numFrames;

%calculate overall probability of merging/splitting
%note that this considers that each merge/split involves two particles
probFeatMerge = 2 * aveMergePerFrame / ( aveFeatPerFrame * (aveFeatPerFrame-1) );
probFeatSplit = 2 * aveSplitPerFrame / ( aveFeatPerFrame * (aveFeatPerFrame-1) );

%% merges/splits vs. type and mode

%get the types/modes of segments participating in merges and splits
if numMergesTot > 0
    listOfMergeTypes = [trackSegmentType(listOfMerges(:,1)) ...
        trackSegmentType(listOfMerges(:,2)) trackSegmentType(listOfMerges(:,3))];
    listOfMergeModes = [trackSegmentMode(listOfMerges(:,1)) ...
        trackSegmentMode(listOfMerges(:,2)) trackSegmentMode(listOfMerges(:,3))];
else
    [listOfMergeTypes,listOfMergeModes] = deal(zeros(0,3));
end
if numSplitsTot > 0
    listOfSplitTypes = [trackSegmentType(listOfSplits(:,1)) ...
        trackSegmentType(listOfSplits(:,2)) trackSegmentType(listOfSplits(:,3))];
    listOfSplitModes = [trackSegmentMode(listOfSplits(:,1)) ...
        trackSegmentMode(listOfSplits(:,2)) trackSegmentMode(listOfSplits(:,3))];
else
    [listOfSplitTypes,listOfSplitModes] = deal(zeros(0,3));
end

%sort the lists so that the more dynamic type/mode among the segments
%before merging or after splitting comes first
listOfMergeTypes(:,2:3) = sort(listOfMergeTypes(:,2:3),2,'descend');
listOfSplitTypes(:,2:3) = sort(listOfSplitTypes(:,2:3),2,'descend');
listOfMergeModes(:,2:3) = sort(listOfMergeModes(:,2:3),2,'descend');
listOfSplitModes(:,2:3) = sort(listOfSplitModes(:,2:3),2,'descend');

%get number of merges/splits based on the types of the three track segments
%involved
if isempty(ignoreUn)
    typeList = [3 2 1 0 -1];
else
    typeList = [3 2 1 0];
end
numTypeIndx = length(typeList);
[numMergesTypeDetail,numSplitsTypeDetail] = deal(NaN(numTypeIndx*ones(1,3)));
for kType = 1 : numTypeIndx %segment after merging
    for iType = 1 : numTypeIndx %more dynamic segment before merging
        for jType = iType : numTypeIndx %less dynamic segment before merging
            numMergesTypeDetail(iType,jType,kType) = length(find( ...
                listOfMergeTypes(:,1) == typeList(kType) & ...
                listOfMergeTypes(:,2) == typeList(iType) & ...
                listOfMergeTypes(:,3) == typeList(jType) ));
        end
    end
end
for kType = 1 : numTypeIndx %segment before splitting
    for iType = 1 : numTypeIndx %more dynamic segment after splitting
        for jType = iType : numTypeIndx %less dynamic segment after splitting
            numSplitsTypeDetail(iType,jType,kType) = length(find( ...
                listOfSplitTypes(:,1) == typeList(kType) & ...
                listOfSplitTypes(:,2) == typeList(iType) & ...
                listOfSplitTypes(:,3) == typeList(jType) ));
        end
    end
end

%get number of merges/splits based on the modes of the three track segments
%involved
if isempty(ignoreUn)
    modeList = [maxMode:-1:1, -1];
else
    modeList = [maxMode:-1:1];
end
numModeIndx = length(modeList);
[numMergesModeDetail,numSplitsModeDetail] = deal(NaN(numModeIndx*ones(1,3)));
for kType = 1 : numModeIndx %segment after merging
    for iType = 1 : numModeIndx %more dynamic segment before merging
        for jType = iType : numModeIndx %less dynamic segment before merging
            numMergesModeDetail(iType,jType,kType) = length(find( ...
                listOfMergeModes(:,1) == modeList(kType) & ...
                listOfMergeModes(:,2) == modeList(iType) & ...
                listOfMergeModes(:,3) == modeList(jType) ));
        end
    end
end
for kType = 1 : numModeIndx %segment before splitting
    for iType = 1 : numModeIndx %more dynamic segment after splitting
        for jType = iType : numModeIndx %less dynamic segment after splitting
            numSplitsModeDetail(iType,jType,kType) = length(find( ...
                listOfSplitModes(:,1) == modeList(kType) & ...
                listOfSplitModes(:,2) == modeList(iType) & ...
                listOfSplitModes(:,3) == modeList(jType) ));
        end
    end
end

%calculate average merges/splits per frame
aveMergesTypeDetail = numMergesTypeDetail / numFrames;
aveSplitsTypeDetail = numSplitsTypeDetail / numFrames;
aveMergesModeDetail = numMergesModeDetail / numFrames;
aveSplitsModeDetail = numSplitsModeDetail / numFrames;

%% probability of merging/splitting vs. type/mode of segments before merging/after splitting

%collapse merge/split details to retain only type/mode information of
%segments before merging/after splitting
aveMergesTypeSeparate = sum(aveMergesTypeDetail,3);
aveSplitsTypeSeparate = sum(aveSplitsTypeDetail,3);
aveMergesModeSeparate = sum(aveMergesModeDetail,3);
aveSplitsModeSeparate = sum(aveSplitsModeDetail,3);

%calculate probability of merging/splitting per motion type pair,
%equivalent to calculation above for overall probability
[probMergeTypeSeparate,probSplitTypeSeparate] = deal(NaN(numTypeIndx));
for iType = 1 : numTypeIndx
    probMergeTypeSeparate(iType,iType) = 2 * aveMergesTypeSeparate(iType,iType) / ( aveFeatTypePerFrame(iType) * (aveFeatTypePerFrame(iType)-1) );
    for jType = iType + 1 : numTypeIndx
        probMergeTypeSeparate(iType,jType) = aveMergesTypeSeparate(iType,jType) / ( aveFeatTypePerFrame(iType) * aveFeatTypePerFrame(jType) );
    end
end
for iType = 1 : numTypeIndx
    probSplitTypeSeparate(iType,iType) = 2 * aveSplitsTypeSeparate(iType,iType) / ( aveFeatTypePerFrame(iType) * (aveFeatTypePerFrame(iType)-1) );
    for jType = iType + 1 : numTypeIndx
        probSplitTypeSeparate(iType,jType) = aveSplitsTypeSeparate(iType,jType) / ( aveFeatTypePerFrame(iType) * aveFeatTypePerFrame(jType) );
    end
end

%calculate probability of merging/splitting per motion mode pair,
%equivalent to calculation above for overall probability
[probMergeModeSeparate,probSplitModeSeparate] = deal(NaN(numModeIndx));
for iMode = 1 : numModeIndx
    probMergeModeSeparate(iMode,iMode) = 2 * aveMergesModeSeparate(iMode,iMode) / ( aveFeatModePerFrame(iMode) * (aveFeatModePerFrame(iMode)-1) );
    for jMode = iMode + 1 : numModeIndx
        probMergeModeSeparate(iMode,jMode) = aveMergesModeSeparate(iMode,jMode) / ( aveFeatModePerFrame(iMode) * aveFeatModePerFrame(jMode) );
    end
end
for iMode = 1 : numModeIndx
    probSplitModeSeparate(iMode,iMode) = 2 * aveSplitsModeSeparate(iMode,iMode) / ( aveFeatModePerFrame(iMode) * (aveFeatModePerFrame(iMode)-1) );
    for jMode = iMode + 1 : numModeIndx
        probSplitModeSeparate(iMode,jMode) = aveSplitsModeSeparate(iMode,jMode) / ( aveFeatModePerFrame(iMode) * aveFeatModePerFrame(jMode) );
    end
end

%% probability of merging/splitting vs. type/mode of segment after merging/before splitting

%collapse merge/split details to retain only type/mode information of
%segment after merging/before splitting
aveMergesTypeTogether = squeeze(nansum(nansum(aveMergesTypeDetail,1),2));
aveSplitsTypeTogether = squeeze(nansum(nansum(aveSplitsTypeDetail,1),2));
aveMergesModeTogether = squeeze(nansum(nansum(aveMergesModeDetail,1),2));
aveSplitsModeTogether = squeeze(nansum(nansum(aveSplitsModeDetail,1),2));

%calculate probability of merging/splitting per motion type/mode after
%merging/before splitting
probMergeTypeTogether = aveMergesTypeTogether ./ aveFeatTypePerFrame;
probSplitTypeTogether = aveSplitsTypeTogether ./ aveFeatTypePerFrame;
probMergeModeTogether = aveMergesModeTogether ./ aveFeatModePerFrame;
probSplitModeTogether = aveSplitsModeTogether ./ aveFeatModePerFrame;

%% probability of type/mode after merging/before splitting vs. types/modes before merging/after splitting and vice versa

probMergeTypeTogVsSep = aveMergesTypeDetail ./ repmat(aveMergesTypeSeparate,1,1,numTypeIndx);
probSplitTypeTogVsSep = aveSplitsTypeDetail ./ repmat(aveSplitsTypeSeparate,1,1,numTypeIndx);
probMergeModeTogVsSep = aveMergesModeDetail ./ repmat(aveMergesModeSeparate,1,1,numModeIndx);
probSplitModeTogVsSep = aveSplitsModeDetail ./ repmat(aveSplitsModeSeparate,1,1,numModeIndx);

probMergeTypeSepVsTog = aveMergesTypeDetail ./ repmat(reshape(aveMergesTypeTogether,1,1,numTypeIndx),numTypeIndx,numTypeIndx,1);
probSplitTypeSepVsTog = aveSplitsTypeDetail ./ repmat(reshape(aveSplitsTypeTogether,1,1,numTypeIndx),numTypeIndx,numTypeIndx,1);
probMergeModeSepVsTog = aveMergesModeDetail ./ repmat(reshape(aveMergesModeTogether,1,1,numModeIndx),numModeIndx,numModeIndx,1);
probSplitModeSepVsTog = aveSplitsModeDetail ./ repmat(reshape(aveSplitsModeTogether,1,1,numModeIndx),numModeIndx,numModeIndx,1);

%% output

%general statistics
statsGeneral = struct('numTracks',numTracks,'numTrackSegments',numTrackSegments,...
    'numFeatPerFrame',aveFeatPerFrame);

%motion statistics
statsType = struct('fracSegments',fracSegmentsType,...
    'prob',probFeatType);
statsMode = struct('fracSegments',fracSegmentsMode,...
    'prob',probFeatMode);
statsMotion = struct('type',statsType,'mode',statsMode);

%merging statistics
statsType = struct('numTotDetail',numMergesTypeDetail,...
    'numPerFrameDetail',aveMergesTypeDetail,...
    'probVsSepSegments',probMergeTypeSeparate,...
    'probVsTogSegment',probMergeTypeTogether,...
    'probMotionTogVsSep',probMergeTypeTogVsSep,...
    'probMotionSepVsTog',probMergeTypeSepVsTog);
statsMode = struct('numTotDetail',numMergesModeDetail,...
    'numPerFrameDetail',aveMergesModeDetail,...
    'probVsSepSegments',probMergeModeSeparate,...
    'probVsTogSegment',probMergeModeTogether,...
    'probMotionTogVsSep',probMergeModeTogVsSep,...
    'probMotionSepVsTog',probMergeModeSepVsTog);
statsMerging = struct('numTot',numMergesTot,'numPerFrame',aveMergePerFrame,...
    'probOverall',probFeatMerge,'type',statsType,'mode',statsMode);

%splitting statistics
statsType = struct('numTotDetail',numSplitsTypeDetail,...
    'numPerFrameDetail',aveSplitsTypeDetail,...
    'probVsSepSegments',probSplitTypeSeparate,...
    'probVsTogSegment',probSplitTypeTogether,...
    'probMotionTogVsSep',probSplitTypeTogVsSep,...
    'probMotionSepVsTog',probSplitTypeSepVsTog);
statsMode = struct('numTotDetail',numSplitsModeDetail,...
    'numPerFrameDetail',aveSplitsModeDetail,...
    'probVsSepSegments',probSplitModeSeparate,...
    'probVsTogSegment',probSplitModeTogether,...
    'probMotionTogVsSep',probSplitModeTogVsSep,...
    'probMotionSepVsTog',probSplitModeSepVsTog);
statsSplitting = struct('numTot',numSplitsTot,'numPerFrame',aveSplitPerFrame,...
    'probOverall',probFeatSplit,'type',statsType,'mode',statsMode);



%% ~~~ the end ~~~

