% %  Steps for reforming tracks:

% Step 1: reform tracks with artifacts coming from simultaneous merges and splits
tracksReformSimul = removeSimultaneousMergeSplit(tracksFinal);

% Step 2: reshape the seqOfEvents when a segment split and merge with the same segment again. If the oldest segment merge with the newest this function will exchange segment index to maintain one long segment and a short that split and merge. 

tracksReformat = reformatSplitToMergeEvents( tracksFinal );

% Step 3: reform tracks with artifacts coming from born to merge, split to end or split to merge. This function take care of the possible artifacts in a chronological way.

tracksReformChronol = removeSplitMergeArtifactsChronological(tracksReformat,thresholdTime);
 
% Step4: calculate aggregState, clustHistory and probability to be an artifact, these inputs are needed to reform merges to splits 

%  calculate aggregate state
 functionCalcMultipleAggregStateUtrack(saveRoot,saveRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,intensityQuantum,flagTrackName);

%  calculate clust history
 functionCalcClustHistoryUtrack(saveRoot,rDDir,aPDir,dRDir,outDirNum,lRDir);
 

% calculate probability to be an artifact
 functionCalcProbArtifactsMerge2SplitTimes(saveRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime)

% Step 5:  reform tracks with artifacts coming from finite merges and splits
 tracksReformFinal = removeMerge2SplitClustHist(clustHistoryAll,tracksIn,probArtifactMerge2split);

