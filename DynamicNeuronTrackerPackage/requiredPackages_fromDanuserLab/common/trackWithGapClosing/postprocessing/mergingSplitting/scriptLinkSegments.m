 % first check the intensity
                    % Intensity before the event
                    
                    meanIntensitySeg1Bef=mean(tracksFeatIndxCG(uniqueSegment,1:timeRepetitionTemp-1));
                    meanIntensitySeg2Bef=mean(tracksFeatIndxCG(segmentMerge,1:timeRepetitionTemp-1));
                    
                    % Intensity after the event
                    
                    meanIntensitySeg1Aft=mean(tracksFeatIndxCG(uniqueSegment,timeRepetitionTemp:end));
                    meanIntensitySeg2Aft=mean(tracksFeatIndxCG(segmentSplit,timeRepetitionTemp:end));
                    
                    % Calculate the min distance
                    
                    distInten1=abs(meanIntensitySeg1Bef-[meanIntensitySeg1Aft,meanIntensitySeg2Aft]);
                    distInten2=abs(meanIntensitySeg2Bef-[meanIntensitySeg1Aft,meanIntensitySeg2Aft]);
                    
                    
                    % check the displacements
                                       
                    % before the event
                    
                    meanDispSeg1xBef=nanmean(tracksCoordAmpCG(uniqueSegment,1:8:8*timeRepetitionTemp-1));
                    meanDispSeg1yBef=nanmean(tracksCoordAmpCG(uniqueSegment,2:8:8*timeRepetitionTemp-1));
                    meanDispSeg2xBef=nanmean(tracksCoordAmpCG(segmentMerge,1:8:8*timeRepetitionTemp-1));
                    meanDispSeg2yBef=nanmean(tracksCoordAmpCG(segmentMerge,2:8:8*timeRepetitionTemp-1));
                    
                    % after the event
                    meanDispSeg1xAft=nanmean(tracksCoordAmpCG(uniqueSegment,8*(timeRepetitionTemp)+1:8:end));
                    meanDispSeg1yAft=nanmean(tracksCoordAmpCG(uniqueSegment,8*(timeRepetitionTemp)+2:8:end));
                    
                    meanDispSeg1xAft=nanmean(tracksCoordAmpCG(segmentSplit,8*(timeRepetitionTemp)+1:8:end));
                    meanDispSeg1yAft=nanmean(tracksCoordAmpCG(segmentSplit,8*(timeRepetitionTemp)+2:8:end));

meanDisple1Bef=sqrt(meanDispSeg1xBef^2+meanDispSeg1yBef^2);
meanDisple2Bef=sqrt(meanDispSeg2xBef^2+meanDispSeg2yBef^2);
meanDisple3=sqrt(meanDispSeg3x^2+meanDispSeg3y^2);
                    % determine which should be the continuation of the
                    % segment
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%% pair of segments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   % calculate mean intensity and displacements for all
                        % merges and splits
                        
                        %intensity
                        meanIntSegBefore=mean(tracksFeatIndxCG(segmentNumber,1:timeMerge-1));
                        meanIntSegAfter =mean(tracksFeatIndxCG(segmentNumber,timeSplit:end));
                        meanIntMerge=mean(tracksFeatIndxCG(segmentNumber,1:timeMerge-1));
                        meanIntSplit=mean(tracksFeatIndxCG(segmentNumber,timeSplit:end));
                        
                        % take the x and y values
                        
                        
                        mergeXvals=tracksCoordAmpCG(segmentMerge,1:8:8*timeRepetitionTemp-1);
                        mergeYvals=tracksCoordAmpCG(segmentMerge,2:8:8*timeRepetitionTemp-1);
                        
                        splitXvals=tracksCoordAmpCG(segmentSplit,8*(timeRepetitionTemp)+1:8:end);
                        splitYvals=tracksCoordAmpCG(segmentSplit,8*(timeRepetitionTemp)+2:8:end);
                        
                        %calculate displacement
                        
                        dispMergeX=mergeXvals(:,2:end)-mergeXvals(:,1:end-1);
                        dispMergeY=mergeYvals(:,2:end)-mergeYvals(:,1:end-1);
                        
                        dispSplitX=splitXvals(:,2:end)-splitXvals(:,1:end-1);
                        dispSplitY=splitYvals(:,2:end)-splitYvals(:,1:end-1);
                        
                        % calculate mean displacement
                        
                        meanDispMerge=nanmean(sqrt(dispMergeX.^2+dispMergeY.^2),2);
                        meanDispSplit=nanmean(sqrt(dispSplitX.^2+dispSplitY.^2),2);
                        
                        % calculate ratio, always doing larger/smaller
                        if length(meanDispMerge)>length(meanDispSplit)
                            dispLarger=meanDispMerge;
                            dispSmaller=meanDispSplit;
                            intLarger=meanIntMerge;
                            intSmaller=meanIntSplit;
                            segSmall=1;
                        else
                            dispLarger=meanDispSplit;
                            dispSmaller=meanDispMerge;
                            intLarger=meanIntSplit;
                            intSmaller=meanIntMerge;
                            segSmall=2;
                        end
                        
                        
                        for indexSeg=1:length(dispSmaller)
                            %calculate ratio
                            ratioIntensity=dispLarger/dispSmaller(indexSeg);
                            ratioMeanDisp=intLarger/intSmaller(indexSeg);
                            
                            ratioFinal=ratioIntensity.*ratioMeanDisp;
                            
                            [~,positionMinDistance]=min(abs(ratioFinal-1));
                            
                            % pair segments
                            if segSmall==1
                                % define which should be the pair of merge and split
                                
                                segMergeTemp=segmentMerge(indexSeg);
                                segSplitTemp=segmentSplit(positionMinDistance);
                                
                                % Remove segments, here will define the pairs based in the
                                % merge and split index
                                removeIndex=[possibleArtifactsTime(mergeIndex(indexSeg)),possibleArtifactsTime(splitIndex(positionMinDistance))];
                                
                                
                                
                            elseif segSmall==2
                                
                                segMergeTemp=segmentMerge(positionMinDistance);
                                segSplitTemp=segmentSplit(indexSeg);
                                
                                % Remove segments, here will define the pairs based in the
                                % merge and split index
                                removeIndex=[possibleArtifactsTime(mergeIndex(positionMinDistance)),possibleArtifactsTime(splitIndex(indexSeg))];
                                
                            end
                            
                            %remove the merge and the split
                            seqOfEvents(removeIndex,:)=[];
                            
                            % replace the segment number to have only one continuos
                            % segment
                            seqOfEvents(seqOfEvents(:,3)==segSplitTemp,3)=segMergeTemp;
                            seqOfEvents(seqOfEvents(:,4)==segSplitTemp,4)=segMergeTemp;
                            
                            % for all segments with index larger than removed
                            % segment, reduce their index by 1
                            
                            seqOfEventsIndex= seqOfEvents(:,3:4);
                            segmentNumbersOld=unique(seqOfEventsIndex);
                            segmentNumbersOld(isnan(segmentNumbersOld))=[];
                            for i=1:length(segmentNumbersOld)
                                if segmentNumbersOld(i)~=i
                                    seqOfEventsIndex(seqOfEventsIndex==segmentNumbersOld(i))=i;
                                end
                                
                                %update seqOfEvents
                                seqOfEvents(:,3:4)=seqOfEventsIndex;
                                
                                
                            end
                            
                        end
        