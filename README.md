# DynamicNeuronTracker (DyNT)

1. ***DynamicNeuronTracker*** extracts dynamic Region-Of-Interests (ROIs) of firing neurons in 3D calcium images of deforming live tissues. 
2. The package is built on Matlab (ver. 2020a). 
3. The algorithm utilizes local patch-matching via spatial correlation to track jittering and flickering neurons in a 3D Ca2+ movie. 

> Toy example of tracked and segmented neurons
 
 https://user-images.githubusercontent.com/38955889/128269347-8013f98b-85e6-4180-a563-4adb23962847.mp4

<img src="/DynamicNeuronTrackerPackage/Doc/HCL-indexedActivityMap.png" width = 49.5%>


<p>&nbsp;</p>

## User Manual

- See a User Manual ([/DynamicNeuronTrackerPackage/Doc/UserManualDynamicNeuronTracker.pdf](/DynamicNeuronTrackerPackage/Doc/UserManualDynamicNeuronTracker.pdf)).
- Or directly go to [a single script](/DynamicNeuronTrackerPackage/DynamicNeuronTracker/Pipelines/masterScript_toExtractDynamicROIsOf3DCaImaging.m) to run the whole pipeline.
- Before running the pipeline, an imaging dataset needs to be registered to a 'movieData' object using 'movieSelectorGUI'. For a detailed guide, see the first half of [a document](DynamicNeuronTrackerPackage/Doc/TFMPackage.pdf) from [Danuser Lab repositories](https://github.com/DanuserLab?tab=repositories).


## Contact

Jungsik Noh (jungsik.noh@utsouthwestern.edu)


