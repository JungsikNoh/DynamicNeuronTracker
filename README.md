# DynamicNeuronTracker (DyNT)

1. ***DynamicNeuronTracker*** extracts dynamic Region-Of-Interests (ROIs) of firing neurons in 3D calcium images of deforming live tissues. 
2. The package is built on Matlab (ver. 2020a). 
3. The algorithm utilizes local patch-matching via spatial correlation to track jittering and flickering neurons in a 3D Ca2+ movie. 

## Toy example of tracked and segmented neurons

> Maximum intensity projection video and overlaid dynamic ROIs

https://github.com/JungsikNoh/DynamicNeuronTracker/assets/38955889/b3b84c60-8b68-43b3-b31d-52c7a397e6f8

> Normalized neuron activity time courses
<img src="Doc/HCL-indexedActivityMap.png" width="49.5%"/>




## User Manual

- See a User Manual ([/DynamicNeuronTrackerPackage/Doc/UserManualDynamicNeuronTracker.pdf](/DynamicNeuronTrackerPackage/Doc/UserManualDynamicNeuronTracker.pdf)).
- Or directly go to [a single script](/DynamicNeuronTrackerPackage/DynamicNeuronTracker/Pipelines/masterScript_toExtractDynamicROIsOf3DCaImaging.m) to run the whole pipeline.
- Before running the pipeline, an imaging dataset needs to be registered to a 'movieData' object using 'movieSelectorGUI'. For a detailed guide, see the first half of [a document](DynamicNeuronTrackerPackage/Doc/TFMPackage.pdf) from [Danuser Lab repositories](https://github.com/DanuserLab?tab=repositories).


## Contact

Jungsik Noh (jungsik.noh@utsouthwestern.edu)


