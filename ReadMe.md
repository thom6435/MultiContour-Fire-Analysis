# MultiContour-Fire-Analysis

This project uses a contour analysis approach to process RGB and thermal (IR) imagery such that the active burning area can be detected and analyzed. By extended the analysis over a time sequence of imagery, the rate of spread, and major direction can be determed.

The current version is built using OpenCV 4.5, and requires the following OpenCV libraries:

* opencv_highgui 
* opencv_core 
* opencv_imgcodecs 
* opencv_imgproc

The _Config_ folder contains a sample configuration file. The configuration file defines the processing parameters, the input imagery, and the results path. To create a config file use the _BatchCreator_ bash script.

To run an analysis:

`MultiContour <<Config-File-Name>>`

The sample config file is setup for the included sample imagery. Note, depending on your install path, you may need to update the paths in the config.

The _Background-Docs_ folder contains imagery examples and short documents that define the processing approach. 

