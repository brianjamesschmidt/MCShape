MCShape
Multiparticle Tracking algorithm for cell motion and morphology
=============

This is a multiparticle version of the CMorph algorithm previously developed. The algorithm takes the gradient of the image, and then applies decreasing thresholds and fill operations to attempt to identify objects. One it has found the objects, it attempts to pair them up with objects in the previous frame. It was designed to track cells and microbeads in flow cell assays, so it is designed to be biased to look for a flow direction.

Instructions:
Download the tracking algorithm and sample parameter file. 
The parameter file is set up for the sample movie, 
or use your own uncompressed avi and adjust the parameter file.
Place all 3 in your Matlab root directory.
Type: MCShape03('Partwobeads')

The program will output two files: 
1. A color-coded *.avi with object numbers so it is easy to check the results. 
2. A *.mat file with the tracking results stored as a structure. The field id's will match the object numbers in the video.

I hope this tool will help you in your research. If you use this to help with a publication, please cite the work in your methods section.

Reference:
Schmidt, B.J., J.A. Papin, and M.B. Lawrence, Nano-motion Dynamics are Determined by Surface-Tethered Selectin Mechanokinetics and Bond Formation. PLoS Comput Biol, 2009. 5(12): p. e1000612.

Reference for original, single particle algorithm (CMorph), conference paper available at IEEXplore: 
Schmidt, B.J., C.D. Paschall, W.H. Guilford, and M.B. Lawrence, "High-Resolution Optical Tracking to Identify Adhesive Events in Vitro." 
Conference Proceedings of the Asilomar Conference on Signals, Systems, and Computers, Pacific Grove, CA, November 4-7, 2007.  p. 1856-1860.

Requirements:
MATLAB release		developed with MATLAB 7.5 (R2007b)
Other requirements	Image Processing Toolbox
Other requirements	Optional for reading compressed formats not supported by mmreader function: videoIO LIBRARY (Gerald Dalley)

One day, when there is time, I might update this and add a GUI...

Release history
Release		Date           Notes
-----------------------------------------------------------------------------------------
03.0		12.21.2009     -Last update to package as distributed on MATLAB Central.
