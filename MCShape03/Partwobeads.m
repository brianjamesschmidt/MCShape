% Parameter file for MCShape03.m
% Filename of video stack to analyze, including extension
inputvideo='Twobeads.avi';

% If there is header data, indicate whether it is appended to the top of
% the images or bottom. " headertop = true " or " headertop = false "
headertop=true;
topheaderheight=67;
headerbottom=false;
bottomheaderheight=0;

% Enter the approximate diameter of the bead in micrometers.
beadsize=15;

% Enter the coefficient of variation in the diameter.  Program will search 
% for objects +/- 5 standard deviations in diameter.  If left as NaN, 
% program will default to object size +50 / -80% unless special sizing
% parameters are entered.
beadsizecv=NaN;
% Enter any special sizing guidlines for tracking.  If left as NaN, program
% will use +50 / -20%
specialmaxsize=NaN;
specialminsize=NaN;

% Minimum threshold to search for objects in terms of fraction of maximum
% intensity (0-1).  Enter 0 to search all the way, but if the image is 
% noisy this may result in false object identification.
absoluteminthresh=.1;

% Image resolution in pixels/µm.  If the video is interlaced, enter the
% camera pixel resolution and the program will perform the conversion on 
% the deitnerlaced image stack.
xres=5.4;
yres=5.4;

% Enter record speed in fps
recordfps=250;

% Enter whether to deinterlace the video.  Boolean true
% or false.
deinterlace=false;

% Enter any additional identifiers you want appended to the name of the
% output avi's and centroid matrix files
outputappend='tracked';

% Set the image feature desired to track based on.  By default, the program
% thresholds the gradient of the image to detect edges.  Alternatively, it 
% will invert and threshold the image itself.
% Two options: trackgradient = true
%              trackgradient = false
trackgradient=true;

% Specify whether the program should perform deconvolution while tracking
deconvolve=false;
psffile='measure_your_psf_and_specify_it_here.tif';

% Specify whether to correct for gain errors.  Each should be in a *.mat
% file saved as Slope and Intercept.  Program will append *.mat.
gaincorrect=false;
slopefile='Photronslope071102';
interceptfile='Photronintercept071102';

% Enter the beginning frame in the avi to use for tracking
% Input NaN, and the program will assume 1.
startframe=NaN;

% Enter the ending frame to use for image tracking
% Input NaN, and the program will stop at the end 
endframe=NaN;

% Enter the integer height of the image data (not including headers) in 
% pixels.  The program assumes there may be a header in the top of the
% vertical direction but no borders around the left and right edges
% Enter NaN if the entire image is data.
%dheight=380;

% Compression setting to use for output.  Pick 'indeo5' or 'none'.
compressionset='None';

% Indicate whether you would like to show the morphological line axes
% on the generated movie file
showline=false;

% Boolean variable indicating whether to apply a spatial mean "blurring" 
% filter while picking out edges.  Sometimes helpful for weak signal 
% situations or with ragged edges.
% Blursize is the nXn mask used for blurring.  Choose an odd number if
% possible to keep the tracked area in the center of the object.  If an
% even number is chosen, however, the tracking will still be good since
% frame-to-frame differences are important.
meanfilter=false;
meanfiltersize=5;

% Boolean variable indicating whether to apply a spatial median filter
medianfilter=false;
medianfiltersize=0;

% The acceleration factor will speed the search by reducing the spacing of
% threshold levels to search the binary image for objects.  Set equal to 
% an integer.  A factor of 1 is the default and provides the most sensitive
% search.
accelerationfactor=1;
