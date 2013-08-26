function []=MCShape03(parameterfilename)
%  MCShape03(parameter filename as 'parameterfilename'.m)
%  NOTE: USE OF THIS M-FILE REQUIRES:
%  1. MATLAB 2007b OR LATER 
%  2. THE PROPER CODECS INSTALLED (TRY K-LITE, FREE, IF NEEDED)
%  3. OPTIONAL: FREE 3RD PARTY videoIO LIBRARY (Gerald Dalley)
%     AVAILABLE THROUGH MATLAB CENTRAL
%  Cell tracking program originally coded by Brian Schmidt at the 
%  University of Virginia, beginning with GradTrackXX.m in May 2005.
%  Two operation modes.  In gradient mode, This program uses the 
%  contrast in surrounding objects (Becke rings if defocused) to attempt 
%  to locate the edge.  It then uses the centroid of the generated shape to
%  track position and records shape properties.
%  Results are output into two files: one is a *.mat file with the tracking
%  results and the other is a video file that shows the tracked region.
%  *.mat file formating:
%  Structure, each field is a tracked object named 'ObjectN' with:
%  Column 1: Apparent object radius
%  Column 2: Time
%  Column 3: Centroid, X-Coordinate
%  Column 4: Centroid, Y-Coordinate
%  Column 5: Morphology data: bounding box x-width
%  Column 6: Morphology data: bounding box y-height
%  Column 7: Morphology data: orientation angle (degrees)
%  Column 8: Morphology data: object length along orientation angle axis
%  Column 9: Morphology data: Object width along line perpendicular to
%           orientation angle axis through centroid
%  The Shape Strain Index used to color the object is the ratio of
%  column 8 : column 9. 
%  The video should be watched to verify the algorithm correctly identified
%  the objects of interest.
%  -----REVISION HISTORY-----
%  V03:   2/05/09:    Added error messages to help trouble-shoot
%                     problems reading in compressed video formats.
%  V02:   1/19/09:    Incorporated allowance for data header sections at 
%                     the bottom of the data in addition to the top.
%  V01:   5/29/08:    MCShape.m: Major revision to previous algorithm to 
%                     incorporate multiparticle tracking.
%  ---:   11/02/07:   Cmorph.m: Major revision to previous GradTrackXX.m
%                     algorithms to incorporate shape-based parameter
%                     measurement.

tic
input_file=parameterfilename;

% Read the input parameters into the workspace
eval(input_file);

accelerationfactor=uint8(accelerationfactor);

if deinterlace==true
    yres=yres/2;
    nsubframes=2;
    outputfps=recordfps*2;
else
    nsubframes=1;
    outputfps=recordfps*1;
end

maxres=max(xres,yres);
minres=min(xres,yres);
resratio=minres/maxres;

% Find the root of the input video file without the extension
stringlength=length(inputvideo);
inputvideoroot=inputvideo(1:stringlength-4);

% AVI output file to write the overlay to.
avioutfile=strcat(inputvideoroot,outputappend,'overlay.avi');

% The filename of the centroid location is automatically appended with .mat
% by matlab.
centertrackerfile=strcat(inputvideoroot,outputappend);

% Variable for remembering the last user notification of tracking status
lastinformeduser=uint8(0);

% Create binary image for individual numbers to label the objects on screen
Numpic(:,:,1)=[1, 1, 1;
                1, 0, 1;
                1, 0, 1;
                1, 0, 1;
                1, 1, 1];
            
Numpic(:,:,2)=[0, 1, 0;
               1, 1, 0;
               0, 1, 0;
               0, 1, 0;
               1, 1, 1];
     
Numpic(:,:,3)=[1, 1, 0;
               0, 0, 1;
               0, 1, 0;
               1, 0, 0;
               1, 1, 1];

Numpic(:,:,4)=[1, 1, 0;
               0, 0, 1;
               1, 1, 0;
               0, 0, 1;
               1, 1, 0];

Numpic(:,:,5)=[1, 0, 1;
               1, 0, 1;
               1, 1, 1;
               0, 0, 1;
               0, 0, 1];
           
Numpic(:,:,6)=[1, 1, 1;
               1, 0, 0;
               1, 1, 0;
               0, 0, 1;
               1, 1, 0];     
           
Numpic(:,:,7)=[0, 1, 1;
               1, 0, 0;
               1, 1, 1;
               1, 0, 1;
               1, 1, 1];     
           
Numpic(:,:,8)=[1, 1, 1;
               0, 0, 1;
               0, 1, 0;
               1, 0, 0;
               1, 0, 0];            

Numpic(:,:,9)=[1, 1, 1;
               1, 0, 1;
               1, 1, 1;
               1, 0, 1;
               1, 1, 1];
           
Numpic(:,:,10)=[1, 1, 1;
               1, 0, 1;
               1, 1, 1;
               0, 0, 1;
               1, 1, 1];

% Assume the beginning frame for tracking is 1 if the user does not
% specify.  
if isnan(startframe) == 1
    startframe = 1;
end

% Pick the best reader for the data.  mmreader will return
% an empty number of frames if it cannot read the data well.
% Create the video reader object for getting frames.
vrvideo = mmreader(inputvideo);
if isempty(vrvideo.NumberOfFrames) == true
    fprintf('Could not read %s using Matlab mmreader.\n', inputvideo)
    fprintf('Attempting to open the video using the videoIO library by Gerald Dailey, available for free at Matlab Central.\n')
    fprintf('If the read fails, videoIO must be downloaded and Matlab set up properly, or the poper codec must be installed on this computer.\n')
    fprintf('To use videoIO, the current directory must be set to to the location videoIO has been downloaded to, and the video file to be tracked must be copied to the same directory.\n')
    fprintf('Do not put videoIO in a directory that the OS has protected, e.g. do not place in a subfolder of C:\Program Files with UAC active in Vista.\n')
    fprintf('Note: if the installed Matlab version is outdated, mmreader may have added support for the video format in the most recent release.\n');
    videoReaderin=true;
    % Use the third party videoIO library if mmreader is having difficulty.
    vrvideo = videoReader(inputvideo);
    % Need to call the next function at least once for the videoreader 
    % object.  "Fast forward" to the desired starting frame.
    for inputframecounter = 1 : startframe
        next(vrvideo);
    end
    fprintf('Successfully read the video format for %s using the videoIO library.\n', inputvideo)
    % Also create a mmreader object for additional file information  
    mmvideo=mmreader(inputvideo);
    % Bugs in mmreader for Matlab 2007b prevented the file info for the
    % number of frames from working directly.  Hence this indirect 
    % calculation.
    numberofinputframes=round(mmvideo.FrameRate*mmvideo.Duration)-1;
else
    fprintf('Successfully read the video format for %s using Matlab mmreader.\n', inputvideo)
    videoReaderin=false;
    vrvideo = mmreader(inputvideo);
    numberofinputframes=vrvideo.NumberOfFrames;
end

% If the user does not specify an ending frame, the program
% will continue until all of the frames are read.
if isnan(endframe) == 1
    endframe = numberofinputframes;
end
ninputmovieframes=endframe-startframe+1;
if deinterlace == true
    noutputmovieframes=ninputmovieframes*2;
else
    noutputmovieframes=ninputmovieframes;
end

% Calculate the input and output image size.
if videoReaderin==true
    [height,width,ncolorindexes]=size(getframe(vrvideo));
else
    [height,width,ncolorindexes]=size(read(vrvideo, 1));
end

% Calculate the height of a movie matrix image & output.
if headertop ~= false && headerbottom ~= false
    moviedataheight=height-topheaderheight-bottomheaderheight;
    movieimagetoppixel = topheaderheight + 1;
    movieimagebottompixel = height - bottomheaderheight;
    
elseif headertop ~= false
    moviedataheight=height-topheaderheight;
    movieimagetoppixel = topheaderheight + 1;
    movieimagebottompixel = height;    
elseif headerbottom ~= false
    moviedataheight=height-bottomheaderheight;
    movieimagetoppixel = 1;
    movieimagebottompixel = height - bottomheaderheight - 1;      
else
    moviedataheight=height;
    movieimagetoppixel = 1;
    movieimagebottompixel = height;
end
dheight=movieimagebottompixel-movieimagetoppixel+1;
if deinterlace == true
    outputdheight=dheight/2;
else
    outputdheight=dheight;
end

% Initialize the Centertracker matrix
centertrackerdatasize=uint8(9);
Centertracker=zeros(1,centertrackerdatasize);

% Initialize the avi used to verify the centroid tracker.
aviobj = avifile(avioutfile);
aviobj.compression=compressionset;
aviobj.fps=30;
aviobj.keyframe=30;
aviobj.quality=100;

% Set a threshold size to make sure the entire background isn't being 
% picked up as one object.  Also eliminate the detection of small objects 
% to avoid tracking debris.
if isnan(beadsizecv) == false
    beadsizemax=beadsize+5*beadsizecv;
    beadsizemin=max(1/max(xres,yres),beadsize-5*beadsizecv);
elseif isnan(specialmaxsize)==false || isnan(specialminsize)==false 
    beadsizemax=specialmaxsize;
    beadsizemin=specialminsize;
else
    beadsizemax=1.5*beadsize;
    beadsizemin=.8*beadsize;
end
pixelsmax=min(dheight*width*.5,pi*beadsizemax^2*.25*xres*yres);
% Set the minimum area hased on the expected object size.
pixelsmin=pi*beadsizemin^2*.25*xres*yres;

% Initialize the area check
lastarea=round((beadsize/2*maxres)^2*pi*resratio);

% Set up the filter for mean filtering.
if meanfilter == 1
    meanfilterimage = ones(meanfiltersize,meanfiltersize) / meanfiltersize^2;
end

if deconvolve == true
    % Read in the point spread function if it has been made available.
    PSF=imread(psffile);
    PSF=double(PSF);
    % Scale the PSF so it has a maximum of 1
    PSF=PSF./(max(max(PSF)));
end

if gaincorrect==true
    % Read in the slope and intercept files.  Reference the 2007 Asilomar
    % conference paper by Schmidt et al. for a description of how it can 
    % be used to correct the image.
    Temp=open(strcat(interceptfile,'.mat'));
    Intercept=Temp.Intercept;
    clear Temp
    Temp=open(strcat(slopefile,'.mat'));
    Slope=Temp.Slope;
    clear Temp
end

% Zero the flow direction vector
Flowdir=zeros(2,1);
% Zero the counter for total objects found
totalobjectsfound=uint32(0);
% Reset the inputframe counter to 1.
inputframecounter=uint32(1);
outputframecounter=uint32(1);
lastminthresh=uint8(0);
startthresh=uint8(255);
absoluteminthreshdetect=uint8(255);
maxdeltathreshobjectdetect=uint8(0);
% Iterate until all of the frames have been processed
while inputframecounter <= ninputmovieframes
    if videoReaderin==true
        Tempframe = getframe(vrvideo);
    else
        Tempframe = read(vrvideo, inputframecounter + startframe - 1);
    end

    if ncolorindexes==3
        Tempframe = rgb2gray(Tempframe);
    end

     % Pick the image data out of the matrix without the header data
     Movmatrix=Tempframe(movieimagetoppixel:movieimagebottompixel,1:width);
    
    % Apply differential pixel gain correction if instructed to.
    if gaincorrect==true
        Movmatrix=double(Movmatrix);
        Movmatrix=uint8(Slope.*Movmatrix+Intercept);
    end
    
    for subframecounter= 1:nsubframes

    % First create subframe matrices to correct for interlacing if
    % necessary.
    if nsubframes > 1
        % First zero movmatrixs to reset it to 1 z-element color depth
        % rather than 3.
        Movmatrixs=zeros(outputdheight,width);
        Movmatrixs(1:outputdheight,1:width)=Movmatrix(subframecounter:nsubframes:moviedataheight,1:width);
    else
        Movmatrixs=Movmatrix;      
    end
    % Also get the current header data
    if headertop ~= false
        Headertop=Tempframe(subframecounter:nsubframes:topheaderheight,1:width);
    end
    if headerbottom ~= false
        Headerbottom=Tempframe(movieimagebottompixel+subframecounter:nsubframes:height,1:width);
    end      

    % The area is used as a flag to check that a cell has been detected.
    % Initialize it to zero for the current frame. 
    AREASTATS.Area=0;

    Tempmatrixs=Movmatrixs;
    Tempmatrixs=double(Tempmatrixs);
    
    % Use the Lucy-Richardson deconvolution algorithm deconvolution if
    % deconvolution is specified.  
    if deconvolve == true
        J = deconvlucy(Tempmatrixs,PSF);
        maxj=max(max(J));
        maxtemp=max(max(Tempmatrixs));
        % Scale this image to the same dynamic range as the original image.
        Tempmatrixs=(J.*(maxtemp/maxj));
        Movmatrixs=uint8(Tempmatrixs);
     end

    % Denoising the image with a mean filter may help in some processing 
    % situations.  Helps to blur cell boundaries, resulting in less
    % frame-to-frame noise in object boundaries and hence shape and
    % orientation.
    if meanfilter == 1
        Tempmatrixs = imfilter(Tempmatrixs,meanfilterimage,'replicate');
    end
        
    % Clear the image to track from the previous frame.
    if trackgradient == false
        if outputframecounter < 2
            tempmatrixmaxval=max(max(Tempmatrixs));
            tempmatrixminval=min(min(Tempmatrixs));
        end
        % If the intensity will be used to track, invert the image to
        % take advantage of the dark rings at the beads edges when
        % thresholding.  Otherwise, use regular for tracking Phase or 
        % DIC image
        if invertimage==true
            FT=tempmatrixminval+tempmatrixmaxval-Tempmatrixs;
        else
            FT=-tempmatrixminval-tempmatrixmaxval+Tempmatrixs;
        end
    else
    % Use the sum of the absolute value of the gradient in the x and
    % y-directions as an indication where image intensity is changing to 
    % detect edges.  Note that this value is arbitrarily assigned to the 
    % upper-left hand corner of the 3 pixels involved in the subtraction, 
    % but this is consistently applied so the changes in position should 
    % still be accurate.
        [FX,FY] = gradient(Tempmatrixs);
        FT = sqrt((FX).^2+(FY).^2);
        
    end 
    if medianfilter == 1
        FT = medfilt2(FT,[medianfiltersize medianfiltersize]);
    end    
    % Scale the gradient image based on the maximum intensity of the
    % gradient in the first frame.
    if outputframecounter < 2
        ftmaxval=max(max(FT));
        ftminval=min(min(FT));
    end
    FT=(255/(ftmaxval-ftminval).*(FT-ftminval));    
    FT=uint8(FT);
    
    % Remember the current FT, before detected objects are subtracted out.
    % It will be used to help deduce the flow direction.
    FTcurrent=FT;

    % Find the threshold value for which there could potentially be enough
    % pixels to form a ring with the given diameter.
    Pixvals=sort(FT(:),1,'descend');
    % Remove degenerate entries
    Pixvals=unique(Pixvals);
    Pixvals=sort(Pixvals,1,'descend');
    [pixvalssize,temp]=size(Pixvals);
    % Intelligently select the index to start trying the threshold at based
    % on the last success.
    if outputframecounter < 2        
        threshvalind=1;
    else
        threshvalind = max(max(find(Pixvals>=startthresh)),1);
    end

    Currentbinaryoutput=im2bw((zeros(outputdheight,width)),0);
    clear Objecttracker
    objectsinframecounter=uint32(0);
    % Verify the minimum search level is not below the minimum restriction
    % set in the parameter file.
    if lastminthresh < absoluteminthresh*255
        lastminthresh = absoluteminthresh*255;
    end
    % Apply the current threshold and search for suitably-sized objects.    
    while ((threshvalind <= pixvalssize)) && (Pixvals(threshvalind)>=lastminthresh)
        objectsatthresh=0;
        % If this is the first time through the loop and the program is
        % starting with the default threshvalind of 1
        if threshvalind < 2
            threshval=Pixvals(threshvalind);
        % Otherwise check if this is the first time thresholding the image
        % but with the adaptive threshold calculated based on previous
        % frames.
        elseif threshvalind == max(max(find(Pixvals>=startthresh)),1)
            threshval = max(find(Pixvals>=startthresh));
            threshval = Pixvals(threshval);
            threshvalind = find(Pixvals==threshval);      
        % If not, decrease the threshold if it results in an acceptable
        % values
        elseif threshval - accelerationfactor > lastminthresh
            threshval = threshval-accelerationfactor;
            threshvalind = find(Pixvals<=threshval);
            threshvalind = threshvalind(1);
            threshval = Pixvals(threshvalind);
        % The only remaining option is this must be the last time through
        % the loop.
        else
            threshval=lastminthresh;
            threshvalind=find(Pixvals<=threshval);
            threshvalind=threshvalind(1);            
        end
        
        % Threshold the FT to look for objects at the current search
        % threshold.
        BW = (im2bw(FT,double(threshval-1)/255));
        % First fill any holes in the objects.
        BW1 = imfill(BW,'holes');
        % Now detect all of the objects at the current threshold value.
        BW2=bwlabel(BW1);
        nobjectsatthresh=max(max(BW2));
        
        % Scan each detected object to see if it meets acceptance criterion
        % for size and location
        for objectsatthreshcounter = 1 : nobjectsatthresh
            [objectrow,objectcolumn]=find(BW2==objectsatthreshcounter);
            [AREASTATS.Area,temp]=size(objectrow);
            % Check to see if the detected object is about the right size.
            if (AREASTATS.Area <= pixelsmax && AREASTATS.Area >= pixelsmin)
                BW3=bwselect(BW1,objectcolumn(1),objectrow(1));                
                % If area criterion is met, check whether the object does
                % not overlap an existing detection
                if max(max(uint8(BW3)+uint8(Currentbinaryoutput)))<=1;
                    % Verify the morphology of the object.  First reject
                    % if there are any lines going into the object.
                    BOUNDINGSTATS = regionprops(uint8(BW3),'BoundingBox');                    
                    BWbound=BW3(round(BOUNDINGSTATS.BoundingBox(2)):round(BOUNDINGSTATS.BoundingBox(2))+BOUNDINGSTATS.BoundingBox(4)-1,round(BOUNDINGSTATS.BoundingBox(1)):round(BOUNDINGSTATS.BoundingBox(1))+BOUNDINGSTATS.BoundingBox(3)-1);
                    BWboundinv=~BWbound;
                    BWneighbors=imfilter(uint8(BWboundinv),[0,1,0;1,0,1;0,1,0]);
                    BWneighbors=uint8(BWboundinv).*BWneighbors;
                    if max(max(BWneighbors(2:BOUNDINGSTATS.BoundingBox(4)-1,2:BOUNDINGSTATS.BoundingBox(3)-1) == 1)) == 0
                        % Use erosion to test the goodness of the
                        % detected objects shape.  Good shapes generally 
                        % erode to a single point.
                        Erodedbw=imerode(BWbound,[0,1,0;1,1,1;0,1,0]);
                        stoperode=false;
                        shapegood=false;
                        while stoperode == false
                            if max(max(bwlabel(Erodedbw))) > 1
                                stoperode=true;
                            elseif max(max(Erodedbw)) == 0
                                stoperode=true;
                                shapegood=true;
                            else
                                Erodedbw=imerode(Erodedbw,[0,1,0;1,1,1;0,1,0]);
                            end
                        end                      
                        if shapegood == true
                            % Check for shapes that might have
                            % interfered with the eroding shape check
                            if max(max(imfilter(uint8(bwperim(BWbound)),[1,1;1,1]))) < 4
                                Currentbinaryoutput=Currentbinaryoutput+BW3;
                                objectsatthresh=objectsatthresh+1;
                                objectsinframecounter=objectsinframecounter+1;
                                % Zero out these pixels in the original FT to avoid
                                % repeated object detections
                                FT=FT.*uint8(~BW3);
                                % Size in column 1
                                Objecttracker(objectsinframecounter,1)=sqrt(AREASTATS.Area/(pi*xres*yres));
                                % Current time in column 2
                                Objecttracker(objectsinframecounter,2)=double(outputframecounter-1)/outputfps;
                                CENTERSTATS = regionprops(uint8(BW3),'Centroid');
                                % X-Coordinate in column 3
                                Objecttracker(objectsinframecounter,3)=CENTERSTATS.Centroid(1)/xres;
                                % Y-Coordinate in column 4
                                Objecttracker(objectsinframecounter,4)=CENTERSTATS.Centroid(2)/yres;
                                % Object width in column 5
                                Objecttracker(objectsinframecounter,5)=BOUNDINGSTATS.BoundingBox(3)/xres;
                                % Object height in column 6
                                Objecttracker(objectsinframecounter,6)=BOUNDINGSTATS.BoundingBox(4)/yres;       
                                % Search for the object orientation in degrees for column 7
                                % First get the perimeter image
                                BWperimeterimage = bwperim(BW3,8);
                                [Perimeterpixelsrow,Perimeterpixelscolumn,Perimeterpixelsvalue]=find(BWperimeterimage);
                                Perimeterpixelscolumn=Perimeterpixelscolumn;
                                Perimeterpixelsrow=Perimeterpixelsrow;
                                % Now find the orientation angle & pixels
                                [majoraxisangle(objectsinframecounter),majoredgepixelx1(objectsinframecounter),majoredgepixely1(objectsinframecounter),majoredgepixelx2(objectsinframecounter),majoredgepixely2(objectsinframecounter),majoraxislength(objectsinframecounter)] = findorientation(Perimeterpixelsrow,Perimeterpixelscolumn,CENTERSTATS.Centroid(1),CENTERSTATS.Centroid(2),xres,yres);
                                Objecttracker(objectsinframecounter,7)=majoraxisangle(objectsinframecounter);
                                Objecttracker(objectsinframecounter,8)=majoraxislength(objectsinframecounter);
                                % Find the distances along the minor axis
                                % Define the minor axis to be perpendicular
                                % to the major axis
                                [minoredgepixelx1(objectsinframecounter),minoredgepixely1(objectsinframecounter),minoredgepixelx2(objectsinframecounter),minoredgepixely2(objectsinframecounter),minoraxislength(objectsinframecounter)] = findedgepointsonline(Perimeterpixelsrow,Perimeterpixelscolumn,CENTERSTATS.Centroid(1),CENTERSTATS.Centroid(2),xres,yres,majoraxisangle(objectsinframecounter)+90);
                                Objecttracker(objectsinframecounter,9) = minoraxislength(objectsinframecounter);
                                % Base the starting threshold in the next image on the current
                                % threshold where the first successful
                                % image detection occurs
                                if objectsinframecounter == 1
                                    startthresh=min(max(Pixvals),max(threshval+5,threshval+accelerationfactor));
                                else
                                    if minthresh-threshval>maxdeltathreshobjectdetect
                                        maxdeltathreshobjectdetect = minthresh-threshval;
                                    end
                                	if threshval < absoluteminthreshdetect
                                        absoluteminthreshdetect = threshval;
                                    end
                                end
                                minthresh=threshval;
                            end
                        end
                    end
                end
            end
        end
        threshvalind=threshvalind+1;
        % Go to the next threshval
    end
    lastminthresh=uint8(double(outputframecounter)/double(outputframecounter+1)*double(lastminthresh)+1/double(outputframecounter+1)*double(max(absoluteminthreshdetect-maxdeltathreshobjectdetect,0)));
    
    % Use the binary frames to determine the direction of flow based on
    % the shift in the centroid of images based on moving objects.
    if outputframecounter > 1
        Deltaforward = FTcurrent-FTlast;
        Deltabackward = FTlast-FTcurrent;
        Correlationmatrix=xcorr2(double(im2bw(Deltabackward)),double(im2bw(Deltaforward)));
        [Temp(2,:),Temp(1,:),temp]=find(Correlationmatrix==max(max(Correlationmatrix)));
        Temp=[Temp(1,1);Temp(2,1)];
        % Figure out the X-component of average motion 
        Temp(1)=(width-Temp(1))/xres;
        % Figure out the Y-component of average motion 
        Temp(2)=(outputdheight-Temp(2))/yres;
        % Make it a distance-unit normal vector
        Temp=Temp./(Temp(1)^2+Temp(2)^2)^.5;
        % Find the rolling average
        Flowdir=(1/(double(outputframecounter))).*Temp+(double(outputframecounter-1)/double(outputframecounter)).*Flowdir;
        Flowdir=Flowdir./(Flowdir(1)^2+Flowdir(2)^2)^.5;
        clear Temp
    end 
    
    % Special initialization operations for the first frame
    if outputframecounter < 2
        % Sort the objects in the first frame based on X-Coordinate
        Temp=[Objecttracker,majoraxisangle',majoredgepixelx1',majoredgepixely1',majoredgepixelx2',majoredgepixely2',majoraxislength',minoredgepixelx1',minoredgepixely1',minoredgepixelx2',minoredgepixely2',minoraxislength'];
        Temp=sortrows(Temp,3);
        Objecttracker=Temp(:,1:centertrackerdatasize);
        majoraxisangle=Temp(:,centertrackerdatasize+1)';
        majoredgepixelx1=Temp(:,centertrackerdatasize+2)';
        majoredgepixely1=Temp(:,centertrackerdatasize+3)';
        majoredgepixelx2=Temp(:,centertrackerdatasize+4)';
        majoredgepixely2=Temp(:,centertrackerdatasize+5)';
        majoraxislength=Temp(:,centertrackerdatasize+6)';
        minoredgepixelx1=Temp(:,centertrackerdatasize+7)';
        minoredgepixely1=Temp(:,centertrackerdatasize+8)';
        minoredgepixelx2=Temp(:,centertrackerdatasize+9)';
        minoredgepixely2=Temp(:,centertrackerdatasize+10)';
        minoraxislength=Temp(:,centertrackerdatasize+11)';                
        Objecttracker(1,centertrackerdatasize+1)=1;
        [nobjectsinfirstframe,temp]=size(Objecttracker);
        Centertracker=struct('Object1',Objecttracker(1,1:centertrackerdatasize));
        if nobjectsinfirstframe > 1
            for objectsinfirstframecounter = 2 : nobjectsinfirstframe
                fieldname=strcat('Object',num2str(objectsinfirstframecounter));
                Centertracker = setfield(Centertracker, fieldname, Objecttracker(objectsinfirstframecounter,1:centertrackerdatasize));
                Objecttracker(objectsinfirstframecounter,centertrackerdatasize+1)=objectsinfirstframecounter;
            end
        end
        totalobjectsfound=nobjectsinfirstframe;
        clear nobjectsinfirstframe temp objectsinfirstframecounter Temp
    else
    % If this is the second frame or greater then attempt to 
    % establish the object's identity based on the previous frame.
    % First find an object with a centroid in the previous frame.
        previoustestobjectcounter=uint32(0);
        Objectsnotpairedtracker=Objecttracker;
        for objectcounter = 1 : totalobjectsfound
            fieldname=strcat('Object',num2str(objectcounter));
            Temp=getfield(Centertracker,fieldname);
            % Check to see if the last entry for the object was for the
            % previous frame
            [temp1,temp2]=size(Temp);
            if Temp(temp1, 2) >= double(outputframecounter-2)/outputfps
                % Now that data for an object in the previous frame has
                % been identified, load its last position data for object 
                % matching.
                previoustestobjectcounter=previoustestobjectcounter+1;
                Centertemp(previoustestobjectcounter,1:centertrackerdatasize)=Temp(temp1,:);
                % Add a pointer to remember where this was in the
                % original Centertracker.
                Centertemp(previoustestobjectcounter,centertrackerdatasize+1)=objectcounter;
            end
            clear temp1 temp2 Temp
        end
        npreviousobjects=previoustestobjectcounter;
        clear previousobjectcounter
        % Now try to match all of the identified objects with ones in the
        % current frame.  Begin by calculating pair-wise distance metrics
        [ncurrentframeobjects,temp]=size(Objecttracker);
        clear temp
        for centertrackercounter = 1 : npreviousobjects
            for objectcounter = 1 : ncurrentframeobjects
                Delta(1,1)=Objectsnotpairedtracker(objectcounter,3)-Centertemp(centertrackercounter,3);
                Delta(2,1)=Objectsnotpairedtracker(objectcounter,4)-Centertemp(centertrackercounter,4);
                % Calculate the angle between the Delta vector and the flow
                % direction
                theta=acos(dot(Delta,Flowdir)./((sum(Delta.^2))^.5*(sum(Flowdir.^2))^.5));
                % Now compute the component of the Delta vector
                % perpendicular to the flow direction
                perpendicular=dot(Delta,Flowdir)*tan(theta);
                % Verify the second object is in approximately the line of
                % flow path of that from the previous frame, if not dismiss 
                % it as a possibility
                if abs(perpendicular) > beadsize/2 
                    Comparematrix(centertrackercounter,objectcounter,1)=NaN;
                    Comparematrix(centertrackercounter,objectcounter,2)=NaN;
                else                    
                    % Compute the total distance
                    Comparematrix(centertrackercounter,objectcounter,2)=(Delta(1,1)^2+Delta(2,1)^2)^.5;
                    % Verify the second object is downstream or
                    % shifted less than 1 radius upstream (allow for some 
                    % noise)
                    Upstreampoint(1,1)=Centertemp(centertrackercounter,3)-beadsize/2*Flowdir(1);
                    Upstreampoint(2,1)=Centertemp(centertrackercounter,4)-beadsize/2*Flowdir(2);
                    Delta(1,1)=Objecttracker(objectcounter,3)-Upstreampoint(1,1);
                    Delta(2,1)=Objecttracker(objectcounter,4)-Upstreampoint(2,1);
                    theta=acos(dot(Delta,Flowdir)./((sum(Delta.^2))^.5*(sum(Flowdir.^2))^.5));
                    if theta > pi/2
                        Comparematrix(centertrackercounter,objectcounter,1)=NaN;
                        Comparematrix(centertrackercounter,objectcounter,2)=NaN;
                    else
                        Comparematrix(centertrackercounter,objectcounter,1)=perpendicular;
                    end
                end
            end
        end
        clear Delta Upstreampoint perpendicular        
        % Now that all of the metrics are established begin pairing objects
        % with previous center tracker entries.
        Objecttracker(:,centertrackerdatasize+1)=NaN;
        while min(min(isnan(Comparematrix(:,:,1))))<1
            [mindrow,mindcolumn,mindv]=find(Comparematrix(:,:,2)==min(min(Comparematrix(:,:,2))));
            field=strcat('Object',num2str(Centertemp(mindrow(1),centertrackerdatasize+1)));
            Centertracker = setfield(Centertracker,field,cat(1,getfield(Centertracker,field),Objecttracker(mindcolumn(1),1:centertrackerdatasize)));
            % Also add a reference in the Objecttracker to the object
            % ID, needed for image processing
            Objecttracker(mindcolumn(1),centertrackerdatasize+1)=Centertemp(mindrow(1),centertrackerdatasize+1);
            % Eliminate the entry from the Objecttracker and
            % Comparematrix
            Objectsnotpairedtracker(mindcolumn(1),:)=NaN;
            Comparematrix(mindrow(1),:,1)=NaN;
            Comparematrix(mindrow(1),:,2)=NaN;
            Comparematrix(:,mindcolumn(1),1)=NaN;
            Comparematrix(:,mindcolumn(1),2)=NaN;
        end
        % Now add any detected objects that could not be paired as new
        % entries to the Centertracker.
        while min(min(isnan(Objectsnotpairedtracker)))<1
            newobjectrow=find(Objectsnotpairedtracker(:,3)==min(Objectsnotpairedtracker(:,3)));
            totalobjectsfound=totalobjectsfound+1;
            field=strcat('Object',num2str(totalobjectsfound));
            Centertracker=setfield(Centertracker,field,Objectsnotpairedtracker(newobjectrow,:));
            Objecttracker(newobjectrow,centertrackerdatasize+1)=totalobjectsfound;
            Objectsnotpairedtracker(newobjectrow,:)=NaN;
        end
        clear Centertemp Comparematrix mindrow mindcolumn mindv Objectsnotpairedtracker newobjectrow prevrow prevcolumn prevv
    end
        
    % Perform output frame and file writing operations
    % Movmatrixs: Output frame image without tracking modifications
    % Currentbinaryoutput: Binary image of objects in Movmatrixs
    % Currentobjectbinary: Binary image of selected object in Movmatrixs
    % Currentobjectcolor: Color image of selected object with data from
    % tracker
    % Currentobjectinv: Inverted image of binary object for selecting the
    % background region
    % Movmatrixrgb: Final modified version of Movmatrixs.
    cmap=gray(256);
    Movmatrixs=ind2rgb(Movmatrixs,cmap);     
    [ncurrentframeobjects,temp]=size(Objecttracker);
    clear temp
    % First write the streaklines to the current image
    for objectcounter = 1 : ncurrentframeobjects
        % Check to see whether the object is in more than 1 frame
        fieldname=strcat('Object',num2str(Objecttracker(objectcounter,centertrackerdatasize+1)));
        Centertemp=getfield(Centertracker,fieldname);
        [npreviousframes,temp]=size(Centertemp);
        if npreviousframes > 1
           for lineframecounter = 2 : npreviousframes
              slope=(Centertemp(lineframecounter,4)-Centertemp(lineframecounter-1,4))*yres/((Centertemp(lineframecounter,3)-Centertemp(lineframecounter-1,3))*xres);
              intercept=Centertemp(lineframecounter,4)*yres-slope*Centertemp(lineframecounter,3)*xres;
              if slope < 1 && slope > -1
                  if round(Centertemp(lineframecounter,3)*xres) >= round(Centertemp(lineframecounter-1,3)*xres)
                    X = (round(Centertemp(lineframecounter-1,3)*xres) : round(Centertemp(lineframecounter,3)*xres))';
                  else
                      X= (round(Centertemp(lineframecounter,3)*xres) : round(Centertemp(lineframecounter-1,3)*xres))';
                  end
                Y = (round(X.*slope+intercept));
              else
                  if round(Centertemp(lineframecounter,4)*yres) >= round(Centertemp(lineframecounter-1,4)*yres)
                    Y = (round(Centertemp(lineframecounter-1,4)*yres) : round(Centertemp(lineframecounter,4)*yres))';
                    if isnan(slope)==1 || isinf(slope) == 1
                        X = ones(length(Y),1)*round(Centertemp(lineframecounter-1,3)*xres);
                    else
                        X = (round((Y-intercept)./slope));
                    end
                  else
                      Y = (round(Centertemp(lineframecounter,4)*yres) : round(Centertemp(lineframecounter-1,4)*yres))';
                      if isnan(slope)==1 || isinf(slope) == 1
                        X = ones(length(Y),1)*round(Centertemp(lineframecounter-1,3)*xres);  
                      else
                        X = (round((Y-intercept)./slope));
                      end
                 end
              end             
              Z = ones(size(X))*2;
              Movmatrixs(sub2ind(size(Movmatrixs),Y,X,Z))=1;
           end
        end
    end
    Movmatrixrgb=Movmatrixs;
    clear npreviousframes temp slope intercept X Y Z Centertemp
    
    % Now color the objects to be written to file
    for objectcounter = 1 : ncurrentframeobjects
        Currentobjectbinary=bwselect(Currentbinaryoutput,majoredgepixelx1(objectcounter),majoredgepixely1(objectcounter));
        % Find the left, right, bottom, and top pixels of the binary image
        % to reference for the line writing.
        BOUNDINGSTATS = regionprops(uint8(Currentobjectbinary),'BoundingBox');
        leftpixel=round(BOUNDINGSTATS.BoundingBox(1));
        rightpixel=round(BOUNDINGSTATS.BoundingBox(1))+BOUNDINGSTATS.BoundingBox(3)-1;
        toppixel=round(BOUNDINGSTATS.BoundingBox(2));
        bottompixel=round(BOUNDINGSTATS.BoundingBox(2))+BOUNDINGSTATS.BoundingBox(4)-1;
        % Find the pixel locations of the centroid of the current
        % object.
        nverpointo=Objecttracker(objectcounter,4)*yres;
        nhorpointo=Objecttracker(objectcounter,3)*xres;
        % Get the object ID number of the current object.
        objectid=Objecttracker(objectcounter,centertrackerdatasize+1);
        % Create a binary object id number image
        decades=round((log10(objectid))+.5);
        Objectidbinary=Numpic(:,:,mod(objectid,10)+1);
        [objectidbinaryrows,objectidbinarycolumns,objectidbinarypages]=size(Objectidbinary);
        if decades > 1
            for decadecounter = 2: decades
                Objectidbinary=cat(2,zeros(objectidbinaryrows,1),Objectidbinary);
                Objectidbinary=cat(2,Numpic(:,:,ceil(mod(objectid/10^(decadecounter-1),10)+1/(10^(decades)))),Objectidbinary);
            end
        end
        [objectidbinaryrows,objectidbinarycolumns,objectidbinarypages]=size(Objectidbinary);
        % Add the label onto the current objectid binary.
        if leftpixel >= width/2
            if bottompixel <= outputdheight/2
                % Write to lower left corner
                objectidrightindex=leftpixel-1;
                objectidtopindex=bottompixel+1;
                objectidleftindex=objectidrightindex-objectidbinarycolumns+1;
                objectidbottomindex=objectidtopindex+objectidbinaryrows-1;
            else
                % Write to upper left corner
                objectidrightindex=leftpixel-1;
                objectidbottomindex=toppixel-1;
                objectidleftindex=objectidrightindex-objectidbinarycolumns+1;
                objectidtopindex=objectidbottomindex-objectidbinaryrows+1;                
            end
        else
            if bottompixel <= outputdheight/2
                % Write to lower left corner
                objectidleftindex=rightpixel+1;
                objectidtopindex=bottompixel+1;
                objectidrightindex=objectidleftindex+objectidbinarycolumns-1;
                objectidbottomindex=objectidtopindex+objectidbinaryrows-1;                   
            else
                % Write to upper left corner
                objectidleftindex=rightpixel+1;
                objectidbottomindex=toppixel-1;
                objectidrightindex=objectidleftindex+objectidbinarycolumns-1;
                objectidtopindex=objectidbottomindex-objectidbinaryrows+1;                 
            end
        end
        if objectidleftindex > 0 && objectidtopindex > 0 && objectidbottomindex <= outputdheight && objectidrightindex <= width
            Currentobjectbinary(objectidtopindex:objectidbottomindex,objectidleftindex:objectidrightindex)=(Objectidbinary);
        end
        Currentobjectinv(:,:,1)=double(~(Currentobjectbinary));
        Currentobjectinv(:,:,2)=double(~(Currentobjectbinary));
        Currentobjectinv(:,:,3)=double(~(Currentobjectbinary));
        % Reference strain index colors to a jet colormap
        cmapmax=256;
        cmapj=jet(cmapmax);
        cmapg=gray(256);
        % Calculate the shape strain index for the current object
        shapestrainindex=Objecttracker(objectcounter,8)/Objecttracker(objectcounter,9);
        % Index color according to morphological strain index
        if shapestrainindex >= 1.5
            colorindex=cmapmax;
        else
            colorindex=round((shapestrainindex-1)*(cmapmax-1)/(1.5-1));
        end
        Currentobjectcolor=uint8(Currentobjectbinary);
        Currentobjectcolor=colorindex*Currentobjectcolor;
        Currentobjectcolor=ind2rgb((Currentobjectcolor),cmapj);
        % Color the object from the original image and brighten if desired        
        Currentobjectcolor=imadjust(Movmatrixs,[0;1],[0;1],0.5).*Currentobjectcolor.*cat(3,cat(3,double(Currentobjectbinary),double(Currentobjectbinary)),double(Currentobjectbinary));     
        % Add in lines to show the strain axes.
        if showline == true
            for strainaxiscounter = 1 : 2
                if strainaxiscounter == 1
                    edgepixely1=majoredgepixely1(objectcounter);
                    edgepixelx1=majoredgepixelx1(objectcounter);
                    edgepixely2=majoredgepixely2(objectcounter);
                    edgepixelx2=majoredgepixelx2(objectcounter);
                    Linecolor=[.8,0,0];
                else
                    edgepixely1=minoredgepixely1(objectcounter);
                    edgepixelx1=minoredgepixelx1(objectcounter);
                    edgepixely2=minoredgepixely2(objectcounter);
                    edgepixelx2=minoredgepixelx2(objectcounter);
                    Linecolor=[0,0,.8];                    
                end
                [temprows,tempcolumns] = size(Currentbinaryoutput);
                slopeaxis=(edgepixely1-nverpointo)/(edgepixelx1-nhorpointo);
                interceptaxis=(nverpointo)-slopeaxis*(nhorpointo);   
                columndecimal=(nhorpointo)-round(nhorpointo);
                rowdecimal=(nverpointo)-round(nverpointo);
                if slopeaxis < 1 && slopeaxis > -1
                    % Draw a more horizonal line
                    if edgepixelx1 <= edgepixelx2
                        spacing = 1;
                    else
                        spacing = -1;
                    end
                    X = edgepixelx1 : spacing : edgepixelx2;
                    X = X'+columndecimal;
                    X = sort(X);
                    linelength=length(X);
                    if round(X(1))<1
                        X = X(2:linelength);
                    end
                    linelength=length(X);
                	if round(X(linelength)) > width
                        X = X(1:linelength-1);
                    end
                    Y = X.*slopeaxis + interceptaxis;
                    linelength = length(Y);
                    if round(Y(1))<1 || round(Y(1))>width
                        X = X(2:linelength);
                        Y = Y(2:linelength);
                    end
                    linelength = length(Y);
                	if round(Y(linelength)) > width || round(Y(linelength)) < 1
                        X = X(1:linelength-1);
                        Y = Y(1:linelength-1);
                    end
                    linelength = length(Y);
                    Z = ones(linelength,1);
                else
                % Draw a more vertical line
                    if edgepixely1 <= edgepixely2
                        spacing = 1;
                    else
                        spacing = -1;
                    end                    
                    Y = edgepixely1 : spacing : edgepixely2;
                    Y = Y'+rowdecimal;
                    Y = sort(Y);                    
                    linelength=length(Y);
                    if round(Y(1))<1
                        Y = Y(2:linelength);
                    end
                    linelength=length(Y);
                	if round(Y(linelength)) > width
                        Y = Y(1:linelength-1);
                    end
                    linelength = length(Y);
                    % Make sure line is not perfectly vertical
                    if isnan(slopeaxis)==1 || isinf(slopeaxis) == 1
                        X = ones(linelength,1)*edgepixelx1;
                    else
                        X = (Y - interceptaxis)./slopeaxis;
                    end
                    linelength = length(X);
                    if round(X(1))<1 || round(X(1))>width
                        X = X(2:linelength);
                        Y = Y(2:linelength);
                    end
                    linelength = length(Y);
                	if round(X(linelength)) > outputdheight || round(X(linelength)) < 1
                        X = X(1:linelength-1);
                        Y = Y(1:linelength-1);
                    end
                    linelength = length(Y);
                end
                Z = ones(linelength,1);                    
                Currentobjectcolor(sub2ind(size(Currentobjectcolor),round(Y),round(X),Z))=Linecolor(1);
                Currentobjectcolor(sub2ind(size(Currentobjectcolor),round(Y),round(X),Z*2))=Linecolor(2);
                Currentobjectcolor(sub2ind(size(Currentobjectcolor),round(Y),round(X),Z*3))=Linecolor(3);
                % Be sure the axis points are marked
                Currentobjectcolor(edgepixely1,edgepixelx1,1)=Linecolor(1);
                Currentobjectcolor(edgepixely1,edgepixelx1,2)=Linecolor(2);
                Currentobjectcolor(edgepixely1,edgepixelx1,3)=Linecolor(3);
                Currentobjectcolor(edgepixely2,edgepixelx2,1)=Linecolor(1);
                Currentobjectcolor(edgepixely2,edgepixelx2,2)=Linecolor(2);
                Currentobjectcolor(edgepixely2,edgepixelx2,3)=Linecolor(3);                
            end
        end    
        % Highlight the approximate location of the centroid in white.
        Currentobjectcolor(round(nverpointo),round(nhorpointo),:)=[1,1,1];
        % Zero out any pixels that may have been included from the
        % line calculation but are on the outside of the object
        Currentobjectcolor=Currentobjectcolor.*(-1*Currentobjectinv+1);
        % Add the colored current object to the Movmatrixrgb.
        Movmatrixrgb=(Movmatrixrgb.*Currentobjectinv)+Currentobjectcolor;
        clear X Y Z
    end
    % Add the header data back on
    if headertop~=false
        Movmatrixrgb=cat(1,ind2rgb(Headertop,cmapg),Movmatrixrgb);
    end
    if headerbottom~=false
        Movmatrixrgb=cat(1,Movmatrixrgb,ind2rgb(Headerbottom,cmapg));
    end    

    frame = im2frame(Movmatrixrgb);
    aviobj = addframe(aviobj,frame);

    % Let the user know every 1% increase in the processing status.
    if round(double(outputframecounter)/double(noutputmovieframes)*100-.5)>lastinformeduser;
        lastinformeduser = round(double(outputframecounter)/double(noutputmovieframes)*100-.5);
        fprintf('Currently processing output frame %g of %g, %g%% complete.  ', outputframecounter, nsubframes*ninputmovieframes, lastinformeduser)
        temp=toc;
        fprintf('Elapsed time is %u hr %u min %u s. \n', round(temp/3600-.5), round(mod(round(temp/60-.5),60)), round(mod(temp,60)))
        % Save the center tracker every 10 frames in case execution is 
        % terminated early. 
        save(centertrackerfile,'Centertracker')
    end
                
    % Update variables for the to recall for the output frame
    Lastbinaryoutput=Currentbinaryoutput;
    outputframecounter=outputframecounter+1;
    FTlast=FTcurrent;
    % Go to the next input "subframe" and output frame
    end
    % Update variables to recall for the input movie frame
    if videoReaderin==true
        next(vrvideo);
    end
    inputframecounter=inputframecounter+1;
    % Go to the next input movie frame
end

aviobj = close(aviobj);
% Now save the output of the centroid location
save(centertrackerfile,'Centertracker')
% This can be reopened with the load command  
end
  
function [majoraxisangle,majoredgepixelx1,majoredgepixely1,majoredgepixelx2,majoredgepixely2,majoraxislength] = findorientation(Perimeterpixelsrow,Perimeterpixelscolumn,xcentroidpixel,ycentroidpixel,xres,yres)
    % This function tests the centroid against the perimeter pixels to find 
    % the one that is located the furthest.
    % Convert pixels to distances
    Xcoordinates=Perimeterpixelscolumn./xres;
    Ycoordinates=Perimeterpixelsrow./yres;
    xcentroid=xcentroidpixel/xres;
    ycentroid=ycentroidpixel/yres;
    % Calculate the angles each pixel makes with the centroid 
    Angles=180/pi*atan2((Ycoordinates-ycentroid),(Xcoordinates-xcentroid));
    % Initialize
    longestlength=0;
    % Test the pixel values

    [npixels,temp]=size(Perimeterpixelsrow);
    for counter = 1 : npixels
        angledifferencemin=180;        
        xcurrent=Xcoordinates(counter);
        ycurrent=Ycoordinates(counter);
        anglecurrent=Angles(counter);
        if anglecurrent >= 0
            oppositeangle = anglecurrent-180;
        elseif anglecurrent < 0
            oppositeangle = anglecurrent+180;
        end
        angledifferencemin=180;
        for anglecounter = 1 : npixels
            angledifference=abs(Angles(anglecounter)-oppositeangle);
            if angledifference < angledifferencemin
                oppositeindex=anglecounter;
                angledifferencemin=angledifference;
                xopposite=Xcoordinates(anglecounter);
                yopposite=Ycoordinates(anglecounter);
            end
        end         
        length=((xcurrent-xopposite)^2+(ycurrent-yopposite)^2)^.5;
        if length > longestlength
            longestlength=length;
            bestindex=counter;
            bestoppositeindex=oppositeindex;
        end        
    end
    % Calculate the angle 
     majoraxisangle=Angles(bestindex);    
    % Restrict this angle to between +/- 90º
    if majoraxisangle > 90
        majoraxisangle=majoraxisangle-180;
    elseif majoraxisangle<90
        majoraxisangle=majoraxisangle+180;
    end
    % Record the values for reporting back to the program
    majoredgepixelx1=Perimeterpixelscolumn(bestindex);
    majoredgepixely1=Perimeterpixelsrow(bestindex);
    majoredgepixelx2=Perimeterpixelscolumn(bestoppositeindex);
    majoredgepixely2=Perimeterpixelsrow(bestoppositeindex);
    majoraxislength=(((majoredgepixelx1-majoredgepixelx2)/xres)^2+((majoredgepixely1-majoredgepixely2)/yres)^2)^.5;
end
    
function [edgepixelx1,edgepixely1,edgepixelx2,edgepixely2,axislength] = findedgepointsonline(Perimeterpixelsrow,Perimeterpixelscolumn,xcentroidpixel,ycentroidpixel,xres,yres,angle)
    % Restrict the angle of search to between -180º and +180º
    if angle > 180
        angle=angle-180;
    elseif angle <= -180
        angle=angle+180;
    end
    % Convert pixels to distances
    Xcoordinates=Perimeterpixelscolumn./xres;
    Ycoordinates=Perimeterpixelsrow./yres;
    xcentroid=xcentroidpixel/xres;
    ycentroid=ycentroidpixel/yres;
    [npixels,temp]=size(Perimeterpixelsrow);
    % Calculate the angles each pixel makes with the centroid 
    Angles=180/pi*atan2((Ycoordinates-ycentroid),(Xcoordinates-xcentroid));
    % Find the corresponding edge pixels along the major axis on the less
    % stretched side
    angledifferencemin=180;
    for anglecounter = 1 : npixels
        angledifference=abs(Angles(anglecounter)-angle);
        if angledifference < angledifferencemin
            bestcounter=anglecounter;
            angledifferencemin=angledifference;
        end
    end 
    % Now look on the other side for the nearest edge pixel.
    if angle > 0
        oppositeangle=angle-180;
    elseif angle <= 0
        oppositeangle=angle+180;
    end
    angledifferencemin=180;
    for anglecounter = 1 : npixels
        angledifference=abs(Angles(anglecounter)-oppositeangle);
        if angledifference < angledifferencemin
            oppositecounter=anglecounter;
            angledifferencemin=angledifference;
        end
    end     
    % Record the values for reporting back to the program
    edgepixelx1=Perimeterpixelscolumn(bestcounter);
    edgepixely1=Perimeterpixelsrow(bestcounter);
    edgepixelx2=Perimeterpixelscolumn(oppositecounter);
    edgepixely2=Perimeterpixelsrow(oppositecounter);
    axislength=(((edgepixelx1-edgepixelx2)/xres)^2+((edgepixely1-edgepixely2)/yres)^2)^.5;
end
