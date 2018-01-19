%SingleObjectTracker by Carolyn Silverman
%
%Tracks a single object in real time. Acquires a list of the brightest pixels
%in each frame and outputs a video stream with the annotated object.

%% Initialization
vidDevice = imaq.VideoDevice('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode', ... % Acquire input video stream
                    'ROI', [1 1 512 512], ...
                    'ReturnedColorSpace', 'rgb');
vidInfo = imaqhwinfo(vidDevice); % Acquire input video property
htextins = vision.TextInserter('Text', 'No objects ', ... % Set text for when no objects exist
                                    'Location',  [7 2], ... // upper left corner
                                    'Color', [1 0 0], ... // red color
                                    'FontSize', 12);
hVideoIn = vision.VideoPlayer('Name', 'Final Video', ... % Output video player
                                'Position', [100 100 vidInfo.MaxWidth+20 vidInfo.MaxHeight+30]);
nFrame = 1; % Frame number initialization
keepRunning = true;

%% Processing Loop
disp('Press Ctrl-C to exit');
while keepRunning
    
    I = step(vidDevice); % Acquire single frame
    
    %create a binary image of blobs using filters
    thresholdValue = 0.15; %brightness threshold
    binaryImage = im2bw(I, thresholdValue); %convert to binary
    binaryImage = imfill(binaryImage, 'holes'); %fill in the holes
    seD = strel('diamond', 1); %structural element for further filtering
                               %may want to change to disk
    BWfinal = imerode(binaryImage,seD); %erode the image based on the strutural element
    BWfinal = imerode(BWfinal,seD); %repeat
    BW2 = bwareaopen(BWfinal,600); %remove blobs smaller than 600 pixels

    blobMeasurements = regionprops(BW2, 'BoundingBox','Centroid','PixelList');
    %Acquire region properties for the blobs. Can add more properties or change to 'all'
    
    numberOfBlobs = size(blobMeasurements, 1); %determine the number of distinct objects
    colors = hsv(numberOfBlobs); %create an M-by-3 matrix color map for the output viedo stream
    

    for k=1 : numberOfBlobs
        blobID = k; %Set blob identity based on place on screen
        thisBlobsPixels = blobMeasurements(k).PixelList; %acquire pixels for blob k
        blobCentroid = blobMeasurements(k).Centroid; %determine centroid for blob k
        blobBBox = blobMeasurements(k).BoundingBox; %determine bounding box for blob k
        
        %Uncomment this section of code to display info to the screen
        %fprintf('Blob Number %d\n', k);
        %fprintf('Centroid: %.2f\n', blobCentroid);
        %fprintf('Pixel List\n');
        %disp(thisBlobsPixels);
        
        blobCentroid = uint16(blobCentroid); %convert to an integer for display
        centX = blobCentroid(1,1); centY = blobCentroid(1,2); %separate X and Y values
        colors(k,:); %Change color based on blobID
        label = ['Blob Number: ' num2str(blobID) ' Centroid: ' ...
                num2str(centX) ' ' num2str(centY)]; %Create label for output screen
        vidIn = insertObjectAnnotation(I, 'rectangle', blobBBox, label, ...
                'Color', colors(blobID,:)); %Annotate objects with bounding box and label
        if k==1
            vidIn2 = vidIn;
        else
            vidIn2 = vidIn2 + vidIn; %Add new annotations to old annotations
        end
    end
    
    if numberOfBlobs==0
        vid = step(htextins, I); %display "No Objects"
        step(hVideoIn, vid); %Output video stream 
    else
        step(hVideoIn, vidIn2); %Output video stream
    end
    nFrame = nFrame+1; %Increment the frame number
end

%% Clearing Memory
release(hVideoIn); % Release all memory and buffer used
release(vidDevice);
% clear all;
clc;
clear;