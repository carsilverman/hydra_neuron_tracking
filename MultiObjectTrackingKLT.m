% MultiObjectTrackingKLT by Carolyn Silverman
%
% Tracks multiple objects in real time based on MATLAB's
% built-in KLT point tracker. The output video stream places a
% bounding box around all of the objects based on a brightness
% threshold and a '+' marker on the approximate centroids of 
% the objects according to the KLT tracker. 
% NOTE: As in stands, this program cannot associate the correct pixel  
% list with the tracked centroids

%% Initialization
vidDevice = imaq.VideoDevice('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode', ... % Acquire input video stream
                    'ROI', [1 1 512 512], ...
                    'ReturnedColorSpace', 'rgb');
vidInfo = imaqhwinfo(vidDevice); % Acquire input video property
htextins = vision.TextInserter('Text', 'No Objects', ... % Set text for when no objects exist
                                    'Location',  [7 2], ...
                                    'Color', [1 0 0], ... // red color
                                    'FontSize', 12);
hVideoIn = vision.VideoPlayer('Name', 'Final Video', ... % Output video player
                                'Position', [100 100 vidInfo.MaxWidth+20 vidInfo.MaxHeight+30]);
nFrame = 1; % Frame number initialization
last = 0; % Variable that defines the number of objects in the last frame; initialize to 0
pointTracker = vision.PointTracker; %Create tracker object
keepRunning = true;

%% 
disp('Press Ctrl-C to exit');
while(keepRunning)
    
    I = step(vidDevice); % Acquire single frame

    thresholdValue = 0.15; %brightness threshold
    binaryImage = im2bw(I, thresholdValue); %convert to binary
    binaryImage = imfill(binaryImage, 'holes'); %fill in the holes
    seD = strel('diamond', 1); %structural element for further filtering
    BWfinal = imerode(binaryImage,seD); %erode the image based on the strutural element
    BWfinal = imerode(BWfinal,seD); %repeat
    BW2 = bwareaopen(BWfinal,600); %remove blobs smaller than 600 pixels

    %determine region properties
    blobMeasurements = regionprops(BW2, 'BoundingBox','Centroid','PixelList');
    %Acquire region properties for the blobs. Can add more properties or change to 'all'
    
    numberOfBlobs = size(blobMeasurements, 1); %determine the number of distinct objects

    %for each blob determine the centroid; create an M-by-2 matrix array of centroids
    centroids = zeros(numberOfBlobs,2); %new array of centroids
    for k=1 : numberOfBlobs
        blobCentroid = blobMeasurements(k).Centroid; %determine centroid for blob k
        blobBBox = blobMeasurements(k).BoundingBox; %determine bounding box for blob k
        centroids(k,:) = blobCentroid; %add the new centroid value   
        label = ' ';
        out2 = insertObjectAnnotation(I, 'rectangle', blobBBox, label, ...
                          'Color', [1 0 0]); %Annotate each object with a red bounding box
        if k==1
            out = out2;
        else
            out = out+out2; %Add new annotations to old annotations
        end
    end
    
    if numberOfBlobs ~= 0 %skip tracking if there are no objects in the field of view
        if nFrame == 1
            initialize(pointTracker, centroids, I); %initialize the tracker based on centroids
        else
            if last ~= numberOfBlobs %number of points has changed
                setPoints(pointTracker, centroids); %reset the points to be tracked
            end
            [points,point_validity] = step(pointTracker,I); %track the points;
            %points is an M-by-2 array that corresponds to the new location
            %of the points; point_validity provides an M-by-1 logical array, 
            %indicating whether or not each point has been reliably tracked
            
            colors = hsv(numberOfBlobs); %create an M-by-3 matrix color map for the output viedo stream
            
            for k=1:numberOfBlobs %iterate through all of the objects
                blobID = ['blob: ', num2str(k)];
                blobCentroid = blobMeasurements(k).Centroid; %determine centroid for blob k
                blobCentroid = uint16(blobCentroid); %convert to an integer to display
                blobBBox = blobMeasurements(k).BoundingBox; %determine bounding box for blob k
                blobPixels = blobMeasurements(k).PixelList; %determine the pixels
                centX = blobCentroid(1,1); centY = blobCentroid(1,2); %separate X and Y values
                label = ['Blob Number: ' num2str(k) ' Centroid: ' ...
                    num2str(points(k,1)) ' ' num2str(points(k,2))]; %Create label for output screen
                out3 = insertText(I, points(k,:), blobID);
                out4 = insertMarker(I, points(k,:), '+' ,'Color', colors(k,:), 'Size', 10); %add marker label
                out = out + out3 + out4;
            end  
        end
    end    
   
    
    if nFrame > 1
         if numberOfBlobs==0
            vidIn = step(htextins, I); % display "No Objects"
            step(hVideoIn, vidIn); % Output video stream
         else
             step(hVideoIn, out); % Output video stream with annotations
         end
    end
    last = size(blobMeasurements, 1); % make a copy of the number of objects in
    %the most recent frame to compare to future frames
    nFrame = nFrame+1; %Increment the frame number    
end

%% Clearing Memory
release(hVideoIn); % Release all memory and buffer used
release(vidDevice);
% clear all;
clc;
clear;