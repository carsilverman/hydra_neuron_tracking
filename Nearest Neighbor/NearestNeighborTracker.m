% NearestNeighborTracker by Carolyn Silverman
%
% Tracks multiple objects in real time based on the 
% nearestneigborlinker function written by Jean-Yves Tinevez

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

    blobMeasurements = regionprops(BW2, 'BoundingBox','Centroid','PixelIdxList');
    %Acquire region properties for the blobs. Can add more properties or change to 'all'
    
    numberOfBlobs = size(blobMeasurements, 1); %determine the number of distinct objects
    colors = hsv(numberOfBlobs); %create an M-by-3 matrix color map

    centroids = NaN(numberOfBlobs,2); %new array of centroids
    pixels = cell(numberOfBlobs, 1); %new 3D array of pixels; 
         %may increase or decrease 2nd dimension depending on size of object being tracked
    bboxes = NaN(numberOfBlobs,4); %new 2D array of bounding boxes
    
    for k=1 : numberOfBlobs %iterate through all the objets
        centroids(k,:) = blobMeasurements(k).Centroid; %add the centroid for blob k
        pixels{k} = blobMeasurements(k).PixelIdxList; %add the pixels for blob k
        bboxes(k,:) = blobMeasurements(k).BoundingBox; %add the bounding box for blob k
    end
     
    if numberOfBlobs ~= 0 %skip tracking if there are no objects in the field of view
        
        if nFrame == 1 || last == 0 %this is the first frame or there were 
                                    %no objects detected in the previous frame
            lastCentroids = zeros(numberOfBlobs,2); %create empty array of centroids from the last frame
            for k=1 : numberOfBlobs
                blobIndex = k; %set blob identity based on place on screen
                blobCentroid = centroids(k,:); %determine centroid for blob k
                blobPixels = pixels{k}; %acquire pixels for blob k
                blobBBox = bboxes(k,:); %determine bounding box for blob k
                blobCentroid = uint16(blobCentroid); %convert the centroid to an integer for display
                centX = blobCentroid(1,1); centY = blobCentroid(1,2); %separate X and Y values
                lastCentroids(k,:) = centroids(k,:);
                label = ['Blob Number: ' num2str(blobIndex) ' Centroid: ' ...
                        num2str(centX) ' ' num2str(centY)]; %Create label for output screen
                out2 = insertObjectAnnotation(I, 'rectangle', blobBBox, label, ...
                        'Color', colors(blobIndex,:)); %Annotate objects with bounding box and label
                if k==1
                    out = out2;
                else
                    out = out + out2; %Add new annotations to old annotations
                end
            end
            
        else
            %Use nearestneighborlinker to link two sets of centroids. 
            %target_indices is an M-by-2 array of centroids that have been
            %successfully tracked. unassigned_targets is an M-by-2 array of
            %centroids that have been detected in the current frame but 
            %cannot be linked to the previous frame. target_distances
            %returns the distance to the matched target point.
            
            [target_indices, target_distances, unmatched_targets] = ...
                 nearestneighborlinker(lastCentroids, centroids);
            numTargets = numel(target_indices); %number of targets, including missing elements
            numAssigned = sum(target_indices ~= -1); %number of linked centroids
            numUnassigned = numel(unmatched_targets); %number of unlinked centroids
            numElements = numAssigned + numUnassigned; %total number of centroids in frame
            k=1;
            while k<=numElements
                blobID = k; %set the index for identifying pixels and bounding boxes      
                if k <= numTargets %blob has been detected in previous frame
                    blobIndex = target_indices(k);
                    while blobIndex == -1
                        k = k+1;
                        blobIndex = target_indices(k);      
                    end
                    objCentroid(1,:) = centroids(blobIndex,:); %find centroid in array  
                else %blob has not been detected in previous frames
                    j = k-numTargets; 
                    blobIndex = unmatched_targets(j); 
                    objCentroid(1,:) = centroids(blobIndex,:); %find centroid in array 
                end         
                blobPixels = pixels{blobIndex}; %get pixels corresponding to blobID
                blobBBox = bboxes(blobIndex,:); %get bounding box corresponding to blobID
                objCentroid = uint16(objCentroid); %convert the centroid to an integer for display
                centX = objCentroid(1,1); centY = objCentroid(1,2); %separate X and Y values
                label = ['Blob Number: ' num2str(blobID) ' Centroid: ' ...
                    num2str(centX) ' ' num2str(centY)]; %Create label for output screen
                out2 = insertObjectAnnotation(I, 'rectangle', blobBBox, label, ...
                        'Color', colors(blobID,:)); %Annotate objects with bounding box and label
                if k==1
                    out = out2;
                else
                    out = out + out2; %Add new annotations to old annotations
                end
                lastCentroids(blobID,:) = centroids(blobIndex,:); %fill the array of last 
                                        %centroids so objects are in order of blobID
                k = k+1;
            end
        end    
    end
    
    if nFrame > 1
         if numberOfBlobs==0 %there are no bright objects detected
            vidOut = step(htextins, I); % display "No Objects"
            step(hVideoIn, vidOut); % Output video stream
         else
             step(hVideoIn, out); % Output video stream with annotations
         end
    end
    
    % make a copy of the number of objects in the most recent frame to 
    % compare to future frames
    last = size(blobMeasurements, 1); 
    nFrame = nFrame+1; %Increment the frame number 
    
end

%% Clearing Memory
release(hVideoIn); % Release all memory and buffer used
release(vidDevice);
% clear all;
clc;
clear;