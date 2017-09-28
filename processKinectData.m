function [normals,center,centroids,meanDistanceToGreen,finalPtCloud] = processKinectData(txtname,bPlot)
% PROCESSKINECTDATA Takes a '*.txt' file of the form (points) X (6) [where 
% the first three columns represent xyz coordinates and the last three
% represent colors (in RGB), which is the usual output of a Kinect sensor
% when used with libfreenect2] and returns: 
% 1. NORMALS - Normal vectors of each of the three planes (again, RGB 
% planes)
% 2. CENTER - Center of the cube
% 3. CENTROIDS - Centroids of each of the planes
% 4. MEANDISTANCETOGREEN - The mean distance of all the points in the green 
% plane to the centroid of that same plane
%
% Tip: NORMALS and CENTROIDS are 3x3 matrices where each row is a vector
% corresponding to the proper characteristic.

% Run the following code from the terminal:
% for idx = 1:4
% generateImage(sprintf('mk2cap%1.gcubo.txt',idx),true)
% end

% If no argument for plotting is given -> assume user WANTS to plot.
fig = length(findobj('type','figure')) + 1;
if nargin < 2, bPlot = true; end

%% Load data
data = load(txtname);
% Convert to uint8 for colors
color = uint8(data(:,4:6));
locs = data(:,1:3);
% Generate singleton dimension
color = permute(color,[1 3 2]);
locs = permute(locs,[1 3 2]);
% Reshape for actual image
r = 512; c = 424;
img = reshape(color,r,c,3);
locs = reshape(locs,r,c,3);
img = rot90(img,-1);
locs = rot90(locs,-1);
% Show image
if bPlot
    figure(fig)
    subplot(2,3,1), imshow(img), title 'Original image'
end
% Crop the image (rectangle = [xmin ymin width height])
xmin = round(c/3);
ymin = round(r/5);
rectangle = [xmin,ymin,2*xmin,2*ymin];
cimg = imcrop(img,rectangle);
croppedlocs = imcrop(locs,rectangle);
dists = sum(croppedlocs.^2,3);

%% Apply distance threshold
distanceThr = 2;
where = find(dists > distanceThr | isnan(dists));
[i,j] = ind2sub(size(dists),where);
for idx = 1:length(i)
    cimg(i(idx),j(idx),:) = 0;
end

%% Apply color threshold
colorThr = 80;% Blacks!
where = find((cimg(:,:,1) < colorThr) & (cimg(:,:,2) < colorThr) & ...
    (cimg(:,:,3) < colorThr));
[i,j] = ind2sub(size(dists),where);
for idx = 1:length(i)
    cimg(i(idx),j(idx),:) = 0;
end

if bPlot
    figure(fig)
    subplot(2,3,2), imshow(cimg), title 'Image after processing'
end

%% Create point cloud
% Find location where all three channels are NOT set to 0
where = find(sum(cimg,3) ~= 0);
[i,j] = ind2sub(size(dists),where);
xyz = zeros(length(i),3);
cols = zeros(length(i),3);
for idx = 1:length(i)
    xyz(idx,:) = croppedlocs(i(idx),j(idx),:);
    cols(idx,:) = cimg(i(idx),j(idx),:);
end
totalCloud = pointCloud(xyz,'Color',uint8(cols));
if bPlot
    figure(fig)
    subplot(3,3,3), pcshow(totalCloud)
    title 'Point cloud of relevant points'
end

%% Perform segmentation with k-means
lbls = kmeans(cols,3,'Start',255*eye(3));

ptCloud = cell(1,3);
planesPtCloud = cell(1,3);
planes = cell(1,3);
maxDistance = .001;

finalData = [];
finalColors = [];
for idx = 1:3
    where = lbls == idx;
    % Define color purely:
    color = [0 0 0];
    color(idx) = 255;
    ptCloud{idx} = pointCloud(xyz(where,:),...
                    'Color',uint8(repmat(color,nnz(where),1)));
    [planes{idx},inlierIdxs,~] = pcfitplane(ptCloud{idx},maxDistance);
    planesPtCloud{idx} = select(ptCloud{idx},inlierIdxs);
    finalData = [finalData; planesPtCloud{idx}.Location];
    finalColors = [finalColors; planesPtCloud{idx}.Color];
end

finalPtCloud = pointCloud(finalData,'Color',uint8(finalColors));

%% Plot the planes
if bPlot
    figure(fig)
    for idx = 1:3
        subplot(2,3,6+idx), pcshow(planesPtCloud{idx})
        title(sprintf('Plane number %g',idx))
    end
end

%% Getting the crossing point of three planes
ptCorner = -[planes{1}.Normal;planes{2}.Normal;planes{3}.Normal]\...
    [planes{1}.Parameters(4);planes{2}.Parameters(4);planes{3}.Parameters(4)];

% Plot the corner
if bPlot
    figure(fig), subplot(2,3,3)
    hold on, scatter3(ptCorner(1),ptCorner(2),ptCorner(3),100)
end
% msgbox(sprintf(['Corner was found at coordinates: %5.3f,%5.3f,%5.3f\n',...
%    'Look at subplot (2,3,3)'],ptCorner(1),ptCorner(2),ptCorner(3)))

%% Plotting normal vectors to each surface
centers = zeros(3);
normals = zeros(3);
for idx = 1:3
    centers(idx,:) = mean(ptCloud{idx}.Location);
    normal = planes{idx}.Normal;
    if normal(2) > 0, normal = -normal; end
    normals(idx,:) = normal;
end

if bPlot
    figure(fig)
    subplot(2,3,3)
    hold on, quiver3(centers(:,1),centers(:,2),centers(:,3),...
        normals(:,1),normals(:,2),normals(:,3))
end

%% Get the center of the box
boxSize = .2;% 20 cm!
center = ptCorner' - boxSize/2*(normals(1,:) + normals(2,:) + normals(3,:));
if bPlot
    figure(fig)
    subplot(2,3,3)
    hold on, scatter3(center(1),center(2),center(3),100)
end

%% Get centroid of each of the planes
centroids = zeros(3);
for idx = 1:3
    centroids(idx,:)    = mean(planesPtCloud{idx}.Location);
end

%% Get mean distance from green point cloud to the center of the plane
% Remember that, from kmeans, R = 1, G = 2 and B = 3
meanDistanceToGreen     = mean(pdist2(planesPtCloud{2}.Location,...
    centroids(2,:)));

