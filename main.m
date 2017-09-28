%% Compute normals, centers, centroids and mean-distances
nKinects = 4;
bPlot = true;
% Preallocate memory
normals = zeros(3,3,nKinects);
centers = zeros(nKinects,3);
centroids = zeros(3,3,nKinects);
meanDistances = zeros(1,nKinects);
ptClouds = cell(1,nKinects);

% Compute things
for idx = 1:4
    [normals(:,:,idx),centers(idx,:),centroids(:,:,idx),...
        meanDistances(idx),ptClouds{idx}] = ...
        processKinectData(sprintf('KinectData3/mk2cap%1g.txt',idx),true);
end

%% Plot original data
fig = length(findobj('type','figure')) + 1;
if bPlot 
    figure(fig)
    for idx = 1:4
        subplot(2,4,idx), pcshow(ptClouds{idx})
        hold on, scatter3(centers(idx,1),centers(idx,2),centers(idx,3),100)
        hold on, quiver3(repmat(centers(idx,1),3,1),...
                         repmat(centers(idx,2),3,1),...
                         repmat(centers(idx,3),3,1),...
                         normals(:,1,idx),...
                         normals(:,2,idx),...
                         normals(:,3,idx))
        title 'Original clouds and vectors'
    end
end

%% 'Decide' locations
% Check which are the shortest distances
[~,b] = sort(meanDistances);
% Extract the first two
indices = b(1:2);
% Iterate over all the point clouds
for idx = 1:4
    % Decide if RED normal vector should be put upside down, that is, check
    % if RED plane is on the top or on the bottom
    zCoords = centroids(:,3,idx);
    if zCoords(1) < zCoords(2) && zCoords(1) < zCoords(3)% Flip vector
        normals(1,:,idx) = -normals(1,:,idx);
    end
    
    % Decide if GREEN normal vector should be put upside down, that is,
    % check if it doesn't have a hole in it by checking if the mean 
    % distances are the smallest ones.
    % If it doesn't have a hole, the mean distances are smaller!
    if ismember(idx,indices)% This means that it doesn't have a hole
        normals(2,:,idx) = -normals(2,:,idx);
    end
    
    % Decide if BLUE normal vector should be put upside down, that is,
    % check if:
    % 1. GREEN plane is to the left of BLUE plane && that GREEN doesn't
    % have a hole
    % OR
    % 2. GREEN plane is to the right of BLUE plane && that GREEN does have
    % a hole
    xCoords = centroids(:,1,idx);
    if xCoords(2) < xCoords(3)% GREEN plane is to the left of BLUE plane
        if ismember(idx,indices)% This means that it doesn't have a hole
            normals(3,:,idx) = -normals(3,:,idx);
        end
    else% GREEN plane is to right of BLUE plane
        if ~ismember(idx,indices)% This means that it does have a hole
            normals(3,:,idx) = -normals(3,:,idx);
        end
    end
end

%% Plot to see if everything worked out fine
if bPlot
    figure(fig)
    for idx = 1:4
        subplot(2,4,4 + idx), pcshow(ptClouds{idx})
        hold on, scatter3(centers(idx,1),centers(idx,2),centers(idx,3),100)
        hold on, quiver3(repmat(centers(idx,1),3,1),...
            repmat(centers(idx,2),3,1),repmat(centers(idx,3),3,1),...
            normals(:,1,idx),normals(:,2,idx),normals(:,3,idx))
        title 'Fixed clouds and vectors'
    end
end

%% Construct transformation matrices and apply transformation
% These matrices convert from the cube coordinates to the coordinates of
% the respective Kinect
transformMatrices = zeros(4,4,4);

% Inverse transformation matrices convert from the coordinates of the
% corresponding Kinect to the coordinates of the cube
inverseTransformMatrices = zeros(4,4,4);

% % Define transformed data (xyz coordinates) and color information
transformed = [];

% Define new ptClouds that are relatively close to each other
newPtClouds = cell(1,nKinects);

for idx = 1:4
    % Transformation matrix
    temporalMat = [normals(:,:,idx)',centers(idx,:)'; 0 0 0 1];
    transformMatrices(:,:,idx) = temporalMat;
    % Inverse transformation matrix
    invTemporalMat = inv(temporalMat);
    inverseTransformMatrices(:,:,idx) = invTemporalMat;
    % Extract xyz data
    data = ptClouds{idx}.Location;
    % Number of points
    [ndata,~] = size(data);
    % Save in matrix with appropriate dimensions
    dataMat = [data';ones(1,ndata)];
    % Apply transformation by using '\' operator -> more accurate
    temporalTransformed = (temporalMat\dataMat)';
    transformed = [transformed;...
        temporalTransformed(:,1:3),double(ptClouds{idx}.Color)];
    newPtClouds{idx} = pointCloud(temporalTransformed(:,1:3),...
        'Color',ptClouds{idx}.Color);
end

%% Create final point cloud with transformed data
finalPtCloud = pointCloud(transformed(:,1:3),...
    'Color',uint8(transformed(:,4:6)));
figure(fig+1)
pcshow(finalPtCloud)
% fig = length(findobj('type','figure')) + 1;
% figure(fig)
% for idx = 1:4
%     subplot(2,2,idx), pcshow(newPtClouds{idx})
% end
% 
% %% Apply ICP 
% % BUT only to planes that have the same colors and then use that
% % transformation for the whole cube
% 
% % Get new centroids (in these new coordinates): each colum is a vector
% centroids = zeros(3,3,nKinects);
% 
% % Iterate over the Kinects (same as newPtClouds)
% for idx = 1:4
%     colorData = newPtClouds{idx}.Color;
%     xyzData = newPtClouds{idx}.Location;
%     % Iterate over the colors
%     for jdx = 1:3
%         where = colorData(:,jdx) == 255;
%         % Compute centroid
%         thisCentroid = mean(xyzData(where,:));
%         centroids(:,jdx,idx) = thisCentroid';
%     end
% end
% 
% % Find the color of the plane that the point cloud shares with the 1st
% % point cloud
% final = newPtClouds{1};
% firstPlaneCentroids = centroids(:,:,1);
% colorsInChar = 'RGB';
% figure
% for idx = 2:4
%     thisCentroid = centroids(:,:,idx);
%     distances = sum((firstPlaneCentroids - thisCentroid).^2);
%     [~,which] = min(distances);
%     % 'which' can have three possible values (1, 2 and 3) that correspond
%     % in RGB to the plane they share with the 'basic' point cloud
%     
%     % Extract the points that have that color
%     basic = find(newPtClouds{1}.Color(:,which) == 255);
%     temporalBasic = select(newPtClouds{1},basic);
%     
%     where = find(newPtClouds{idx}.Color(:,which) == 255);
%     temporal = select(newPtClouds{idx},where);
%     
%     % Compute transformation with ICP
%     tform = pcregrigid(temporal,temporalBasic,...
%         'Metric','pointToPlane','Extrapolate',true);
%     
%     % Apply transformation
%     ptCloudAligned = pctransform(newPtClouds{idx},tform);
%     
%     % Merge the point clouds
%     mergeSize = 0.001;
%     final = pcmerge(final, ptCloudAligned, mergeSize);
%     endgg
% 
% 
% 
