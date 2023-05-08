%% 

% MA 493 PROJECT 1
% NAME: Srinivas Ganesan
% DATE: 03/07/2023


% Clear the workspace and close all figure windows
clear all
close all

% Load the provided data vectors into the variable XData
load Q1data.mat

% Get the dimensions(shape) of the set fo data vectors using 'size' command
% Set the number of data vectors (n) and the dimension of the data space (m)
shape = size(XData);
n = shape(1);
m = shape(2);

% Set the number of clusters (k)
k=5;

% PART I)a) Random initiatialization

% Create data structures to store the current cluster representative vectors, and the 
% cluster representative vectors from the previous iteration (cPrev)
[c] = Part1_a_random_init(k,m);
cPrev = c; 


% PART I)b) Kmeans++ initiatialization

% Create a structures to store the current cluster representative vectors, and the 
% cluster representative vectors from the previous iteration (cPrev)
%[c] = Part1_b_kplusplus_init(k,n,m,XData);
%cPrev = c; 


% Create a data structure to store intial closest cluster representative vectors for each data
% point
IndexSet=zeros(n,1);
% Assign each data vector to its closest cluster vector
for d=1:n
    % Set the minimum distance tracker to be a very large number
    sqDistMin=1e16;

    xd = XData(d,:);
    % Find the closest cluster representative vector to the current data
    % vector
    for i=1:k
        sqDist = norm(c(i,:)-xd,2);

        % If the distance is less than the current min, assign the
        % current data vector to this cluster
        if sqDist<sqDistMin
            IndexSet(d)=i;
            sqDistMin=sqDist;
        end

    end
end

% Plot the data
scatter(XData(:,1),XData(:,2),64,IndexSet,'filled');
hold on


%The Alternating Minimization Scheme
doneFlag=0;

%Keep alternating updates to weight vectors and cluster assignments for 10
%iterations

while (~doneFlag)
   
   % Update the weight vectors in each cluster via the centroid formula
    for i=1:k 

        % Find the indices for all data vectors currently in cluster i
        ClusterIndices = find(IndexSet==i);

        % Find the number of data vectors currently in cluster i
        NumVecsInCluster = size(ClusterIndices,1);

        % Create a data structure to store weight vector for the current
        % cluster
        c(i,:)=0; 

        % Update cluster vector using the centroid formula
        for j=1:NumVecsInCluster
            for l=1:m
                c(i,l) = c(i,l) + XData(ClusterIndices(j,1),l)/NumVecsInCluster;
            end
        end

    end

    % Plot the updated weight vectors for each cluster
    scatter(XData(:,1),XData(:,2),64,IndexSet,'filled')
    hold on
    scatter(c(:,1),c(:,2),200,linspace(1,k,k),'filled')
    pause
    
    % Now reassign all data vectors to the closest weight vector (cluster)

    % Create a data structure to store closest weight vector for each data
    % point
    closestCluster=zeros(n,1);

    % Reassign each data vector to the new, closest cluster
    for d=1:n

        % Store the coordinates of the current data vector
        xD = XData(d,:);

        % Set the minimum distance tracker to be a very large number
        sqDistMin=1e16;

        % Find the closest weight vector (cluster) to the current data
        % vector
        for i=1:k
            sqDist = norm(c(i,:)-xD,2);

            % If the distance is less than the current min, assign the
            % current data vector to this cluster
            if sqDist<sqDistMin
                closestCluster(d)=i;
                sqDistMin=sqDist;
            end

        end
    end

    % Update the assignments of the data vectors to their new clusters
    IndexSet = closestCluster;
    
    % Plot the data and the updated weight vectors
    hold off
    scatter(XData(:,1),XData(:,2),64,IndexSet,'filled')
    hold on
    scatter(c(:,1),c(:,2),200,linspace(1,k,k),'filled')
    hold off
    pause
    
    % Terminate the alternating scheme if the weight vectors are unaltered
    % relative to the previous iteration
    if c==cPrev
        doneFlag=1;
    else
        cPrev=c;
    end
end