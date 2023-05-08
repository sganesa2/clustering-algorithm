%% MA 493 - Project#1 - Part 3 - Intro Script for visualizing and resizing
%% images of handwritten digits from MNIST data set. Each image is a 20x20 
%% matrix of greyscale values
%% Data downloaded from: http://yann.lecun.com/exdb/mnist/
%%
%% by Dr. Mansoor Haider - NCSU Mathematics
%%

%clear all
close all

% Set the number of images to extract from the data set of test images
NImages = 100;

% Uses the script by S Hegde to extract the images ('testImages') and
% correct labels ('testLabels') from the two files in the active path
% file "readMNIST.m" should be in an active path

[imgs,labels] = readMNIST('testImages','testLabels', NImages, 0);

% Example of how to images 6 through 10 in separate figures

%for i=6:10
 %  figure(i) 
  % imshow(imgs(:,:,i),'InitialMagnification',1000) 
%end


% Example of how to convert the first 100 images into vectors for input into
% the clustering algorithm

m = 20*20;
v = zeros(1,m);
XData = zeros(NImages,m);

for i=1:NImages
    XData(i,:) = reshape(imgs(:,:,i),[1,m]);
end

% CLUSTERING PART WITH ELBOW PLOT

% Get the dimensions(shape) of the set fo data vectors using 'size' command
% Set the number of data vectors (n) and the dimension of the data space (m)
shape = size(XData);
n = shape(1);
m = shape(2);

% Get elbow plots for k values from 1,2...,8
% Initialize a data structure to hold overall coherence for each k value. 
overall_coherence=zeros(8,1);
% Since k=6 is the ideal number of clusters for this image set, we want to
% create a data structure to store the cluster indices of each image in a
% data structure k6_Clusterindices
k6_Clusterindices = zeros(n,1);
for k = 3:10


    % Kmeans++ initiatialization

    % Create data structures to store the current cluster representative vectors, and the 
    % cluster representative vectors from the previous iteration (cPrev)
    % Initialize the first cluster representative vector to be the vector
    % of the 42nd image.
    
    c = zeros(k,m);
    %Assign the ith row of XData to our initial cluster c1
    c(1,:) = XData(42,:);
    
    % use nested for loops to obtain the intital clusters c2, c3,...,ck
    for i = 2:k
        % Create an nx(i-1) matrix to store the euclidean distance between each data
        % point and each previously chosen cluster.
        dist = zeros(n,i-1);
        min_dist = zeros(n,1); % min_dist is used to find the nearest cluster representative for each data point
        for j = 1:n % denotes the index of the data vector
            for l = 1:(i-1) % denotes the index of the previously chosen cluster representatives
                dist(j,l) = norm(XData(j,:)- c(l,:),2); %euclidian distance between data point and cluster representative
            end
            min_dist(j) = min(dist(j,:));
        end
        [val,idx] = max(min_dist); % idx is used to obtain the index of the data vector that is farthest away from the nearest cluster representatives
        c(i,:) = XData(idx,:);
    end
    cPrev = c; 



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


    %The Alternating Minimization Scheme
    doneFlag=0;

    %Keep alternating updates to weight vectors and cluster assignments for 5
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

        % Get the cluster indices when k=6(ideal number of clusters)
        if k==6
            k6_Clusterindices = closestCluster;
        end


       % Terminate the alternating scheme if the weight vectors are unaltered
        % relative to the previous iteration
        if c==cPrev
            doneFlag=1;
        else
            cPrev=c;
        end
            
    end
    % Store the coherence for every cluster in a data structure q
    q = zeros(k,1);
    % Get the coherence value for each cluster
    for i=1:k 
        % Find the indices for all data vectors currently in cluster i
        ClusterIndices = find(IndexSet==i);
    
        % Find the number of data vectors currently in cluster i
        NumVecsInCluster = size(ClusterIndices,1);
    
        for j=1:NumVecsInCluster
            q(i) = q(i) + (norm(XData(ClusterIndices(j),:)-c(i,:),2)).^2;
        end
    end
    % Calculate overall coherance by performing summation of the cluster
    % coherence values
    for s = 1:k
        overall_coherence(k-2) = overall_coherence(k-2) + q(s);
    end

end
% Let k_set denote successive k values from 1 to 8
k_set = linspace(3,10,8);
% Now perform the elbow plot of 'overall coherence' vs 'k values'
figure(50)
plot(k_set,overall_coherence,'k-','LineWidth',2)

% STEP 4 of Part 3 ( Evaluating the success score for k=6)
% initialize a data structure to get the most repeated element in each
% cluster
mostrepeated = zeros(6,1);
%initialize a data structure to get the frequency of the most repeated
%elements in each cluster
freq_mostrepeated = zeros(6,1);
for i = 1:6
    ClusterIndices = find(k6_Clusterindices==i);
    NumVecsInCluster = size(ClusterIndices,1);
    % Initialize a data structure to store all the digit labels present in a given cluster i
    cluster_digits = zeros(NumVecsInCluster,1);
    for j = 1:NumVecsInCluster
        cluster_digits(j) = labels(ClusterIndices(j));
    end
    %use mode function to obtain the most repeated element, and its frequency in cluster_digits
    [A,F] = mode(cluster_digits);
    mostrepeated(i) = A;
    freq_mostrepeated(i) = F;

    % The following 'if' statements have been written after observing the
    % most frequent element in each cluster an their respective
    % frequencies.
    
    % CLUSTER 1 with the representative digit label as '1' gives the worst
    % clustering
    if i==1
        best_cluster_1 = find(k6_Clusterindices==i);
        best_vector_count_1 = size(best_cluster_1,1);
         %Define the number of rows and columns in the subplot grid
        num_rows = 13; % 13 rows
        num_cols = 4; % 4 columns
        % Loop through the image files and display each one in a different subplot
        for h = 1:best_vector_count_1       
            % Compute the row and column index of the subplot
            row = ceil(h/num_cols);
            col = mod(h-1, num_cols) + 1;
            
            % Create the subplot and display the image
            figure(i)
            subplot(num_rows, num_cols, (row-1)*num_cols + col);
            currImg = reshape(XData(best_cluster_1(h),:),[20,20]);
            imshow(currImg,'InitialMagnification',1000)
            title(sprintf('Image %d', h));
        end
    end

    % CLUSTER 5 with the representative digit label as '4' gives a good
    % clustering
    if i==5
        best_cluster_5 = find(k6_Clusterindices==i);
        best_vector_count_5 = size(best_cluster_5,1);  
       %Define the number of rows and columns in the subplot grid
        num_rows = 2; % 2 rows
        num_cols = 9; % 9 columns
        %Loop through the image files and display each one in a different subplot
        for h = 1:best_vector_count_5       
            %Compute the row and column index of the subplot
            row = ceil(h/num_cols);
            col = mod(h-1, num_cols) + 1;
            
            % Create the subplot and display the image
            figure(i)
            subplot(num_rows, num_cols, (row-1)*num_cols + col);
            currImg = reshape(XData(best_cluster_5(h),:),[20,20]);
            imshow(currImg,'InitialMagnification',1000)
            title(sprintf('Image %d', h));
        end
    end

    % CLUSTER 2 with the representative digit label as '0' also gives a good
    % clustering
    if i==2
        best_cluster_2 = find(k6_Clusterindices==i);
        best_vector_count_2 = size(best_cluster_2,1);
        %Define the number of rows and columns in the subplot grid
        num_rows = 2; % 2 rows
        num_cols = 2; % 2 columns
        %Loop through the image files and display each one in a different subplot
        for h = 1:best_vector_count_2       
            %Compute the row and column index of the subplot
            row = ceil(h/num_cols);
            col = mod(h-1, num_cols) + 1;
            
            % Create the subplot and display the image
            figure(i)
            subplot(num_rows, num_cols, (row-1)*num_cols + col);
            currImg = reshape(XData(best_cluster_2(h),:),[20,20]);
            imshow(currImg,'InitialMagnification',1000)
            title(sprintf('Image %d', h));
        end
    end

end
% Initialize a variable to store the Success score S for the clustering.
% Note: We used k=6 because this was the best value of k given to us by our
% elbow plot
% We calculate success score by adding the frequencies of the most
% repeated element in each cluster for all the clusters. This is because,
% if an digit label in the cluster is not identical to the Cluster
% representative a.k.a 'The most repeated element', then the clustering
% isn't perfect.
% The perfect elements in each cluster woukld only be the ones that are
% identical to the cluster representative. For that reason, we calculate our
% success score value using the method outlined above.

S=0;
for i =1:6
    S = S + freq_mostrepeated(i);
end
mostrepeated
freq_mostrepeated



% Last, lets plot the correct labels for each image and the histogram of
% labels in images extracted from the data set
figure(101)
plot(labels)
figure(102)
hist(labels)
