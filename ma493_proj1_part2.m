% Clear the workspace and close all figure windows
%clear all
%close all

% Load the provided data vectors into the variable XData
load Q1data.mat

% Get the dimensions(shape) of the set fo data vectors using 'size' command
% Set the number of data vectors (n) and the dimension of the data space (m)
shape = size(XData);
n = shape(1);
m = shape(2);

% Get elbow plots for k values from 1,2...,8
% Initialize a data structure to hold overall coherence for each k value. 
overall_coherence=zeros(8,1);
for k = 1:8
    % PART I)a) Random initiatialization

    % Create data structures to store the current cluster representative vectors, and the 
    % cluster representative vectors from the previous iteration (cPrev)
    [c] = Part1_a_random_init(k,m);
    cPrev = c; 


    % PART I)b) Kmeans++ initiatialization

    % Create data structures to store the current cluster representative vectors, and the 
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
        overall_coherence(k) = overall_coherence(k) + q(s);
    end

end
% Let k_set denote successive k values from 1 to 8
k_set = linspace(1,8,8);
% Now perform the elbow plot of 'overall coherence' vs 'k values'
plot(k_set,overall_coherence,'y-','LineWidth',2)
xlabel("k-values")
ylabel("Overall Coherence")




