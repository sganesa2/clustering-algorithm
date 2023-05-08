function [c] = Part1_b_kplusplus_init(k,n,m,XData)

%Get a random index between 1 and n.
i = randi(n);

c = zeros(k,m);
%Assign the ith row of XData to our initial cluster c1
c(1,:) = XData(i,:);

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
end

