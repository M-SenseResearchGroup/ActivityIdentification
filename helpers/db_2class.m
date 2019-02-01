function db_rank = db_2class(dataTrain,labelsTrain)
    
    % Define logical index array that identifies one of the classes
    classes = unique(labelsTrain);
    ind = labelsTrain == classes(1);
    
    % Compute davies-bouldin index for each feature
    db_feat = zeros(1,size(dataTrain,2));
    for feat_ind = 1:size(dataTrain,2)
        feat = dataTrain(:,feat_ind);
        c1 = median(feat(ind)); % centroid of class 1
        c2 = median(feat(~ind)); % centroid of class 2
        s1 = rms(feat(ind)-c1); % average distance from centroid for class 1
        s2 = rms(feat(~ind)-c2); % average distance from centroid for class 2
        m12 = norm(c1-c2); % distance between centroids of class 1 and class 2
        db_feat(feat_ind) = (s1+s2)/(m12+0.1); %adding 1 in denom to prevent Inf.
    end
    
    % Rank davies-bouldin index for features (low to high, low=good separation)
    [~,db_rank] = sort(db_feat,'ascend');
end