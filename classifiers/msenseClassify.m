function [ labels, confidence ] = msenseClassify(classifier,features)
%% Reed Gurchiek, 2018
%   description
%
%----------------------------------INPUTS----------------------------------
%
%   in:
%       data type, description
%
%---------------------------------OUTPUTS----------------------------------
%
%   out:
%       data type, description
%
%--------------------------------------------------------------------------
%% msenseClassify

% for each feature vector
n = size(features,2);
model = classifier.modelDetails.model;
labels = zeros(1,n);
confidence = zeros(1,n);
for f = 1:n

    % predict
    [~,score] = predict(model.mdl,features(:,f)');
    
    % assign label
    labels(f) = -1;
    if score(2) > model.scoreThreshold
        labels(f) = 1;
    end

    % get column associated with prediction
    column = model.mdl.ClassNames == labels(f);

    % save confidence
    confidence(f) = score(column);
    
end
        
    
end