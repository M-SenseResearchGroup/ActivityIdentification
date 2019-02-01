function [ labels,confidence ] = msenseClassify(classifier,features)
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
modelLabels = zeros(length(model),n);
modelConfidence = zeros(length(model),n);
for f = 1:n
    
    % for each model
    for m = 1:length(model)

        % predict
        [modelLabels(m,f),score] = predict(model(m).mdl,features(:,f)');
        
        % get column associated with prediction
        column = model(m).mdl.ClassNames == modelLabels(m,f);
        
        % save confidence
        modelConfidence(m,f) = score(column);
        
    end
    
end

% if single classifier
if strcmp(classifier.modelDetails.type,'singleClassifier')

    % then done
    labels = modelLabels;
    confidence = modelConfidence;

% otherwise democratic
else
    % if weighted combination
    if strcmp(classifier.modelDetails.combination,'weightedConfidence')

        % get weights
        weights = zeros(1,length(model));
        for m = 1:length(model)
            weights(m) = model(m).crossValidationAccuracy;
        end
        
    % otherwise unweighted
    else
        
        % weights all same to result in normal average
        weights = ones(1,length(model));
        
    end
    
    % for each label
    labels = zeros(1,n);
    confidence = zeros(1,n);
    for l = 1:n
        
        % combine
        confidence(l) = sum(weights.*modelLabels(:,l)'.*modelConfidence(:,l)')/sum(weights);
        labels(l) = sign(confidence(l));
        confidence(l) = abs(confidence(l));
        
    end
    
end
        
    
end