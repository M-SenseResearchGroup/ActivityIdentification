function [ manipulatorInfo ] = featureManipulator_LassoDB(manipulatorInfo)
%Reed Gurchiek, 2018
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
%% featureManipulator_DaviesBouldin
not ready dont use
% if initialize
if strcmpi(manipulatorInfo.action,'train')
    
    % get mu, sigma, and normalize features
    [features,mu,sigma] = zscore(manipulatorInfo.features,0,2);
    
    % Define logical index array that identifies one of the classes
    labels = manipulatorInfo.labels;
    classes = unique(labels);
    ind = labels == classes(1);
    
    % fit linear model with lasso regularization
    linearModel = fitclinear(features',labels,'Regularization','lasso','Learner','logistic');
    linearModel = lassoglm(features',labels,'NumLambda',2);
    
    % Compute davies-bouldin index for each feature
    dbIndex = zeros(1,size(features,1));
    for f = 1:length(dbIndex)
        featureVector = features(f,:);
        c1 = median(featureVector(ind)); % centroid of class 1
        c2 = median(featureVector(~ind)); % centroid of class 2
        s1 = rms(featureVector(ind)-c1); % average distance from centroid for class 1
        s2 = rms(featureVector(~ind)-c2); % average distance from centroid for class 2
        m12 = norm(c1-c2); % distance between centroids of class 1 and class 2
        dbIndex(f) = (s1+s2)/m12; %adding 1 in denom to prevent Inf.
    end
    
    % Rank davies-bouldin index for features (low to high, low=good separation)
    [dbIndex,dbRank] = sort(dbIndex,'ascend');
    
    % if no keep criteria
    if ~isfield(manipulatorInfo,'keepCriteria')
        
        % get criteria
        keepCriteria = questdlg('Discrete number of features or Davies-Bouldin Index based criteria for keeping features?','Keep Criteria','Number of Features','Davies-Bouldin Index','Davies-Bouldin Index');
        
        % display features, db index and instruct user
        instruction = sprintf('Observe the list of Davies-Bouldin Indices and determine the threshold %s',keepCriteria);
        questdlg(instruction,'Threshold Determination','OK','OK');
        dbDescription = cell(1,length(dbRank));
        for i = 1:length(dbRank)
            if isfield(manipulatorInfo,'originalFeatureNames')
                dbDescription{i} = sprintf('(%d) DB Index = %4.2f, Feature: %s',i,dbIndex(i),manipulatorInfo.originalFeatureNames{dbRank(i)});
            else
                dbDescription{i} = sprintf('(%d) DB Index = %4.2f',i,dbIndex(i));
            end
        end
        i = listdlg('ListString',dbDescription,'PromptString','Select OK or Cancel when ready.','SelectionMode','single','ListSize',[300 200]);
        
        % confirm or alter
        if strcmp(keepCriteria,'Number of Features') 
            threshold = i;
            manipulatorInfo.keepCriteria = 'nFeatures';
        else
            threshold = dbIndex(i); 
            manipulatorInfo.keepCriteria = 'dbIndex';
        end
        
        prompt = sprintf('Confirm or alter %s threshold = %4.2f',keepCriteria,threshold);
        threshold = inputdlg(prompt,'Keep Threshold',[1 50],{num2str(threshold)});

        % save
        manipulatorInfo.keepThreshold = str2double(threshold);
        
        
        
    end
    
    % get keep indices
    if strcmp(manipulatorInfo.keepCriteria,'nFeatures')
        keepIndices = 1:manipulatorInfo.keepThreshold;
    else
        keepIndices = find(dbIndex <= manipulatorInfo.keepThreshold);
    end
    
    % save featureIndices, davies bouldin indices, and mu/sigma
    manipulatorInfo.featureIndices = dbRank(keepIndices);
    manipulatorInfo.dbIndices = dbIndex(keepIndices);
    manipulatorInfo.mu = mu(keepIndices);
    manipulatorInfo.sigma = sigma(keepIndices);
    manipulatorInfo.featureNames = manipulatorInfo.originalFeatureNames(dbRank(keepIndices));
    
    % no need to save features or labels
    manipulatorInfo = rmfield(manipulatorInfo,{'features','labels'});
    
% if manipulate
elseif strcmp(manipulatorInfo.action,'manipulate')
    
    % get mu, sigma, feature indices
    mu = manipulatorInfo.mu;
    sigma = manipulatorInfo.sigma;
    i = manipulatorInfo.featureIndices;
    
    % manipulate
    manipulatorInfo.features = diag(1./sigma)*(manipulatorInfo.features(i,:) - mu);
    
end