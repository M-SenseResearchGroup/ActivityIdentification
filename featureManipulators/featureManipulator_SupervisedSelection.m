function [ manipulatorInfo ] = featureManipulator_SupervisedSelection(manipulatorInfo)
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
%% featureManipulator_SupervisedSelection

% if initialize
if strcmpi(manipulatorInfo.action,'train')
    
    % get mu, sigma, and normalize features
    [features,mu,sigma] = zscore(manipulatorInfo.features,0,2);
    
    % if no keep criteria
    if ~isfield(manipulatorInfo,'keepCriteria')
    
        % get first principal component
        pc = pca(features');
        pc = pc(:,1);

        % Define logical index array that identifies one of the classes
        labels = manipulatorInfo.labels;
        classes = unique(labels);
        ind = labels == classes(1);

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

        % sort feature names and pc weights
        featureNames = manipulatorInfo.originalFeatureNames(dbRank);
        pc = pc(dbRank);
        
        % display features and info and instruct user
        keepCriteria = 'supervisedSelection';
        questdlg('Observe the list of features and select the ones to keep','Feature Selection','OK','OK');
        description = cell(1,length(dbRank));
        for i = 1:length(featureNames)
            description{i} = sprintf('(%d) Feature: %s, DB = %4.2f, PC weight = %10.6f',i,featureNames{i},dbIndex(i),pc(i));
        end
        i = listdlg('ListString',description,'PromptString','Select OK or Cancel when ready.','SelectionMode','multiple','ListSize',[300 300]);

        % save
        manipulatorInfo.keepCriteria = keepCriteria; % this essentially confirms whether or not the supervised selection has been done
        manipulatorInfo.featureIndices = dbRank(i);
        manipulatorInfo.dbIndices = dbIndex(i);
        manipulatorInfo.featureNames = featureNames(i);
        
    end
    
    % save mu and sigma
    manipulatorInfo.mu = mu(manipulatorInfo.featureIndices);
    manipulatorInfo.sigma = sigma(manipulatorInfo.featureIndices);
    
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