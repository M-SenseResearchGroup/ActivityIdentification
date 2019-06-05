function [ manipulatorInfo ] = featureManipulator_DB2PCA(manipulatorInfo)
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
%% featureManipulator_DB2PCA

% if initialize
if strcmpi(manipulatorInfo.action,'train')
    
    % get mu, sigma, and normalize features
    [features,mu,sigma] = zscore(manipulatorInfo.features,0,2);
    
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
    
    % if no DB keep criteria
    if ~isfield(manipulatorInfo,'dbKeepCriteria')
        
        % get criteria
        dbKeepCriteria = questdlg('Discrete number of features or Davies-Bouldin Index based criteria for keeping features?','Keep Criteria','Number of Features','Davies-Bouldin Index','Davies-Bouldin Index');
        
        % display features, db index and instruct user
        instruction = sprintf('Observe the list of Davies-Bouldin Indices and determine the threshold %s',dbKeepCriteria);
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
        if strcmp(dbKeepCriteria,'Number of Features') 
            dbThreshold = i;
            manipulatorInfo.dbKeepCriteria = 'nFeatures';
        else
            dbThreshold = dbIndex(i); 
            manipulatorInfo.dbKeepCriteria = 'dbIndex';
        end
        
        prompt = sprintf('Confirm or alter %s threshold = %4.2f',dbKeepCriteria,dbThreshold);
        dbThreshold = inputdlg(prompt,'Keep Threshold',[1 50],{num2str(dbThreshold)});

        % save
        manipulatorInfo.dbKeepThreshold = str2double(dbThreshold);
        
    end
    
    % get db keep indices
    if strcmp(manipulatorInfo.dbKeepCriteria,'nFeatures')
        dbKeepIndices = 1:manipulatorInfo.dbKeepThreshold;
    else
        dbKeepIndices = find(dbIndex <= manipulatorInfo.dbKeepThreshold);
    end
    
    % save featureIndices, davies bouldin indices, and mu/sigma
    manipulatorInfo.dbFeatureIndices = dbRank(dbKeepIndices);
    manipulatorInfo.dbIndices = dbIndex(dbKeepIndices);
    manipulatorInfo.mu = mu(dbKeepIndices);
    manipulatorInfo.sigma = sigma(dbKeepIndices);
    manipulatorInfo.dbFeatureNames = manipulatorInfo.originalFeatureNames(dbRank(dbKeepIndices));
    
    % get pca of db selected features
    [pc,~,pcvar] = pca(features(dbKeepIndices,:)');
    pcvar = 100*pcvar/sum(pcvar);
    
    % if no PCA keep criteria
    if ~isfield(manipulatorInfo,'pcaKeepCriteria')
        
        % get criteria
        pcaKeepCriteria = questdlg('Discrete number of features or component variance based criteria for keeping features?','Keep Criteria','Number of Features','Minimum PC Variance','Minimum PC Variance');
        
        % display pcs and variance and instruct user
        instruction = sprintf('Observe the list of principal components and determine the threshold %s',pcaKeepCriteria);
        questdlg(instruction,'Threshold Determination','OK','OK');
        pcDescription = cell(1,length(pcvar));
        for ipc = 1:length(pcvar)
            pcDescription{ipc} = sprintf('(%d) Variance = %4.2f',ipc,pcvar(ipc));
        end
        ipc = listdlg('ListString',pcDescription,'PromptString','Select OK or Cancel when ready.');
        
        % confirm or alter
        if strcmp(pcaKeepCriteria,'Number of Features') 
            pcaThreshold = ipc;
            manipulatorInfo.pcaKeepCriteria = 'nFeatures';
        else
            pcaThreshold = pcvar(ipc); 
            manipulatorInfo.pcaKeepCriteria = 'componentVariance';
        end
        prompt = sprintf('Confirm or alter %s threshold = %4.2f',pcaKeepCriteria,pcaThreshold);
        pcaThreshold = inputdlg(prompt,'Keep Threshold',[1 50],{num2str(pcaThreshold)});
        
        % save
        manipulatorInfo.pcaKeepThreshold = str2double(pcaThreshold);
        
    end
    
    % get pca keep indices
    if strcmp(manipulatorInfo.pcaKeepCriteria,'nFeatures')
        pcaKeepIndices = 1:manipulatorInfo.pcaKeepThreshold;
    else
        pcaKeepIndices = find(pcvar >= manipulatorInfo.pcaKeepThreshold);
    end
    
    % save pcs and variance
    manipulatorInfo.principalComponents = pc(:,pcaKeepIndices);
    manipulatorInfo.pcVariance = pcvar(pcaKeepIndices);
    
    % no need to save features or labels
    manipulatorInfo = rmfield(manipulatorInfo,{'features','labels'});
    
% if manipulate
elseif strcmp(manipulatorInfo.action,'manipulate')
    
    % get mu, sigma, db feature indices
    mu = manipulatorInfo.mu;
    sigma = manipulatorInfo.sigma;
    i = manipulatorInfo.dbFeatureIndices;
    
    % project db selected features onto pcs
    manipulatorInfo.features = manipulatorInfo.principalComponents'*diag(1./sigma)*(manipulatorInfo.features(i,:) - mu);
    
end