function [ manipulatorInfo ] = featureManipulator_PCA(manipulatorInfo)
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
%% featureManipulator_PCA

% if train
if strcmpi(manipulatorInfo.action,'train')
    
    % get mu, sigma, and normalize features
    [allFeatures,mu,sigma] = zscore(manipulatorInfo.features,0,2);
    manipulatorInfo.mu = mu;
    manipulatorInfo.sigma = sigma;
    [pc,~,~,~,pcvar] = pca(allFeatures');
    
    % if no keep criteria
    if ~isfield(manipulatorInfo,'keepCriteria')
        
        % get criteria
        keepCriteria = questdlg('Discrete number of features or component variance based criteria for keeping features?','Keep Criteria','Number of Features','Minimum PC Variance','Minimum PC Variance');
        
        % display pcs and variance and instruct user
        instruction = sprintf('Observe the list of principal components and determine the threshold %s',keepCriteria);
        questdlg(instruction,'Threshold Determination','OK','OK');
        pcDescription = cell(1,length(pcvar));
        for ipc = 1:length(pcvar)
            pcDescription{ipc} = sprintf('(%d) Variance = %4.2f',ipc,pcvar(ipc));
        end
        ipc = listdlg('ListString',pcDescription,'PromptString','Select OK or Cancel when ready.');
        
        % confirm or alter
        if strcmp(keepCriteria,'Number of Features') 
            threshold = ipc;
            manipulatorInfo.keepCriteria = 'nFeatures';
        else
            threshold = pcvar(ipc); 
            manipulatorInfo.keepCriteria = 'componentVariance';
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
        keepIndices = find(pcvar >= manipulatorInfo.keepThreshold);
    end
    
    % save pcs and variance
    manipulatorInfo.principalComponents = pc(:,keepIndices);
    manipulatorInfo.pcVariance = pcvar(keepIndices);
    
    % no need to save features or labels
    manipulatorInfo = rmfield(manipulatorInfo,{'features','labels'});
    
% if manipulate
elseif strcmp(manipulatorInfo.action,'manipulate')
    
    % get mu, sigma, pcs
    mu = manipulatorInfo.mu;
    sigma = manipulatorInfo.sigma;
    pc = manipulatorInfo.principalComponents;
    
    % manipulate
    manipulatorInfo.features = pc'*diag(1./sigma)*(manipulatorInfo.features - mu);
    
end
    
    
            



end