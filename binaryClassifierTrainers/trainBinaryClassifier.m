function [ classifier ] = trainBinaryClassifier(trainer)
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
%% trainBinaryClassifier


% if no trainer
if nargin == 0
        
    if isfield(trainer,'msensePath')
        msense = trainer.msensePath;
    else
        % get msenseresearchgroup path
        ok = questdlg('Select the msenseresearchgroup folder','MSENSE Path','OK',{'OK'});
        if isempty(ok); error('Initialization terminated.'); end
        msense = uigetdir;
    end

    % get possible trainer structs
    trainDir = dir(fullfile(msense,'Project Data','ActivityIdentification','ClassifierTrainers'));
    trainerNames = cell(0);
    i = 1;
    while i <= numel(trainDir)

        % delete if hidden
        if trainDir(i).name(1) == '.'; trainDir(i) = [];

        % if not .mat extension
        elseif ~strcmpi(trainDir(i).name(end-3:end),'.mat'); trainDir(i) = [];

        % otherwise save
        else
            trainerNames{i} = trainDir(i).name;
            i = i + 1;
        end

    end
    
    % if no trainer
    if isempty(trainerNames)
        error('No trainer structs exist in ActivityIdentification/ClassifierTrainers.  Use initializeBinaryClassifierTrainer first.');
    end

    % select trainer
    itrainerNames = listdlg('ListString',trainerNames,'PromptString','Select trainer:','SelectionMode','single','ListSize',[300,100]);
    trainDir = trainDir(itrainerNames);
    trainerNames = trainerNames(itrainerNames);

    % load trainer
    trainer = load(fullfile(msense,'Project Data','ActivityIdentification','ClassifierTrainers',trainerNames{1}),'trainer');
    trainer = trainer.trainer;
    
end

% get features and labels
classNames = trainer.featureSet.classNames;
features = horzcat(trainer.featureSet.class.(classNames{1}).features,trainer.featureSet.class.(classNames{2}).features);
labels = horzcat(repmat(trainer.featureSet.class.(classNames{1}).label,[1 trainer.featureSet.class.(classNames{1}).nObservations]),...
                 repmat(trainer.featureSet.class.(classNames{2}).label,[1 trainer.featureSet.class.(classNames{2}).nObservations]));
             
% manipulate features
manipulator = str2func(trainer.featureManipulator.name);
manipulatorInfo = trainer.featureManipulator.manipulatorInfo;
manipulatorInfo.action = 'manipulate';
manipulatorInfo.features = features;
manipulatorInfo = manipulator(manipulatorInfo);
features = manipulatorInfo.features;

% for each model
modelDetails = trainer.modelDetails;
for m = 1:numel(modelDetails.model)
    
    % get fitter
    fitter = str2func(modelDetails.model(m).fitter);
    
    % train
    modelDetails.model(m).mdl = fitter(features',labels,'OptimizeHyperparameters',modelDetails.model(m).optimize,...
                                                         modelDetails.model(m).parameterName1,modelDetails.model(m).parameterValue1,...
                                                         modelDetails.model(m).parameterName2,modelDetails.model(m).parameterValue2);
    
    % if svm
    if strcmp(modelDetails.model(m).fitter,'fitcsvm')
        
        % get score transform
        modelDetails.model(m).mdl = fitSVMPosterior(modelDetails.model(m).mdl);
        
    end
    
    % get 10 fold cross validation accuracy
     cv = fitter(features',labels,'CrossVal','on',...
                                  modelDetails.model(m).parameterName1,modelDetails.model(m).parameterValue1,...
                                  modelDetails.model(m).parameterName2,modelDetails.model(m).parameterValue2);
                              
    % save
    modelDetails.model(m).crossValidationAccuracy = 1-kfoldLoss(cv);
    
end

trainer.modelDetails = modelDetails;

%% save to classifier

trainer = rmfield(trainer,'dateInitialized');
classifier = trainer;
classifier.trainerName = classifier.name;
classifier.dateTrained = date;

% if democratic
keepFeatures = 0;
if strcmp(modelDetails.type,'democratic')
    
    % and loops > 0
    if modelDetails.nLoops > 0
        
        % then keep features
        keepFeatures = 1;
        
    end
    
end

% keep features?
if ~keepFeatures
    
    classNames = classifier.featureSet.classNames;
    classifier.featureSet.class.(classNames{1}).features = [];
    classifier.featureSet.class.(classNames{2}).features = [];
    
end
                                  
end