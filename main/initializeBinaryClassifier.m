function [ classifier ] = initializeBinaryClassifier()
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
%% initializeBinaryClassifier

clear
close all
clc

% get ActID data path
done = 0;
while ~done
    done = 1;
    ActIDCodePath = fxndir('initializeBinaryClassifier');
    storedPath = fullfile(ActIDCodePath,'storedActIDCodePath');
    storedActID = fullfile(storedPath,'ActivityIdentificationDataPath.txt');
    if ~exist(storedPath,'dir')
        ok = questdlg('Select the ActivityIdentification data folder (contains the FeatureSets and Classifier folders)','ActID Path','OK',{'OK'});
        if isempty(ok); error('Initialization terminated.'); end
        ActID = uigetdir;
        if ActID == 0; error('Initialization terminated.'); end

        % save?
        ok = questdlg('Use this directory everytime?','ActID Path','Yes','No','Yes');
        if strcmp(ok,'Yes')
            mkdir(storedPath);
            fActID = fopen(storedActID,'w');
            fprintf(fActID,'%s',ActID);
            fclose(fActID);
        end
    else
        fActID = fopen(storedActID,'r');
        ActID = fscanf(fActID,'%c');
        fclose(fActID);
        if ~exist(ActID,'dir')
            error('The directory %s does not exist. Quit and try again or delete the folder ''storedActIDCodePath'' in %s and define new ActID data directory.',ActID,ActIDCodePath);
        end
    end
end

% initialize trainer
trainer.dateTrained = date;

%% model selection

% select model
models = {'Support Vector Machine' 'K-Nearest Neighbors' 'Naive Bayes' 'Decision Tree' 'Linear' 'Discriminant Analysis' 'Random Forest'};
imodel = listdlg('ListString',models,'PromptString','Select classification model:','SelectionMode','single','ListSize',[350 300]);
if isempty(imodel); error('Initialization terminated.'); end
model.name = models{imodel};

%% model specification

% specify model
model = specifyModel(model);

% save models
modelDetails.model = model;
trainer.modelDetails = modelDetails;

%% feature set selection

% get possible feature sets to train with
fsetdir = dir(fullfile(ActID,'FeatureSets'));
fsets = specifyFeatureSet(fsetdir);

% save feature set
trainer.featureSet = fsets;

%% specify feature manipulation

% get manipulator details
manipulatorInfo = specifyFeatureManipulator(trainer);

% save to trainer
featureManipulator = manipulatorInfo;
trainer.featureManipulator = featureManipulator;

%% specify scoreThresholdSelector

% 0.5 or use selector
selectors = {'Use Standard 0.5','none';'Youden''s Index','youden'; 'Distance to Corner', 'cornerDistance'; 'Positive Likelihood Ratio', 'PLR'; 'Negative Likelihood Ratio', 'NLR'; 'Diagnostic Odds Ratio', 'DOR'};
iselect = listdlg('ListString',selectors(:,1),'PromptString','Choose the selector function to determine the score cutoff (or select ''Use Standard 0.5'')','SelectionMode','single','ListSize',[400,200],'OKString','Select','CancelString','Use Standard 0.5');
if isempty(iselect); iselect = 1; end
trainer.modelDetails.scoreThresholdSelector = selectors{iselect,2};

%% name classifier

% get/verify name
unverified = 1;
while unverified
    incorrectStart = 1;
    while incorrectStart
        classifierName = inputdlg('Name the classifier (must start with Classifier_):','Classifier Name',[1 100],{['Classifier_' trainer.featureSet.classNames{1} '_' trainer.featureSet.classNames{2} '_']});
        if length(classifierName{1}) < 11
            questdlg('Classifier name must start with Classifier_. Try again.','Warning','OK','OK');
        elseif ~strcmp(classifierName{1}(1:11),'Classifier_')
            questdlg('Classifier name must start with Classifier_. Try again.','Warning','OK','OK');
        else
            incorrectStart = 0;
        end
    end
    classifierPath = fullfile(ActID,'Classifiers');
    if isfile(fullfile(classifierPath,[classifierName{1} '.mat']))
        warned = questdlg('A classifier with this name already exists.  Overwrite or Rename?','Warning','Overwrite','Rename','Rename');
        if strcmp(warned,'Overwrite')
            unverified = 0;
        end
    else
        unverified = 0;
    end
end
trainer.name = classifierName{1};

%% train/validate classifier

% loso validation?
performLoso = questdlg('Perform LOSO validation after training?','LOSO','Yes','No','No');
    
% train classifier
classifier = trainBinaryClassifier(trainer);

% loso
if strcmp(performLoso,'Yes')
    classifier.validation = validateClassifier_loso(classifier);
    disp(classifier.validation.loso)
end

% remove features
classifier.featureSet.class.(classifier.featureSet.classNames{1}).features = [];
classifier.featureSet.class.(classifier.featureSet.classNames{2}).features = [];

% save
save(fullfile(classifierPath,classifierName{1}),'classifier');

end