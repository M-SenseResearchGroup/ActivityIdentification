function [ featureManipulator ] = specifyFeatureManipulator(trainer)
%Reed Gurchiek, 2019
%   select feature manipulator to use in initializeBinaryClassifier
%
%----------------------------------INPUTS----------------------------------
%
%   trainer:
%       trainer struct from initializeBinaryClassifier
%
%---------------------------------OUTPUTS----------------------------------
%
%   featureManipulator:
%       struct with details about how to manipulate/reduce features when
%       training
%
%--------------------------------------------------------------------------
%% specifyFeatureManipulator

% get featureManipulators
manipulatorDir = dir(fxndir('featureManipulator_PCA'));
manipulatorFxns = cell(0);
i = 1;
while i <= numel(manipulatorDir)
    
    % delete if hidden
    if manipulatorDir(i).name(1) == '.'; manipulatorDir(i) = [];
    
    % delete if not at least 21 characters (.mat extension)
    elseif length(manipulatorDir(i).name) <= 21; manipulatorDir(i) = [];
        
    % otherwise save
    else
        manipulatorFxns{i} = manipulatorDir(i).name;
        i = i + 1;
    end
    
end

% select manipulator fxn
ifxn = listdlg('ListString',manipulatorFxns,'PromptString','Select the feature manipulator to use:','SelectionMode','Single','ListSize',[200 100]);
featureManipulator.name = replace(manipulatorFxns{ifxn},'.m','');

% set manipulatorInfo
manipulatorInfo.action = 'train';
manipulatorInfo.features = horzcat(trainer.featureSet.class.(trainer.featureSet.classNames{1}).features,trainer.featureSet.class.(trainer.featureSet.classNames{2}).features);
manipulatorInfo.labels = horzcat(repmat(trainer.featureSet.class.(trainer.featureSet.classNames{1}).label,[1 trainer.featureSet.class.(trainer.featureSet.classNames{1}).nObservations]),...
                                 repmat(trainer.featureSet.class.(trainer.featureSet.classNames{2}).label,[1 trainer.featureSet.class.(trainer.featureSet.classNames{2}).nObservations]));
                             
% get feature names
extractor = str2func(trainer.featureSet.featureDetails.featureExtractor);
[~,manipulatorInfo.originalFeatureNames] = extractor();

% train manipulator
manipulator = str2func(featureManipulator.name);
manipulatorInfo = manipulator(manipulatorInfo);
featureManipulator.manipulatorInfo = manipulatorInfo;

end