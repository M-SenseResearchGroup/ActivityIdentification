function [ classifier ] = trainBinaryClassifier(trainer)
%Reed Gurchiek, 2018
%   trains binary classifier based on user specifications determined in
%   initializeBinaryClassifier
%
%----------------------------------INPUTS----------------------------------
%
%   trainer:
%       trainer struct from intializeBinaryClassifer, contains training
%       details, model training params, feature details and
%       manipulation/reduction details
%
%---------------------------------OUTPUTS----------------------------------
%
%   classifier:
%       same as trainer except with fitted model in modelDetails.model.mdl
%       and without the features
%
%--------------------------------------------------------------------------
%% trainBinaryClassifier

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
    
% get fitter
model = trainer.modelDetails.model;
fitter = str2func(model.fitter);

% get params
paramNames = fieldnames(model.fitterParameters);
param = model.fitterParameters;

% train
% 0 params
if isempty(paramNames)
    mdl = fitter(features',labels);
% 1 param
elseif length(paramNames) == 1
    mdl = fitter(features',labels,paramNames{1},param.(paramNames{1}));
% 2 param
elseif length(paramNames) == 2
    mdl = fitter(features',labels,paramNames{1},param.(paramNames{1}),paramNames{2},param.(paramNames{2}));
% 3 param
elseif length(paramNames) == 3
    mdl = fitter(features',labels,paramNames{1},param.(paramNames{1}),paramNames{2},param.(paramNames{2}),paramNames{3},param.(paramNames{3}));
% 4 param
elseif length(paramNames) == 4
    mdl = fitter(features',labels,paramNames{1},param.(paramNames{1}),paramNames{2},param.(paramNames{2}),paramNames{3},param.(paramNames{3}),paramNames{4},param.(paramNames{4}));
% 5 param
elseif length(paramNames) == 5
    mdl = fitter(features',labels,paramNames{1},param.(paramNames{1}),paramNames{2},param.(paramNames{2}),paramNames{3},param.(paramNames{3}),paramNames{4},param.(paramNames{4}),paramNames{5},param.(paramNames{5}));
% 6 param
elseif length(paramNames) == 6
    mdl = fitter(features',labels,paramNames{1},param.(paramNames{1}),paramNames{2},param.(paramNames{2}),paramNames{3},param.(paramNames{3}),paramNames{4},param.(paramNames{4}),paramNames{5},param.(paramNames{5}),paramNames{6},param.(paramNames{6}));
% 7 param
elseif length(paramNames) == 7
    mdl = fitter(features',labels,paramNames{1},param.(paramNames{1}),paramNames{2},param.(paramNames{2}),paramNames{3},param.(paramNames{3}),paramNames{4},param.(paramNames{4}),paramNames{5},param.(paramNames{5}),paramNames{6},param.(paramNames{6}),paramNames{7},param.(paramNames{7}));
end

% get score transform if svm
if strcmp(model.fitter,'fitcsvm')
    mdl = fitSVMPosterior(mdl);
elseif strcmp(model.fitter,'fitclinear')
    if strcmp(model.fitterParameters.Learner,'svm')
        mdl = fitSVMPosterior(mdl);
    end
end

% score threshold
if strcmp(trainer.modelDetails.scoreThresholdSelector,'none'); trainer.modelDetails.model.scoreThreshold = 0.5;
else
    [~,predictionScores] = predict(mdl,features');
    [xroc,yroc,thresholds] = perfcurve(labels',predictionScores(:,2),1);
    trainer.modelDetails.model.scoreThreshold = scoreThresholdSelector(trainer.modelDetails.scoreThresholdSelector,xroc,yroc,thresholds);
end

trainer.modelDetails.model.mdl = mdl;

% save to classifier
classifier = trainer;
                                  
end