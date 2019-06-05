function [ validation ] = validateClassifier_loso(classifier)
%Reed Gurchiek, 2019
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
%% validateClassifier_loso

% save original feature set
class = classifier.featureSet.class;
cnames = fieldnames(class);
totalObservations = class.(cnames{1}).nObservations + class.(cnames{2}).nObservations;

% duplicate classifier as trainer for loso
trainer = classifier;

% feature manipulator
manipulator = str2func(classifier.featureManipulator.name);
manipulatorInfo = classifier.featureManipulator.manipulatorInfo;

% initialize true and predicted label arrays
trueLabels = zeros(1,totalObservations);
predictedLabels = zeros(1,totalObservations);
predictionConfidences = zeros(1,totalObservations);
originalClass = cell(1,totalObservations);

% for each subject
istart = 1;
snames = fieldnames(class.(cnames{1}).subject);
nsubjects = length(snames);
for s = 1:nsubjects
    
    % leave out subject
    out = snames{s};
    fprintf('-Leaving out %s (%d of %d)\n',out,s,nsubjects);
    
    % get indices for subject s to leave out
    indices = struct();
    for c = 1:2
        indices.(cnames{c}) = class.(cnames{c}).subject.(out).indices;
    end
    
    % save test features
    testFeatures = horzcat(class.(cnames{1}).features(:,indices.(cnames{1})),class.(cnames{2}).features(:,indices.(cnames{2})));
    
    % save test labels and original class names
    iend = istart + size(testFeatures,2) - 1;
    trueLabels(istart:iend) = horzcat(repmat(class.(cnames{1}).label,[1 length(indices.(cnames{1}))]),repmat(class.(cnames{2}).label,[1 length(indices.(cnames{2}))]));
    originalClass(istart:iend) = horzcat(class.(cnames{1}).originalClass(indices.(cnames{1})),class.(cnames{2}).originalClass(indices.(cnames{2})));
    subject.(out).trueLabels = trueLabels(istart:iend);
    subject.(out).originalClass = originalClass(istart:iend);
    
    % create trainer set to send to trainBinaryClassifier with this subject
    % features removed
    trainer.featureSet.class = class;
    for c = 1:2
        trainer.featureSet.class.(cnames{c}).features(:,indices.(cnames{c})) = [];
        trainer.featureSet.class.(cnames{c}).originalClass(indices.(cnames{c})) = [];
        trainer.featureSet.class.(cnames{c}).nObservations = size(trainer.featureSet.class.(cnames{c}).features,2);
    end
    
    % train feature manipulator
    manipulatorInfo.action = 'train';
    manipulatorInfo.features = horzcat(trainer.featureSet.class.(cnames{1}).features,trainer.featureSet.class.(cnames{2}).features);
    manipulatorInfo.labels = [repmat(trainer.featureSet.class.(cnames{1}).label,[1 trainer.featureSet.class.(cnames{1}).nObservations])...
                              repmat(trainer.featureSet.class.(cnames{2}).label,[1 trainer.featureSet.class.(cnames{2}).nObservations])];
    trainer.feautureManipulator.manipulatorInfo = manipulator(manipulatorInfo);
    
    % train classifier
    cf = trainBinaryClassifier(trainer);
    
    % manipulate test features
    manipulatorInfo = cf.featureManipulator.manipulatorInfo;
    manipulatorInfo.action = 'manipulate';
    manipulatorInfo.features = testFeatures;
    manipulatorInfo = manipulator(manipulatorInfo);
    
    % classify
    [predictedLabels(istart:iend),predictionConfidences(istart:iend)] = msenseClassify(cf,manipulatorInfo.features);
    subject.(out).predictedLabels = predictedLabels(istart:iend);
    subject.(out).predictionConfidence = predictionConfidences(istart:iend);
    
    % keep feature names if has them
    if isfield(manipulatorInfo,'featureNames')
        subject.(out).featureNames = manipulatorInfo.featureNames;
    end
    
    % update istart
    istart = iend + 1;
    
end

% evaluate
[acc,sens,spec,prec,auc,roc,err] = evalbc(trueLabels,predictedLabels,predictionConfidences,originalClass);

% save results
validation.loso.accuracy = acc;
validation.loso.sensitivity = sens;
validation.loso.specificity = spec;
validation.loso.precision = prec;
validation.loso.auc = auc;
validation.loso.roc = roc;
validation.loso.error = err;
validation.loso.subject = subject;

end

