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

% get msense path if not given
if ~isfield(classifier,'msensePath')
    
    % get msenseresearchgroup path
    ok = questdlg('Select the msenseresearchgroup folder','MSENSE Path','OK',{'OK'});
    if isempty(ok); error('Initialization terminated.'); end
    msense = uigetdir;
    
else
    msense = classifier.msensePath;
end

% load trainer
trainer = load(fullfile(msense,'Project Data','ActivityIdentification','ClassifierTrainers',classifier.trainerName),'trainer');
trainer = trainer.trainer;

% save original feature set
class = trainer.featureSet.class;
cnames = fieldnames(class);
totalObservations = class.(cnames{1}).nObservations + class.(cnames{2}).nObservations;

% initialize true and predicted label arrays
trueLabels = zeros(1,totalObservations);
predictedLabels = zeros(1,totalObservations);
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
    
    % create trainer set to send to trainBinaryClassifier with this subject
    % features removed
    trainer.featureSet.class = class;
    for c = 1:2
        trainer.featureSet.class.(cnames{c}).features(:,indices.(cnames{c})) = [];
        trainer.featureSet.class.(cnames{c}).originalClass(indices.(cnames{c})) = [];
        trainer.featureSet.class.(cnames{c}).nObservations = size(trainer.featureSet.class.(cnames{c}).features,2);
    end
    
    % train feature manipulator
    manipulator = str2func(trainer.featureManipulator.name);
    manipulatorInfo = trainer.featureManipulator.manipulatorInfo;
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
    predictedLabels(istart:iend) = msenseClassify(cf,manipulatorInfo.features);
    
    % update istart
    istart = iend + 1;
    
end

% evaluate
[acc,sens,spec,prec,err] = evalbc(trueLabels,predictedLabels,originalClass);

% save results
validation.loso.accuracy = acc;
validation.loso.sensitivity = sens;
validation.loso.specificity = spec;
validation.loso.precision = prec;
validation.loso.error = err;

end

