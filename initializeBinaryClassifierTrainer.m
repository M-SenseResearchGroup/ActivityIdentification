function [ trainer, classifier ] = initializeBinaryClassifierTrainer()
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
%% initializeBinaryClassifierTrainer

clear
close all
clc

% get msenseresearchgroup path
ok = questdlg('Select the msenseresearchgroup folder','MSENSE Path','OK',{'OK'});
if isempty(ok); error('Initialization terminated.'); end
msense = uigetdir;

% initialize trainer
trainer.dateInitialized = date;

%% model selection

% select model
models = {'Support Vector Machine' 'K-Nearest Neighbors' 'Naive Bayes' 'Decision Tree' 'Linear'};
imodel = listdlg('ListString',models,'PromptString','Select classification model (select multiple for democratic co-learning):','SelectionMode','multiple','ListSize',[350 300]);
if isempty(imodel); error('Initialization terminated.'); end
models = models(imodel);

% if multiple models
if length(models) > 1
    
    % then democratic
    modelDetails.type = 'democratic';
    
    % combination method
    combination = questdlg('Select combination method:','Combination Method','Unweighted Average Confidence','Accuracy Weighted Confidence','Accuracy Weighted Confidence');
    if isempty(combination); error('Initialization terminated.'); end
    if strcmp(combination,'Unweighted Average Confidence')
        modelDetails.combination = 'unweightedConfidence';
    else
        modelDetails.combination = 'weightedConfidence';
    end
    
    % how many loops
    loops = inputdlg('Input number of correction loops (0 to use democratic combination only)','Democratic Loops',[1 75],{'0'});
    if isempty(loops); error('Initialization terminated'); end
    modelDetails.nLoops = str2double(loops);
    
    % if loops selected
    if modelDetails.nLoops > 0
        
        % get criteria and metric for keeping unlabeled data
        keepCriteria = questdlg('Percentage or minimum confidence based criteria for keeping unlabeled data?','Keep Criteria','Percentage','Minimum Confidence','Minimum Confidence');
        if isempty(keepCriteria); error('Initialization terminated.'); end
        modelDetails.keepCriteria = keepCriteria;
        keepMetric = questdlg('Rank strength of unlabeled observation by disagreement strength or weighted confidence?','Keep Metric' ,'Disagreement Strength','Accuracy Weighted Confidence','Accuracy Weighted Confidence');
        if isempty(keepMetric); error('Initialization terminated.'); end
        modelDetails.keepMetric = keepMetric;
        
        % get threshold
        if strcmp(keepCriteria,'Percentage')
            modelDetails.keepCriteria = 'percentage';
            prompt = sprintf('Keep top what percent of observations with highest %s (0.01 to 1.0)',keepMetric);
            keepThreshold = inputdlg(prompt,'Keep Threshold',[1 100],{'0.1'});
            if isempty(keepThreshold); error('Initilization terminated.'); end
            keepThreshold = str2double(keepThreshold);
            if keepThreshold < 0.01 || keepThreshold > 1; error('Keep Threshold must be between 0.01 and 1.0'); end
            modelDetails.keepThreshold = keepThreshold;
        else
            if strcmp(keepMetric,'Disagreement Strength')
                minimum = 0.5;
                maximum = 1.0;
            else
                minimum = 0.0;
                maximum = 1.0;
            end
            prompt = sprintf('Input the minimum %s required to keep an unlabeled observation (%2.1f to %2.1f)',keepMetric,minimum,maximum);
            keepThreshold = inputdlg(prompt,'Keep Threshold',[1 100],{'0.5'});
            if isempty(keepThreshold); error('Initialization terminated.'); end
            keepThreshold = str2double(keepThreshold);
            if keepThreshold > maximum || keepThreshold < minimum; error('Threshold for %s must be between %2.1f and %2.1f.',keepMetric,minimum,maximum); end
            modelDetails.keepThreshold = keepThreshold;
        end
    end
else
    
    % single classifier
    modelDetails.type = 'singleClassifier';
    
end

%% model specification

% specify models
model(length(models)) = struct();
for m = 1:length(models)
    
    % name
    model(m).name = models{m};
    
    % svm
    if models{m}(1) == 'S'
        
        % fitter function
        model(m).fitter = 'fitcsvm';
        
        % kernel
        kernel = {'Gaussian (RBF)' 'Cubic Polynomial' 'Linear'; 'gaussian'       'polynomial'       'linear'; 'rbf' 'cubic' 'linear'};
        ikernel = listdlg('ListString',kernel(1,:),'PromptString','Select the Kernel Function:','ListSize',[200 100],'SelectionMode','Single');
        model(m).parameterName1 = 'KernelFunction';
        model(m).parameterValue1 = kernel{2,ikernel};
        
        % abbreviation
        model(m).abbreviation = ['svm_' kernel{3,ikernel}];
        
        % solver
        solver = {'Iterative Single Data Algorithm (ISDA)' 'Sequential Minimal Optimization (SMO)'; 'ISDA' 'SMO'};
        isolver = listdlg('ListString',solver(1,:),'PromptString','Select the solver algorithm:','ListSize',[250 100],'SelectionMode','Single');
        model(m).parameterName2 = 'Solver';
        model(m).parameterValue2 = solver{2,isolver};
        
        % optimize box constraint
        optimize = questdlg('Optimize Box Constraint (cost)?','Optimize Box Constraint','Yes','No','No');
        model(m).optimize = 'none';
        if strcmp(optimize,'Yes'); model(m).optimize = {'BoxConstraint'}; end
            
            
    % knn
    elseif models{m}(1) == 'K'
        
        % fitter function
        model(m).fitter = 'fitcknn';
        
        % num neighbors
        nNeighbors = inputdlg('How many neighbors?','K',[1 50],{'3'});
        model(m).parameterName1 = 'NumNeighbors';
        model(m).parameterValue1 = str2double(nNeighbors{1});
        
        % abbreviation
        model(m).abbreviation = ['nn' nNeighbors{1}];
        
        % default prior probs
        model(m).parameterName2 = 'Prior';
        model(m).parameterValue2 = 'empirical';
        
        % no optimization
        model(m).optimize = 'none';
        
    % nb
    elseif models{m}(1) == 'N'
        
        % abbreviation and fitter function
        model(m).abbreviation = 'nb';
        model(m).fitter = 'fitcnb';
        
        % no score transform and empirical priors (defaults)
        model(m).parameterName1 = 'ScoreTransform';
        model(m).parameterValue1 = 'none';
        model(m).parameterName2 = 'Prior';
        model(m).parameterValue2 = 'empirical';
        
        % no optimization
        model(m).optimize = 'none';
        
    % tree
    elseif models{m}(1) == 'D'
        
        % abbreviation and fitter function
        model(m).abbreviation = 'tree';
        model(m).fitter = 'fitctree';
        
        % Minimum leaf size
        mls = inputdlg('Input minimum leaf size:','Minimimum Leaf Size',[1 50],{'10'});
        model(m).parameterName1 = 'MinLeafSize';
        model(m).parameterValue1 = str2double(mls{1});
        
        % merge leaves (default)
        model(m).parameterName2 = 'MergeLeaves';
        model(m).parameterValue2 = 'on';
        
        % no optimization
        model(m).optimize = 'none';
    
    % linear
    elseif models{m}(1) == 'L'
        
        % fitter function
        model(m).fitter = 'fitclinear';
        
        % regularization
        model(m).parameterName1 = 'Regularization';
        reg = questdlg('Lasso or ridge regularization?','Regularization','Ridge','Lasso','Ridge');
        if strcmp(reg,'Ridge'); model(m).parameterValue1 = 'ridge';
        elseif strcmp(reg,'Lasso'); model(m).parameterValue1 = 'lasso'; end
        
        % learner
        model(m).parameterName2 = 'Learner';
        learner = questdlg('SVM or Logistic linear model?','Learner','SVM','Logistic','Logistic');
        if strcmp(learner,'SVM'); model(m).parameterValue2 = 'svm';
        elseif strcmp(learner,'Logistic'); model(m).parameterValue2 = 'logistic'; end
        
        % abbreviation
        model(m).abbreviation = ['lin_' model(m).parameterValue2];
        
        % no optimization
        model(m).optimize = 'none';
        
    end
    
end

% save models
modelDetails.model = model;
trainer.modelDetails = modelDetails;

%% feature set selection

% get possible feature sets to train with
fsetdir = dir(fullfile(msense,'Project Data','ActivityIdentification','FeatureSets'));
fsetNames = cell(0);
i = 1;
while i <= numel(fsetdir)
    
    % delete if hidden
    if fsetdir(i).name(1) == '.'; fsetdir(i) = [];
    
    % delete if not at least 4 characters (.mat extension)
    elseif length(fsetdir(i).name) <= 4; fsetdir(i) = [];
        
    % if not .mat extension
    elseif ~strcmpi(fsetdir(i).name(end-3:end),'.mat'); fsetdir(i) = [];
        
    % otherwise save
    else
        fsetNames{i} = fsetdir(i).name;
        i = i + 1;
    end
    
end

% select feature sets
ifsetNames = listdlg('ListString',fsetNames,'PromptString','Select feature compatible feature sets:','SelectionMode','multiple','ListSize',[300,100]);
fsetdir = fsetdir(ifsetNames);

% for each feature set
fsets = cell(1,length(fsetdir));
originalClasses = cell(0);
for f = 1:length(fsetdir)
    
    % load
    fset = load(fullfile(fsetdir(f).folder,fsetdir(f).name));
    fsets{f} = fset.fset;
    clear fset
    
    % get classNames
    originalClasses = horzcat(originalClasses,fsets{f}.classNames);
    
end

% concatenate feature sets and thereby confirm compatibility
originalClasses = unique(originalClasses);
newClasses = originalClasses;
fsets = concatenateFeatureSets(originalClasses,newClasses,fsets);

% name the class with label 1
class1 = inputdlg('Name the class with label = 1:','Label 1 Class',[1 50],{'activity1'});

% get originalClass names associated with label 1 class
prompt = sprintf('Select classes belonging to the %s class:',class1{1});
inames = listdlg('ListString',originalClasses,'PromptString',prompt,'SelectionMode','multiple','ListSize',[300,300]);
class1names = originalClasses(inames);
originalClasses(inames) = [];

% name the class with label -1
notClass1 = inputdlg('Name the class with label = -1:','Label -1 Class',[1 50],{['not' cap(class1{1})]});
% get originalClass names associated with label 1 class
prompt = sprintf('Select classes belonging to the %s class:',notClass1{1});
inames = listdlg('ListString',originalClasses,'PromptString',prompt,'SelectionMode','multiple','ListSize',[300,300]);
notClass1names = originalClasses(inames);

% unpack the feature sets with new naming convention
originalClasses = horzcat(class1names,notClass1names);
newClasses = horzcat(repmat(class1,[1 length(class1names)]),repmat(notClass1,[1 length(notClass1names)]));
fsets = unpackFeatureSet(fsets,originalClasses,newClasses);

% force classNames order so that 1 corresponds to class for label 1 and 2
% correspond to class for label -1
fsets.classNames{1} = class1{1};
fsets.classNames{2} = notClass1{1};

% report class distribution and verify ok
% n observations
n(1) = fsets.class.(fsets.classNames{1}).nObservations;
n(2) = fsets.class.(fsets.classNames{2}).nObservations;
ntotal = n(1) + n(2);

% report and give option to make 50/50
question = sprintf('Total Observations: %d, %s observations: %d (%2.2f%%), %s observations: %d (%2.2f%%)',ntotal,fsets.classNames{1},n(1),n(1)/ntotal*100,fsets.classNames{2},n(2),n(2)/ntotal*100);
answer = questdlg(question,'Class Distribution','Make 50/50','Leave As Is','Leave As Is');

% if make 50/50
if strcmp(answer,'Make 50/50')
    
    % get smaller set
    nmin = min(n);

    % for each class
    for c = 1:2

        % if this class is the larger
        if n(c) > nmin

            % get n to remove
            nremove = n(c) - nmin;

            % randomly permutate indices
            iremove = randperm(n(c));

            % keep the first nremove
            iremove = sort(iremove(1:nremove),'ascend');

            % remove these features
            fsets.class.(fsets.classNames{c}).features(:,iremove) = [];
            
            % remove corresponding original class names
            fsets.class.(fsets.classNames{c}).originalClass(iremove) = [];

            % update observationsPerTrial
            indicesPerTrial = cell(length(fsets.class.(fsets.classNames{c}).observationsPerTrial),1);
            i = 0;
            % for each trial
            for t = 1:length(fsets.class.(fsets.classNames{c}).observationsPerTrial)

                % for each observation
                indicesPerTrial{t} =[];
                k = 0;
                for o = 1:fsets.class.(fsets.classNames{c}).observationsPerTrial(t)
                    i = i+1;
                    % if not removing index i
                    if ~any(i == iremove)

                        % increment k and save index
                        k = k + 1;
                        indicesPerTrial{t}(k) = i;

                    end
                end
            end
            
            % update subject indices
            snames = fieldnames(fsets.class.(fsets.classNames{c}).subject);
            
            % for each subject
            for s = 1:length(snames)
                
                % remove any indices apart of iremove
                for i = 1:length(iremove)
                    check = fsets.class.(fsets.classNames{c}).subject.(snames{s}).indices == iremove(i);
                    fsets.class.(fsets.classNames{c}).subject.(snames{s}).indices(check) = [];
                end
                
                % for each previous index to keep
                for i = 1:length(fsets.class.(fsets.classNames{c}).subject.(snames{s}).indices)
                    
                    % get number of indices removed before this one
                    nbefore = sum(iremove < fsets.class.(fsets.classNames{c}).subject.(snames{s}).indices(i));
                    
                    % subtract this from this index to get new index
                    fsets.class.(fsets.classNames{c}).subject.(snames{s}).indices(i) = fsets.class.(fsets.classNames{c}).subject.(snames{s}).indices(i) - nbefore;
                    
                end
                
            end

            % get observations per trial after removal
            fsets.class.(fsets.classNames{c}).originalObservationsPerTrial = fsets.class.(fsets.classNames{c}).observationsPerTrial;
            fsets.class.(fsets.classNames{c}).observationsPerTrial = zeros(1,length(fsets.class.(fsets.classNames{c}).originalObservationsPerTrial));
            for t = 1:length(fsets.class.(fsets.classNames{c}).observationsPerTrial)
                fsets.class.(fsets.classNames{c}).observationsPerTrial(t) = length(indicesPerTrial{t});
            end
            
            % update nObservations and report
            fsets.class.(fsets.classNames{c}).nObservations = sum(fsets.class.(fsets.classNames{c}).observationsPerTrial);
            n(1) = fsets.class.(fsets.classNames{1}).nObservations;
            n(2) = fsets.class.(fsets.classNames{2}).nObservations;
            ntotal = n(1) + n(2);
            statement = sprintf('Total Observations: %d, %s observations: %d (%2.2f%%), %s observations: %d (%2.2f%%)',ntotal,fsets.classNames{1},n(1),n(1)/ntotal*100,fsets.classNames{2},n(2),n(2)/ntotal*100);
            questdlg(statement,'New Distribution','OK','OK');

            break;

        end
    end
end

% set label
fsets.class.(class1{1}).label = 1;
fsets.class.(notClass1{1}).label = -1;

% save feature set
trainer.featureSet = fsets;

%% specify feature manipulation

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
manipulatorInfo.features = horzcat(trainer.featureSet.class.(class1{1}).features,trainer.featureSet.class.(notClass1{1}).features);
manipulatorInfo.labels = horzcat(repmat(trainer.featureSet.class.(class1{1}).label,[1 trainer.featureSet.class.(class1{1}).nObservations]),...
                                 repmat(trainer.featureSet.class.(notClass1{1}).label,[1 trainer.featureSet.class.(notClass1{1}).nObservations]));
                             
% get feature names
extractor = str2func(fsets.featureDetails.featureExtractor);
[~,manipulatorInfo.originalFeatureNames] = extractor();

% train manipulator
manipulator = str2func(featureManipulator.name);
manipulatorInfo = manipulator(manipulatorInfo);

% save to trainer
featureManipulator.manipulatorInfo = manipulatorInfo;
trainer.featureManipulator = featureManipulator;

%% save trainer

trainerName = inputdlg('Name the trainer (recommended to start with Trainer_):','Trainer Name',[1 100],{['Trainer_' class1{1} '_' notClass1{1} '_']});
trainerPath = fullfile(msense,'Project Data','ActivityIdentification','ClassifierTrainers');
trainer.name = trainerName{1};
save(fullfile(trainerPath,trainerName{1}),'trainer');

%% train now?
    
% train classifier
trainer.msensePath = msense;
classifier = trainBinaryClassifier(trainer);

% save
classifierName = replace(trainerName,'Trainer_','Classifier_');
classifierPath = fullfile(msense,'Project Data','ActivityIdentification','Classifiers');
classifier.name = classifierName{1};
classifier = rmfield(classifier,'msensePath');
save(fullfile(classifierPath,classifierName{1}),'classifier');

% loso validation?
valNow = questdlg('Perform LOSO validation?','LOSO','Yes','No','No');
if strcmp(valNow,'Yes')
    classifier.msensePath = msense;
    classifier.validation = validateClassifier_loso(classifier);
    classifier = rmfield(classifier,'msensePath');
    disp(classifier.validation.loso)

    % resave
    save(fullfile(classifierPath,classifierName{1}),'classifier');
end

end