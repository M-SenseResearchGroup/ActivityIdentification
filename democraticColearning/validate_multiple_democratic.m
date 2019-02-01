%% validateMultipleDemocraticClassifiers
% Reed Gurchiek, 2018
% compare different democratic classifier approaches for classifying walk
% bouts using accelerometer data.  Classifier details set here in trainer
% (after loading) and then loso implemented in losoDemocratic()

clear
clc

% msense open?
ok = questdlg('Make sure MSENSE shared drive open','Trainer','OK','OK');
if isempty(ok); error('Validation Terminated'); end

% get msense dir
ok = questdlg('Select MSENSE shared drive','LOSO','OK','OK');
if isempty(ok); error('Validation Terminated'); end
msensedir = uigetdir();

% load trainer
ok = questdlg('Select trainer file','Trainer','OK','OK');
if isempty(ok); error('Validation Terminated'); end
[ftrainer,ptrainer] = uigetfile();
trainer = load(fullfile(ptrainer,ftrainer));

% directory to save loso results to
ok = questdlg('Select directory to save loso results to','LOSO','OK','OK');
if isempty(ok); error('Validation Terminated'); end
ploso = uigetdir(ptrainer);
floso = inputdlg('Save LOSO results as:','LOSO',[1 100],{['DCL_losoValidation_' date]});
floso = floso{1};

% for each subject
for k = 1:length(trainer.subject)
    % get right msense dir
    fld = trainer.subject(k).folder;
    fld = fullfile(msensedir,fld(regexp(fld,'Raw Data'):end));
    
    % compensate for windows
    if filesep == '\'; replace(fld,'/','\'); end
    
    % get home data folder and posture calibration name
    if trainer.subject(k).dataset(1) == 'H'
        trainer.subject(k).homeMC10Folder = fullfile(replace(fld,'Lab','Home'),'MC10');
        trainer.subject(k).calibrationName = 'Posture Calibration';
        subnum = str2double(trainer.subject(k).ID(2:end));
        if any(subnum == [3 7 8])
            trainer.subject(k).calibrationName = 'Standing';
        end
    else
        trainer.subject(k).homeMC10Folder = fullfile(fld,'MC10');
        trainer.subject(k).calibrationName = 'StandingCal';
    end
end
        

%% democratic settings

% feature details
trainer.featureManipulator = 'pca'; % must be pca for this setup
trainer.pcVarianceThreshold = 2.5;

% democratic training
trainer.democraticLoops = 20;
% others set in loop below

% classifier 1
trainer.classifier(1).name = 'Gaussian SVM';
trainer.classifier(1).abb = 'svm';
trainer.classifier(1).fitter = str2func('fitcsvm');
trainer.classifier(1).optimize = 'none';
trainer.classifier(1).parameter(1).name = 'KernelFunction';
trainer.classifier(1).parameter(1).value = 'gaussian';
trainer.classifier(1).parameter(2).name = 'Solver';
trainer.classifier(1).parameter(2).value = 'SMO';

% classifier 2
trainer.classifier(2).name = 'Naive Bayes';
trainer.classifier(2).abb = 'nb';
trainer.classifier(2).fitter = str2func('fitcnb');
trainer.classifier(2).optimize = 'none';
trainer.classifier(2).parameter(1).name = 'ClassNames';
trainer.classifier(2).parameter(1).value = [-1 1]; % default
trainer.classifier(2).parameter(2).name = 'Prior';
trainer.classifier(2).parameter(2).value = 'empirical'; % default

% classifier 3
trainer.classifier(3).name = 'Decision Tree';
trainer.classifier(3).abb = 'tree';
trainer.classifier(3).fitter = str2func('fitctree');
trainer.classifier(3).optimize = 'none';
trainer.classifier(3).parameter(1).name = 'MinLeafSize';
trainer.classifier(3).parameter(1).value = 10;
trainer.classifier(3).parameter(2).name = 'MergeLeaves';
trainer.classifier(3).parameter(2).value = 'on'; % default

% classifier 4
trainer.classifier(4).name = '3-Nearest Neighbor';
trainer.classifier(4).abb = '3nn';
trainer.classifier(4).fitter = str2func('fitcknn');
trainer.classifier(4).optimize = 'none';
trainer.classifier(4).parameter(1).name = 'NumNeighbors';
trainer.classifier(4).parameter(1).value = 3;
trainer.classifier(4).parameter(2).name = 'Prior';
trainer.classifier(4).parameter(2).value = 'empirical'; % default

% classifier 5
trainer.classifier(5).name = 'Logistic';
trainer.classifier(5).abb = 'log';
trainer.classifier(5).fitter = str2func('fitclinear');
trainer.classifier(5).optimize = 'none';
trainer.classifier(5).parameter(1).name = 'Learner';
trainer.classifier(5).parameter(1).value = 'logistic';
trainer.classifier(5).parameter(2).name = 'Regularization';
trainer.classifier(5).parameter(2).value = 'ridge'; % default

trainer.nClassifiers = length(trainer.classifier);

% save trainer
validation.date = date;
validation.trainer = fullfile(ptrainer,ftrainer);

%% Loso Loop

criteria = {'percent'; 'percent'; 'confidence'; 'confidence'};
metric = {'disagreementStrength'; 'democraticConfidence'; 'disagreementStrength'; 'democraticConfidence'};
threshold = { 0.5;...
             [0.1 0.2 0.3 0.4 0.5];...
             [0.95 0.9 0.85 0.8 0.75];...
             [0.95 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1]};
    
% for each criteria (columns of threshold)
ival = 0;
for k = 1:length(criteria)
    
    % for each metric (rows of threshold)
    for j = 1:length(threshold{k})
        
        fprintf('\n\nCriteria = %s, Metric = %s, threshold = %3.2f\n\n',criteria{k},metric{k},threshold{k}(j));

        % set criteria and threshold
        trainer.newDataCriteria = criteria{k};
        trainer.newDataMetric = metric{k};
        trainer.newDataThreshold = threshold{k}(j);

        % run loso
        ival = ival + 1;
        validation.loop(ival).newDataCriteria = criteria{k};
        validation.loop(ival).newDataMetric = metric{k};
        validation.loop(ival).newDataThreshold = threshold{k}(j);
        validation.loop(ival).loso = losoDemocratic(trainer);
        
    end
    
end

% save
save(fullfile(ploso,floso),'validation');

