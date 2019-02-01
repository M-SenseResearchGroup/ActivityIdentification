%% validateDemocraticWalkClassifier
%   Reed Gurchiek, reed.gurchiek@uvm.edu, 2018
%   this code follows nearly the same path as train2StageWalkClassifier_v1
%   except the user input is automated so that only the principal
%   components which explain more than x`% of the variance are kept.  The
%   loso loop pulls out one subject, trains a model on the others, reads in
%   the at home data for the out subject, modifies the classifier using the
%   democratic method, tests the classifier on the in lab data, determines
%   accuracy, and then loops through that process 10 times, and then again
%   for each healthy subject.
%   
%--------------------------------------------------------------------------
%% trainer initialization

clear
close all
clc

% get new data
getNewData = 0;

% this is a two stage classifier.  first decides whether window is class 1
% or not then decides if class 1 is subclass 1 or other class1
classA = 'locomotion';
notClassA = 'notLocomotion';
subClassA = 'walk';
otherSubClassA = 'otherLocomotion';

% keywords in trial names associated with class1 (not case sensitive)
classAwords = {'Walk','Stair','Crutches','run'};

% keywords in trial names associated with other class1 (not subclass1)
otherSubClassAwords = {'laterally','backwards','Stair','Crutches','run'};

% feature settings
trainer.featureManipulator = 'pca'; % must be pca for this setup
trainer.pcVarianceThreshold = 2.5;
trainer.featureExtractor = 'extractFeatures_Acc2';
extractor = str2func(trainer.featureExtractor);
[~,featNames] = extractor();
trainer.nInitialFeatures = length(featNames);

% democratic setting
trainer.democraticLoops = 10;
trainer.percentNewData = 0.5; % percent of new data added to a classifier if it was minority

% window settings
trainer.sf = 31.25; % samples per second
trainer.windowSize = 4; % seconds
trainer.overlap = 0; % decimal between 0 and 1

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
trainer.classifier(2).parameter(1).name = 'CrossVal';
trainer.classifier(2).parameter(1).value = 'off'; % default
trainer.classifier(2).parameter(2).name = 'Prior';
trainer.classifier(2).parameter(2).value = 'empirical'; % default

% classifier 3
trainer.classifier(3).name = 'Decision Tree';
trainer.classifier(3).abb = 'tree';
trainer.classifier(3).fitter = str2func('fitctree');
trainer.classifier(3).optimize = 'none';
trainer.classifier(3).parameter(1).name = 'MinLeafSize';
trainer.classifier(3).parameter(1).value = 10;
trainer.classifier(3).parameter(2).name = 'CrossVal';
trainer.classifier(3).parameter(2).value = 'off'; % default

trainer.nClassifiers = length(trainer.classifier);

% msense datasets
trainer.datasets = {'Healthy Subjects' 'Healthy Subjects v2' 'Healthy Subjects v3'};

% notes
trainer.notes = {'This classifier uses the extractFeatures_Acc2 feature extractor',...
                 'The first input (a1) should be the right leg',...
                 'The second input (a2) should be the left leg'};
             
fprintf('Datsets: ')
for k = 1:length(trainer.datasets)
    fprintf('%s, ',trainer.datasets{k});
end
fprintf('\b\b\n')
fprintf('Data details: sf = %5.2f Hz, window size = %2.1f seconds, overlap = %2.2f\n',trainer.sf,trainer.windowSize,trainer.overlap);
fprintf('Feature details: feature extractor: %s, featureManipulator: %s\n',trainer.featureExtractor,trainer.featureManipulator);
fprintf('Classifier details:\n')
for k = 1:length(trainer.classifier)
    fprintf('\t(%d) %s\n',k,trainer.classifier(k).name);
end

%% dataset initialization

% open shared drive
msense = mntsmb('fs1.cems.uvm.edu','msensegroup');
pause(3);

if getNewData
    
    % for each dataset
    isub = 0;
    subject = struct();
    for k = 1:length(trainer.datasets)

        % dataset directory
        subjectDirectory = fullfile(msense,'Raw Data',trainer.datasets{k});

        % get subjects/folders
        sdir = dir(subjectDirectory);
        idir = 1;
        while idir <= length(sdir)
            if ~isfolder(fullfile(sdir(idir).folder,sdir(idir).name)); sdir(idir) = [];
            elseif sdir(idir).name(1) == '.' || sdir(idir).name(1) == 'p'; sdir(idir) = [];
            else
                isub = isub + 1;
                subject(isub).ID = sdir(idir).name;
                subject(isub).folder = fullfile(subjectDirectory,subject(isub).ID);
                subject(isub).dataset = trainer.datasets{k};
                idir = idir + 1;
            end
        end

    end
    
    % ACLR subject
    isub = isub + 1;
    subject(isub).ID = 'S009';
    subject(isub).folder = fullfile(msense,'Raw Data','Toth_Study','S009');
    subject(isub).dataset = 'Toth_Study';
    
end

%% import and process data

if getNewData
    
    % for each subject
    All.(subClassA).features = zeros(trainer.nInitialFeatures,1);
    All.(otherSubClassA).features = zeros(trainer.nInitialFeatures,1);
    All.(notClassA).features = zeros(trainer.nInitialFeatures,1);
    All.(subClassA).nObservations = 0;
    All.(otherSubClassA).nObservations = 0;
    All.(notClassA).nObservations = 0;
    for k = 1:numel(subject)
        %status
        fprintf('-Processing subject %d of %d\n',k,numel(subject))
        
        % specify the body segment names (row 1) and associated mc10 sensor
        % location names (row 2)
        segmentNames = {'right_thigh','left_thigh';...                 
                        'anterior_thigh_right','anterior_thigh_left'};
            
        %get mc10 data
        if subject(k).dataset(1) == 'M'
            subject(k).folder = fullfile(subject(k).folder,'Session_1','Lab');
            % TUG/T25W not included here since we want to ignore transitions to
            % walking.  See notes in study folder
            trialNames = {'ADL: Normal Walking','ADL: Normal Standing','ADL: Normal Sitting','ADL: Slouch sitting','ADL: Lying on back','30s Chair Stand Test'};
            subject(k).calibration.trialName = 'ADL_Normal_Standing';
            
        elseif subject(k).dataset(1) == 'T'
            subject(k).folder = fullfile(subject(k).folder,'Session_2','Home');
            trialNames = {'Sitting','Lie Supine','StandingCal','Walking'};
            subject(k).calibration.trialName = 'StandingCal';
            
            % specify the body segment names (row 1) and associated mc10 sensor
            % location names (row 2)
            segmentNames = {'right_thigh','left_thigh';...                 
                            'rectus_femoris_right','rectus_femoris_left'};
        else
            subject(k).folder = fullfile(subject(k).folder,'Session_1','Lab');
            trialNames = {'Static calibration','Stand to sit','Sit','Sit to stand','Squat','Treadmill walk slow','Treadmill walk normal','Treadmill walk fast','Treadmill walk laterally left','Treadmill walk laterally right','Treadmill walk backwards','Treadmill run slow','Treadmill run comfortable','Treadmill run fast','Lie prone','Lie supine','Lie left','Lie right','Stairs ascent','Stairs descent','Crutches left','Crutches right','Hallway walk slow','Hallway walk normal','Hallway walk fast'};
            subject(k).calibration.trialName = 'Static_calibration';

            % names changed for S0007 and later
            if str2double(subject(k).ID(2:end)) > 6 || subject(k).dataset(end-1) == 'v'
                trialNames = {'Static Calibration','Stand to Sit','Sit','Sit to Stand','Squat','Treadmill Walk Slow','Treadmill Walk Normal','Treadmill Walk Fast','Treadmill Walk Laterally Left','Treadmill Walk Laterally Right','Treadmill Walk Backwards','Treadmill Run Slow','Treadmill Run Comfortable','Treadmill Run Fast','Lie Prone','Lie Supine','Lie Left','Lie Right','Stair Ascent','Stair Descent','Crutches Left','Crutches Right','Hallway Walk Slow','Hallway Walk Normal','Hallway Walk Fast'};
                subject(k).calibration.trialName = 'Static_Calibration';

                % if v2 or v3, then include calf raises, air squats
                % dont include multispeed overground walking trial since not
                % consistent enough within 4 seconds
                if subject(k).dataset(end-1) == 'v'
                    trialNames{length(trialNames)+1} = 'Calf Raises';
                    trialNames{length(trialNames)+1} = 'Air Squats';
                end
            end
        end
        sensorLocations = segmentNames(2,:);
        sensors = {'acc'};
        subject(k).trial = importMC10(fullfile(subject(k).folder,'MC10'),'trialNames',trialNames,'locationNames',sensorLocations,'sensorNames',{'acc'},'resample',trainer.sf,'reportStatus',1,'storeSameTrials','last');
        subject(k).calibration.static = subject(k).trial.(subject(k).calibration.trialName);
        subject(k).trial = rmfield(subject(k).trial,'move');
        trialNames = fieldnames(subject(k).trial);

        % build model
        if subject(k).dataset(1) == 'T'
            subject(k).calibration.static.time = subject(k).calibration.static.rectus_femoris_right.acc.t;
        else
            subject(k).calibration.static.time = subject(k).calibration.static.anterior_thigh_right.acc.t;
        end
        subject(k).calibration.static.segments = struct('right_thigh',[],'left_thigh',[]);
        for j = 1:2

            %characterize bias/uncertainty
            istill = staticndx([subject(k).calibration.static.(segmentNames{2,j}).acc.a],30);
            subject(k).calibration.static.(segmentNames{2,j}).acc.g = mean(vecnorm(subject(k).calibration.static.(segmentNames{2,j}).acc.a(:,istill(1):istill(2))));

            %segment acceleration in sensor frame
            subject(k).calibration.static.segments.(segmentNames{1,j}).acceleration = subject(k).calibration.static.(segmentNames{2,j}).acc.a/subject(k).calibration.static.(segmentNames{2,j}).acc.g;

            %get sensor to segment orientation
            subject(k).calibration.static.segments.(segmentNames{1,j}).orientation.s2b = ...
                getrot(normalize(mean(subject(k).calibration.static.segments.(segmentNames{1,j}).acceleration(:,istill(1):istill(2)),2)),[0;0;1],'dcm');

        end

        % for each trial
        fprintf('-Extracting features\n');
        subject(k).class.(subClassA).features = zeros(trainer.nInitialFeatures,1);
        subject(k).class.(subClassA).trialNames = {''; 1};
        isubClassA = 0;

        subject(k).class.(otherSubClassA).features = zeros(trainer.nInitialFeatures,1);
        subject(k).class.(otherSubClassA).trialNames = {''; 1};
        iotherSubClassA = 0;

        subject(k).class.(notClassA).features = zeros(trainer.nInitialFeatures,1);
        subject(k).class.(notClassA).trialNames = {''; 1};
        inotClassA = 0;
        for j = 1:length(trialNames)

            % for each sub trial
            for jj = 1:numel(subject(k).trial.(trialNames{j}))

                % for each segment
                for i = 1:2

                    % normalize acceleration, get in segment frame
                    subject(k).trial.(trialNames{j})(jj).segments.(segmentNames{1,i}).acceleration = subject(k).calibration.static.segments.(segmentNames{1,i}).orientation.s2b*...
                                                                                                            subject(k).trial.(trialNames{j})(jj).(segmentNames{2,i}).acc.a/...
                                                                                                            subject(k).calibration.static.(segmentNames{2,i}).acc.g;
                end

                % supervised labeling
                trimStart = round(trainer.sf*3);
                trimEnd = round(trainer.sf*3);
                if contains(trialNames{j},classAwords,'IgnoreCase',1)
                    if contains(trialNames{j},otherSubClassAwords,'IgnoreCase',1)

                        % set className and save
                        className = otherSubClassA;
                        iotherSubClassA = iotherSubClassA + 1;
                        ind = iotherSubClassA;
                        subject(k).class.(otherSubClassA).trialNames{1,ind} = trialNames{j};

                        % trim 2 seconds for these to get at least 1 window
                        % for stairs
                        trimStart = round(trainer.sf*2);
                        trimEnd = trimStart;

                    else

                        % set className and save
                        className = subClassA;
                        isubClassA = isubClassA + 1;
                        ind = isubClassA;
                        subject(k).class.(subClassA).trialNames{1,ind} = trialNames{j};

                    end
                else

                    % set className and save
                    className = notClassA;
                    inotClassA = inotClassA + 1;
                    ind = inotClassA;
                    subject(k).class.(notClassA).trialNames{1,ind} = trialNames{j};

                    % use all data
                    trimStart = 0;
                    trimEnd = 0;

                end

                % extract features
                feat = extractor(subject(k).trial.(trialNames{j})(jj).segments.right_thigh.acceleration(:,1+trimStart:end-trimEnd),...
                                                          subject(k).trial.(trialNames{j})(jj).segments.left_thigh.acceleration(:,1+trimStart:end-trimEnd),...
                                                          trainer.sf,trainer.windowSize,trainer.overlap);
                subject(k).class.(className).features = [subject(k).class.(className).features feat];

                % record number of features for this trial
                subject(k).class.(className).trialNames{2,ind} = size(feat,2);

            end
        end
        % remove dummy column, save for all and report
        subject(k).class.(subClassA).features(:,1) = [];
        subject(k).class.(subClassA).nObservations = size(subject(k).class.(subClassA).features,2);
        All.(subClassA).nObservations = All.(subClassA).nObservations + subject(k).class.(subClassA).nObservations;

        subject(k).class.(otherSubClassA).features(:,1) = [];
        subject(k).class.(otherSubClassA).nObservations = size(subject(k).class.(otherSubClassA).features,2);
        All.(otherSubClassA).nObservations = All.(otherSubClassA).nObservations + subject(k).class.(otherSubClassA).nObservations;

        subject(k).class.(notClassA).features(:,1) = [];
        subject(k).class.(notClassA).nObservations = size(subject(k).class.(notClassA).features,2);
        All.(notClassA).nObservations = All.(notClassA).nObservations + subject(k).class.(notClassA).nObservations;

        All.(subClassA).features = horzcat(All.(subClassA).features,subject(k).class.(subClassA).features);
        All.(otherSubClassA).features = horzcat(All.(otherSubClassA).features,subject(k).class.(otherSubClassA).features);
        All.(notClassA).features = horzcat(All.(notClassA).features,subject(k).class.(notClassA).features);

        fprintf('-Extracted %d %s observations (%d %s), %d %s observations\n',(subject(k).class.(subClassA).nObservations + subject(k).class.(otherSubClassA).nObservations),classA,...
                                                                                   subject(k).class.(subClassA).nObservations,subClassA,subject(k).class.(notClassA).nObservations,notClassA);

    end

    % delete trial and calibration data
    subject = rmfield(subject,'trial');
    subject = rmfield(subject,'calibration');

    % all data structure
    All.(subClassA).features(:,1) = [];
    All.(subClassA).nObservations = size(All.(subClassA).features,2);

    All.(otherSubClassA).features(:,1) = [];
    All.(otherSubClassA).nObservations = size(All.(otherSubClassA).features,2);

    All.(notClassA).features(:,1) = [];
    All.(notClassA).nObservations = size(All.(notClassA).features,2);

    All.features = [All.(subClassA).features All.(otherSubClassA).features All.(notClassA).features];
    All.nObservations = size(All.features,2);
    All.labels = [ones(1,All.(subClassA).nObservations) ones(1,All.(otherSubClassA).nObservations) -ones(1,All.(notClassA).nObservations)];

    All.(classA).features = [All.(subClassA).features All.(otherSubClassA).features];
    All.(classA).nObservations = size(All.(classA).features,2);
    All.(classA).labels = [ones(1,All.(subClassA).nObservations) -ones(1,All.(otherSubClassA).nObservations)];

    fprintf('-Import Finished\n')
    fprintf('\t-Totals: %d total observations: %d %s (%2.1f%%), %d %s (%2.1f%%)\n',All.nObservations,All.(classA).nObservations,classA,All.(classA).nObservations/All.nObservations*100,...
                                                                                                     All.(notClassA).nObservations,notClassA,All.(notClassA).nObservations/All.nObservations*100);
    fprintf('\t-%s: %d %s observations (%2.1f%%), %d %s observations (%2.1f%%)\n',cap(classA),All.(subClassA).nObservations,subClassA,...
                                                                                                        All.(subClassA).nObservations/All.(classA).nObservations*100,...
                                                                                                        All.(otherSubClassA).nObservations,otherSubClassA,...
                                                                                                        All.(otherSubClassA).nObservations/All.(classA).nObservations*100);

    % save to trainer
    trainer.(['n' cap(subClassA) 'Observations']) = All.(subClassA).nObservations;
    trainer.(['n' cap(otherSubClassA) 'Observations']) = All.(otherSubClassA).nObservations;
    trainer.(['n' cap(notClassA) 'Observations']) = All.(notClassA).nObservations;
    trainer.trainingData.subject = subject;
    
end
                                                                                       
%% LOSO validation for democratic method

fprintf('\n============================================\n LOSO\n============================================\n');

if ~getNewData
    subject = load(fullfile('/Users/reedgurchiek/Google Drive/Code/matlab_rg/projects/UVM/activityIdentification/democraticMethod/trainerStruct_HSn8_ACLRn1_7Dec2018'),'trainingData');
    subject = subject.trainingData.subject;
end

% get classifiers
model = trainer.classifier;

% for each subject
extractor = str2func(trainer.featureExtractor);
loso(numel(subject)).(classA) = struct('train',[],'test',[]);
class = {classA subClassA};
for k = 1:numel(subject)
    
    % load home data
    fprintf('-Leaving out subject %d of %d\n',k,numel(subject))
    fprintf('-Importing home MC10 data\n');
    if subject(k).dataset(1) == 'T'
        homeData = importMC10(fullfile(subject(k).folder,'MC10'),'trialNames','StandingCal','locationNames',{'rectus_femoris_right' 'rectus_femoris_left'},'sensorNames',{'acc'},'resample',trainer.sf,'reportStatus',1,'storeSameTrials','last');
        calName = 'StandingCal';
    else
        homeData = importMC10(fullfile(replace(subject(k).folder,'Lab','Home'),'MC10'),'trialNames','Posture Calibration','locationNames',{'rectus_femoris_right' 'rectus_femoris_left'},'sensorNames',{'acc'},'resample',trainer.sf,'reportStatus',1,'storeSameTrials','last');
        calName = 'Posture_Calibration';
    end
    % if no data for right/left leg
    continueFlag = 1;
    if ~isfield(homeData.move,'rectus_femoris_right') || ~isfield(homeData.move,'rectus_femoris_left')
        % flip continue flag
        continueFlag = 0;
    end
    
    if continueFlag
        
        fprintf('-Processing Data\n');
        
        % RIGHT ACC
        % get bias
        istill = staticndx(homeData.(calName).rectus_femoris_right.acc.a,round(trainer.sf));
        g1 = mean(vecnorm(homeData.(calName).rectus_femoris_right.acc.a(:,istill(1):istill(2))));

        % normalize posture acc data
        acc1 = homeData.(calName).rectus_femoris_right.acc.a/g1;

        % get sensor to segment orientation
        s2b1 = getrot(normalize(mean(acc1(:,istill(1):istill(2)),2)),[0;0;1],'dcm');

        % normalize home acceleration, get in segment frame
        acc1 = s2b1*homeData.move.rectus_femoris_right.acc.a/g1;

        % LEFT ACC
        % get bias
        istill = staticndx(homeData.(calName).rectus_femoris_left.acc.a,30);
        g2 = mean(vecnorm(homeData.(calName).rectus_femoris_left.acc.a(:,istill(1):istill(2))));

        % normalize posture acc data
        acc2 = homeData.(calName).rectus_femoris_left.acc.a/g1;

        % get sensor to segment orientation
        s2b2 = getrot(normalize(mean(acc2(:,istill(1):istill(2)),2)),[0;0;1],'dcm');

        % normalize home acceleration, get in segment frame
        acc2 = s2b2*homeData.move.rectus_femoris_left.acc.a/g2;

        clear homeData
        
        % get features
        [feat,~,windowIndices] = extractor(acc1,acc2,trainer.sf,trainer.windowSize,0,1);
        nWindows = size(windowIndices,2);
        fprintf('-Extracted %d windows\n',nWindows);

        % build train and test datasets (leave out subject k)
        loso(k).subjectOut = subject(k).ID;
        loso(k).(classA).train.features = zeros(trainer.nInitialFeatures,1);
        loso(k).(classA).train.labels = 1;
        loso(k).(classA).test.features = zeros(trainer.nInitialFeatures,1);
        loso(k).(classA).test.labels = 1;
        loso(k).(subClassA).train.features = zeros(trainer.nInitialFeatures,1);
        loso(k).(subClassA).train.labels = 1;
        loso(k).(subClassA).test.features = zeros(trainer.nInitialFeatures,1);
        loso(k).(subClassA).test.labels = 1;
        for j = 1:numel(subject)
            group = 'train';
            if j == k; group = 'test'; end
            
            % classA
            loso(k).(classA).(group).features = [loso(k).(classA).(group).features, subject(j).class.(subClassA).features, subject(j).class.(otherSubClassA).features, subject(j).class.(notClassA).features];
            loso(k).(classA).(group).labels = [loso(k).(classA).(group).labels, ones(1,subject(j).class.(subClassA).nObservations), ones(1,subject(j).class.(otherSubClassA).nObservations), -ones(1,subject(j).class.(notClassA).nObservations)];
        
            % subClassA
            loso(k).(subClassA).(group).features = [loso(k).(subClassA).(group).features, subject(j).class.(subClassA).features, subject(j).class.(otherSubClassA).features];
            loso(k).(subClassA).(group).labels = [loso(k).(subClassA).(group).labels, ones(1,subject(j).class.(subClassA).nObservations) -ones(1,subject(j).class.(otherSubClassA).nObservations)];
        end
        loso(k).(classA).train.features(:,1) = [];
        loso(k).(classA).train.labels(1) = [];
        loso(k).(classA).test.features(:,1) = [];
        loso(k).(classA).test.labels(1) = [];
        loso(k).(subClassA).train.features(:,1) = [];
        loso(k).(subClassA).train.labels(1) = [];
        loso(k).(subClassA).test.features(:,1) = [];
        loso(k).(subClassA).test.labels(1) = [];

        % for each class
        for c = 1:2
            
            % democratic corrections loop
            for j = 1:trainer.democraticLoops + 1 % + 1 since first loop is based on first extraction
                fprintf('\t-Class %s: Democratic Loop %d\n',class{c},j);

                % for each classifier
                continueFlag = 1;
                predictedLabels = zeros(length(loso(k).(class{c}).test.labels),length(model));
                for m = 1:length(model) + 1

                    % if first loop, then feat/labels are those from first extraction
                    if j == 1 && m < length(model) + 1

                        trainFeat = loso(k).(class{c}).train.features;
                        labels = loso(k).(class{c}).train.labels;
                    
                    % if last iteration
                    elseif m == length(model) + 1
                        
                        continueFlag = 0;

                    % if new data
                    elseif newData.model(m).n > 0

                        % concatenate with new features and labels
                        trainFeat = [loso(k).(class{c}).train.features newData.model(m).features];
                        labels = [loso(k).(class{c}).train.labels newData.model(m).labels];

                    % otherwise, no new data, dont continue
                    else
                        continueFlag = 0;
                    end

                    % if continue
                    if continueFlag

                        % normalize features for this training set
                        [trainFeat,mu,sigma] = zscore(trainFeat,0,2);
                        loso(k).(class{c}).loop(j).model(m).mu = mu;
                        loso(k).(class{c}).loop(j).model(m).sigma = sigma;

                        % pca for feature reduction
                        [pc,~,~,~,pcvar] = pca(trainFeat');
                        loso(k).(class{c}).loop(j).model(m).pcvar = pcvar;
                        pc = pc';

                        % get nKeep
                        pcvar(pcvar < trainer.pcVarianceThreshold) = [];
                        nKeep = length(pcvar);
                        loso(k).(class{c}).loop(j).model(m).nKeep = nKeep;
                        pc(nKeep+1:end,:) = [];
                        loso(k).(class{c}).loop(j).model(m).pc = pc;

                        % project features
                        trainFeat = pc*trainFeat;
                        testFeat = pc*diag(1./sigma)*(loso(k).(class{c}).test.features - mu);
                        
                        % train
                        loso(k).(class{c}).loop(j).model(m).mdl = ...
                            model(m).fitter(trainFeat',labels,model(m).parameter(1).name,model(m).parameter(1).value,model(m).parameter(2).name,model(m).parameter(2).value,'OptimizeHyperparameters',model(m).optimize);
                        close all
                        
                        % if svm
                        if strcmp(model(m).abb,'svm')
                            % get score transform
                            loso(k).(class{c}).loop(j).model(m).mdl = fitSVMPosterior(loso(k).(class{c}).loop(j).model(m).mdl);
                        end
                        
                        % predict test labels using test feat
                        [predictedLabels(:,m),score] = predict(loso(k).(class{c}).loop(j).model(m).mdl,testFeat');
                        loso(k).(class{c}).loop(j).model(m).predictedLabels = predictedLabels(:,m)';

                        % performance
                        [accuracy,sensitivity,specificity,precision] = evalbc(loso(k).(class{c}).test.labels,predictedLabels(:,m)');
                        loso(k).(class{c}).loop(j).model(m).accuracy = accuracy;
                        loso(k).(class{c}).loop(j).model(m).sensitivity = sensitivity;
                        loso(k).(class{c}).loop(j).model(m).specificity = specificity;
                        loso(k).(class{c}).loop(j).model(m).precision = precision;

                    % else use previous model, reset continueFlag and go to next kernel
                    else
                        
                        % if last iteration
                        if m == length(model) + 1
                            
                            % get majority
                            loso(k).(class{c}).loop(j).democratic.predictedLabels = zeros(1,size(predictedLabels,1));
                            for i = 1:size(predictedLabels,1)
                                loso(k).(class{c}).loop(j).democratic.predictedLabels(i) = sign(mean(predictedLabels(i,:)));
                            end
                            
                            % performance
                            [accuracy,sensitivity,specificity,precision] = evalbc(loso(k).(class{c}).test.labels,loso(k).(class{c}).loop(j).democratic.predictedLabels);
                            loso(k).(class{c}).loop(j).democratic.accuracy = accuracy;
                            loso(k).(class{c}).loop(j).democratic.sensitivity = sensitivity;
                            loso(k).(class{c}).loop(j).democratic.specificity = specificity;
                            loso(k).(class{c}).loop(j).democratic.precision = precision;
                            
                            
                        else

                            % copy old model to new
                            loso(k).(class{c}).loop(j).model(m).mu = loso(k).(class{c}).loop(j-1).model(m).mu;
                            loso(k).(class{c}).loop(j).model(m).sigma = loso(k).(class{c}).loop(j-1).model(m).sigma;
                            loso(k).(class{c}).loop(j).model(m).pcvar = loso(k).(class{c}).loop(j-1).model(m).pcvar;
                            loso(k).(class{c}).loop(j).model(m).nKeep = loso(k).(class{c}).loop(j-1).model(m).nKeep;
                            loso(k).(class{c}).loop(j).model(m).pc = loso(k).(class{c}).loop(j-1).model(m).pc;
                            loso(k).(class{c}).loop(j).model(m).mdl = loso(k).(class{c}).loop(j-1).model(m).mdl;
                            loso(k).(class{c}).loop(j).model(m).predictedLabels = loso(k).(class{c}).loop(j-1).model(m).predictedLabels;
                            loso(k).(class{c}).loop(j).model(m).accuracy = loso(k).(class{c}).loop(j-1).model(m).accuracy;
                            loso(k).(class{c}).loop(j).model(m).sensitivity = loso(k).(class{c}).loop(j-1).model(m).sensitivity;
                            loso(k).(class{c}).loop(j).model(m).specificity = loso(k).(class{c}).loop(j-1).model(m).specificity;
                            loso(k).(class{c}).loop(j).model(m).precision = loso(k).(class{c}).loop(j-1).model(m).precision;
                            
                            predictedLabels(:,m) = loso(k).(class{c}).loop(j).model(m).predictedLabels';
                            
                        end

                        % reset flag
                        continueFlag = 1;

                    end

                    % report
                    if j > 1 && m < length(model) + 1
                        fprintf('\t\t-%s classifier specs: nKeep = %d, accuracy = %5.4f (prev %5.4f), sensitivity = %5.4f (prev %5.4f), specificity = %5.4f (prev %5.4f), precision = %5.4f (prev %5.4f)\n',...
                            model(m).name,nKeep,accuracy,loso(k).(class{c}).loop(j-1).model(m).accuracy,...
                            sensitivity,loso(k).(class{c}).loop(j-1).model(m).sensitivity,...
                            specificity,loso(k).(class{c}).loop(j-1).model(m).specificity,...
                            precision,loso(k).(class{c}).loop(j-1).model(m).precision);
                    elseif j == 1 && m < length(model) + 1
                        fprintf('\t\t-%s classifier specs: nKeep = %d, accuracy = %5.4f, sensitivity = %5.4f, specificity = %5.4f, precision = %5.4f\n',...
                        model(m).name,nKeep,accuracy,sensitivity,specificity,precision);
                    elseif j > 1 && m == length(model) + 1
                        fprintf('\t\t-Democratic classifier specs: accuracy = %5.4f (prev %5.4f), sensitivity = %5.4f (prev %5.4f), specificity = %5.4f (prev %5.4f), precision = %5.4f (prev %5.4f)\n',...
                            accuracy,loso(k).(class{c}).loop(j-1).democratic.accuracy,...
                            sensitivity,loso(k).(class{c}).loop(j-1).democratic.sensitivity,...
                            specificity,loso(k).(class{c}).loop(j-1).democratic.specificity,...
                            precision,loso(k).(class{c}).loop(j-1).democratic.precision);
                    elseif j == 1 && m == length(model) + 1
                        fprintf('\t\t-Democratic classifier specs: accuracy = %5.4f, sensitivity = %5.4f, specificity = %5.4f, precision = %5.4f\n',...
                        accuracy,sensitivity,specificity,precision);
                    end

                end

                % initialize newData struct
                newData = struct();
                newData.model(trainer.nClassifiers) = struct('n',0,'features',zeros(trainer.nInitialFeatures,1),'labels',0,'indices',zeros(2,1));
                for m = 1:trainer.nClassifiers
                    newData.model(m).features = zeros(trainer.nInitialFeatures,1);
                    newData.model(m).labels = 0;
                    newData.model(m).indices = zeros(2,1);
                    newData.model(m).n = 0;
                    newData.model(m).majorityConfidence = 0;
                    newData.model(m).minorityConfidence = 0;
                    newData.model(m).disagreementStrength = 0;
                end

                % initialize loop params
                fprintf('\t\t-Classifying Data\n');
                wb = waitbar(0,'','Name','Classifying Data');
                nBouts = 0; % number of bouts
                activityIndices = zeros(2,1);
                activityWindows = zeros(1,1);
                activityProb = zeros(1,1);
                for w = 1:nWindows
                    
                    waitbar(w/nWindows,wb,sprintf('Classifying window %d of %d',w,nWindows))
                    
                    % if subClassA
                    continueFlag2 = 1;
                    if c == 2
                        
                        % if this window is not classA
                        if all(w ~= loso(k).(class{1}).loop(end).activityWindows)
                            
                            continueFlag2 = 0;
                            
                        end
                        
                    end
                    
                    % continue?
                    if continueFlag2

                        % manipulate, predict and score
                        lbl = zeros(length(model),1);
                        prob = zeros(length(model),1);
                        for m = 1:length(model)
                            locoFeat = loso(k).(class{c}).loop(j).model(m).pc*diag(1./loso(k).(class{c}).loop(j).model(m).sigma)*(feat(:,w) - loso(k).(class{c}).loop(j).model(m).mu);
                            [lbl(m),prob0] = predict(loso(k).(class{c}).loop(j).model(m).mdl,locoFeat');
                            if loso(k).(class{c}).loop(j).model(m).mdl.ClassNames(1) == -1
                                prob(m) = prob0(2);
                            else
                                prob(m) = prob0(1);
                            end
                        end

                        % all agree is class{c}?
                        if all(lbl == 1)

                            % store
                            nBouts = nBouts + 1;
                            activityIndices(:,nBouts) = windowIndices(:,w);
                            activityWindows(:,nBouts) = w;
                            activityProb(nBouts) = max(prob);

                        % if disagreement
                        elseif any(lbl == 1)

                            % go with majority
                            majority = sign(mean(lbl));
                            if majority == 1
                                nBouts = nBouts + 1;
                                activityIndices(:,nBouts) = windowIndices(:,w);
                                activityProb(nBouts) = max(prob);
                            end

                            % find disagreer and strength of their vote
                            dis = find(lbl ~= majority);
                            newData.model(dis).n = newData.model(dis).n + 1;
                            newData.model(dis).minorityConfidence(newData.model(dis).n) = prob(dis);
                            prob(dis) = [];
                            newData.model(dis).majorityConfidence(newData.model(dis).n) = mean(prob);

                            % add feature set, label and details
                            newData.model(dis).features(:,newData.model(dis).n) = feat(:,w);
                            newData.model(dis).labels(newData.model(dis).n) = majority;
                            newData.model(dis).indices(:,newData.model(dis).n) = windowIndices(:,w);
                            newData.model(dis).disagreementStrength(newData.model(dis).n) = abs(newData.model(dis).majorityConfidence(end) - newData.model(dis).minorityConfidence(end));
                        end
                    end
                end
                close(wb)
                
                % save indices
                loso(k).(class{c}).loop(j).nBouts = nBouts;
                loso(k).(class{c}).loop(j).activityIndices = activityIndices;
                loso(k).(class{c}).loop(j).activityWindows = activityWindows;
                loso(k).(class{c}).loop(j).activityProb = activityProb;
                
                % sort new data according to strength and report
                nDisagreements = 0;
                for m = 1:trainer.nClassifiers
                    nDisagreements = nDisagreements + newData.model(m).n;
                end
                fprintf('\t\t-Disagreements (Total = %d):\n',nDisagreements);
                nKeepDisagreements = 0;
                for m = 1:trainer.nClassifiers
                    fprintf('\t\t\t-%s: %d total disagreements, ',model(m).name,newData.model(m).n);
                    newData.model(m).n = round(trainer.percentNewData*newData.model(m).n);
                    fprintf('keeping %d, ',newData.model(m).n);
                    nKeepDisagreements = nKeepDisagreements + newData.model(m).n;
                    if newData.model(m).n > 0
                        [reorder,ireorder] = sort(newData.model(m).disagreementStrength,'descend');
                        reorder = reorder(1:newData.model(m).n);
                        fprintf('strength range (%3.2f - %3.2f), ',reorder(end),reorder(1));
                        ireorder = ireorder(1:newData.model(m).n);
                        newData.model(m).disagreementStrength = reorder;
                        newData.model(m).features = newData.model(m).features(:,ireorder);
                        newData.model(m).labels = newData.model(m).labels(ireorder);
                        newData.model(m).indices = newData.model(m).indices(:,ireorder);
                    end
                    fprintf('\b\b\n');
                end
                    
                if nKeepDisagreements == 0
                    fprintf('\t\t\t-No disagreements to keep.  Terminating democratic loop.\n');
                    % copy loop j data to loop j+1, j+2, ...
                    for i = j+1:trainer.democraticLoops + 1
                        loso(k).(class{c}).loop(i) = loso(k).(class{c}).loop(j);
                    end
                    break;

                end
            end

        end
    
    else
        fprintf('-No rectus_femoris_right and rectus_femoris_left data for subject %s.  Skipping\n',subject(k).ID);
    end
    
end
