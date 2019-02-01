function [loso] = losoDemocratic(trainer)

fprintf('\n============================================\n LOSO\n============================================\n');

% get classifiers
model = trainer.classifier;

% subject
subject = trainer.subject;

% classes
classA = trainer.classA;
notClassA = trainer.notClassA;
subClassA = trainer.subClassA;
otherSubClassA = trainer.otherSubClassA;

% for each subject
extractor = str2func(trainer.featureExtractor);
loso(numel(subject)).(classA) = struct('train',[],'test',[]);
class = {classA subClassA};
for k = 1:numel(subject)
    
    % load home data
    fprintf('\n-Leaving out subject %d of %d, Dataset: %s, ID: %s\n',k,numel(subject),subject(k).dataset,subject(k).ID)
    fprintf('\t-Importing home MC10 data\n');
    homeData = importMC10(subject(k).homeMC10Folder,'trialNames',subject(k).calibrationName,'locationNames',{'rectus_femoris_right' 'rectus_femoris_left'},'sensorNames',{'acc'},'resample',trainer.sf,'storeSameTrials','last');
    calName = valfname(subject(k).calibrationName);
    
    % if no data for right/left leg
    continueFlag = 1;
    if ~isfield(homeData.move,'rectus_femoris_right') || ~isfield(homeData.move,'rectus_femoris_left')
        % flip continue flag
        continueFlag = 0;
    end
    
    if continueFlag
        
        fprintf('\t-Processing Data\n');
        
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
        fprintf('\t-Extracted %d windows\n',nWindows);

        % build train and test datasets (leave out subject k)
        loso(k).subjectOut = subject(k).ID;
        loso(k).skip = 0;
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
                fprintf('\n\t\t-Classifying %s: Democratic Loop %d of %d, Leaving out subject %d of %d, criteria: %s, metric: %s, threshold = %3.2f\n',...
                    class{c},j,trainer.democraticLoops + 1,k,numel(subject),trainer.newDataCriteria,trainer.newDataMetric,trainer.newDataThreshold);

                % for each classifier
                continueFlag = 1;
                predictedLabels = zeros(length(loso(k).(class{c}).test.labels),length(model));
                predictedScores = zeros(length(loso(k).(class{c}).test.labels),length(model));
                crossvalAccuracy = zeros(1,length(model));
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
                        
                        % cross validation
                        cv = model(m).fitter(trainFeat',labels,model(m).parameter(1).name,model(m).parameter(1).value,model(m).parameter(2).name,model(m).parameter(2).value,'OptimizeHyperparameters',model(m).optimize,'CrossVal','on');
                        crossvalAccuracy(m) = 1-kfoldLoss(cv);
                        clear cv
                        
                        % performance on train labels
                        [predictedTrainLabels,trainScores] = predict(loso(k).(class{c}).loop(j).model(m).mdl,trainFeat');
                        [~,trainSensitivity,trainSpecificity,trainPrecision] = evalbc(labels,predictedTrainLabels');
                        loso(k).(class{c}).loop(j).model(m).trainClassification.accuracy = crossvalAccuracy(m);
                        loso(k).(class{c}).loop(j).model(m).trainClassification.sensitivity = trainSensitivity;
                        loso(k).(class{c}).loop(j).model(m).trainClassification.specificity = trainSpecificity;
                        loso(k).(class{c}).loop(j).model(m).trainClassification.precision = trainPrecision;
                        
                        % predict test labels using test feat
                        [predictedLabels(:,m),scores0] = predict(loso(k).(class{c}).loop(j).model(m).mdl,testFeat');
                        
                        % for each label
                        for i = 1:size(scores0,1)
                            
                            % get column associated with prediction
                            scoreColumn = loso(k).(class{c}).loop(j).model(m).mdl.ClassNames == predictedLabels(i,m);
                            
                            % get confidence in prediction
                            predictedScores(i,m) = scores0(i,scoreColumn);
                            
                        end
                        loso(k).(class{c}).loop(j).model(m).predictedLabels = predictedLabels(:,m)';
                        loso(k).(class{c}).loop(j).model(m).predictedScores = predictedScores(:,m);

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
                            
                            % get combined classification
                            % this combination considers each classifiers
                            % accuracy during the training stage, the
                            % predicted classification, and the confidence
                            % in the predicted classification
                            loso(k).(class{c}).loop(j).democratic.individualAccuracy = crossvalAccuracy;
                            loso(k).(class{c}).loop(j).democratic.predictedLabels = zeros(1,size(predictedLabels,1));
                            for i = 1:size(predictedLabels,1)
                                democraticCombination = sum(crossvalAccuracy.*predictedLabels(i,:).*predictedScores(i,:))/sum(crossvalAccuracy);
                                loso(k).(class{c}).loop(j).democratic.predictionConfidence(i) = abs(democraticCombination);
                                loso(k).(class{c}).loop(j).democratic.predictedLabels(i) = sign(democraticCombination);
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
                            loso(k).(class{c}).loop(j).model(m).predictedScores = loso(k).(class{c}).loop(j-1).model(m).predictedScores;
                            loso(k).(class{c}).loop(j).model(m).accuracy = loso(k).(class{c}).loop(j-1).model(m).accuracy;
                            loso(k).(class{c}).loop(j).model(m).sensitivity = loso(k).(class{c}).loop(j-1).model(m).sensitivity;
                            loso(k).(class{c}).loop(j).model(m).specificity = loso(k).(class{c}).loop(j-1).model(m).specificity;
                            loso(k).(class{c}).loop(j).model(m).precision = loso(k).(class{c}).loop(j-1).model(m).precision;
                            loso(k).(class{c}).loop(j).model(m).trainClassification = loso(k).(class{c}).loop(j-1).model(m).trainClassification;
                            
                            % update global vars for this loop
                            predictedLabels(:,m) = loso(k).(class{c}).loop(j).model(m).predictedLabels';
                            predictedScores(:,m) = loso(k).(class{c}).loop(j).model(m).predictedScores';
                            crossvalAccuracy(m) = loso(k).(class{c}).loop(j).model(m).trainClassification.accuracy;
                            
                        end

                        % reset flag
                        continueFlag = 1;

                    end

                    % report
                    if j == 1 && m < length(model) + 1
                        fprintf('\t\t\t(%d) %s: nKeep = %d, accuracy = %5.4f (cv %5.4f), sensitivity = %5.4f (train %5.4f), specificity = %5.4f (train %5.4f), precision = %5.4f (train %5.4f)\n',...
                            m,model(m).name,nKeep,loso(k).(class{c}).loop(j).model(m).accuracy,loso(k).(class{c}).loop(j).model(m).trainClassification.accuracy,...
                            loso(k).(class{c}).loop(j).model(m).sensitivity,loso(k).(class{c}).loop(j).model(m).trainClassification.sensitivity,...
                            loso(k).(class{c}).loop(j).model(m).specificity,loso(k).(class{c}).loop(j).model(m).trainClassification.specificity,...
                            loso(k).(class{c}).loop(j).model(m).precision,loso(k).(class{c}).loop(j).model(m).trainClassification.precision);
                    elseif j > 1 && m < length(model) + 1
                        fprintf('\t\t\t(%d) %s: nKeep = %d, accuracy = %5.4f (cv %5.4f, init %5.4f), sensitivity = %5.4f (init %5.4f), specificity = %5.4f (init %5.4f), precision = %5.4f (init %5.4f)\n',...
                            m,model(m).name,nKeep,loso(k).(class{c}).loop(j).model(m).accuracy,loso(k).(class{c}).loop(j).model(m).trainClassification.accuracy,loso(k).(class{c}).loop(1).model(m).accuracy,...
                            loso(k).(class{c}).loop(j).model(m).sensitivity,loso(k).(class{c}).loop(1).model(m).sensitivity,...
                            loso(k).(class{c}).loop(j).model(m).specificity,loso(k).(class{c}).loop(1).model(m).specificity,...
                            loso(k).(class{c}).loop(j).model(m).precision,loso(k).(class{c}).loop(1).model(m).precision);
                    elseif j > 1
                        fprintf('\t\t\t-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                        fprintf('\t\t\t-Democratic classifier specs: accuracy = %5.4f (init %5.4f), sensitivity = %5.4f (init %5.4f), specificity = %5.4f (init %5.4f), precision = %5.4f (init %5.4f)\n',...
                            accuracy,loso(k).(class{c}).loop(1).democratic.accuracy,...
                            sensitivity,loso(k).(class{c}).loop(1).democratic.sensitivity,...
                            specificity,loso(k).(class{c}).loop(1).democratic.specificity,...
                            precision,loso(k).(class{c}).loop(1).democratic.precision);
                        fprintf('\t\t\t-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                    elseif j == 1
                        fprintf('\t\t\t-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                        fprintf('\t\t\t-Democratic classifier specs: accuracy = %5.4f, sensitivity = %5.4f, specificity = %5.4f, precision = %5.4f\n',...
                        accuracy,sensitivity,specificity,precision);
                        fprintf('\t\t\t-----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
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
                    newData.model(m).majorityVote = 0;
                    newData.model(m).majorityConfidence = 0;
                    newData.model(m).democraticVote = 0;
                    newData.model(m).democraticConfidence = 0;
                    newData.model(m).disagreerPositiveConfidence = 0;
                    newData.model(m).disagreementStrength = 0;
                end

                % initialize loop params
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
                        prob1 = zeros(length(model),1);
                        prob = zeros(length(model),1);
                        for m = 1:length(model)
                            locoFeat = loso(k).(class{c}).loop(j).model(m).pc*diag(1./loso(k).(class{c}).loop(j).model(m).sigma)*(feat(:,w) - loso(k).(class{c}).loop(j).model(m).mu);
                            [lbl(m),prob0] = predict(loso(k).(class{c}).loop(j).model(m).mdl,locoFeat');
                            
                            % get column associated with prediction
                            columnScore = loso(k).(class{c}).loop(j).model(m).mdl.ClassNames == lbl(m);
                            
                            % get confidence in prediction
                            prob(m) = prob0(columnScore);
                            
                            % get column associated with class{c}
                            column1 = loso(k).(class{c}).loop(j).model(m).mdl.ClassNames == 1;
                            
                            % get probability is a 1
                            prob1(m) = prob0(column1);
                            
                        end
                            
                        % get democratic vote
                        democraticCombination = sum(crossvalAccuracy'.*lbl.*prob)/sum(crossvalAccuracy);
                        democraticVote = sign(democraticCombination);
                        democraticConfidence = abs(democraticCombination);
                            
                        % get majority vote and confidence
                        majorityVote = sign(mean(lbl));
                        majorityConfidence = mean(prob1(lbl == majorityVote));
                        
                        % if metric is democraticConfidence or last iteration (to get final hypothesis)
                        if strcmp(trainer.newDataMetric,'democraticConfidence') || j == trainer.democraticLoops + 1
                            
                            % then vote is democratic vote
                            vote = democraticVote;
                            voteConfidence = democraticConfidence;
                        
                        % otherwise
                        else
                            
                            % go with majority
                            vote = majorityVote;
                            voteConfidence = majorityConfidence;
                            
                        end

                        % if is class{c}
                        if vote == 1

                            % store
                            nBouts = nBouts + 1;
                            activityIndices(:,nBouts) = windowIndices(:,w);
                            activityWindows(:,nBouts) = w;
                            activityProb(nBouts) = voteConfidence;

                        % if disagreement
                        elseif any(lbl ~= vote)

                            % find disagreers and strength of their vote
                            iminority = find(lbl ~= vote);
                            for i = 1:length(iminority)
                                
                                % increment count and save features/label
                                newData.model(iminority(i)).n = newData.model(iminority(i)).n + 1;
                                newData.model(iminority(i)).features(:,newData.model(iminority(i)).n) = feat(:,w);
                                newData.model(iminority(i)).labels(newData.model(iminority(i)).n) = vote;
                                newData.model(iminority(i)).indices(:,newData.model(iminority(i)).n) = windowIndices(:,w);
                                
                                % save voting statistics
                                newData.model(iminority(i)).disagreerPositiveConfidence(newData.model(iminority(i)).n) = prob1(iminority(i));
                                newData.model(iminority(i)).majorityVote(newData.model(iminority(i)).n) = majorityVote;
                                newData.model(iminority(i)).majorityConfidence(newData.model(iminority(i)).n) = majorityConfidence;
                                newData.model(iminority(i)).democraticVote(newData.model(iminority(i)).n) = democraticVote;
                                newData.model(iminority(i)).democraticConfidence(newData.model(iminority(i)).n) = democraticConfidence;
                                newData.model(iminority(i)).disagreementStrength(newData.model(iminority(i)).n) = ...
                                    abs(newData.model(iminority(i)).majorityConfidence(end) - newData.model(iminority(i)).disagreerPositiveConfidence(end));
                            end
                                
                        end
                    end
                end
                close(wb)
                
                % save indices
                loso(k).(class{c}).loop(j).nBouts = nBouts;
                loso(k).(class{c}).loop(j).activityIndices = activityIndices;
                loso(k).(class{c}).loop(j).activityWindows = activityWindows;
                loso(k).(class{c}).loop(j).activityProb = activityProb;
                
                % get total disagreements
                nDisagreements = 0;
                for m = 1:trainer.nClassifiers
                    nDisagreements = nDisagreements + newData.model(m).n;
                end
                
                % if this is last loop
                if j == trainer.democraticLoops + 1
                     
                    % report
                    fprintf('\t\t\t\t-Total Disagreements %d, Identified %d %s bouts.\n',nDisagreements,nBouts,class{c});
                    
                else
                    
                    %report
                    fprintf('\t\t\t\t-Disagreements (Total = %d):\n',nDisagreements);
                    
                    % get name of metric associated with metric requested by user for selecting observations to keep
                    if strcmp(trainer.newDataMetric,'disagreementStrength')
                        metric = 'disagreementStrength';
                    elseif strcmp(trainer.newDataMetric,'democraticConfidence')
                        metric = 'democraticConfidence';
                    end

                    % determine which observations to keep
                    nKeepDisagreements = 0;
                    for m = 1:trainer.nClassifiers
                        fprintf('\t\t\t\t\t-%s: %d total disagreements, ',model(m).name,newData.model(m).n);

                        % if percentile criteria
                        if strcmp(trainer.newDataCriteria,'percent')

                            % only keep top x% ("top" part is done below with sort)
                            newData.model(m).n = round(trainer.newDataThreshold*newData.model(m).n);
                            if newData.model(m).n > 0
                                % reorder metric in descending order (if value based then doesn't matter, keeping all anyway since those below threshold were deleted before)
                                [~,ikeep] = sort(newData.model(m).(metric),'descend');
                                ikeep = ikeep(1:newData.model(m).n);
                            end

                        % otherwise is confidence related criteria
                        elseif strcmp(trainer.newDataCriteria,'confidence')

                            % only keep observations where strength of metric is 'enough' where enough is set by newDataThreshold
                            ikeep = newData.model(m).(metric) > trainer.newDataThreshold;
                            newData.model(m).n = length(ikeep);
                        end

                        % report what keeping
                        fprintf('keeping %d, ',newData.model(m).n);
                        nKeepDisagreements = nKeepDisagreements + newData.model(m).n;

                        % if keeping any
                        if newData.model(m).n > 0

                            % report and store values for ikeep
                            newData.model(m).(metric) = newData.model(m).(metric)(ikeep);
                            fprintf('%s range (%3.2f - %3.2f), ',trainer.newDataMetric,newData.model(m).(metric)(end),newData.model(m).(metric)(1));
                            newData.model(m).features = newData.model(m).features(:,ikeep);
                            newData.model(m).labels = newData.model(m).labels(ikeep);
                            newData.model(m).indices = newData.model(m).indices(:,ikeep);
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
            
        end
        
    else
        
        fprintf('-No rectus_femoris_right and rectus_femoris_left data for subject %s.  Skipping.\n',subject(k).ID);
        loso(k).subjectOut = subject(k).ID;
        loso(k).skip = 1;
    end
    
end

end

