function [ newFeatureSet ] = concatenateFeatureSets(originalClasses,newClasses,FeatureSetsCellArray)
%Reed Gurchiek, 2018
%   builds a feature set out of multiple feature sets and can assign old 
%   classes to new ones.  For example, if walk, run, sit, lie were the
%   original classes, one could concatenate multiple (or single) feature
%   sets assigning all walk and run observations to locomotion and all sit
%   and lie observations to notLocomotion classes.
%
%----------------------------------INPUTS----------------------------------
%
%   originalClasses, newClasses:
%       cell arrays of strings, each is 1xn, the strings in originalClasses
%       should be classNames of the feature sets inputted and elements in 
%       newClasses with corresponding indices in originalClasses will 
%       replace the className in originalClasses from the original 
%       feature set in the new one.  For example, to assign walk and run 
%       observations to locomotion and sit to notLocmotion use:
%           originalClasses = {'walk' 'run' 'sit'}
%           newClasses = {'locomotion' 'locomotion' 'notLocomotion'}
%
%   FeatureSetsCellArray:
%       feature set structures from buildFeatureSet_, feature sets to be
%       concatenated must have the same feature extractor, sampling
%       frequency, window size, data processor, and data calibrator. The
%       same subject cannot be in more than one feature set.  Otherwise an
%       error will be thrown
%
%---------------------------------OUTPUTS----------------------------------
%
%   newFeatureSet:
%       feature set struct with all subjects from input feature sets with
%       all originalClasses names replaced with newClasses names
%
%--------------------------------------------------------------------------
%% concatenateFeatureSets

% error check
if length(originalClasses) ~= length(newClasses)
    error('originalClasses and newClasses cell arrays must be same length')
end

% originalClasses must be a unique set
if length(unique(originalClasses)) ~= length(originalClasses)
    error('originalClasses must be a unique set (no repeats)')
end

% get feature set characteristics
n = length(FeatureSetsCellArray);
fset0 = FeatureSetsCellArray{n};
dataProcessor = fset0.dataDetails.dataProcessor;
dataCalibrator = fset0.dataDetails.dataCalibrator;
featureExtractor = fset0.featureDetails.featureExtractor;
samplingFrequency = num2str(fset0.featureDetails.extractorInfo.samplingFrequency);
windowSize = num2str(fset0.featureDetails.extractorInfo.windowSize);
[rows,columns] = size(fset0.subject(1).session(1).environment(1).data.dataID);
dataIDs = cell(1,rows*columns);
element = 0;
for r = 1:rows
    for c = 1:columns
        element = element + 1;
        dataIDs{element} = fset0.subject(1).session(1).environment(1).data.dataID{r,c};
    end
end
fsetCharacteristics0 = horzcat(dataProcessor,dataCalibrator,featureExtractor,samplingFrequency,windowSize,dataIDs);

% get feature sets and verify same characteristics
for k = 1:n
    fset(k) = FeatureSetsCellArray{k};
    dataProcessor = fset(k).dataDetails.dataProcessor;
    dataCalibrator = fset(k).dataDetails.dataCalibrator;
    featureExtractor = fset(k).featureDetails.featureExtractor;
    samplingFrequency = num2str(fset(k).featureDetails.extractorInfo.samplingFrequency);
    windowSize = num2str(fset(k).featureDetails.extractorInfo.windowSize);
    [rows,columns] = size(fset(k).subject(1).session(1).environment(1).data.dataID);
    dataIDs = cell(1,rows*columns);
    element = 0;
    for r = 1:rows
        for c = 1:columns
            element = element + 1;
            dataIDs{element} = fset(k).subject(1).session(1).environment(1).data.dataID{r,c};
        end
    end
    fsetCharacteristics = horzcat(dataProcessor,dataCalibrator,featureExtractor,samplingFrequency,windowSize,dataIDs);
    if ~all(strcmp(fsetCharacteristics,fsetCharacteristics0)); error('Characteristics for input FeatureSet %d do not match those for FeatureSet %d',k,n);end
end

% verify no repeat subjects and build subjectController housing unique
% subjects across all input feature sets and unique sessions
subjectController = struct();
for f = 1:n
    nsub = numel(fset(f).subject);
    for s = 1:nsub
        dataset = valfname(fset(f).subject(s).dataset);
        id = fset(f).subject(s).ID;
        session = fset(f).subject(s).session(1).name;
        env = fset(f).subject(s).session(1).environment(1).name;
        if isfield(subjectController,dataset)
            if isfield(subjectController.(dataset),id)
                if isfield(subjectController.(dataset).(id),session)
                    if isfield(subjectController.(dataset).(id).(session),env)
                        error('Data for subject %s from dataset %s during session %s in the %s environement appears in more than one feature set.  There must be no repeat observations.',id,dataset,session,env);
                    end
                end
            end
        end
        
        % if no error, then add 
        subjectController.(dataset).(id).(session).(env).ifset = f;
        subjectController.(dataset).(id).(session).(env).isub = s;
        
    end
end
            
        
% initialize newFeatureSet classes
newFeatureSet.dateUpdated = {date};
uNewClasses = unique(newClasses);
newFeatureSet.classNames = uNewClasses;
for c = 1:length(uNewClasses)
    newFeatureSet.class.(uNewClasses{c}).nObservations = 0;
end
newFeatureSet.featureDetails = fset(1).featureDetails;
newFeatureSet.dataDetails = fset(1).dataDetails;

% for each original dataset
dataset = fieldnames(subjectController);
isub = 0;
for d = 1:length(dataset)
    
    % for each subject
    id = fieldnames(subjectController.(dataset{d}));
    for i = 1:length(id)
        
        % initialize new subject struct
        isub = isub + 1;
        
        % initialize new classes
        for c = 1:length(uNewClasses)
            newFeatureSet.subject(isub).class.(uNewClasses{c}) = struct('features',zeros(fset(1).featureDetails.nFeatures,1),'nObservations',0,'trialNames',{'0'},'session',{'0'},'environment',{'0'},'observationsPerTrial',zeros(1,1));
        end
        
        % for each session
        session = fieldnames(subjectController.(dataset{d}).(id{i}));
        for s = 1:length(session)
            
            % for each environment
            env = fieldnames(subjectController.(dataset{d}).(id{i}).(session{s}));
            for e = 1:length(env)
                
                % get feature set and subject indices
                ifset = subjectController.(dataset{d}).(id{i}).(session{s}).(env{e}).ifset;
                isub0 = subjectController.(dataset{d}).(id{i}).(session{s}).(env{e}).isub;
                
                % dataset and isub
                newFeatureSet.subject(isub).dataset = fset(ifset).subject(isub0).dataset;
                newFeatureSet.subject(isub).ID = fset(ifset).subject(isub0).ID;
                
                % get data
                newFeatureSet.subject(isub).session(s).environment(e).data = fset(ifset).subject(isub0).session(1).environment(1).data;
                
                % save subject index in newFeatureSet
                subjectController.(dataset{d}).(id{i}).(session{s}).(env{e}).isubnew = isub;
                
            end
        
        end
        
    end
    
end

% for each original dataset
dataset = fieldnames(subjectController);
for d = 1:length(dataset)
    
    % for each subject
    id = fieldnames(subjectController.(dataset{d}));
    for i = 1:length(id)
        
        % for each session
        session = fieldnames(subjectController.(dataset{d}).(id{i}));
        for s = 1:length(session)
            
            % for each environment
            env = fieldnames(subjectController.(dataset{d}).(id{i}).(session{s}));
            for e = 1:length(env)
                
                % get indices for original feature set, subnum, and new
                % subnum
                ifset = subjectController.(dataset{d}).(id{i}).(session{s}).(env{e}).ifset;
                isub0 = subjectController.(dataset{d}).(id{i}).(session{s}).(env{e}).isub;
                isub = subjectController.(dataset{d}).(id{i}).(session{s}).(env{e}).isubnew;
                
                % get original class
                class0 = fset(ifset).subject(isub0).class;
                
                % for each original class
                for c = 1:length(originalClasses)
                    
                    % if class available
                    if isfield(class0,originalClasses{c})
                        
                        % add to concatenated feature set under new class
                        newFeatureSet.subject(isub).class.(newClasses{c}).features = horzcat(newFeatureSet.subject(isub).class.(newClasses{c}).features,class0.(originalClasses{c}).features);
                        newFeatureSet.subject(isub).class.(newClasses{c}).trialNames = horzcat(newFeatureSet.subject(isub).class.(newClasses{c}).trialNames,class0.(originalClasses{c}).trialNames);
                        newFeatureSet.subject(isub).class.(newClasses{c}).session = horzcat(newFeatureSet.subject(isub).class.(newClasses{c}).session,class0.(originalClasses{c}).session);
                        newFeatureSet.subject(isub).class.(newClasses{c}).environment = horzcat(newFeatureSet.subject(isub).class.(newClasses{c}).environment,class0.(originalClasses{c}).environment);
                        newFeatureSet.subject(isub).class.(newClasses{c}).observationsPerTrial = horzcat(newFeatureSet.subject(isub).class.(newClasses{c}).observationsPerTrial,class0.(originalClasses{c}).observationsPerTrial);
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

% for each new subject
for s = 1:length(newFeatureSet.subject)
    
    % for each new class
    for c = 1:length(uNewClasses)
        
        % delete dummy first element and get total observations
        newFeatureSet.subject(s).class.(uNewClasses{c}).features(:,1) = [];
        newFeatureSet.subject(s).class.(uNewClasses{c}).trialNames(1) = [];
        newFeatureSet.subject(s).class.(uNewClasses{c}).session(1) = [];
        newFeatureSet.subject(s).class.(uNewClasses{c}).environment(1) = [];
        newFeatureSet.subject(s).class.(uNewClasses{c}).observationsPerTrial(1) = [];
        newFeatureSet.subject(s).class.(uNewClasses{c}).nObservations = sum(newFeatureSet.subject(s).class.(uNewClasses{c}).observationsPerTrial);
        
        % increment grand total of class observations
        newFeatureSet.class.(uNewClasses{c}).nObservations = newFeatureSet.class.(uNewClasses{c}).nObservations + newFeatureSet.subject(s).class.(uNewClasses{c}).nObservations;
        
    end
    
end


end