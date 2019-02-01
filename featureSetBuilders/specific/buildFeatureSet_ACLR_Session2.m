function [ fset ] = buildFeatureSet_ACLR_Session2(msense)
% Reed Gurchiek, 2018
%
%--------------------------------------------------------------------------
%% initialization

if nargin == 0
    % get msenseresearchgroup path
    ok = questdlg('Select the msenseresearchgroup folder','MSENSE Path','OK',{'OK'});
    if isempty(ok); error('Building stopped'); end
    msense = uigetdir;
end

% feature set path
fsetp = fullfile(msense,'Project Data','ActivityIdentification','FeatureSets');

% feature set structure name
fsetf = 'FeatureSet_ACLR_Session2';

% if feature set already exists
if exist(fullfile(fsetp,[fsetf '.mat']),'file')
    
    % set update flag
    updateExisting = 1;
    
    % load
    load(fullfile(fsetp,[fsetf '.mat']));
    subject0 = fset.subject;
    fset = rmfield(fset,'subject');
    
    % new date
    fset.dateUpdated{end+1} = date;
    
% otherwise create new
else
    
    % set update flag
    updateExisting = 0;
    
    % class
    fset.classNames = {'walk' 'crutchRight' 'crutchLeft' 'sit' 'stand' 'lie'};
    for c = 1:length(fset.classNames)
        fset.class.(fset.classNames{c}).nObservations = 0;
    end

    % feature details
    fset.featureDetails.featureExtractor = 'extractFeatures_Acc2';
    fset.featureDetails.nFeatures = 106;
    fset.featureDetails.extractorInfo.samplingFrequency = 31.25; % hertz
    fset.featureDetails.extractorInfo.windowSize = 4; % seconds
    fset.featureDetails.extractorInfo.overlap = 0; % percent, 0 - 1

    % data details
    fset.dataDetails.dataProcessor = 'dataProcessor_normRotAcc';
    fset.dataDetails.processorInfo = struct();
    fset.dataDetails.dataCalibrator = 'dataCalibrator_staticAccAttitude';
    fset.dataDetails.calibratorInfo.nSamples = 32;
    
    % date
    fset.dateUpdated = {date};
    
end

% save msense path for current session to feature set
fset.msensePath = msense;
    
%% subject data

% read master data sheet
[numdat,txtdat] = xlsread(fullfile(msense,'Raw Data','Toth_Study','master.xlsx'));

% version for session 1
version = numdat(:,7);

% only version 2 has supervised labeling
iv2 = find(version == 2) + 1;

% keep these subjects
subid = txtdat(iv2,1);

% injured side
injuredSide = txtdat(iv2,4);

% for each subject
subject = struct();
for k = 1:length(subid)
    subject(k).dataset = 'Toth_Study';
    subject(k).ID = subid{k};
    subject(k).session(1).name = 'Session_2';
    subject(k).session(1).environment(1).name = 'Home';
    subject(k).session(1).environment(1).data.sensors = {'acc'};
    subject(k).session(1).environment(1).data.locations = {'rectus_femoris_right' 'rectus_femoris_left'};
    subject(k).session(1).environment(1).data.sensorDataName = {'a' 'a'}; % row is sensors, column is locations
    subject(k).session(1).environment(1).data.dataName = {'a1' 'a2'}; % name for data, match element for element in sensorDataName
    
    % give general name to data specifying data type, body
    % segment/muscle, side of the body, etc.  This will be used to
    % make sure that if two feature sets are combined these names
    % must match to ensure features are grabbed from the exact same
    % type of data.  For example, if acc and emg data is collected
    % from the rectus_femoris_right and from another study acc is
    % collected from anterior_thigh_right, these seem incompatible
    % when really they are compatible.  Thus they would be given
    % the same dataID.  Each element in dataID is the dataID for
    % each element in dataName.
    subject(k).session(1).environment(1).data.dataID = {'right_mid_thigh_acceleration' 'left_mid_thigh_acceleration'};

    % for each class
    for c = 1:length(fset.classNames)
        subject(k).class.(fset.classNames{c}).features = [];
        subject(k).class.(fset.classNames{c}).nObservations = 0;
        subject(k).class.(fset.classNames{c}).trialNames = cell(0);
        subject(k).class.(fset.classNames{c}).session = cell(0);
        subject(k).class.(fset.classNames{c}).environment = cell(0);
        subject(k).class.(fset.classNames{c}).observationsPerTrial = [];
    end
    
    % calibration
    subject(k).session(1).environment(1).data.calibrationTrial = 'StandingCal';
    crutch = 'crutchRight';
    if injuredSide{k} == 'L'; crutch = 'crutchLeft'; end
                
    % supervised labeling
    %              mc10 trial name               class names to assign              nMiddleSeconds trimStartFromStart   trimEndFromStart trimStartFromEnd   trimEndFromEnd
    trialNames0 = {'Lie Supine'                 {'lie'}                             {''}            {0}                 {''}             {''}               {0};...
                   'Sitting'                    {'sit'}                             {''}            {0}                 {''}             {''}               {0};...
                   'StandingCal'                {'stand'}                           {''}            {0}                 {''}             {''}               {0};...
                   'Walking'                    {'walk'}                            {''}            {2}                 {''}             {''}               {2};...
                   'Crutching'                  {crutch}                            {''}            {2}                 {''}             {''}               {2}};
    subject(k).session(1).environment(1).data.trialNames = trialNames0(:,1)';
    subject(k).session(1).environment(1).data.classNames = trialNames0(:,2)';
    subject(k).session(1).environment(1).data.nMiddleSeconds =   trialNames0(:,3)';
    subject(k).session(1).environment(1).data.trimStartFromStart = trialNames0(:,4)';
    subject(k).session(1).environment(1).data.trimEndFromStart =   trialNames0(:,5)';
    subject(k).session(1).environment(1).data.trimStartFromEnd =   trialNames0(:,6)';
    subject(k).session(1).environment(1).data.trimEndFromEnd =   trialNames0(:,7)';

end

% if updating existing fset
if updateExisting
    
    % loop through all subjects
    k = 1;
    while k <= numel(subject)
        
        % get identifiers
        dataset = subject(k).dataset;
        id = subject(k).ID;
        
        % increment index assuming this is new subject
        % if it's not, index will be decremented in for loop below
        k = k + 1;
        
        % for each already identified subject
        for j = 1:length(subject0)
            
            % if match
            if strcmp(dataset,subject0(j).dataset) && strcmp(id,subject0(j).ID)
                
                % delete subject k
                k = k - 1;
                subject(k) = [];
                break;
                
            end
            
        end
        
    end
    
end

% nothing to extract if nothing new
extractNew = 1;
if isempty(subject)
    extractNew  = 0;
    
    % save old
    fset.subject = subject0;
end

% save to fset
if extractNew
    fset.subject = subject;
end
                
%% build feature set
if extractNew
    fset = buildMC10FeatureSet(fset);

    % if updating existing fset
    if updateExisting

        % new subjects
        nNew = length(fset.subject);

        % concatenate with old
        subject0(end+1:end+nNew) = fset.subject;

        % save concatenation
        fset.subject = subject0;

    end
end

%% save feature set

summarizeFeatureSet(fset);

if extractNew
    fset = rmfield(fset,'msensePath');
    save(fullfile(fsetp,fsetf),'fset');
else
    fprintf('-Nothing new to extract\n');
end

end