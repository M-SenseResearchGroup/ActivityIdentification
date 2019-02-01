function [ fset ] = buildMC10FeatureSet(fset)
%Reed Gurchiek, 2018
%   builds a feature set using MC10 sensor data
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
%--------------------------------REFERENCES--------------------------------
%
%   LastName et al. (year), Title.
%
%--------------------------------------------------------------------------
%% initialization

% pull vars;
extractor = str2func(fset.featureDetails.featureExtractor);
extractorInfo = fset.featureDetails.extractorInfo;
sf = extractorInfo.samplingFrequency;
subject = fset.subject;
msense = fset.msensePath;

% calibrate?
calibrate = 0;
if isfield(fset.dataDetails,'dataCalibrator')
    if ~isempty(fset.dataDetails.dataCalibrator)
        calibrate = 1;
        calibrator = str2func(fset.dataDetails.dataCalibrator);
        calibratorInfo = fset.dataDetails.calibratorInfo;
    end
end

% process?
process = 0;
if isfield(fset.dataDetails,'dataProcessor')
    if ~isempty(fset.dataDetails.dataProcessor)
        process = 1;
        processor = str2func(fset.dataDetails.dataProcessor);
        processorInfo = fset.dataDetails.processorInfo;
    end
end

%% import and process data

% for each subject
for k = 1:numel(subject)
    
    %status
    fprintf('\n-Processing subject %d of %d\n',k,numel(subject))
    
    % sensors and locations
    sensors = subject(k).session(1).environment(1).data.sensors;
    locations = subject(k).session(1).environment(1).data.locations;
    sensorDataName = subject(k).session(1).environment(1).data.sensorDataName;
    
    % get trialNames to read and associated classNames and trim start/end
    trialNames = subject(k).session(1).environment(1).data.trialNames;
    classNames = subject(k).session(1).environment(1).data.classNames;
    nMiddleSeconds = subject(k).session(1).environment(1).data.nMiddleSeconds;
    trimStartFromStart = subject(k).session(1).environment(1).data.trimStartFromStart;
    trimEndFromStart = subject(k).session(1).environment(1).data.trimEndFromStart;
    trimStartFromEnd = subject(k).session(1).environment(1).data.trimStartFromEnd;
    trimEndFromEnd = subject(k).session(1).environment(1).data.trimEndFromEnd;
    
    % get mc10 data
    fprintf('-Importing MC10 data\n');
    dataFolder = fullfile(msense,'Raw Data',subject(k).dataset,subject(k).ID,subject(k).session(1).name,subject(k).session(1).environment(1).name,'MC10');
    trial = importMC10(dataFolder,'trialNames',horzcat(trialNames,subject(k).session(1).environment(1).data.calibrationTrial),'locations',locations,'sensors',sensors,'resample',sf,'reportStatus',0,'storeSameTrials','last');
    trial = rmfield(trial,'move');
    
    % calibrate?
    if calibrate
        calibrationTrial = trial.(valfname(subject(k).session(1).environment(1).data.calibrationTrial));
        for l = 1:length(locations)
            for s = 1:length(sensors)
                calibratorInfo.data = calibrationTrial.(locations{l}).(sensors{s}).(sensorDataName{s});
                processorInfo.calibration.(subject(k).session(1).environment(1).data.dataName{s,l}) = calibrator(calibratorInfo);
            end
        end
    end

    % for each trial
    for j = 1:length(trialNames)
        
        % get field valid trial name to match with mc10
        trialName = valfname(trialNames{j});
        
        % if no trial with this name
        if ~isfield(trial,trialName)

            % report
            fprintf('\t-Warning: no %s trial.\n',trialNames{j});
            
        else
        
            % for each division
            for idiv = 1:length(classNames{j})

                % start/finish indices to isolate class idiv from trial j
                finish = '';
                nsamp = length(trial.(trialName).(locations{1}).(sensors{1}).(sensorDataName{1}));
                if ~isempty(nMiddleSeconds{j}{idiv})
                    start = round(0.5*nsamp-0.5*sf*nMiddleSeconds{j}{idiv});
                    finish = start + round(sf*nMiddleSeconds{j}{idiv});
                elseif ~isempty(trimStartFromStart{j}{idiv})
                    start = round(sf*trimStartFromStart{j}{idiv})+1;
                elseif ~isempty(trimStartFromEnd{j}{idiv})
                    start = nsamp-round(sf*trimStartFromEnd{j}{idiv})+1;
                else
                    error('nMiddleSeconds, trimStartFromStart, and trimStartFromEnd cannot all be empty.');
                end
                if isempty(finish)
                    if ~isempty(trimEndFromStart{j}{idiv})
                        finish = round(sf*trimEndFromStart{j}{idiv})+1;
                    elseif ~isempty(trimEndFromEnd{j}{idiv})
                        finish = nsamp-round(sf*trimEndFromEnd{j}{idiv});
                    else
                        error('trimEndFromStart and trimEndFromEnd cannot both be empty.');
                    end
                end
                if start == 0; start = 1; end
                if finish > nsamp; finish = nsamp; end
                className = classNames{j}{idiv};

                % report
                fprintf('-Trial: %s (%d of %d), class: %s',trialName,j,length(trialNames),className);

                % get data
                for l = 1:length(locations)
                    for s = 1:length(sensors)
                        data.(subject(k).session(1).environment(1).data.dataName{s,l}) = trial.(trialName).(locations{l}).(sensors{s}).(sensorDataName{s})(:,start:finish);
                    end
                end

                % process?
                if process
                    processorInfo.data = data;
                    data = processor(processorInfo);
                end

                % save to extractor
                extractorInfo.data = data;

                % extract features
                feat = extractor(extractorInfo);
                subject(k).class.(className).features = [subject(k).class.(className).features feat];

                % record number of features for this trial/class
                itrial = length(subject(k).class.(className).observationsPerTrial) + 1;
                subject(k).class.(className).observationsPerTrial(itrial) = size(feat,2);
                subject(k).class.(className).trialNames{itrial} = trialName;
                subject(k).class.(className).session{itrial} = subject(k).session(1).name;
                subject(k).class.(className).environment{itrial} = subject(k).session(1).environment(1).name;

                % report
                fprintf(', extracted %d observations\n',size(feat,2));

            end
            
        end
            
    end
    
    % get nObservations for each class
    totalObservations = 0;
    classNamesUnique = cell(0);
    for c = 1:length(classNames)
        for idiv = 1:length(classNames{c})
            if all(~strcmp(classNames{c}{idiv},classNamesUnique))
                classNamesUnique = horzcat(classNamesUnique,classNames{c}(idiv));
            end
        end
    end
    classNames = classNamesUnique;
    for c = 1:length(classNames)
        subject(k).class.(classNames{c}).nObservations = sum(subject(k).class.(classNames{c}).observationsPerTrial);
        totalObservations = totalObservations + subject(k).class.(classNames{c}).nObservations;
        fset.class.(classNames{c}).nObservations = fset.class.(classNames{c}).nObservations + subject(k).class.(classNames{c}).nObservations;
    end
    
    %report
    fprintf('-Total Observations = %d:\n\t',totalObservations);
    classNames = fieldnames(subject(k).class);
    lineCount = 0;
    maxLineCount = 5;
    for c = 1:length(classNames)
        lineCount = lineCount + 1;
        if lineCount > maxLineCount
            lineCount = 0;
            fprintf('\n\t');
        end
        fprintf('%d %s (%2.1f%%),',subject(k).class.(classNames{c}).nObservations,classNames{c},subject(k).class.(classNames{c}).nObservations/totalObservations*100);
    end
    fprintf('\b\n');

end

% report
classNames = fset.classNames;
totalObservations = 0;
for c = 1:length(classNames)
    totalObservations = totalObservations + fset.class.(classNames{c}).nObservations;
end
fprintf('\n-Total Observations = %d',totalObservations);
for c = 1:length(classNames)
    fprintf(', %d %s (%2.1f%%)',fset.class.(classNames{c}).nObservations,classNames{c},fset.class.(classNames{c}).nObservations/totalObservations*100);
end
fprintf('\n');

% save to feature set
fset.subject = subject;

end