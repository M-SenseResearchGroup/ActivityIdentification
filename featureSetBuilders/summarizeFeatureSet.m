function [ fset ] = summarizeFeatureSet(FeatureSet)
%Reed Gurchiek, 2018
%   description
%
%----------------------------------INPUTS----------------------------------
%
%   FeatureSet:
%       struct, feature set output from buildFeatureSet_
%
%---------------------------------OUTPUTS----------------------------------
%
%   none
%
%--------------------------------------------------------------------------
%% summarizeFeatureSet

% if no feature set
if nargin == 0
    
    % select and load
    ok = questdlg('Select the feature set to load','Feature Set Selection','OK','OK');
    if isempty(ok); error('Summary ended'); end
    [fsetf,fsetp] = uigetfile;
    load(fullfile(fsetp,fsetf),'fset');
    
else
    fset = FeatureSet;
end
   
fprintf('-Last Updated: %s\n',fset.dateUpdated{end});

% feature details
fDetails = fset.featureDetails;
exInfo = fDetails.extractorInfo;
fprintf('\n-Feature Details:\n');
fprintf('\t-Feature Extractor: %s\n',fDetails.featureExtractor);
fprintf('\t-N Features: %d\n',fDetails.nFeatures);
fprintf('\t-Sampling Frequency: %6.2f Hz\n\t-Window Size: %4.2f s\n\t-Overlap: %3.2f\n',exInfo.samplingFrequency,exInfo.windowSize,exInfo.overlap);

% data details
fprintf('\n-Data Details:\n')
dDetails = fset.dataDetails;
if isfield(dDetails,'dataProcessor')
    if ~isempty(dDetails.dataProcessor)
        fprintf('\t-Data Processor: %s\n',dDetails.dataProcessor);
    end
end
if isfield(dDetails,'dataCalibrator')
    if ~isempty(dDetails.dataCalibrator)
        fprintf('\t-Data Calibrator: %s\n',dDetails.dataCalibrator);
    end
end

% subject details            
fprintf('\nSubject Details:\n')
sub = fset.subject;
fprintf('\t-N Subjects: %d\n',numel(sub));
datasets = cell(1,numel(sub));
for k = 1:numel(sub)
    datasets{k} = sub(k).dataset;
end
datasets = unique(datasets);
fprintf('\t-Datasets: %s',datasets{1});
for k = 2:numel(datasets)
    fprintf(', %s',datasets{k});
end
fprintf('\n');

% Class Representation
classNames = fset.classNames;
totalObservations = 0;
for c = 1:numel(classNames)
    totalObservations = totalObservations + fset.class.(classNames{c}).nObservations;
end
fprintf('\n-Total Observations: %d\n',totalObservations);
for c = 1:numel(classNames)
    fprintf('\t(%d) %s: N Observations = %d (%4.2f%%)\n',c,classNames{c},fset.class.(classNames{c}).nObservations,fset.class.(classNames{c}).nObservations/totalObservations*100);
end

end