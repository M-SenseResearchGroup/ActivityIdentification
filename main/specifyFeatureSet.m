function [ fsets ] = specifyFeatureSet(featureSetsDirectory)
%Reed Gurchiek, 2019
%   specify feature sets to use in initializeBinaryClassifier. Loads,
%   unpacks, concatenates, and evenly distributes multiple compatible
%   feature set structs for use within ActivityIdentification
%
%----------------------------------INPUTS----------------------------------
%
%   featureSetsDirectory:
%       path to feature sets
%
%---------------------------------OUTPUTS----------------------------------
%
%   fset:
%       feature set struct
%
%--------------------------------------------------------------------------
%% featureSetsDirectory

% get possible feature sets to train with
fsetdir = featureSetsDirectory;
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

end