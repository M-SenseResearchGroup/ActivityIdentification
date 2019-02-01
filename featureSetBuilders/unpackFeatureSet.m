function [ unpackedFeatureSet ] = unpackFeatureSet(FeatureSet,originalClasses,newClasses)
%Reed Gurchiek, 2018
%   takes a feature set structure (output from buildFeatureSet_) with
%   features for each class packed within the subject for whom they
%   represent and unpacks them returning a feature set with exactly the
%   same structure as the input FeatureSet but with no subject structure
%   and the features matrix for all subjects are concatenated in to a
%   single matrix in:
%
%   unpackedFeatureSet.class.(className).features
%
%   Only the classes listed in originalClasses are grabbed for each subject
%   and they are reassigned the class name according to the corresponding
%   element in newClasses
%
%----------------------------------INPUTS----------------------------------
%
%   originalClasses, newClasses:
%       cell arrays of strings, each is 1xn, the strings in originalClasses
%       should be classNames of the FeatureSet inputted and elements in 
%       newClasses with corresponding indices in originalClasses will 
%       replace the className in originalClasses from the original 
%       feature set in the new one.  For example, to assign walk and run 
%       observations to locomotion and sit to notLocmotion use:
%           originalClasses = {'walk' 'run' 'sit'}
%           newClasses = {'locomotion' 'locomotion' 'notLocomotion'}
%
%       -if originalClasses is empty then all originalClasses will be kept
%       (ie those in FeatureSet.classNames
%
%       -if newClasses is empty then all classes listed in originalClasses
%       will be kept
%
%   FeatureSet:
%       struct, output from buildFeatureSet_
%
%---------------------------------OUTPUTS----------------------------------
%
%   unpackedFeatureSet:
%       struct, exact same structure as FeatureSet except with all features
%       within the subject.class.(className).features matrix moved to the
%       class.(className).features matrix with newClass names and the
%       subject structure deleted
%
%--------------------------------------------------------------------------
%% unpackFeatureSet

% check inputs
if nargin > 1 
    if isempty(originalClasses)
        originalClasses = FeatureSet.classNames;
        newClasses = originalClasses;
    end
else
    originalClasses = FeatureSet.classNames;
    newClasses = originalClasses;
end

if nargin > 2
    if isempty(newClasses)
        newClasses = originalClasses;
    end
else
    newClasses = originalClasses;
end

% error check
if length(originalClasses) ~= length(newClasses)
    error('originalClasses and newClasses cell arrays must be same length')
end

% originalClasses must be a unique set
if length(unique(originalClasses)) ~= length(originalClasses)
    error('originalClasses must be a unique set (no repeats)')
end

% initialize unpacked feature set
subject = FeatureSet.subject;
unpackedFeatureSet = rmfield(FeatureSet,{'subject','class'});
unpackedFeatureSet.dateUpdated = {date};
unpackedFeatureSet.classNames = unique(newClasses);

% initialize class specific features matrices
for c = 1:length(unpackedFeatureSet.classNames)
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).features = zeros(unpackedFeatureSet.featureDetails.nFeatures,1);
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).nObservations = 0;
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).trialNames = {'0'};
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).observationsPerTrial = zeros(1,1);
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).originalClass = {'0'};
end

% for each subject
for s = 1:length(subject)
    
    % initialize subject indices
    ID = valfname([subject(s).dataset '_' subject(s).ID]);
    for c = 1:length(unpackedFeatureSet.classNames)
        unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).subject.(ID).indices = [];
    end
    
end
    
% for each subject
for s = 1:length(subject)
    
    
    % for each originalClass
    ID = valfname([subject(s).dataset '_' subject(s).ID]);
    for c = 1:length(originalClasses)
        
        % get indices of features being added
        firstIndex = size(unpackedFeatureSet.class.(newClasses{c}).features,2); % would need to add 1, but since the first element is a dummy zero to be removed later, no worries
        lastIndex = size(subject(s).class.(originalClasses{c}).features,2) + firstIndex - 1;
        
        % save indices corresponding to this subject
        unpackedFeatureSet.class.(newClasses{c}).subject.(ID).indices = horzcat(unpackedFeatureSet.class.(newClasses{c}).subject.(ID).indices,firstIndex:lastIndex);
        
        % concatenate to unpacked set
        unpackedFeatureSet.class.(newClasses{c}).features = horzcat(unpackedFeatureSet.class.(newClasses{c}).features,subject(s).class.(originalClasses{c}).features);
        unpackedFeatureSet.class.(newClasses{c}).trialNames = horzcat(unpackedFeatureSet.class.(newClasses{c}).trialNames,subject(s).class.(originalClasses{c}).trialNames);
        unpackedFeatureSet.class.(newClasses{c}).observationsPerTrial = horzcat(unpackedFeatureSet.class.(newClasses{c}).observationsPerTrial,subject(s).class.(originalClasses{c}).observationsPerTrial);
        unpackedFeatureSet.class.(newClasses{c}).originalClass = horzcat(unpackedFeatureSet.class.(newClasses{c}).originalClass,repmat(originalClasses(c),[1 sum(subject(s).class.(originalClasses{c}).observationsPerTrial)]));
    end
    
end

% for each new class
for c = 1:length(unpackedFeatureSet.classNames)
    
    % remove dummy element
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).features(:,1) = [];
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).trialNames(:,1) = [];
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).observationsPerTrial(:,1) = [];
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).originalClass(1) = [];
    
    % sum to get total observations
    unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).nObservations = sum(unpackedFeatureSet.class.(unpackedFeatureSet.classNames{c}).observationsPerTrial);
    
end
        

end