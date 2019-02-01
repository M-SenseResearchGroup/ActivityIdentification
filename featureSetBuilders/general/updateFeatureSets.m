%% updateFeatureSets
%   Reed Gurchiek, 2018
%   
%   updateFeatureSets runs all specific featureSetBuilders and thus updates
%   all feature sets they correspond to if new data is available
%
%--------------------------------------------------------------------------

%% initialization

clear
clc

% get msenseresearchgroup path
ok = questdlg('Select the msenseresearchgroup folder','MSENSE Path','OK',{'OK'});
if isempty(ok); error('Building stopped'); end
msense = uigetdir;

% get builders
builders = dir(fullfile(replace(fxndir('updateFeatureSets'),'general','specific'),'buildFeatureSet*.m'));

% for each
for b = 1:numel(builders)
    
    % update
    builder = builders(b).name(1:end-2);
    fprintf('\n-Calling: %s\n\n',builder);
    builder = str2func(builder);
    builder(msense);
    
end