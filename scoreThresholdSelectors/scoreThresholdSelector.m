function [ scoreThreshold ] = scoreThresholdSelector(criteria,X,Y,T)
%Reed Gurchiek, 2018
%   description
%
%----------------------------------INPUTS----------------------------------
%
%   criteria:
%       string, name of criterion, acceptable values:
%           (1) youden, Youden's index
%           (2) distance, distance to corner
%           (3) PLR, positive likelihood ratio
%           (4) NLR, negative likelihood ratio
%           (5) DOR, diagnostic odds ratio
%
%   X,Y,T:
%       outputs from perfcurve, x axis and y axis values on ROC curve and
%       corresponding thresholds
%
%---------------------------------OUTPUTS----------------------------------
%
%   scoreThreshold:
%       double, between 0 and 1, optimal score according to selected
%       criterion above which predictions should be labeled true
%
%--------------------------------------------------------------------------
%% scoreThresholdSelector

% distance
if strcmp(criteria,'cornerDistance')
    
    % Compute distance from (0,1) for each point
    scoreThreshold = sqrt(sum([X, 1-Y].^2,2));

    % get min distance location
    [~,imin]=min(scoreThreshold);

    % Extract threshold of min distance
    scoreThreshold = T(imin);
    
% youden    
elseif strcmp(criteria,'youden')
    
    % youden value for each point
    scoreThreshold = Y-X;

    % Find max youdin value and its location
    [~,imax] = max(scoreThreshold);

    % Extract threshold of max youdin value
    scoreThreshold = T(imax);
    
% otherwise DOR, PLR, or NLR
elseif strcmp(criteria,'DOR')
    
    % get false negative rate and specificity
    T(X == 0 | X == 1) = [];
    Y(X == 0 | X == 1) = [];
    X(X == 0 | X == 1) = [];
    fnr = 1-Y;
    spe = 1-X;

    % get positive and negative likelihood ratio and diagnostic odds ratio
    plr = Y./X;
    nlr = fnr./spe;
    dor = plr./nlr;
    
    % select which to use
    if strcmp(criteria,'PLR'); scoreThreshold = plr;
    elseif strcmp(criteria,'NLR'); scoreThreshold = nlr;
    elseif strcmp(criteria,'DOR'); scoreThreshold = dor;
    end

    % get max
    [~,imax] = max(scoreThreshold);

    % get corresponding threshold
    scoreThreshold = T(imax);
    
end