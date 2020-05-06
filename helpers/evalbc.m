function [ accuracy, sensitivity, specificity, precision, auc, roc, err ] = evalbc(trueLabels,predictedLabels,predictionConfidence,trueOriginalClasses)
%Reed Gurchiek, 2018
%   evalbc evaluates a binary classifier in terms of accuracy, sensitivity
%   (recall or true positive rate), specificity (true negative rate), and
%   precision (positive predictive value).
%
%---------------------------INPUTS-----------------------------------------
%
%   trueLabels:
%       1xn array of 1 or -1, 1 is true and -1 is false.  trueLabels
%       contains the true labels
%
%   predictedLabels:
%       1xn array of 1 or -1, labels predicted from binary classifier
%
%   predictionConfidence (optional):
%       1xn array of 0 to 1, confidence observation belongs to class 1. If
%       empty ( = []), then auc and roc will be empty
%
%   trueOriginalClasses (optional):
%       1xn cell array whose elements are strings describing the original
%       class name corresponding to the same element in trueLabel.  If
%       these are given then the error output will be populated.
%
%--------------------------OUTPUTS-----------------------------------------
%
%   accuracy:
%       0 to 1, percent of correct predictions
%
%   sensitivity:
%       0 to 1, percent of actual positives that are correctly identified
%
%   specificity:
%       0 to 1, percent of actual negatives that are correctly identified
%
%   precision:
%       0 to 1, percent of labeled positives that are actually positive
%
%   auc:
%       area under roc curve
%
%   roc:
%       3xn, row 1: 1-specificity, row 2: sensitivity. row 3: thresholds.
%       >> plot(roc(1:2,:)') for roc curve
%
%   error:
%       if trueOriginalClasses is given then this will be populated as
%       structure with a field for each unique trueOriginalClass name.
%       These fields are 3x2 cell arrays.
%           -Row 1 is a header file = {'falsePositives','falseNegatives','total'}
%           -Row 2 contains the number of errors
%           -Row 3 contains the percent of all errors accounted for
%           -Row 4 contains the percent of class specific errors accounted
%           for
%
%--------------------------------------------------------------------------
%% evalbc

% error check
nLabels = length(trueLabels);
if nLabels ~= length(predictedLabels)
    error('predictedLabels and trueLabels must be same length.')
end

% accuracy
correctIndices = trueLabels == predictedLabels;
accuracy = sum(correctIndices)/nLabels;

% sensitivity
truePositives = sum((predictedLabels == 1) & correctIndices);
falseNegatives = sum((predictedLabels == -1) & ~correctIndices);
sensitivity = truePositives/(truePositives + falseNegatives);

% specificity
trueNegatives = sum((predictedLabels == -1) & correctIndices);
falsePositives = sum((predictedLabels == 1) & ~correctIndices);
specificity = trueNegatives/(trueNegatives + falsePositives);

% precision
precision = truePositives/(truePositives + falsePositives);

% initialize err
err = cell(0);

% if confidences given
if nargin > 2
    if ~isempty(predictionConfidence)
        % prediction confidence from msenseClassify is the confidence the
        % predicted label belongs to that class. perfcurve expect the
        % scores. so subtract predictionConfidence from one where all
        % predicted labels are -1.
        scores = predictionConfidence;
        scores(predictedLabels == -1) = 1 - scores(predictedLabels == -1);
        [xroc,yroc,thresh,auc] = perfcurve(trueLabels',scores',1);
        roc = [xroc';yroc';thresh'];
    else
        auc = [];
        roc = [];
    end
else
    auc = [];
    roc = [];
end
        

% if original classes given
if nargin == 4
    
    % error check
    if length(trueOriginalClasses) ~= nLabels
        warning('trueOriginalClasses must have same number of elements as labels.  Skipping class based error analysis.')
    else
        
        % initialize error matrix
        class = unique(trueOriginalClasses);
        err = cell(length(class)+1,6);
        err(1,:) = {'ClassName','CorrectLabel','nObservations','nErrors','WithinClassPercentError','AllPercentError'};
        for c = 1:length(class)
            trueLabel = trueLabels(strcmp(class{c},trueOriginalClasses));
            err(c+1,:) = {class{c},trueLabel(1),length(trueLabel),0,0,0};
        end
        
        % get indices of errors
        errorIndices = find(~correctIndices);
        nerrors = length(errorIndices);
        
        % count class specific errors
        for i = 1:nerrors
            row = find(strcmp(trueOriginalClasses{errorIndices(i)},class));
            err{row+1,4} = err{row+1,4} + 1;
        end
        
        % get percentages and total
        for c = 1:length(class)
            err{c+1,5} = err{c+1,4}/err{c+1,3}*100;
            err{c+1,6} = err{c+1,4}/nerrors*100;
        end
        
    end
    
end

%% references

% Halilaj et al. (2018) Machine learning in human movement biomechanics:
% best practices, common pitfalls and new opportunities

end