function [ model ] = specifyModel(model)
%Reed Gurchiek, 2019
%   specify model details (fitting params) for use in
%   initializeBinaryClassifier
%
%----------------------------------INPUTS----------------------------------
%
%   model:
%       model struct, one field with model name: model.name = 'name';
%
%---------------------------------OUTPUTS----------------------------------
%
%   model:
%       model struct, includes parameters to be used during fitting
%
%--------------------------------------------------------------------------
%% specifyModel

% SVM
if strcmp(model.name,'Support Vector Machine')

    % fitter function
    model.fitter = 'fitcsvm';
    model.fitterParameters.OptimizeHyperparameters = {};
    
    % set defaults
    model.fitterParameters.BoxConstraint = 1;
    model.fitterParameters.KernelScale = 1;
    model.fitterParameters.KernelFunction = 'linear';
    model.fitterParameters.Solver = 'ISDA';
    
    % parameters
    parameters = {'BoxConstraint (1)' 'KernelScale (1)' 'KernelFunction (linear)' 'Solver (ISDA)'};
    modify = listdlg('ListString',parameters,'PromptString','Select the parameters to modify (or else Use Defaults):','SelectionMode','multiple',...
                     'OKString','Select','CancelString','Use Defaults','ListSize',[300 200]);
    % if modifications
    if ~isempty(modify)
        
        % for each modification
        for imod = 1:length(modify)
            
            parameter = parameters{modify(imod)};
    
            % box constraint
            if strcmp(parameter,'BoxConstraint (1)')
                
                % optimize or manual
                model.fitterParameters.BoxConstraint = questdlg('Select BoxConstraint (Misclassification Cost)','BoxConstraint','Manual Input','Optimize','Manual Input');
                if strcmp(model.fitterParameters.BoxConstraint,'Manual Input')
                    model.fitterParameters.BoxConstraint = str2double(inputdlg('Input BoxConstraint','BoxConstraint',[1 50],{'1'}));
                else
                    model.fitterParameters = rmfield(model.fitterParameters,'BoxConstraint');
                    model.fitterParameters.OptimizeHyperparameters = {'BoxConstraint'};
                end
               
            % kernel scale
            elseif strcmp(parameter,'KernelScale (1)')
                
                % optimize or manual
                model.fitterParameters.KernelScale = questdlg('Select KernelScale','KernelScale','Manual Input','Optimize','Manual Input');
                if strcmp(model.fitterParameters.KernelScale,'Manual Input')
                    model.fitterParameters.KernelScale = str2double(inputdlg('Input KernelScale','KernelScale',[1 50],{'1'}));
                else
                    model.fitterParameters = rmfield(model.fitterParameters,'KernelScale');
                    model.fitterParameters.OptimizeHyperparameters{length(model.fitterParameters.OptimizeHyperparameters)+1} = 'KernelScale';
                end
              
            % kernel function
            elseif strcmp(parameter,'KernelFunction (linear)')
                
                % get kernel
                kernels = {'linear','gaussian','polynomial'};
                model.fitterParameters.KernelFunction = listdlg('ListString',{'linear','gaussian','polynomial'},'PromptString','Select the kernel:','SelectionMode','single','ListSize',[200 100]);
                model.fitterParameters.KernelFunction = kernels{model.fitterParameters.KernelFunction};
                
                % if polynomial, get order
                if strcmp(model.fitterParameters.KernelFunction,'polynomial')
                    model.fitterParameters.PolynomialOrder = str2double(inputdlg('Input polynomial order','Order',[1 50],{'3'}));
                end
                
            % solver
            elseif strcmp(parameter,'Solver (ISDA)')
                
                % get solver
                model.fitterParameters.Solver = questdlg('Select solver:','Solver','Iterative Single Data Algorithm (ISDA)','Sequential Minimal Optimization (SMO)','Iterative Single Data Algorithm (ISDA)');
                if strcmp(model.fitterParameters.Solver,'Iterative Single Data Algorithm (ISDA)')
                    model.fitterParameters.Solver = 'ISDA';
                else
                    model.fitterParameters.Solver = 'SMO';
                end
                
            end
            
        end
    
    end

% K NEAREST NEIGHBORS
elseif strcmp(model.name,'K-Nearest Neighbors')
    
    % ensemble or single model
    ens = questdlg('Ensemble of weak learners or single model?','Ensemble','Ensemble','Single Model','Single Model');
    
    % ensemble
    if strcmp(ens,'Ensemble')
        model.fitter = 'fitcensemble';
        ens = 1;
    % single model
    else
        model.fitter = 'fitcknn';
        ens = 0;
    end
    
    % set defaults
    model.fitterParameters.OptimizeHyperparameters = {};
    model.fitterParameters.NumNeighbors = 1;
    model.fitterParameters.BreakTies = 'smallest';
    model.fitterParameters.Distance = 'euclidean';
    model.fitterParameters.DistanceWeight = 'equal';
    model.fitterParameters.BucketSize = 50;
    
    % parameters
    parameters = {'NumNeighbors (1)' 'BreakTies (smallest)' 'Distance (euclidean)' 'DistanceWeight (equal)','BucketSize (50)'};
    modify = listdlg('ListString',parameters,'PromptString','Select the parameters to modify (or else Use Defaults):','SelectionMode','multiple',...
                     'OKString','Select','CancelString','Use Defaults','ListSize',[300 200]);
    % if modifications
    if ~isempty(modify)
        
        % for each modification
        for imod = 1:length(modify)
            
            parameter = parameters{modify(imod)};
    
            % num neighbors
            if strcmp(parameter,'NumNeighbors (1)')
                
                if ~ens
                    % optimize or manual
                    model.fitterParameters.NumNeighbors = questdlg('Select number of neighbors','NumNeighbors','Manual Input','Optimize','Manual Input');
                else
                    model.fitterParameters.NumNeighbors = 'Manual Input';
                end
                
                if strcmp(model.fitterParameters.NumNeighbors,'Manual Input')
                    model.fitterParameters.NumNeighbors = str2double(inputdlg('Input NumNeighbors','NumNeighbors',[1 50],{'1'}));
                else
                    model.fitterParameters = rmfield(model.fitterParameters,'NumNeighbors');
                    model.fitterParameters.OptimizeHyperparameters = {'NumNeighbors'};
                end
                
            % break ties
            elseif strcmp(parameter,'BreakTies (smallest)')
                
                % select break tie criteria
                model.fitterParameters.BreakTies = questdlg('Select BreakTies criteria','BreakTies','smallest','nearest','random','smallest');
                
            % distance metric
            elseif strcmp(parameter,'Distance (euclidean)')
                
                % get distance metric
                distances = {'euclidean','cityblock','chebychev','correlation','cosine','hamming','jaccard','mahalanobis','minkowski','seuclidean','spearman'};
                model.fitterParameters.Distance = listdlg('ListString',distances,'PromptString','Select the distance metric:','SelectionMode','single','ListSize',[200 200]);
                model.fitterParameters.Distance = distances{model.fitterParameters.Distance};
                
                % if minkowski, get exponent
                if strcmp(model.fitterParameters.Distance,'minkowski')
                    model.fitterParameters.Exponent = str2double(inputdlg('Input minkowski exponent','Exponent',[1 50],{'3'}));
                end
                
            % distance weighting
            elseif strcmp(parameter,'DistanceWeight (equal)')
                
                % get distance weight
                model.fitterParameters.DistanceWeight = questdlg('Select distance weight:','DistanceWeight','equal','inverse','squaredinverse','equal');
                
            % bucket size
            elseif strcmp(parameter,'BucketSize (50)')
                
                % get bucket size
                model.fitterParameters.BucketSize = str2double(inputdlg('Input bucket size','Bucket Size',[1 50],{'50'}));
                
            end
            
        end
    
    end
    
    % create template if ensemble
    if ens
        fp0 = model.fitterParameters;
        model = rmfield(model,'fitterParameters');
        if ~strcmp(fp0.Distance,'minkowski')
            model.fitterParameters.Learners = templateKNN('NumNeighbors',fp0.NumNeighbors,'BreakTies',fp0.BreakTies,'Distance',fp0.Distance,'DistanceWeight',fp0.DistanceWeight,'BucketSize',fp0.BucketSize);
        else
            model.fitterParameters.Learners = templateKNN('NumNeighbors',fp0.NumNeighbors,'BreakTies',fp0.BreakTies,'Distance',fp0.Distance,'Exponent',fp0.Exponent,'DistanceWeight',fp0.DistanceWeight,'BucketSize',fp0.BucketSize);
        end
        model.fitterParameters.OptimizeHyperparameters = {};
        
        % method has to be subspace
        model.fitterParameters.Method = 'Subspace';
        
        % Resample?
        model.fitterParameters.Resample = questdlg('Resample in training?','Resample','on','off','off');
        
        if strcmp(model.fitterParameters.Resample,'on')
        
            % sample with replacement?
            model.fitterParameters.Replace = questdlg('Sample with replacement?','Replace','on','off','on');

            % fraction to resample
            model.fitterParameters.FResample = str2double(inputdlg('Input fraction to resample (between 0 and 1):','Resample Fraction',[1 50],{'1'}));
            
        end
        
        % N learners
        nlearners = inputdlg('Input the number of weak learners to use (default = 100)','N Learners',[1 100],{'100'});
        model.fitterParameters.NumLearningCycles = str2double(nlearners{1});
    end

% NAIVE BAYES
elseif strcmp(model.name,'Naive Bayes')

    % fitter function
    model.fitter = 'fitcnb';
    model.fitterParameters.OptimizeHyperparameters = {};
    
    % set defaults
    model.fitterParameters.DistributionNames = 'normal';
    
    % parameters
    parameters = {'normal (Gaussian)','kernel (select kernel smoother)','mn (multinomial)','mvmn (multivariate multinormal)'};
    modify = listdlg('ListString',parameters,'PromptString','Select the distribution model (MATLAB default is normal):','SelectionMode','single',...
                     'OKString','Select','CancelString','Use Default','ListSize',[300 200]);
                 
    % get name
    if modify == 3
        
        model.fitterParameters.DistributionNames = 'mn';
        
    elseif modify == 4
        
        model.fitterParameters.DistributionNames = 'mvmn';
     
    % if kernel
    elseif modify == 2
        
        % get kernel type
        kernelTypes = {'normal','box','epanechnikov','triangle'};
        ikernelType = listdlg('ListString',kernelTypes,'PromptString','Select the kernel smoother (MATLAB default is normal):','SelectionMode','single',...
                         'OKString','Select','CancelString','Use Default','ListSize',[300 200]);
        model.fitterParameters.Kernel = kernelTypes{ikernelType};
        
    end

% DECISION TREE
elseif strcmp(model.name,'Decision Tree')

    % fitter function
    model.fitter = 'fitctree';
    model.fitterParameters.OptimizeHyperparameters = {};
    
    % set defaults
    model.fitterParameters.MergeLeaves = 'on';
    model.fitterParameters.MaxNumSplits = 1;
    model.fitterParameters.SplitCriterion = 'gdi';
    model.fitterParameters.MinLeafSize = 1;
    
    % parameters
    parameters = {'MergeLeaves (on)','MaxNumSplits (1)','SplitCriterion (gdi)','MinimumLeafSize (1)'};
    modify = listdlg('ListString',parameters,'PromptString','Select the parameters to modify (or else Use Defaults):','SelectionMode','multiple',...
                     'OKString','Select','CancelString','Use Defaults','ListSize',[300 200]);
                 
    % if modifications
    if ~isempty(modify)
        
        % for each modification
        for imod = 1:length(modify)
            
            parameter = parameters{modify(imod)};
    
            % merge leaves
            if strcmp(parameter,'MergeLeaves (on)')
                
                % on or off
                model.fitterParameters.MergeLeaves = questdlg('Merge Leaves?','Merge','on','off','on');
                
            % MaxNumSplits
            elseif strcmp(parameter,'MaxNumSplits (1)')
                
                % N splits
                model.fitterParameters.MaxNumSplits = questdlg('Manually set or optimize maximum number of splits?','N Splits','Manual Input','Optimize','Manual Input');
                if strcmp(model.fitterParameters.MaxNumSplits,'Manual Input')
                    model.fitterParameters.MaxNumSplits = str2double(inputdlg('Input MaxNumSplits','N Splits',[1 50],{'1'}));
                else
                    model.fitterParameters = rmfield(model.fitterParameters,'MaxNumSplits');
                    model.fitterParameters.OptimizeHyperparameters = {'MaxNumSplits'};
                end
                
            % split criterion
            elseif strcmp(parameter,'SplitCriterion (gdi)')
                
                % get split criterion
                model.fitterParameters.SplitCriterion = questdlg('Select split criterion','Split Criterion','Gini''s diversity index (gdi)','Twoing rule (twoing)','Maximum deviance reduction (deviance)','Gini''s diversity index (gdi)');
                if model.fitterParameters.SplitCriterion(1) == 'G'; model.fitterParameters.SplitCriterion = 'gdi';
                elseif model.fitterParameters.SplitCriterion(1) == 'T'; model.fitterParameters.SplitCriterion = 'twoing';
                else; model.fitterParameters.SplitCriterion = 'deviance';
                end
                
            % min leaf size
            elseif strcmp(parameter,'MinimumLeafSize (1)')
                
                % optimize or manual
                model.fitterParameters.MinLeafSize = questdlg('Optimize or manually input minimum leaf size?','MinLeafSize','Manual Input','Optimize','Manual Input');
                if strcmp(model.fitterParameters.MinLeafSize,'Manual Input')
                    model.fitterParameters.MinLeafSize = str2double(inputdlg('Input MinLeafSize','MinLeafSize',[1 50],{'1'}));
                else
                    model.fitterParameters = rmfield(model.fitterParameters,'MinLeafSize');
                    model.fitterParameters.OptimizeHyperparameters{length(model.fitterParameters.OptimizeHyperparameters)+1} = 'MinLeafSize';
                end
                
            end
            
        end
    
    end

% LINEAR
elseif strcmp(model.name,'Linear')

    % fitter function
    model.fitter = 'fitclinear';
    model.fitterParameters.OptimizeHyperparameters = {};
    
    % set defaults
    model.fitterParameters.Learner = 'svm';
    model.fitterParameters.Regularization = 'ridge';
    model.fitterParameters.Lambda = 'auto';
    
    % parameters
    parameters = {'Learner (svm)','Regularization (ridge)','Lambda (auto)'};
    modify = listdlg('ListString',parameters,'PromptString','Select the parameters to modify (or else Use Defaults):','SelectionMode','multiple',...
                     'OKString','Select','CancelString','Use Defaults','ListSize',[300 200]);
    % if modifications
    if ~isempty(modify)
        
        % for each modification
        for imod = 1:length(modify)
            
            parameter = parameters{modify(imod)};
    
            % learner
            if strcmp(parameter,'Learner (svm)')
                
                % svm or logistic
                model.fitterParameters.Learner = questdlg('Select learner (loss function)','Learner','svm','logistic','svm');
                
            % regularization
            elseif strcmp(parameter,'Regularization (ridge)')
                
                % ridge or lasso
                model.fitterParameters.Regularization = questdlg('Select regularization (penalty)','Regularization','ridge','lasso','ridge');
                
                % solvers
                solvers = {'sgd (stochastic gradient descent)','asgd (average stochastic gradient descent)','dual (dual SGD for SVM)','lbfgs (limited-memory BFGS)','sparsa (Sparse Reconstruction by Separable Approximation)'};
                if model.fitterParameters.Regularization(1) == 'r'
                    solvers(end) = [];
                    if model.fitterParameters.Learner(1) == 'l'
                        solvers(3) = [];
                    end
                elseif model.fitterParameters.Regularization(1) == 'l'
                    solvers([3 4]) = [];
                end
                isolver = listdlg('ListString',solvers,'PromptString','Select the minimization technique:','SelectionMode','single',...
                                    'OKString','Select','CancelString','Use Defaults','ListSize',[300 200]);
                % if not use defaults
                if ~isempty(isolver)
                    ispace = strfind(solvers{isolver},' ');
                    model.fitterParameters.Solver = solvers{isolver}(1:ispace(1)-1);
                end
                
            % lambda
            elseif strcmp(parameter,'Lambda (auto)')
                
                % optimize or auto
                model.fitterParameters.Lambda = questdlg('Optimize or automatically set Lambda?','Lambda','optimize','auto','auto');
                if model.fitterParameters.Lambda(1) == 'o'
                    model.fitterParameters = rmfield(model.fitterParameters,'Lambda');
                    model.fitterParameters.OptimizeHyperparameters = {'Lambda'};
                end
                
            end
            
        end
    
    end
    
% DISCRIMINANT ANALYSIS
elseif strcmp(model.name,'Discriminant Analysis')
    
     % fitter function
    model.fitter = 'fitcdiscr';
    model.fitterParameters.OptimizeHyperparameters = {};
    
    % get discriminator
    discrimTypes = {'linear','diagLinear','pseudoLinear','quadratic','diagquadratic','pseudoquadratic'};
    idiscrim = listdlg('ListString',discrimTypes,'PromptString','Select the discriminant type (default is linear):','SelectionMode','multiple',...
                     'OKString','Select','CancelString','Use Default','ListSize',[300 200]);
                 
    model.fitterParameters.DiscrimType = discrimTypes{idiscrim};
    
    % set default delta
    model.fitterParameters.Delta = 0;
    
    % if linear model then give option to modify
    if idiscrim <= 3
        model.fitterParameters.Delta = questdlg('Manually set, optimize, or use default (0) linear coefficient threshold (delta)?','Delta','Manual Input','Optimize','Default (0)','Default (0)');
        if strcmp(model.fitterParameters.Delta,'Manual Input')
                model.fitterParameters.Delta = str2double(inputdlg('Input delta (non-negative scalar):','Delta',[1 50],{'0'}));
        elseif strcmp(model.fitterParameters.Delta,'Optimize')
                model.fitterParameters = rmfield(model.fitterParameters,'Delta');
                model.fitterParameters.OptimizeHyperparameters = {'Delta'};
        end
    end
    
    % auto or optimize gamma
    discrGamma = questdlg('Optimize or default gamma (amount of regularization)?','Gamma','Default','Optimize','Default');
    if discrGamma(1) == 'O'
        model.fitterParameters.OptimizeHyperparameters{length(model.fitterParameters.OptimizeHyperparameters)+1} = 'Gamma';
    end
    
% RANDOM FOREST
elseif strcmp(model.name,'Random Forest')
    
     % fitter function
    model.fitter = 'fitcensemble';
    model.fitterParameters.OptimizeHyperparameters = {};
    
    % set params
    model.fitterParameters.Method = 'Bag';
    model.fitterParameters.Learners = 'tree';
    
    % sample with replacement?
    model.fitterParameters.Replace = questdlg('Sample with replacement?','Replace','on','off','on');
    
    % fraction to resample
    model.fitterParameters.FResample = str2double(inputdlg('Input fraction to resample (between 0 and 1):','Resample Fraction',[1 50],{'1'}));
    
    % N learners
    nlearners = inputdlg('Input the number of weak learners to use (default = 100)','N Learners',[1 100],{'100'});
    model.fitterParameters.NumLearningCycles = str2double(nlearners{1});

end
        
% if optimizing
if ~isempty(model.fitterParameters.OptimizeHyperparameters)

    % defaults
    model.fitterParameters.HyperparameterOptimizationOptions.ShowPlots = 0;
    model.fitterParameters.HyperparameterOptimizationOptions.Verbose = 0;
    model.fitterParameters.HyperparameterOptimizationOptions.MaxObjectiveEvaluations = 30;


    % plot display
    answer = questdlg('Show plots during optimization?','Plots','Yes','No','No');
    if strcmp(answer,'Yes')
        model.fitterParameters.HyperparameterOptimizationOptions.ShowPlots = 1;
    end

    % verbosity
    answer = questdlg('Print iterative display during optimization?','Verbosity','Yes','No','No');
    if strcmp(answer,'Yes')
        model.fitterParameters.HyperparameterOptimizationOptions.Verbose = 1;
    end

    % time or iteration stopping criteria
    stopCriteria = questdlg('Select optimization stopping criteria','Stop Criteria','Max Number of Iterations','Max Time','Max Number of Iterations');
    if strcmp(stopCriteria,'Max Number of Iterations')
        model.fitterParameters.HyperparameterOptimizationOptions.MaxObjectiveEvaluations = str2double(inputdlg('Max number of iterations:','Max Iterations',[1 50],{'30'}));
        model.fitterParameters.HyperparameterOptimizationOptions.MaxTime = Inf;
    else
        model.fitterParameters.HyperparameterOptimizationOptions.MaxTime = str2double(inputdlg('Maximum time to spend optimizing (seconds):','Max Time',[1 50],{'Inf'}));
        model.fitterParameters.HyperparameterOptimizationOptions.MaxObjectiveEvaluations = Inf;
    end

else
    
    model.fitterParameters.OptimizeHyperparameters = 'none';
    
end

end