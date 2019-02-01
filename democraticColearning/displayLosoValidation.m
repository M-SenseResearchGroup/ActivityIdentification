% load DCL_losoValidation_17Dec2018
loop = validation.loop;
clear validation

% for each set
task = {'locomotion' 'walk'};
performance.locomotion.true = [];
performance.walk.true = [];

for set = 1:23
    
    for l = 1:21
        
        for t = 1:2
    
            % criteria, metric, threshold
            c = loop(set).newDataCriteria;
            m = loop(set).newDataMetric;
            th = loop(set).newDataThreshold;
            ths = valfname(num2str(th));
            performance.(task{t}).(c).(m).(ths).democratic.acc = zeros(1,21);
            performance.(task{t}).(c).(m).(ths).democratic.spec = zeros(1,21);
            performance.(task{t}).(c).(m).(ths).democratic.iter = 1:21;
            performance.(task{t}).(c).(m).(ths).democratic.threshold = th*ones(1,21);
            performance.(task{t}).(c).(m).(ths).loop(l).democratic.labels = [];
            
            for mod = 1:5
                performance.(task{t}).(c).(m).(ths).loop(l).model(mod).labels = [];
                performance.(task{t}).(c).(m).(ths).model(mod).threshold = th*ones(1,21);
                performance.(task{t}).(c).(m).(ths).model(mod).iter = 1:21;
                performance.(task{t}).(c).(m).(ths).model(mod).acc = zeros(1,21);
                performance.(task{t}).(c).(m).(ths).model(mod).spec = zeros(1,21);
            end
                
        
        end
    
    end
    
end

for set = 1:23
    disp(set)
    % for each subject out
    for s = 1:9
        
        % if not skipping
        if ~loop(set).loso(s).skip
        
            % for each task
            for t = 1:2

                % get predicted labels
                if set == 1
                    performance.(task{t}).true = horzcat(performance.(task{t}).true,loop(1).loso(s).(task{t}).test.labels);
                end

                % for each loop
                for l = 1:21
                    
                    % criteria, metric, threshold
                    c = loop(set).newDataCriteria;
                    m = loop(set).newDataMetric;
                    th = loop(set).newDataThreshold;
                    ths = valfname(num2str(th));

                    % save point
                    performance.(task{t}).(c).(m).(ths).loop(l).democratic.labels = horzcat(performance.(task{t}).(c).(m).(ths).loop(l).democratic.labels,loop(set).loso(s).(task{t}).loop(l).democratic.predictedLabels);
                    
                    % for each model
                    for mod = 1:5
                        
                        % save point
                        performance.(task{t}).(c).(m).(ths).loop(l).model(mod).labels = horzcat(performance.(task{t}).(c).(m).(ths).loop(l).model(mod).labels,loop(set).loso(s).(task{t}).loop(l).model(mod).predictedLabels);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end     
    
% criteria
criteria = {'confidence' 'percent'};

% metric
metric = {'democraticConfidence' 'disagreementStrength'};

% for each task
for t = 1:2
    
    % for each criteria
    for c = 1:2
        
        % for each metric
        for m = 1:2
            
            % get thresholds
            ths = fieldnames(performance.(task{t}).(criteria{c}).(metric{m}));
            figure
            
            % for each threshold
            for th = 1:length(ths)
                
                % for each loop
                for l = 1:21
                    
                    % evaluate
                    [acc,~,spec] = evalbc(performance.(task{t}).true,performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).loop(l).democratic.labels);
                    
                    % store point
                    performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).democratic.acc(l) = acc;
                    performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).democratic.spec(l) = spec;
                    
                    % for each model
                    for mod = 1:5
                        
                        % evaluate
                        [acc,~,spec] = evalbc(performance.(task{t}).true,performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).loop(l).model(mod).labels);

                        % store point
                        performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).model(mod).acc(l) = acc;
                        performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).model(mod).spec(l) = spec;
                        
                    end
                    
                end
                
                
            
                % plot
                plot3(performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).democratic.threshold,...
                      performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).democratic.iter,...
                      performance.(task{t}).(criteria{c}).(metric{m}).(ths{th}).democratic.acc);
                hold on
                
            end
            
            % label
            title(strcat('Task: ',task{t},', Criteria: ',criteria{c},', Metric: ',metric{m}))
            xlabel('Threshold')
            ylabel('Iteration')
            zlabel('Accuracy')
            
        end
        
    end
    
end

% for each task
wacc = zeros(1,6);
wacc(1) = performance.walk.percent.democraticConfidence.x01.democratic.acc(1);
wspec = zeros(1,6);
wspec(1) = performance.walk.percent.democraticConfidence.x01.democratic.spec(1);
lacc = zeros(1,6);
lacc(1) = performance.locomotion.percent.democraticConfidence.x01.democratic.acc(1);
lspec = zeros(1,6);
lspec(1) = performance.locomotion.percent.democraticConfidence.x01.democratic.spec(1);
for m = 1:5
    wacc(m+1) = performance.walk.percent.democraticConfidence.x01.model(m).acc(1);
    wspec(m+1) = performance.walk.percent.democraticConfidence.x01.model(m).spec(1);
    lacc(m+1) = performance.locomotion.percent.democraticConfidence.x01.model(m).acc(1);
    lspec(m+1) = performance.locomotion.percent.democraticConfidence.x01.model(m).spec(1);
end
figure
scatter(1:6,wacc,'ro')
hold on
scatter(1:6,wspec,'b*')
scatter(1:6,lacc,'g+')
scatter(1:6,lspec,'k.')
legend('Walk Acc','Walk Spec','Loco Acc','Loco Spec');
