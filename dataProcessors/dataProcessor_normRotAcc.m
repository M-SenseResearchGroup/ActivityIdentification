function [ processedData ] = dataProcessor_normRotAcc(processorInfo)
%Reed Gurchiek, 2018
%   
%
%----------------------------------INPUTS----------------------------------
%
%   processorInfo:
%       struct.  Fields:
%       (1) data.(dataName(x)): 3xn, accelerometer data during static trial
%               -fieldnames should match fieldnames in calibration
%       (2) calibration.(dataName(x)): struct with fields:
%               (2a) gravityOffset: if perfect sensor then would be 1
%               (2b) Rs2b: matrix to rotate sensor to body such that row 3
%                           is z-axis data
%---------------------------------OUTPUTS----------------------------------
%
%   processedData:
%       -struct.  Same fields as processorInfo.data, but processed
%
%--------------------------------------------------------------------------
%% normRotAcc

% get dataNames
dataNames = fieldnames(processorInfo.data);

% for each sensor
for k = 1:length(dataNames)
    
    % normalize and rotate
    processedData.(dataNames{k}) = processorInfo.calibration.(dataNames{k}).Rs2b*processorInfo.data.(dataNames{k})/processorInfo.calibration.(dataNames{k}).gravityOffset;

end