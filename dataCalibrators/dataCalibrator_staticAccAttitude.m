function [ calibration ] = dataCalibrator_staticAccAttitude(calibratorInfo)
%Reed Gurchiek, 2018
%   takes accelerometer data from static trial and returns the gravity
%   offset and constant dcm to rotate acc data to segment.  Third axis is
%   aligned with segment long axis which is assumed vertical during trial.
%   other two axes (transverse plane) are not necessarily aligned with x =
%   forward and y = lateral.  Only the transverse plan magnitude can be
%   known
%
%----------------------------------INPUTS----------------------------------
%
%   calibratorInfo:
%       struct.  Fields:
%       (1) data: 3xn, accelerometer data during static trial
%       (2) nSamples: number of samples to use for calibration
%
%---------------------------------OUTPUTS----------------------------------
%
%   calibration:
%       calibration.gravityOffset: scalar to normalize acc data
%                  .Rs2b: 3x3 matrix to express sensor referenced acc in 
%                         segment frame (see description)
%
%--------------------------------REFERENCES--------------------------------
%
%   LastName et al. (year), Title.
%
%--------------------------------------------------------------------------
%% staticAcc2Seg

% vars
a = calibratorInfo.data;
n = calibratorInfo.nSamples;

% get most still period
i = staticndx(a,n);

% get gravity offset
g = mean(vecnorm(a(:,i(1):i(2))));

% get sensor to body dcm
r = getrot(normalize(mean(a,2)),[0;0;1],'dcm');

% save
calibration.gravityOffset = g;
calibration.Rs2b = r;

end