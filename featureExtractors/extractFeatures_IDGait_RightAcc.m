function [ feat, featNames ] = extractFeatures_IDGait_RightAcc(rightAcc,sf,windowSize,overlap)
%Reed Gurchiek, 2018
%   extract features from left and right anterior thigh accelerometer data.
%
%---------------------------INPUTS-----------------------------------------
%
%   rightAcc, leftAcc:
%       3xn accelerometer data from left and right anterior thigh in thigh
%       frame. acceleration along segment long axis should be row 3
%
%   sf:
%       sampling frequency
%
%   windowSize:
%       size of the window in seconds
%
%   overlap:
%       0 < x < 1.  percent overlap of windows (e.g. 0 means data is not
%       shared between windows, 1 means 1 sample of data is not shared
%       between windows, 0.5 means 50% of data is shared between windows
%
%--------------------------OUTPUTS-----------------------------------------
%
%   feat:
%       pxm, m windows and p features characterizing each window
%
%   featNames:
%       px1, name of feature in corresponding column of feat
%
%--------------------------------------------------------------------------
%% extractFeatures_IDGait

% data loop
i = 1; %window index
feat = zeros(22,1);
featNames = {'rzrm' 'rhrm' 'rzrv' 'rhrv' 'rzrs' 'rhrs' 'rzrk' 'rhrk' 'rzlm'  'rzlv' 'rzls' 'rzlk'...
             'rzrp1' 'rzrw1' 'rzrp2' 'rzrw2' 'rhrp1' 'rhrw1' 'rhrp2' 'rhrw2'...
             'xrzrrhr' 'lagrzrrhr'};
         
ns = round(sf*windowSize); %samples per window
s1 = 1; %window starting index
s2 = ns; %window ending index
inc = ns - overlap*ns; %samples to slide between windows
if inc < 1; inc = 1; end
nsamp = length(rightAcc); %total samples
w = fftfreq(sf,ns); %frequencies for fft
w([1 round(ns/2):end]) = []; %remove negative and DC frequency
nw = length(w);
while s2 <= nsamp
    
    % window data and extract signals
    rz = rightAcc(3,s1:s2);
    rzl = filtmat_class(1/sf,1,rz',1,2)';
    rh = resultant(rightAcc(1:2,s1:s2));
    
    % raw z/horizontal acc mean, variance, skewness, kurtosis
    feat(1:2,i) = [mean(rz); mean(rh)]; 
    feat(3:4,i) = [var(rz); var(rh)];
    feat(5:6,i) = [skewness(rz); skewness(rh)];
    feat(7:8,i) = [kurtosis(rz); kurtosis(rh)];
    
    % low pass z acc mean, variance, skewness, kurtosis 
    feat(9:12,i) = [mean(rzl); var(rzl); skewness(rzl); kurtosis(rzl)];
    
    % raw z/horizontal in frequency domain (without negative or 0 freq)
    frz = fft(rz); frz(nw+1:end) = []; frz(1) = [];
    frh = fft(rh); frh(nw+1:end) = []; frh(1) = [];
    
    % raw z/horizontal top 2 most dominant frquencies and their percentage power
    [frzm,ifrzm] = extrema(abs(frz)); [frzm,isort] = sort(frzm,'descend'); frzm = frzm(1:3)./sum(abs(frz)); isort(4:end) = []; frzmw = w(ifrzm(isort));
    [frhm,ifrhm] = extrema(abs(frh)); [frhm,isort] = sort(frhm,'descend'); frhm = frhm(1:3)./sum(abs(frh)); isort(4:end) = []; frhmw = w(ifrhm(isort));
    
    % include in feature set
    feat(13:16,i) = [frzm(1); frzmw(1); frzm(2); frzmw(2)]; 
    feat(17:20,i) = [frhm(1); frhmw(1); frhm(2); frhmw(2)];

    % cross correlation
    [x,lag] = xcorr(rz,rh);
    [feat(21,i),imax] = max(x); feat(22,i) = abs(lag(imax));
    
    % next window
    s1 = s1 + inc;
    s2 = s2 + inc;
    i = i + 1;

end