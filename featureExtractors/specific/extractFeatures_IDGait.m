function [ feat, featNames ] = extractFeatures_IDGait(rightAcc,leftAcc,sf,windowSize,overlap)
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
feat = zeros(48,1);
featNames = {'rzrm' 'lzrm' 'rhrm' 'lhrm' 'rzrv' 'lzrv' 'rhrv' 'lhrv' 'rzrs' 'lzrs' 'rhrs' 'lhrs' 'rzrk' 'lzrk' 'rhrk' 'lhrk' 'rzlm' 'lzlm' 'rzlv' 'lzlv' 'rzls' 'lzls' 'rzlk' 'lzlk'...
             'rzrp1' 'rzrw1' 'rzrp2' 'rzrw2' 'lzrp1' 'lzrw1' 'lzrp2' 'lzrw2' 'rhrp1' 'rhrw1' 'rhrp2' 'rhrw2' 'lhrp1' 'lhrw1' 'lhrp2' 'lhrw2'...
             'xrzrlzr' 'lagrzrlzr' 'xrzrlzl' 'lagrzrlzl' 'xrzllzr' 'lagrzllzr' 'xrzllzl' 'lagrzllzl'};
ns = round(sf*windowSize); %samples per window
s1 = 1; %window starting index
s2 = ns; %window ending index
inc = ns - overlap*ns; %samples to slide between windows
if inc < 1; inc = 1; end
nsamp = min([length(rightAcc) length(leftAcc)]); %total samples
w = fftfreq(sf,ns); %frequencies for fft
w([1 round(ns/2):end]) = []; %remove negative and DC frequency
nw = length(w);
while s2 <= nsamp
    
    % window data and extract signals
    rz = rightAcc(3,s1:s2);
    lz = leftAcc(3,s1:s2);
    rzl = filtmat_class(1/sf,1,rz',1,2)';
    lzl = filtmat_class(1/sf,1,lz',1,2)';
    rh = resultant(rightAcc(1:2,s1:s2));
    lh = resultant(leftAcc(1:2,s1:s2));
    
    % raw z/horizontal acc mean, variance, skewness, kurtosis
    feat(1:4,i) = [mean(rz); mean(lz); mean(rh); mean(lh)]; 
    feat(5:8,i) = [var(rz); var(lz); var(rh); var(lh)];
    feat(9:12,i) = [skewness(rz); skewness(lz); skewness(rh); skewness(lh)];
    feat(13:16,i) = [kurtosis(rz); kurtosis(lz); kurtosis(rh); kurtosis(lh)];
    
    % low pass z acc mean, variance, skewness, kurtosis 
    feat(17:18,i) = [mean(rzl); mean(lzl)]; 
    feat(19:20,i) = [var(rzl); var(lzl)];
    feat(21:22,i) = [skewness(rzl); skewness(lzl)];
    feat(23:24,i) = [kurtosis(rzl); kurtosis(lzl)];
    
    % raw z/horizontal in frequency domain (without negative or 0 freq)
    frz = fft(rz); frz(nw+1:end) = []; frz(1) = [];
    flz = fft(lz); flz(nw+1:end) = []; flz(1) = [];
    frh = fft(rh); frh(nw+1:end) = []; frh(1) = [];
    flh = fft(lh); flh(nw+1:end) = []; flh(1) = [];
    
    % raw z/horizontal top 3 most dominant frquencies and their percentage power
    [frzm,ifrzm] = extrema(abs(frz)); [frzm,isort] = sort(frzm,'descend'); frzm = frzm(1:3)./sum(abs(frz)); isort(4:end) = []; frzmw = w(ifrzm(isort));
    [flzm,iflzm] = extrema(abs(flz)); [flzm,isort] = sort(flzm,'descend'); flzm = flzm(1:3)./sum(abs(flz)); isort(4:end) = []; flzmw = w(iflzm(isort));
    [frhm,ifrhm] = extrema(abs(frh)); [frhm,isort] = sort(frhm,'descend'); frhm = frhm(1:3)./sum(abs(frh)); isort(4:end) = []; frhmw = w(ifrhm(isort));
    [flhm,iflhm] = extrema(abs(flh)); [flhm,isort] = sort(flhm,'descend'); flhm = flhm(1:3)./sum(abs(flh)); isort(4:end) = []; flhmw = w(iflhm(isort));
    
    % include in feature set
    feat(25:28,i) = [frzm(1); frzmw(1); frzm(2); frzmw(2)]; 
    feat(29:32,i) = [flzm(1); flzmw(1); frzm(2); flzmw(2)];
    feat(33:36,i) = [frhm(1); frhmw(1); frhm(2); frhmw(2)];
    feat(37:40,i) = [flhm(1); flhmw(1); flhm(2); flhmw(2)];

    % cross correlation
    [x,lag] = xcorr(rz,lz);
    [feat(41,i),imax] = max(x); feat(42,i) = abs(lag(imax));
    [x,lag] = xcorr(rz,lzl);
    [feat(43,i),imax] = max(x); feat(44,i) = abs(lag(imax));
    [x,lag] = xcorr(rzl,lz);
    [feat(45,i),imax] = max(x); feat(46,i) = abs(lag(imax));
    [x,lag] = xcorr(rzl,lzl);
    [feat(47,i),imax] = max(x); feat(48,i) = abs(lag(imax));
    
    % next window
    s1 = s1 + inc;
    s2 = s2 + inc;
    i = i + 1;

end