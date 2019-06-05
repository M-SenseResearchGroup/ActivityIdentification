function [ f, name] = stdfeat(x,sf)
%Reed Gurchiek, 2019
%   extract standard features from signal x. For power spectral density
%   estimate (Welch's method), uses a window of half the signal length
%
%----------------------------------INPUTS----------------------------------
%
%   x:
%       1xn signal in
%
%   sf:
%       sampling frequency
%
%---------------------------------OUTPUTS----------------------------------
%
%   f:
%       nx1, extracted features vector
%
%   name:
%       nx1 cell array, feature names
%
%% stdfeat

% if x isempty then return names only
get = 1;
if isempty(x)
    f = []; 
    get = 0;
else

    % error check
    [nrow,ncol] = size(x);
    if nrow ~= 1 && ncol ~= 1; error('Signal x must be 1D (1xn or nx1).'); end
    
end

% median, 25th/75th quantiles, interquartile range
if get; f = [median(x); quantile(x,0.25); quantile(x,0.75); iqr(x)]; end
name = {'median'; 'quantile25'; 'quantile75'; 'iqr'};

% mean, 2nd, 3rd, 4th central moments
if get; f = vertcat(f,[mean(x); var(x); skewness(x); kurtosis(x)]); end
name = vertcat(name,{'mean'; 'variance'; 'skewness'; 'kurtosis'});

% RMS, absolute mean, max, min, range
if get; f = vertcat(f,[rms(x); mean(abs(x)); max(x); min(x); max(x) - min(x)]); end
name = vertcat(name,{'rms';'absoluteMean';'max';'min';'range'});

% spectral entropy
if get; f = vertcat(f,pentropy(x,sf,'Instantaneous',0)); end
name = vertcat(name,{'spectralEntropy'});

% median frequency
if get; f = vertcat(f,medfreq(x,sf)); end
name = vertcat(name,{'medianFrequency'});

% estimate power spectral density of raw signal
if get; [p,w] = pwelch(x,rectwin(round(length(x)/2)),[],4096,sf); end

% total signal power (autocorrelation at 0 lag), magnitude of low freq (< 0.25 Hz) power & percent power
if get; f = vertcat(f,trapz(w,p),trapz(w(w <= 0.25),p(w <= 0.25)),trapz(w(w <= 0.25),p(w <= 0.25))/trapz(w,p)); end
name = vertcat(name,{'totalSignalPower';'lowFrequencyPower';'lowFrequencyPercentPower'});

% estimate power spectral density of unbiased signal
if get
    [p,w] = pwelch(x-mean(x),rectwin(round(length(x)/2)),[],4096,sf);
    p(w <= 0.25) = [];
    w(w <= 0.25) = [];
end

% get dominant frequency and corresponding magnitude
if get
    [dom,idom] = max(p);
    f = vertcat(f, [w(idom); dom]);
end
name = vertcat(name,{'psdDominantFrequency';'psdMagnitudeAtDominantFrequency'});

% frequency below which contains 25, 50, 75, and 95% of power and total unbiased signal power
if get
    pp = cumtrapz(w,p);
    tp = pp(end);
    pp = pp./tp;
    w25 = w(pp <= 0.25);
    w50 = w(pp <= 0.50);
    w75 = w(pp <= 0.75);
    w95 = w(pp <= 0.95);
    f = vertcat(f,w25(end),w50(end),w75(end),w95(end),tp);
end
name = vertcat(name,{'frequencyBelowContains25';'frequencyBelowContains50';'frequencyBelowContains75';'frequencyBelowContains95';'totalUnbiasedSignalPower'});

% median, 25th/75th quantiles, interquartile range of PSD below w95
if get; f = vertcat(f,[median(p(w <= w95(end))); quantile(p(w <= w95(end)),0.25); quantile(p(w <= w95(end)),0.75); iqr(p(w <= w95(end)))]); end
name = vertcat(name,{'psdBelow95Median'; 'psdBelow95Quantile25'; 'psdBelow95Quantile75'; 'psdBelow95InterquartileRange'});

% mean, 2nd, 3rd, 4th central moments of PSD below w95
if get; f = vertcat(f,mean(p(w <= w95(end))),var(p(w <= w95(end))),skewness(p(w <= w95(end))),kurtosis(p(w <= w95(end)))); end
name = vertcat(name,{'psdBelow95Mean'; 'psdBelow95Variance'; 'psdBelow95Skewness';'psdBelow95Kurtosis'});

end