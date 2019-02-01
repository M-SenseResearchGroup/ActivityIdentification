function [ feat, featNames, indices ] = extractFeatures_Acc2(extractorInfo)
%Reed Gurchiek, 2018
%   extract features from 2 accelerometers with acceleration(3,:) =
%   longitiudinal axis and acceleration(1:2,:) = horizontal plane
%
%---------------------------INPUTS-----------------------------------------
%
%   extractorInfo:
%       struct with following fields:
%
%   data:
%       struct containing a1 and a2 which are accelerometer data from 2 
%       sensors in a segment fixed frame. acceleration along segment long 
%       axis should be row 3.
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
%   reportStatus (optional):
%       if 1 then displays waitbar reporting completion status
%
%--------------------------OUTPUTS-----------------------------------------
%
%   feat:
%       pxm, m windows and p features characterizing each window
%
%   featNames:
%       px1, name of feature in corresponding column of feat
%
%   indices:
%       2xm, row 1 (row 2) column c contains starting (ending) index of 
%       window whose features are in column c of feat
%
%--------------------------------------------------------------------------
%% extractFeaturesAcc2

% initialization
i = 1; %window index
indices = zeros(2,1);
feat = zeros(106,1);
featNames = {'z1m' 'z2m' 'h1m' 'h2m' 'z1v' 'z2v' 'h1v' 'h2v' 'z1s' 'z2s' 'h1s' 'h2s' 'z1k' 'z2k' 'h1k' 'h2k' ...
             'z1lv' 'z2lv' 'h1lv' 'h2lv' 'z1ls' 'z2ls' 'h1ls' 'h2ls' 'z1lk' 'z2lk' 'h1lk' 'h2lk' ...
             'p1z1' 'p1z1w' 'p2z1' 'p2z1w' 'p3z1' 'p3z1w' 'p1z2' 'p1z2w' 'p2z2' 'p2z2w' 'p3z2' 'p3z2w' ...
             'p1h1' 'p1h1w' 'p2h1' 'p2h1w' 'p3h1' 'p3h1w' 'p1h2' 'p1h2w' 'p2h2' 'p2h2w' 'p3h2' 'p3h2w' ...
             'p1z1l' 'p1z1lw' 'p2z1l' 'p2z1lw' 'p3z1l' 'p3z1lw' 'p1z2l' 'p1z2lw' 'p2z2l' 'p2z2lw' 'p3z2l' 'p3z2lw' ...
             'p1h1l' 'p1h1lw' 'p2h1l' 'p2h1lw' 'p3h1l' 'p3h1lw' 'p1h2l' 'p1h2lw' 'p2h2l' 'p2h2lw' 'p3h2l' 'p3h2lw' ...
             'acor0z1' 'acor1z1' 'acor1z1w' 'acor0z2' 'acor1z2' 'acor1z2w' 'acor0h1' 'acor1h1' 'acor1h1w' 'acor0h2' 'acor1h2' 'acor1h2w' ...
             'acor0z1l' 'acor1z1l' 'acor1z1lw' 'acor0z2l' 'acor1z2l' 'acor1z2lw' 'acor0h1l' 'acor1h1l' 'acor1h1lw' 'acor0h2l' 'acor1h2l' 'acor1h2lw' ...
             'xcor0z1l_z2l' 'xcor1z1l_z2l' 'xcor1z1l_z2lw' 'xcor0h1l_h2l' 'xcor1h1l_h2l' 'xcor1h1l_h2lw'};

% if no arguments, just return featNames
if nargin == 0
    feat = [];
    indices = [];
else
    
    % get vars
    a1 = extractorInfo.data.a1;
    a2 = extractorInfo.data.a2;
    sf = extractorInfo.samplingFrequency;
    windowSize = extractorInfo.windowSize;
    overlap = extractorInfo.overlap;

    % report status?
    if isfield(extractorInfo,'reportStatus')
        reportStatus = extractorInfo.reportStatus;
        if ~isa(reportStatus,'numeric'); reportStatus = 0; end
    else
        reportStatus = 0;
    end
    
    % report?
    if reportStatus; wb = waitbar(0,'Progress: 0%','Name','Extracting Features'); end
    
    ns = round(sf*windowSize); %samples per window
    s1 = 1; %window starting index
    s2 = ns; %window ending index
    inc = ns - round(overlap*ns); %samples to slide between windows
    if inc < 1; inc = 1; end
    nsamp = min([length(a1) length(a2)]); %total samples

    % if not enough samples for single window
    if nsamp < ns
        feat = [];
    else
        while s2 <= nsamp
            
            % report?
            if reportStatus; waitbar(s2/nsamp,wb,sprintf('Progress: %3.2f%%',s2/nsamp*100)); end
                
            % save indices
            indices(:,i) = [s1;s2];

            % window data and extract signals
            z1 = a1(3,s1:s2);
            z2 = a2(3,s1:s2);
            z1l = bwfilt(z1,1,sf,'lowpass',4);
            z2l = bwfilt(z2,1,sf,'lowpass',4);
            h1 = vecnorm(a1(1:2,s1:s2));
            h2 = vecnorm(a2(1:2,s1:s2));
            h1l = bwfilt(h1,1,sf,'lowpass',4);
            h2l = bwfilt(h2,1,sf,'lowpass',4);

            % raw z/horizontal acc mean, variance, skewness, kurtosis
            feat(1:4,i) = [mean(z1); mean(z2); mean(h1); mean(h2)];
            feat(5:8,i) = [var(z1); var(z2); var(h1); var(h2)];
            feat(9:12,i) = [skewness(z1); skewness(z2); skewness(h1); skewness(h2)];
            feat(13:16,i) = [kurtosis(z1); kurtosis(z2); kurtosis(h1); kurtosis(h2)];

            % low pass z/horizontal acc variance, skewness, kurtosis (no mean, should = raw)
            feat(17:20,i) = [var(z1l); var(z2l); var(h1l); var(h2l)];
            feat(21:24,i) = [skewness(z1l); skewness(z2l); skewness(h1l); skewness(h2l)];
            feat(25:28,i) = [kurtosis(z1l); kurtosis(z2l); kurtosis(h1l); kurtosis(h2l)];

            % unbiased raw z/horizontal spectral power density estimate
            [pz1,fz1] = pwelch(z1-mean(z1),rectwin(ns),[],4096,sf);
            [pz2,fz2] = pwelch(z2-mean(z2),rectwin(ns),[],4096,sf);
            [ph1,fh1] = pwelch(h1-mean(h1),rectwin(ns),[],4096,sf);
            [ph2,fh2] = pwelch(h2-mean(h2),rectwin(ns),[],4096,sf);

            % unbiased low pass z/horizontal spectral power density estimate
            [pz1l,fz1l] = pwelch(z1l-mean(z1l),rectwin(ns),[],4096,sf);
            [pz2l,fz2l] = pwelch(z2l-mean(z2l),rectwin(ns),[],4096,sf);
            [ph1l,fh1l] = pwelch(h1l-mean(h1l),rectwin(ns),[],4096,sf);
            [ph2l,fh2l] = pwelch(h2l-mean(h2l),rectwin(ns),[],4096,sf);

            % most dominant frequencies for raw and the power density there
            [pz1m,ipz1m] = extrema(pz1); [pz1m,isort] = sort(pz1m,'descend'); pz1mw = fz1(ipz1m(isort));
            [pz2m,ipz2m] = extrema(pz2); [pz2m,isort] = sort(pz2m,'descend'); pz2mw = fz2(ipz2m(isort));
            [ph1m,iph1m] = extrema(ph1); [ph1m,isort] = sort(ph1m,'descend'); ph1mw = fh1(iph1m(isort));
            [ph2m,iph2m] = extrema(ph2); [ph2m,isort] = sort(ph2m,'descend'); ph2mw = fh2(iph2m(isort));

            % most dominant frequencies for lowpass and the power density there
            [pz1lm,ipz1lm] = extrema(pz1l); [pz1lm,isort] = sort(pz1lm,'descend'); pz1lmw = fz1l(ipz1lm(isort));
            [pz2lm,ipz2lm] = extrema(pz2l); [pz2lm,isort] = sort(pz2lm,'descend'); pz2lmw = fz2l(ipz2lm(isort));
            [ph1lm,iph1lm] = extrema(ph1l); [ph1lm,isort] = sort(ph1lm,'descend'); ph1lmw = fh1l(iph1lm(isort));
            [ph2lm,iph2lm] = extrema(ph2l); [ph2lm,isort] = sort(ph2lm,'descend'); ph2lmw = fh2l(iph2lm(isort));

            % include top 3 in feature set for raw
            feat(29:34,i) = [pz1m(1); pz1mw(1); pz1m(2); pz1mw(2); pz1m(3); pz1mw(3)]; 
            feat(35:40,i) = [pz2m(1); pz2mw(1); pz2m(2); pz2mw(2); pz2m(3); pz2mw(3)];
            feat(41:46,i) = [ph1m(1); ph1mw(1); ph1m(2); ph1mw(2); ph1m(3); ph1mw(3)];
            feat(47:52,i) = [ph2m(1); ph2mw(1); ph2m(2); ph2mw(2); ph2m(3); ph2mw(3)];

            % include top 3 in feature set for lowpass
            if length(pz1lm) < 3; pz1lm(end+1:3) = 0; pz1lmw(end+1:3) = 0; end
            if length(pz2lm) < 3; pz2lm(end+1:3) = 0; pz2lmw(end+1:3) = 0; end
            if length(ph1lm) < 3; ph1lm(end+1:3) = 0; ph1lmw(end+1:3) = 0; end
            if length(ph2lm) < 3; ph2lm(end+1:3) = 0; ph2lmw(end+1:3) = 0; end
            
            feat(53:58,i) = [pz1lm(1); pz1lmw(1); pz1lm(2); pz1lmw(2); pz1lm(3); pz1lmw(3)]; 
            feat(59:64,i) = [pz2lm(1); pz2lmw(1); pz2lm(2); pz2lmw(2); pz2lm(3); pz2lmw(3)];
            feat(65:70,i) = [ph1lm(1); ph1lmw(1); ph1lm(2); ph1lmw(2); ph1lm(3); ph1lmw(3)];
            feat(71:76,i) = [ph2lm(1); ph2lmw(1); ph2lm(2); ph2lmw(2); ph2lm(3); ph2lmw(3)];

            % autocorrelation of raw signals
            [acorz1,lz1] = xcorr(z1,'unbiased'); lz1 = lz1/sf; acorz1(1:(end+1)/2-1) = []; lz1(1:(end+1)/2-1) = [];
            [acorz2,lz2] = xcorr(z2,'unbiased'); lz2 = lz2/sf; acorz2(1:(end+1)/2-1) = []; lz2(1:(end+1)/2-1) = [];
            [acorh1,lh1] = xcorr(h1,'unbiased'); lh1 = lh1/sf; acorh1(1:(end+1)/2-1) = []; lh1(1:(end+1)/2-1) = [];
            [acorh2,lh2] = xcorr(h2,'unbiased'); lh2 = lh2/sf; acorh2(1:(end+1)/2-1) = []; lh2(1:(end+1)/2-1) = [];

            % autocorrelation of lowpass signals
            [acorz1l,lz1l] = xcorr(z1l,'unbiased'); lz1l = lz1l/sf; acorz1l(1:(end+1)/2-1) = []; lz1l(1:(end+1)/2-1) = [];
            [acorz2l,lz2l] = xcorr(z2l,'unbiased'); lz2l = lz2l/sf; acorz2l(1:(end+1)/2-1) = []; lz2l(1:(end+1)/2-1) = [];
            [acorh1l,lh1l] = xcorr(h1l,'unbiased'); lh1l = lh1l/sf; acorh1l(1:(end+1)/2-1) = []; lh1l(1:(end+1)/2-1) = [];
            [acorh2l,lh2l] = xcorr(h2l,'unbiased'); lh2l = lh2l/sf; acorh2l(1:(end+1)/2-1) = []; lh2l(1:(end+1)/2-1) = [];

            % get height at 0 lag and height/lag at first peak for raw
            % require first peak be at least 0.3 s past 0
            minSample = floor(0.3*sf);
            [p,ip] = extrema(acorz1(minSample:end)); ip = ip + minSample - 1;
            feat(77:79,i) = [acorz1(1); p(1); lz1(ip(1))];
            [p,ip] = extrema(acorz2(minSample:end)); ip = ip + minSample - 1;
            feat(80:82,i) = [acorz2(1); p(1); lz2(ip(1))];
            [p,ip] = extrema(acorh1(minSample:end)); ip = ip + minSample - 1;
            feat(83:85,i) = [acorh1(1); p(1); lh1(ip(1))];
            [p,ip] = extrema(acorh2(minSample:end)); ip = ip + minSample - 1;
            feat(86:88,i) = [acorh2(1); p(1); lh2(ip(1))];

            % get height at 0 lag and height/lag at first peak for lowpass
            [p,ip] = extrema(acorz1l(minSample:end)); ip = ip + minSample - 1;
            feat(89:91,i) = [acorz1l(1); p(1); lz1l(ip(1))];
            [p,ip] = extrema(acorz2l(minSample:end)); ip = ip + minSample - 1;
            feat(92:94,i) = [acorz2l(1); p(1); lz2l(ip(1))];
            [p,ip] = extrema(acorh1l(minSample:end)); ip = ip + minSample - 1;
            feat(95:97,i) = [acorh1l(1); p(1); lh1l(ip(1))];
            [p,ip] = extrema(acorh2l(minSample:end)); ip = ip + minSample - 1;
            feat(98:100,i) = [acorh2l(1); p(1); lh2l(ip(1))];

            % cross correlation of lowpass sigs z1l-z2l, h2l-h2l
            [xz,lxz] = xcorr(z1l,z2l,'unbiased'); lxz = lxz/sf; xz(1:(end+1)/2-1) = []; lxz(1:(end+1)/2-1) = [];
            [xh,lxh] = xcorr(h1l,h2l,'unbiased'); lxh = lxh/sf; xh(1:(end+1)/2-1) = []; lxh(1:(end+1)/2-1) = [];

            % get correlation at 0 lag and correlation/lag at first peak
            % require first peak be at least 0.3 s past 0
            minSample = floor(0.3*sf);
            [p,ip] = extrema(xz(minSample:end)); ip = ip + minSample - 1;
            feat(101:103,i) = [xz(1); p(1); lxz(ip(1))];
            [p,ip] = extrema(xh(minSample:end)); ip = ip + minSample - 1;
            feat(104:106,i) = [xh(1); p(1); lxh(ip(1))];

            % next window
            s1 = s1 + inc;
            s2 = s2 + inc;
            i = i + 1;
        end
        
        % report?
        if reportStatus; close(wb); end
    end
end


end