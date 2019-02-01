function [ feat, featNames, indices ] = extractFeatures_RLThighAcc(extractorInfo)
%Reed Gurchiek, 2018
%   extract features from 2 accelerometers on right and left thigh
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
%% extractFeatures_RLThighAcc

% initialization
i = 1; %window index
indices = zeros(2,1);
feat = zeros(28,1);
featNames = {'mean_zCC', 'mean_yAP','variance_z','variance_y','variance_x','zVarianceDifference_half1_half2','yVarianceDifference_half1_half2','mean3Peaks_z''mean3Peaks_y',...
             'corr_z1z2','corr_y1y2','corr_zy_ipsilateral','corrDifference_z1y1_z2y2','corr_zy_contralateral','corrDifference_z1y2_z2y1','zPCAWeight','zPCAWeightDifference',...
             'dominantFrequency_z','zPercentPowerBelow5Hz','zDifferencePercentPowerBelow5Hz','zPercentPower5to10Hz','zDifferencePercentPower5to10Hz',...
             'dominantFrequency_y','yPercentPowerBelow5Hz','yDifferencePercentPowerBelow5Hz','yPercentPower5to10Hz','yDifferencePercentPower5to10Hz','differencePercentPowerBelow5Hz_magHalf1_magHalf2'};
             

% if no arguments, just return featNames
if nargin == 0
    feat = [];
    indices = [];
else
    
    % get vars
    dataNames = fieldnames(extractorInfo.data);
    a1 = extractorInfo.data.(dataNames{1});
    a2 = extractorInfo.data.(dataNames{2});
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
        
        % get signals
        z1 = a1(3,:);
        h1 = a1(1:2,:);
        z1l = bwfilt(z1,3,sf,'low',4);
        h1l = bwfilt(h1,3,sf,'low',4);
        z2 = a2(3,:);
        h2 = a2(1:2,:);
        z2l = bwfilt(z2,3,sf,'low',4);
        h2l = bwfilt(h2,3,sf,'low',4);
        
        while s2 <= nsamp
            
            % report?
            if reportStatus; waitbar(s2/nsamp,wb,sprintf('Progress: %3.2f%%',s2/nsamp*100)); end
                
            % save indices
            indices(:,i) = [s1;s2];

            % window data
            iz1 = z1(s1:s2);
            ih1 = h1(:,s1:s2);
            iz1l = z1l(s1:s2);
            ih1l = h1l(:,s1:s2);
            iz2 = z2(s1:s2);
            ih2 = h2(:,s1:s2);
            iz2l = z2l(s1:s2);
            ih2l = h2l(:,s1:s2);
            
            % estimate AP/ML acc as if walking
            pc1 = pca(ih1l');
            pc1(:,2) = [];
            ih1 = [ pc1(2) -pc1(1); pc1(1)  pc1(2)]*ih1; % AP is e2
            ih1l = [ pc1(2) -pc1(1); pc1(1)  pc1(2)]*ih1l;
            scal = sign(mean(ih1(2,:)));
            ih1(2,:) = scal*ih1(2,:); % make anterior direction positive
            ih1l(2,:) = scal*ih1l(2,:);
            pc2 = pca(ih2l');
            pc2(:,2) = [];
            ih2 = [ pc2(2) -pc2(1); pc2(1)  pc2(2)]*ih2; % AP is e2
            ih2l = [ pc2(2) -pc2(1); pc2(1)  pc2(2)]*ih2l;
            scal = sign(mean(ih2(2,:)));
            ih2(2,:) = scal*ih2(2,:);
            ih2l(2,:) = scal*ih2l(2,:);
            
            % y is AP
            iy1 = ih1(2,:);
            iy1l = ih1l(2,:);
            iy2 = ih2(2,:);
            iy2l = ih2l(2,:);
            
            % x is ML, but don't know direction so use magnitude
            ix1 = abs(ih1(1,:));
            ix2 = abs(ih2(1,:));
            
            % mean of z1 and z2 concatenated
            feat(1,i) = mean([iz1 iz2]);
            
            % mean of y1 and y2 concatenated
            feat(2,i) = mean([iy1 iy2]);
            
            % variance of z1 and z2 concatenated
            feat(3,i) = var([iz1 iz2]);
            
            % variance of y1 and y2 concatenated
            feat(4,i) = var([iy1 iy2]);
            
            % variance of x1 and x2 concatenated
            feat(5,i) = var([ix1 ix2]);
            
            % mean difference between variance of first half for z1 & z2
            nhalf = round(length(iz1)/2);
            feat(6,i) = (abs(var(iz1(1:nhalf))-var(iz1(nhalf+1:end))) + abs(var(iz2(1:nhalf))-var(iz2(nhalf+1:end))))/2;
            
            % mean difference between variance of first half for y1 & y2
            nhalf = round(length(iy1)/2);
            feat(7,i) = (abs(var(iy1(1:nhalf))-var(iy1(nhalf+1:end))) + abs(var(iy2(1:nhalf))-var(iy2(nhalf+1:end))))/2;
            
            % mean of 3 biggest peaks of z
            pks1 = findpeaks(iz1,'SortStr','descend');
            pks2 = findpeaks(iz2,'SortStr','descend');
            if isempty(pks1)
                pks1 = 0; 
            elseif length(pks1) < 3
                pks1 = mean(pks1); 
            else
                pks1 = mean(pks1(1:3));
            end
            if isempty(pks2)
                pks2 = 0; 
            elseif length(pks2) < 3
                pks2 = mean(pks2); 
            else
                pks2 = mean(pks2(1:3));
            end
            feat(8,i) = (pks1 + pks2)/2;
            
            % mean of 3 biggest peaks of y
            pks1 = findpeaks(iy1,'SortStr','descend');
            pks2 = findpeaks(iy2,'SortStr','descend');
            if isempty(pks1)
                pks1 = 0; 
            elseif length(pks1) < 3
                pks1 = mean(pks1); 
            else
                pks1 = mean(pks1(1:3));
            end
            if isempty(pks2)
                pks2 = 0; 
            elseif length(pks2) < 3
                pks2 = mean(pks2); 
            else
                pks2 = mean(pks2(1:3));
            end
            feat(9,i) = (pks1 + pks2)/2;
            
            % correlation between low passed z1 and z2
            feat(10,i) = corr(iz1l',iz2l');
            
            % correlation between low passed y1 and y2
            feat(11,i) = corr(iy1l',iy2l');
            
            % mean correlation between ipsilateral low passed y/z
            feat(12,i) = mean([corr(iz1l',iy1l') corr(iz2l',iy2l')]);
            
            % correlation difference between ipsilateral low passed y/z
            feat(13,i) = abs(corr(iz1l',iy1l') - corr(iz2l',iy2l'));
            
            % mean correlation between contralateral legs low passed y/z
            feat(14,i) = mean([corr(iz1l',iy2l') corr(iz2l',iy1l')]);
            
            % correlation difference between contralateral legs low passed y/z
            feat(15,i) = abs(corr(iz1l',iy2l') - corr(iz2l',iy1l'));
            
            % weight of z axis of first principal component
            pcz1 = pca(a1(:,s1:s2)');
            pcz1 = pcz1(3,1);
            pcz2 = pca(a2(:,s1:s2)');
            pcz2 = pcz2(3,1);
            feat(16,i) = (pcz1 + pcz2)/2;
            feat(17,i) = abs(pcz1 - pcz2);
            
            % estimate power spectral density of z1 and z2 using welch's
            % method (2 second window) (unbiased)
            [pz1,w1] = pwelch(iz1 - mean(iz1),rectwin(round(2*sf)),[],4096,sf);
            [pz2,w2] = pwelch(iz2 - mean(iz2),rectwin(round(2*sf)),[],4096,sf);
            pz1(w1 < 0.5) = [];
            w1(w1 < 0.5) = [];
            pz2(w2 < 0.5) = [];
            w2(w2 < 0.5) = [];
            
            % get mean frequency of most power for z1/z2
            [~,ipk1] = max(pz1);
            [~,ipk2] = max(pz2);
            feat(18,:) = (w1(ipk1)+w2(ipk2))/2;
            
            % mean/difference percent power below 5 of z1 and z2
            tpz1 = sum(pz1);
            tpz2 = sum(pz2);
            feat(19,i) = (sum(pz1(w1 <= 5))/tpz1 + sum(pz2(w2 <= 5))/tpz2)/2;
            feat(20,i) = abs(sum(pz1(w1 <= 5))/tpz1 - sum(pz2(w2 <= 5))/tpz2);
            
            % mean/difference percent power between 5 and 10 of z1 and z2
            feat(21,i) = (sum(pz1(w1 > 5 & w1 <= 10))/tpz1 + sum(pz2(w2 > 5 & w2 <= 10))/tpz2)/2;
            feat(22,i) = abs(sum(pz1(w1 > 5 & w1 <= 10))/tpz1 - sum(pz2(w2 > 5 & w2 <= 10))/tpz2);
            
            % estimate power spectral density of y1 and y2 using welch's
            % method (2 second window) (unbiased)
            [py1,w1] = pwelch(iy1 - mean(iy1),rectwin(round(2*sf)),[],4096,sf);
            [py2,w2] = pwelch(iy2 - mean(iy2),rectwin(round(2*sf)),[],4096,sf);
            py1(w1 < 0.5) = [];
            w1(w1 < 0.5) = [];
            py2(w2 < 0.5) = [];
            w2(w2 < 0.5) = [];
            
            % get mean frequency of most power for y1/y2
            [~,ipk1] = max(py1);
            [~,ipk2] = max(py2);
            feat(23,:) = (w1(ipk1)+w2(ipk2))/2;
            
            % mean/difference percent power between 0.5 and 5 of z1 and z2
            tpy1 = sum(py1);
            tpy2 = sum(py2);
            feat(24,i) = (sum(py1(w1 <= 5))/tpy1 + sum(py2(w2 <= 5))/tpy2)/2;
            feat(25,i) = abs(sum(py1(w1 <= 5))/tpy1 - sum(py2(w2 <= 5))/tpy2);
            
            % mean/difference percent power between 5 and 10 of z1 and z2
            feat(26,i) = (sum(py1(w1 > 5 & w1 <= 10))/tpy1 + sum(py2(w2 > 5 & w2 <= 10))/tpy2)/2;
            feat(27,i) = abs(sum(py1(w1 > 5 & w1 <= 10))/tpy1 - sum(py2(w2 > 5 & w2 <= 10))/tpy2);
            
            % estimate power spectral density of magnitude of first/second half
            % using welch's method (2 second window) (unbiased)
            mag11 = vecnorm([ix1;iy1;iz1]);
            mag12 = mag11(nhalf+1:end);
            mag11(nhalf+1:end) = [];
            mag21 = vecnorm([ix2;iy2;iz2]);
            mag22 = mag21(nhalf+1:end);
            mag21(nhalf+1:end) = [];
            [pm11,w1] = pwelch(mag11 - mean(mag11),rectwin(round(sf)),[],4096,sf);
            pm12 = pwelch(mag12 - mean(mag12),rectwin(round(sf)),[],4096,sf);
            [pm21,w2] = pwelch(mag21 - mean(mag21),rectwin(round(sf)),[],4096,sf);
            pm22 = pwelch(mag22 - mean(mag22),rectwin(round(sf)),[],4096,sf);
            pm11(w1 < 0.5) = [];
            pm12(w1 < 0.5) = [];
            w1(w1 < 0.5) = [];
            pm21(w2 < 0.5) = [];
            pm22(w2 < 0.5) = [];
            w2(w2 < 0.5) = [];
            
            % get difference between percent power between 0.5 and 5 hz
            % from first half and second half of window
            tpm11 = sum(pm11);
            tpm12 = sum(pm12);
            tpm21 = sum(pm21);
            tpm22 = sum(pm22);
            feat(28,i) = abs(mean([sum(pm11(w1 < 5))/tpm11 sum(pm21(w2 < 5))/tpm21]) - mean([sum(pm12(w1 < 5))/tpm12 sum(pm22(w2 < 5))/tpm22]));

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