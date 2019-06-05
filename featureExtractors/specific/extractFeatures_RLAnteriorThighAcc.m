function [ feat, featNames, indices ] = extractFeatures_RLAnteriorThighAcc(extractorInfo)
%Reed Gurchiek, 2019
%   extract features from 3-axis accel on R/L anterior thigh. Expects
%   the 3rd acc axis to be aligned with the segment longitudinal axis. Also
%   expects 2nd acc axis close to antero posterior axis.  This fxn
%   optimized for sagittal plane activity (e.g. walking)
%
%---------------------------INPUTS-----------------------------------------
%
%   extractorInfo:
%       struct with following fields:
%
%   data:
%       struct containing a1 and a2 which are accelerometer data from 2 
%       sensors in a segment fixed frame. acceleration along segment long 
%       axis should be row 3. a1 and a2 are 3xn
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
%       2xm, row 1 (row 2) & column c contains starting (ending) index of 
%       window whose features are in column c of feat
%
%--------------------------------------------------------------------------
%% extractFeatures_RLThighAcc

% first few features are mean and difference of standard features from
% same axis of signals from different sensors (see code for details)

% get feat names for standard features
[~,stdnames] = stdfeat([],[]);

% concatenate appropriately for description
featNames = vertcat(append2all(stdnames,'cranialCaudal_',0),append2all(stdnames,'cranialCaudalDifference_',0),append2all(stdnames,'anteriorPosterior_',0),append2all(stdnames,'anteriorPosteriorDifference_',0));
featNames = vertcat(featNames,{'cranialCaudal_PCAWeight';'cranialCaudalDifference_PCAWeight';...
                               'correlation_CC_contralateral';'correlation_AP_contralateral';...
                               'correlation_CC_AP_ipsilateral';'correlationDifference_CC_AP_ipsilateral';...
                               'correlation_CC_AP_contralateral';'correlationDifference_CC_AP_contralateral';...
                               'maxCrossCorrelation_CC_contralateral';'lagMaxCrossCorrelation_CC_contralateral';'maxCrossCorrelation_AP_contralateral';'lagMaxCrossCorrelation_AP_contralateral';...
                               'maxCrossCorrelation_CC_AP_ipsilateral';'maxCrossCorrelationDifference_CC_AP_ipsilateral';...
                               'lagMaxCrossCorrelation_CC_AP_ipsilateral';'lagMaxCrossCorrelationDifference_CC_AP_ipsilateral';...
                               'maxCrossCorrelation_CC_AP_contralateral';'maxCrossCorrelationDifference_CC_AP_contralateral';...
                               'lagMaxCrossCorrelation_CC_AP_contralateral';'lagMaxCrossCorrelationDifference_CC_AP_contralateral'});

% initialization
i = 1; %window index
indices = zeros(2,1);
feat = zeros(length(featNames),1);

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
        
        % low pass filter cutoff = 3 Hz chosen to capture step frequency
        lpfc = 3;
        
        % get signals
        z1 = a1(3,:);
        h1 = a1(1:2,:);
        z1l = bwfilt(z1,lpfc,sf,'low',4); 
        h1l = bwfilt(h1,lpfc,sf,'low',4);
        z2 = a2(3,:);
        h2 = a2(1:2,:);
        z2l = bwfilt(z2,lpfc,sf,'low',4);
        h2l = bwfilt(h2,lpfc,sf,'low',4);
        
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
            
            % estimate AP/ML acc as if sagittal plane activity
            pc1 = pca(ih1l');
            pc1(:,2) = [];
            ih1 = [ pc1(2) -pc1(1); pc1(1)  pc1(2)]*ih1; % AP is e2
            ih1l = [ pc1(2) -pc1(1); pc1(1)  pc1(2)]*ih1l;
            scal = sign(mean(ih1(2,:)));
            ih1(2,:) = scal*ih1(2,:); % make anterior direction positive
            ih1l(2,:) = scal*ih1l(2,:);
            
            % now sensor 2
            pc2 = pca(ih2l');
            pc2(:,2) = [];
            ih2 = [ pc2(2) -pc2(1); pc2(1)  pc2(2)]*ih2;
            ih2l = [ pc2(2) -pc2(1); pc2(1)  pc2(2)]*ih2l;
            scal = sign(mean(ih2(2,:)));
            ih2(2,:) = scal*ih2(2,:);
            ih2l(2,:) = scal*ih2l(2,:);
            
            % y is AP
            iy1 = ih1(2,:);
            iy1l = ih1l(2,:);
            iy2 = ih2(2,:);
            iy2l = ih2l(2,:);
            
            % get standard features for each signal
            z1f = stdfeat(iz1,sf);
            z2f = stdfeat(iz2,sf);
            y1f = stdfeat(iy1,sf);
            y2f = stdfeat(iy2,sf);
            
            % get mean and absolute difference between contralateral sigs
            mz = (z1f + z2f)/2;
            dz = abs(z1f - z2f);
            my = (y1f + y2f)/2;
            dy = abs(y1f - y2f);
            
            % concatenate
            ifeat = vertcat(mz,dz,my,dy);
            
            % weight of z axis of first principal component
            pcz1 = pca(a1(:,s1:s2)');
            pcz1 = pcz1(3,1);
            pcz2 = pca(a2(:,s1:s2)');
            pcz2 = pcz2(3,1);
            ifeat = vertcat(ifeat,(pcz1 + pcz2)/2,abs(pcz1 - pcz2));
            
            % correlation between low passed z1&z2, y1&y2
            ifeat = vertcat(ifeat,corr(iz1l',iz2l'),corr(iy1l',iy2l'));
            
            % mean/difference correlation between z&y of same leg
            ifeat = vertcat(ifeat,mean([corr(iz1l',iy1l') corr(iz2l',iy2l')]),abs(corr(iz1l',iy1l') - corr(iz2l',iy2l')));
            
            % mean/difference correlation between z&y of opposite leg
            ifeat = vertcat(ifeat,mean([corr(iz1l',iy2l') corr(iz2l',iy1l')]),abs(corr(iz1l',iy2l') - corr(iz2l',iy1l')));
            
            % max/lag xcorr between low passed, unbiased z1&z2
            [xc,lags] = xcorr(iz1l-mean(iz1l),iz2l-mean(iz2l),'biased');
            [xc,maxlag] = max(xc); maxlag = lags(maxlag)/sf;
            ifeat = vertcat(ifeat,xc,abs(maxlag));
            
            % max/lag xcorr between low passed y1&y2
            xc = xcorr(iy1l-mean(iy1l),iy2l-mean(iy2l),'biased');
            [xc,maxlag] = max(xc); maxlag = lags(maxlag)/sf;
            ifeat = vertcat(ifeat,xc,abs(maxlag));
            
            % mean/difference max/lag xcorr between z&y of same leg
            xc1 = xcorr(iz1l-mean(iz1l),iy1l-mean(iy1l),'biased');
            [xc1,il1] = max(xc1); il1 = abs(lags(il1)/sf);
            xc2 = xcorr(iz2l-mean(iz2l),iy2l-mean(iy2l),'biased');
            [xc2,il2] = max(xc2); il2 = abs(lags(il2)/sf);
            ifeat = vertcat(ifeat,(xc1+xc2)/2,abs(xc1-xc2),(il1 + il2)/2,abs(il1 - il2));
            
            % mean/difference max/lag xcorr between z&y of opposite legs
            xc1 = xcorr(iz1l-mean(iz1l),iy2l-mean(iy2l),'biased');
            [xc1,il1] = max(xc1); il1 = abs(lags(il1)/sf);
            xc2 = xcorr(iz2l-mean(iz2l),iy1l-mean(iy1l),'biased');
            [xc2,il2] = max(xc2); il2 = abs(lags(il2)/sf);
            ifeat = vertcat(ifeat,(xc1+xc2)/2,abs(xc1-xc2),(il1 + il2)/2,abs(il1 - il2));
            
            % add feature vector
            feat(:,i) = ifeat;
            clear ifeat

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