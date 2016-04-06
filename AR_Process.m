function [B,Pvals,BurstRate,threshVec] = AR_Process(P_Av_Data,Fs,maxThreshold)
%AR_Process.m
% Fit a generalized-linear model with a log link function to the burst 
%  data (which is preprocessed to consider burst activity as a point process,
%  i.e. as a series of ones and zeros.  The firing rate, mu, is modelled as 
%  a function of the activity during 20 seconds of lags up to the present 
%  like ln(mu) = X*b, where X is a matrix of the previous 20 seconds of 
%  activity.  Thus, the log odds of firing a burst at the present moment 
%  depends on a linear combination of the previous 20 seconds of bursting 
%  activity. The output vector b(lag) represents the odds of firing in the 
%  present given a burst lag-1 bins ago.  So, the value of b(2) gives the 
%  log odds of firing a burst given a burst 1 time step ago and no bursts 
%  at any other time
% See Pillow et al. Nature 2008
% See also Gerhard et al. PLoS Comp Bio 2013
%INPUT: P_Av_Data - output of the Ca_im_DataStream code, i.e.
%        Processed_Data.Av_Data
%        This code assumes that, for a given night, the same number of 
%        frames were taken for each video, i.e. each video was recorded for
%        the same amount of time
%       Fs - sampling frequency, Hertz
%       maxThreshold - fraction from 0 to 1, the code will progressively 
%        threshold the burst data taking only those bursts with greater than
%        threshold fraction of the total number of active ROIs, and 
%        maxThreshold tells this threshold procedure when to stop
%OUTPUT: B - list of coefficients for the link function ln(mu) = X*b
%         Each row of B represents a burst threshold (see above), while
%         each column represents the number of lags
%          B(1,1) represents the log odds of firing a burst given no activity
%          in the previous 20 seconds and a burst threshold of only a single
%          ROI, B(1,2) is a lag of one frame, B(1,3) is a lag of two frames
%       Pvals - statistics on the fit of the model
%         equal to stats.p, for the p-value of the t test, which
%         evaluates the null hypothesis that the coefficient b(lag) is equal
%         to zero
%       BurstRate - the baseline rate of bursting at each threshold (Hertz) 
%       threshVec - vector of thresholds up to maxThreshold
%
%Created: 2016/02/19 at 24 Cummington, Boston
%   Byron Price
%Updated: 2016/04/06
%  By: Byron Price

threshVec = 0:0.1:maxThreshold;
%timeLag = 2;
%numLags = round(timeLag*Fs);
numLags = 600;
timeLag = numLags/Fs;

a =15;c=65;
phiVec = 20*pi:pi/2:31.5*pi;
numBases = length(phiVec);
BASIS = zeros(numLags,numBases);
count = 1;
t = 1:numLags;
for phi = phiVec
    for tt=1:numLags
        if a*log(t(tt)+c) > (phi-pi) && a*log(t(tt)+c) < (phi+pi)
            BASIS(tt,count) = 0.5.*cos(a*log(t(tt)+c)-phi)+0.5;
        else
            BASIS(tt,count) = 0;
        end
    end
%     plot(t,b)
%     hold on;
%     title(sprintf('phi=%d',phi));
%     pause(0.2)
    count = count+1;
end

numNights = size(P_Av_Data,2);
numVideos = zeros(numNights,1);
numFrames = zeros(numNights,1);
totalLength = 0;
for ii=1:numNights
    numVideos(ii) = size(P_Av_Data{ii},1);
    numFrames(ii) = size(P_Av_Data{ii},2);
    if numFrames(ii) > numLags
        totalLength = totalLength + numVideos(ii)*(numFrames(ii)-numLags);
    else
        display('Must decrease the variable timeLag')
        return;
    end
end
clear numVideos numFrames;

Time = linspace(1/Fs,timeLag,numLags);
B = zeros(length(threshVec),numBases+1);
Pvals = zeros(length(threshVec),numBases+1);
threshcount = 1;
for threshold=threshVec
    Y = zeros(totalLength,1);
    HIST = zeros(totalLength,numLags);
    count = 1;
    for zz=1:numNights
        Av_Data = P_Av_Data{zz};
        numVideos = size(Av_Data,1);
        numFrames = size(Av_Data,2);
        for jj=1:numVideos
            newData = Av_Data(jj,:)./max(Av_Data(jj,:));
            newData(newData<threshold) = 0; % clear all activity below the threshold
            newData(newData>0) = 1;
            if count == 1
                top = 1;
            else
                top = bottom+1;
            end
            bottom = top+length(newData((numLags+1):end))-1;
            Y(top:bottom,1) = newData((numLags+1):end);
            for kk=1:numLags
                HIST(top:bottom,kk) = newData((numLags-(kk-1)):(numFrames-kk));
            end
            count = count+1;
            clear newData;
        end
    end
    X = HIST*BASIS;
    [b,~,stats] = glmfit(X,Y,'poisson');
    B(threshcount,:) = b';
    Pvals(threshcount,:) = stats.p;
    threshcount = threshcount+1;
end
figure();plot(Time,BASIS*(B(1,2:end)'));xlabel('Time Lag (seconds)');ylabel('Probability of Bursting');
title('Probability a burst will occur time-lag seconds after a previous burst (Threshold for Bursting of 1 ROI');
figure();plot(1:numBases,Pvals(1,2:end));xlabel('Parameter');ylabel('P-value');
title('P-values for each of the model parameters');
figure();imagesc(Time,threshVec,(BASIS*(B(:,2:end)'))');xlabel('Time Lag (seconds)');ylabel('Threshold');
title('Probability of bursting (all thresholds)');
%figure();plot(Time,B(1,2:end));figure();plot(Time,B(end,2:end));
BurstRate = exp(B(:,1)).*Fs;
figure();plot(threshVec,BurstRate);ylabel('Baseline Burst Rate (Hz)');
xlabel('Threshold (fraction of total ROIs active)');
title('Burst rate as a function of burst threshold');
end

