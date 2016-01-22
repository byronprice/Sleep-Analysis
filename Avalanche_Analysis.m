function [Av_Size,mean_Size,Av_IEI,mean_IEI,Time] = Avalanche_Analysis(Av_Data,window,maxThreshold)
%Avalanche_Analysis.m
%   Take "Av_Data" matrix from Ca_im_DataStream output and count 
%    avalanches using different definitions and thresholds.
% Created: 2015/12/15 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/12/15
%  by: Byron Price
% 
% INPUT: Av_Data - a single matrix size videos by frames
%
%        window - for avalanche definition
%        maxThreshold - also for avalanche definition, the code will
%        iterate from a threshold of 1 up to maxThreshold, where each
%        threshold signifies the minimum number of ROIs that must be active
%        in a single frame for an avalanche to be counted 
%
% OUTPUT: Av_Size - a matrix with count data for each of the possible
%         thresholds and avalanche sizes ... if Av_Size(1,1) = 50, then
%         there were 50 instances of an avalanche of size 1
%
%         mean_Av_Size - the mean size of an avalanche at each threshold
%         ... this isn't so interesting as we force the minimum up with
%         each increase in the treshold
%
%         Av_IEI - a matrix with count data for each of the possible
%         avalanche thresholds and inter-event intervals
%
%         mean_Av_IEI - the mean inter-event interval between avalanches of
%         size at least "threshold"
%
%         Time - vector of times corresponding to the bins in Av_IEI,
%         try plotting  plot(Time, Av_IEI(1,:))  for a plot of the number
%         of instances of an inter-event interval at each of the times in
%         Time assuming only those bins with greater than or equal to 1
%         active ROI is counted as an avalanche ... Av_IEI(2,:) would give
%         you the same plot but with IEI's between those avalanches greater
%         than or equal to a size of two ROIs
%
% AVALANCHE DEFINITION
%   If window = 0, then the definition is ... the summed activity across
%   all ROIs in a single time bin, if that sum is greater than the value of 
%   threshold.  So, if 10 ROIs are active in time bin 20,
%   then that's an avalanche of size 10. Each bin is an island and
%   contiguous bins are never counted as part of the same avalanche.  As
%   threshold increases, the minimum number of active ROIs per bin
%   increases.
%      Example: For a threshold of 1, a sequence from Av_Data like 
%      [0,0,1,3,4,0,5,6,9,0,0,0] would have 4 inter-event intervals of 
%      0 seconds (contiguous bins active), 1 IEI of 0.0333 seconds (1 bin 
%      separation), 1 IEI of 0.0666 seconds and at the end 1 of 0.0999 seconds.
%      That would be for threshold 1 ... Av_IEI(1,:). For a threshold of
%      two, stored in Av_IEI(2,:), you would have only 5 avalanches and
%      so the IEIs would change accordingly (the opening 3 bins would be a 
%      single IEI.
%
%   If window = 1, then the definition is ... the summed activity across
%   all ROIs and time bins for which at least "threshold" ROIs are active.
%   So, for threshold = 1, contiguous bins with at least one ROI active are
%   counted as a single avalanche.
%      Example: Given the sequence above [0,0,1,3,4,0,5,6,9,0,0,0] and a
%      threshold of 1, the inter-event inetervals would be the same, except
%      there would be no IEIs of 0 seconds because IEIs of 0 seconds mean
%      two consecutive bins were active and both bins are counted in the
%      same avalanche.  The sizes, however, would change, with 1 avalanche
%      of size 1+3+4=8 ROIs and another avalanche of size 5+6+9=20 ROIs.
%
%   As 'window' increases, the number of bins with activity below threshold
%   between bins with activity above or equal to threshold grows.  Thus,
%   for window = 2 and the sequence above [0,0,1,3,4,0,5,6,9,0,0,0], there
%   would now be only one avalanche of size 28.  The IEIs would be
%   basically the same but with 1 IEI of 2 bins (0.0666 seconds) and one of
%   3 bins (0.0999 seconds).  
%
Fs = 30;
numVideos = size(Av_Data,1);
T = size(Av_Data,2);
Time = 0:1/Fs:T/Fs;
maxROIs = 500;   

Av_Size = zeros(maxThreshold,maxROIs); % avalanche size 
%Av_Size(:,:,2) = ones(maxThreshold,1)*(1:maxROIs);
mean_Size = zeros(maxThreshold,1);

% Av_Len = zeros(maxThreshold,length(Time));
% mean_Len = zeros(maxThreshold,1);

Av_IEI = zeros(maxThreshold,length(Time)); % inter-event interval (seconds), or time 
                % between avalanches
mean_IEI = zeros(maxThreshold,1);



for threshold = 1:maxThreshold
    Sizes = [];
    IEIs = [];
    for i=1:numVideos
            icount = 0;
            scount = 0;
            %lcount = 0;
            if window == 0
                for t=1:T
                    if Av_Data(i,t) >= threshold || t == 1
                        if Av_Data(i,t) >= threshold
                            Av_Size(threshold,Av_Data(i,t),1) = Av_Size(threshold,Av_Data(i,t),1)+1;
                            Sizes = [Sizes,Av_Data(i,t)];
                        end
                        count = 0;
                        k = t+1;
                        while Av_Data(i,k) < threshold && k < T
                            k = k+1;
                            count = count+1;
                        end
                        if count < T-2
                            Av_IEI(threshold,count+1) = Av_IEI(threshold,count+1)+1;
                            IEIs = [IEIs,Time(count+1)];
                        end
                    end
                end

            elseif window ~= 0
                for t=1:T
                    if t >= (T-(window-1))
                        if t == T
                            if scount > 0
                                scount = sum(squeeze(Av_Data(i,t-scount:end)));
                                Av_Size(threshold,scount,1) = Av_Size(threshold,scount,1)+1;
                                Sizes = [Sizes,scount];
                            end
                            if icount < T-1
                                icount = icount+(window-1);
                                Av_IEI(threshold,icount+1,1) = Av_IEI(threshold,icount+1,1)+1;
                                IEIs = [IEIs,Time(icount+1)];
                            end
                        else
                            len = length(Av_Data(i,t:end));
                            if (squeeze(Av_Data(i,t:end)) < threshold) == ones(1,len)
                                icount = icount+1;
                                if scount > 0
                                    scount = sum(squeeze(Av_Data(i,t-scount:t-1)));
                                    Av_Size(threshold,scount,1) = Av_Size(threshold,scount,1)+1;
                                    Sizes = [Sizes,scount];
                                    scount = 0;
                                end
                            else
                                scount = scount+1;
                                if icount > 0
                                    icount = icount+(window-1);
                                    Av_IEI(threshold,icount+1,1) = Av_IEI(threshold,icount+1,1)+1;
                                    IEIs = [IEIs,Time(icount+1)];
                                    icount = 0;
                                end
                            end
                        end
                    else
                        if (squeeze(Av_Data(i,t:(t+(window-1)))) < threshold) == ones(1,window)
                            icount = icount+1;
                            if scount > 0
                                scount = sum(squeeze(Av_Data(i,t-scount:t-1)));
                                Av_Size(threshold,scount,1) = Av_Size(threshold,scount,1)+1;
                                Sizes = [Sizes,scount];
                                scount = 0;
                            end
                        else
                            scount = scount+1;
                            if icount > 0
                                icount = icount+(window-1);
                                Av_IEI(threshold,icount+1,1) = Av_IEI(threshold,icount+1,1)+1;
                                IEIs = [IEIs,Time(icount+1)];
                                icount = 0;
                            end
                        end
                    end
                end
            end
    end
    mean_Size(threshold,1) = mean(Sizes);
    mean_IEI(threshold,1) = mean(IEIs);
end

end


