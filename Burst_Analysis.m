function [Av_Size,Av_IEI,Time] = Burst_Analysis(P_Av_Data,window,maxThreshold,Fs)
%Burst_Analysis.m
%   Take "Av_Data" matrix from Ca_im_DataStream output and count 
%    avalanches using different definitions and thresholds.
% 
% INPUT: P_Av_Data - cell array sized by the number of nights
%         of data processed in the Ca_im_DataStream code and output as
%         Processed_Data.Av_Data
%        window - for avalanche definition (see below)
%        maxThreshold - also for avalanche definition, the code will
%          iterate from a threshold of 0.1 up to maxThreshold, where each
%          threshold signifies the fraction of ROIs that must be active
%          in a single frame for an avalanche to be counted
%          must be between 0 and 1
%        Fs - sampling frequency
%
% OUTPUT: Av_Size - a matrix with count data for each of the possible
%         thresholds and avalanche sizes ... if Av_Size(1,1) = 50, then
%         there were 50 instances of an avalanche of size 1
%
%         Av_IEI - a matrix with count data for each of the possible
%         avalanche thresholds and inter-event intervals, Av_IEI(1,1) would
%         be the number of instances of an avalanche of size at least 1
%         with an inter-event interval of 1/Fs, or one time bin ...
%         Av_IEI(1,2) would be the number of instances of an avalanche of
%         size at least 1 with an IEI of 2/Fs, or two time bins, so the
%         activity would be something like this [1,0,1], the time between
%         the first instance and the third instance is counted as two
%         frames
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
%      1/Fs seconds (contiguous bins active), 1 IEI of 2/Fs seconds (one empty bin
%      jumpted), 1 IEI of 3/Fs seconds, and the beginning and end won't be
%      counted.
%      That would be for threshold 1 ... Av_IEI(1,:). For a threshold of
%      two, stored in Av_IEI(2,:), you would have only 5 avalanches and
%      so the IEIs would change accordingly
%
%   If window = 1, then the definition is ... the summed activity across
%   all ROIs and time bins for which at least "threshold" ROIs are active.
%   So, for threshold = 1, contiguous bins with at least one ROI active are
%   counted as a single avalanche.
%      Example: Given the sequence above, [0,0,1,3,4,0,5,6,9,0,0,0], and a
%      threshold of 1, the inter-event inetervals would be the same, except
%      there would be no IEIs of 1/Fs seconds, because IEIs of 1/Fs seconds mean
%      two consecutive bins were active and, if that is true, those two bins 
%      are counted in the same avalanche.  The sizes, however, would change,
%      with 1 avalanche of size 1+3+4=8 ROIs and another avalanche of size 
%      5+6+9=20 ROIs.
%
%   As 'window' increases, the number of bins with activity below threshold
%   between bins with activity above or equal to threshold grows.  Thus,
%   for window = 2 and the sequence above [0,0,1,3,4,0,5,6,9,0,0,0], there
%   would now be only one avalanche of size 28.  There would now be no IEIs
%   counted because the beginning and end are ignored.
%
% Created: 2015/12/15 at 24 Cummington, Boston
%   Byron Price
% Updated: 2016/02/22
%  By: Byron Price

numNights = size(P_Av_Data,2);
numFrames = zeros(numNights,1);
for ii=1:numNights
    numFrames(ii) = size(P_Av_Data{ii},2);
end

maxFrames = max(numFrames);
Time = 1/Fs:1/Fs:maxFrames/Fs;
maxROIs = 500;   
threshVec = 0.1:0.1:maxThreshold;
Av_Size = zeros(length(threshVec),maxROIs); % avalanche size 

Av_IEI = zeros(length(threshVec),length(Time)); % inter-event interval (seconds), or time 
                % between avalanches

for zz = 1:numNights
    Av_Data = P_Av_Data{zz};
    numVideos = size(Av_Data,1);
    threshcount = 0;
    for threshold=threshVec
        threshcount = threshcount+1;
        for ii=1:numVideos
            newData = Av_Data(ii,:)./max(Av_Data(ii,:));
            newData(newData<threshold) = 0; % clear all activity below the threshold
            burstsAt = find(newData); % find indeces of avalanches, indeces
            % with the number of active ROIs above "threshold"
            if window == 0
                for jj=1:length(burstsAt)
                    Av_Size(threshcount,newData(burstsAt(jj))) = ...
                        Av_Size(threshcount,newData(burstsAt(jj)))+1;
                    if jj > 1
                        Av_IEI(threshcount,burstsAt(jj)-burstsAt(jj-1)) = ...
                            Av_IEI(threshcount,burstsAt(jj)-burstsAt(jj-1))+1;
                    end
                end
                
            elseif window > 0
                extraNewData = zeros(length(newData),1);
                for jj=1:length(burstsAt)
                    count = 0;
                    for kk=jj+1:length(burstsAt)
                        if (burstsAt(kk)-burstsAt(kk-1)) <= window
                            count = count+1;
                        else
                            break;
                        end
                    end
                    if count == 0
                        newcount = count;
                        extraNewData(burstsAt(jj)) = newData(burstsAt(jj));
                    else
                        summedActiv = sum(newData(burstsAt(jj):burstsAt(jj)+count));
                        if mod(count,2) == 0
                            extraNewData(burstsAt(jj)+count/2) = summedActiv;
                            newData(burstsAt(jj)+count/2) = summedActiv;
                            newcount = count/2;
                        elseif mod(count,2) == 1
                            randNum = rand;
                            if randNum <= 0.5
                                newcount = floor(count/2);
                                extraNewData(burstsAt(jj)+floor(count/2)) = summedActiv;
                                newData(burstsAt(jj)+floor(count/2)) = summedActiv;
                            elseif randNum > 0.5
                                newcount = ceil(count/2);
                                extraNewData(burstsAt(jj)+ceil(count/2)) = summedActiv;
                                newData(burstsAt(jj)+ceil(count/2)) = summedActiv;
                            end
                        end
                    end
                    newData([burstsAt(jj):(burstsAt(jj)+newcount-1),(burstsAt(jj)+newcount+1):burstsAt(jj+count)]) = 0;
                end
                clear burstsAt;
                burstsAt = find(extraNewData);
                for ll=1:length(burstsAt)
                    Av_Size(threshcount,extraNewData(burstsAt(ll))) = ...
                        Av_Size(threshcount,extraNewData(burstsAt(ll)))+1;
                    if ll > 1
                        Av_IEI(threshcount,burstsAt(ll)-burstsAt(ll-1)) = ...
                            Av_IEI(threshcount,burstsAt(ll)-burstsAt(ll-1))+1;
                    end
                end
            end
        end
    end
end
end

