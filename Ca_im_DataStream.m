function [Av_Data,Spike_Data,Trace_Data] = Ca_im_DataStream(foldernums,FitType,threshold,Fs,binSize)
%Ca_im_DataStream.m
%Take stored data in different folders and perform Avalanche analysis
%    Folders are named 'mat*' with '*' denoting different numbers
% Created: 2015/10/21 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/11/18
% By: Byron Price
%
% INPUT:   foldernums - cell array with strings for which folders to upload
%          FitType - string with type of fit to perform on data from single
%                       ROIs within individual recordings, this is for 
%                       detrending bleaching and some noise, MATLAB allows
%                       varied input, like 'poly1' for linear or 'fourier8'
%                       for 8-term Fourier series
%          THREE OPTIONAL INPUTS 
%          threshold - z-score threshold (peaks in activity above this value
%                       are counted as spikes, default = 2
%          Fs - sampling frequency (Hz), default = 30
%          binSize - width of bins (secs), spikes are counted as having occurred
%                       within bins of a certain size, default = 1/Fs
% 
% OUTPUT:  Av_Size - binned data for avalanche size, i.e. frequency
%                       of the occurence of avalanches with size x
%          Av_Len - frequency of occurence of avalanches with
%                       length t
%          IEI - frequency of occurence of length-t lulls
%                       in avalanche activity

% ALLOW FOR VARIABLE INPUT AND SET DEFAULTS
% Check number of inputs.
if nargin > 5
    error('Gardner:DataStream:TooManyInputs', ...
        'Requires at most 3 optional inputs');
end

% Fill in unset optional inputs.
switch nargin
    case 2
        threshold = 2;
        Fs = 30;
        binSize = 1/30;
    case 3
        Fs = 30;
        binSize = 1/30;
    case 4
        binSize = 1/30;
end


% LOOP THROUGH THE VARIOUS FOLDERS, each of which has data from a set of 
                   % recordings, each ave_roi.mat is a structure with
                   % a cell inside called 'raw' ...  raw contains an
                   % unknown number of recordings (a set), each ~30 seconds long,
                   % which is 900 frames at 30 Hz, though some are sized
                   % differently, 1800 frames so 60 seconds
maxROIs = 1000;             
Av_Size = zeros(maxROIs+1,2); % avalanche size (number of active ROIs in a 
                % single time bin, unless contiguous bins are active, then
                % its the number of active ROIs summed over continuously
                % active bins)
Av_Size(:,2) = 0:maxROIs;

maxFrames = 1800; % maximum number of frames
times = binSize:binSize:maxFrames/Fs;
Av_IEI = zeros(length(times),2); % inter-event interval (seconds), or time 
                % between avalanches
Av_IEI(:,2) = times;

Av_Len = zeros(length(times),2); % length (seconds) of avalanches
Av_Len(:,2) = times;


for a=1:length(foldernums)
    % ACCESS INFO IN DIFFERENT FOLDERS, ALL OF WHICH HAVE BEEN NAMED
    % FOR SIMPLICITY, 'mat6' or 'mat7' or 'mat11'
    disp(strcat('Currently Running Folder #',char(foldernums(a))))
    folder = strcat('mat',char(foldernums(a)));
    mainfolder = strcat('/Users/gardnerlab/Documents/MATLAB/',folder);
    directory = strcat(mainfolder,'/roi');
    cd(directory)
    
    load('ave_roi.mat')
    
    alldata = roi_ave.raw;
    numVideos = size(alldata,2); 
    numROIs = size(alldata{1},1);
    shift = 1*Fs; % how many frames to cut off from beginning and end
                    % the number is in units of seconds
    numFrames = length(alldata{1,1}(1,shift:end-shift));
    
    % DETREND EACH ROI INDIVIDUALLY FOR EACH VIDEO
    Trace_Data = zeros(numVideos,numROIs,numFrames);
    for i = 1:numVideos
        for j = 1:numROIs
            [alldata{1,i}(j,shift:end-shift)] = Preprocessing(alldata{1,i}(j,shift:end-shift),FitType);
        end
        Trace_Data(i,:,:) = alldata{1,i}(:,shift:end-shift);
    end
    
    
    % COUNT SPIKES AND SUM ACROSS ROIs IN A SINGLE RECORDING
    Spike_Data = zeros(numVideos,numROIs,numFrames);
    Av_Data = zeros(numVideos,numFrames);
    for j=1:numVideos
        binarySpikes = []; % will contain event data for all ROIs in 
                            % a single recording
        for k=1:numROIs
            % detect calcium spikes
            [~,Spikes] = Spike_Detector(Trace_Data(j,k,:),Fs,binSize,threshold);
            binarySpikes = [binarySpikes,Spikes];
        end
        Spike_Data(j,:,:) = binarySpikes';
        sumSpikes = sum(binarySpikes,2); % add up spikes from each ROI
                                            % into a single vector, this
                                            % will be used to deterimine
                                            % avalanche size
                                                         
        Av_Data(j,:) = sumSpikes;
        % SUM SPIKES ACROSS BINS TO COUNT AVALANCHE SIZE, LENGTH, AND ALSO
        % INTER-EVENT INTERVAL
        %scount = 0;
        %lcount = 0;
        %icount = 0;
    
        %for t =2:length(sumSpikes)
        %    [~,~,Av_Size(:,1),Av_Len(:,1),Av_IEI(:,1),scount,lcount,icount] = Av_Counter(sumSpikes,t,Av_Size(:,1),Av_Len(:,1),Av_IEI(:,1),scount,lcount,icount);
        %end
    end
end

% PLOTS
%figure
%subplot(3,2,1)
%bar(Av_Size(2:end,1)) % starting at 2 skips zero-sized bin, for visual
%title('Distribution of Avalanche Size')
%xlabel('Avalanche Size (# ROIs above threshold)')
%ylabel('Frequency')


%subplot(3,2,2)
%scatter(log10(Av_Size(2:600,2)),log10(Av_Size(2:600,1)))
%title('Log-Log Distribution of Avalanche Size')
%xlabel('Log(Avalanche Size)')
%ylabel('Log(Frequency)')

%subplot(3,2,3)
%bar(Av_IEI(1:300,2),Av_IEI(1:300,1))
%title('Distribution of Inter-Event Interval')
%xlabel('Inter-Event Interval (seconds between avalanches)')
%ylabel('Frequency')

%subplot(3,2,4)
%scatter(log10(Av_IEI(:,2)),log10(Av_IEI(:,1)))
%title('Log-Log Distribution of Inter-Event Interval')
%xlabel('Log[Inter-Event Interval (seconds)]')
%ylabel('Log(Frequency)')

%subplot(3,2,5)
%bar(Av_Len(1:150,2),Av_Len(1:150,1))
%title('Distribution of Avalanche Length')
%xlabel('Avalanche Length (seconds of continuous spikes)')
%ylabel('Frequency')

%subplot(3,2,6)
%scatter(log10(Av_Len(:,2)),log10(Av_Len(:,1)))
%title('Log-Log Distribution of Avalanche Length')
%xlabel('Log[Avalanche Length (seconds)]')
%ylabel('Log(Frequency)')
cd('/Users/gardnerlab/Documents/MATLAB/')

end

function [numBins,Spikes] = Spike_Detector(timeSeries,sampleFreq,binSize,threshold)
%Avalanche.m  
%   Detect calcium spikes in imaging data for a single ROI.
% Created: 2015/09/30 at 24 Cummington, Boston
%  Byron Price
% Updated: 2015/12/2
%  By: Byron Price
%
% REFERENCE: Klaus, Plenz 2011 Statistical Analyses Support Power Law ...
%            Ribeiro, Copelli et al. 2010 Spike Avalanches Exhibit ...
%
% INPUT: timeSeries - change in activity over time for a single ROI in
%               units given by T = 1/Fs  ... Fs = sampling frequency
%        sampleFreq - sampling frequency, Hz
%        binSize - width of bin (seconds)
%               Spikes are counted as having occurred within a given 
%               bin if the activity of that bin is above "threshold"
%        threshold - zscore threshold (1.5, 2 etc. standard deviations
%               above the mean) for decision: above threshold = activity 
%               forms part of an avalanche
% OUTPUT: numBins - given the binSize and the length of the recording,
%               numBins = (total # frames) / (# frames / bin)
%         Spikes - a vector containing either a 0 or a 1, 1 being an
%               spike of activity within that bin, 0 otherwise
%         binSize - binSize (sec) should be an integer multiple of the period T
%               (1/Fs), if it is not, it will automatically be adjusted,
%               this is crucial for calculation

%stdActivity = std(timeSeries);
T = 1/sampleFreq; % period for a single frame
framesPerBin = round(binSize/T); 


times = 1:framesPerBin:length(timeSeries);
numBins = length(times); %number of bins, based on the
                         % length of the data and bin size (secs)

Spikes = zeros(numBins,1); 
timeSeries = timeSeries./max(timeSeries);
stdActivity = std(timeSeries);
% SPIKE DETECTION LOOP

% the trace is divided by its maximum value, there are two steps to being
% counted as a calcium spike: 1) have value above 0.2, or be at least 20% 
% of the maximum fluorescence change; and 2) be a local maximum, meaning
% the fluorescence intensity of the current time step should be greater
% than both of its neighbors

for i=2:numBins-1
        if timeSeries(i)/stdActivity > threshold % tentatively counted as spike if above threshold
            if timeSeries(i) > timeSeries(i-1) && timeSeries(i) > timeSeries(i+1)
                Spikes(i) = 1;
            end
        end
end
end  

function [timeSeries] = Preprocessing(timeSeries,FitType)
% Av_Preprocessing.m
%   Detrend a time series and subtract out the mean
% Created:2015/10/21 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/10/21
%   By: Byron Price
% INPUT: timeSeries - calcium imaging pixel intensities over time for a
%                       single ROI
% OUTPUT: timeSeries - detrended data

% LINEAR REGRESSION TO SUBTRACT OUT BLEACHING
% x = 1:1:length(timeSeries);
% X = [ones(length(x),1) x'];
% b = X\timeSeries';
% ycalc = X*b;
% 
% timeSeries = (timeSeries'-ycalc)';

% REGRESSION TO SUBTRACT OUT BLEACHING AND OTHER 
% INHOMOGENEITIES IN THE DATA
x = 1:length(timeSeries);
f = fit(x',timeSeries',FitType);
Y = f(x);

timeSeries = (timeSeries'-Y)';

timeSeries = timeSeries./max(timeSeries);



end

function [Spikes,timeStep,av_size,av_len,iei_len,scount,lcount,icount] = Av_Counter(Spikes,timeStep,av_size,av_len,iei_len,scount,lcount,icount)
% Av_Counter.m
%   Count Avalanches, this is meant to be run as part of a for loop in
%           which it will step through the data in Spikes, counting 
%           avalanche sizes and lengths
% Created:2015/10/21 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/10/21
%   By: Byron Price
% INPUT: Spikes - binned and summed spiking data across all ROIs for a
%                   single recording, output of Avalanche.m for multiple
%                   ROIs summed together
%        timeStep - set initially to 2
%
%               THE NEXT THREE SHOULD ALL INITIALLY BE EMPTY VECTORS
%        av_size - mX1 vector of avalanche sizes, the i-th bin signifies the 
%                    number of avalanches of size i in the recording
%        av_len - tX1 vector of avalanche lengths, the i-th bin signifies the
%                    number of avalanches of length i
%        iei_len - tX1 vector of inter-event intervals, the i-th bin
%                    signifies the number of IEIs of length i
%        
%        scount - lcount - icount ... for coding, set to 0 initially
%
% OUTPUT: same as the inputs, but with av_size, av_len, iei_len filled in,
%           Spikes, scount, lcount, icount, timeStep can be ignored at the end

% WHAT TO DO FOR LAST TIME STEP
%    Check stored values for avalanche size and length, or inter-event 
%    interval, and add them to the total counts 
if timeStep == length(Spikes)
    if scount > 0
        av_size(scount) = av_size(scount)+1;
        av_len(lcount) = av_len(lcount)+1;
    elseif icount > 0
        iei_len(icount) = iei_len(icount)+1;
    end
end

% WHAT TO DO FOR TIME STEP 2, which should be the first step in a for loop
%   Look back at first time step and see if it's part of an avalanche or 
%   just an empty bin.
if timeStep == 2
    if Spikes(1) >= 1
        scount = Spikes(1);
        lcount = 1;
    else
        icount = 1;
    end
end

% FOUR POSSIBILITIES: 
%   1) Avalanche at t and t-1 : add to growing count of avalanche size and 
%       and length, then keep going
%   2) No activity at t and t-1 : add to growing count of inter-event
%       interval and keep going
%   3) No activity at t but activity at t-1 : start counting inter-event 
%       interval, but stop counting avalanches, take counted avalanche size
%       and add 1 to the bin for avalanches of that size, i.e. we added up
%       the total number of ROIs active in the avalanche and that is now
%       one incidence of such an avalanche
%   4) Activity at t but none at t-1 : start counting avalanches and stop
%       counting inter-event interval, take counted inter-event interval
%       and add 1 to the bin for IEIs of that length
if Spikes(timeStep) >= 1 && Spikes(timeStep-1) >=1
    scount = scount+Spikes(timeStep);
    lcount = lcount+1;

elseif Spikes(timeStep) == 0 && Spikes(timeStep-1) == 0
    icount = icount+1;

elseif Spikes(timeStep) == 0 && Spikes(timeStep-1) >= 1
    icount = icount+1;
    av_size(scount) = av_size(scount)+1;
    av_len(lcount) = av_len(lcount)+1;
    scount = 0;
    lcount = 0;

elseif Spikes(timeStep) >= 1 && Spikes(timeStep-1) == 0
    scount = scount+Spikes(timeStep);
    lcount = lcount+1;
    iei_len(icount) = iei_len(icount)+1;
    icount = 0;
end


end
