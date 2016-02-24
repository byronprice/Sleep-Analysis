function [] = Jeffs_Ephys_Analysis(dates,earlyCutOff,lateCutOff)
%Jeff_Ephys_Analysis.m
%   Once you have used the Jeff_Ephys_Conversion.m function to get the data
%    into a usable format, you use this function.  It will perform some
%    fairly simple analyses of the data, including creation spectrograms of
%    the time-series over the night. Multi-taper spectrogram created thanks to 
%    Babadi & Brown 2014, Multi-taper spectral analysis
% INPUT: dates - This is the data for the folder output from the conversion
%         code.
%         Input the dates as a cell array, if you want to analyze
%         data from multiple conversions.  The converted data is
%         automatically save to a folder named 'Converted_Ephys_YYYY-MM-DD'
%         so the code will grab those folders with the dates provided.
%         dates = {'2016-02-16','2016-02-17'};
%        earlyCutOff - amount of time (in hours) to exclude from the 
%         beginning of the night, if the bird was not asleep, for example
%        lateCutOff - amount of time (in hours) to exclude at the end of the
%         night, if the bird woke up, for example
% OUTPUT: pretty pictures

% Created: 2016/01/19 at 24 Cummington, Boston
%   Byron Price
% Updated: 2016/02/23
% By: Byron Price

if nargin < 2
    earlyCutOff = 1; % 1 hour chopped from the beginning and end of the night
    lateCutOff = 1;
end
earlyCutOff = round(earlyCutOff*3600);
lateCutOff = round(lateCutOff*3600);

originalDirectory = pwd;
for ii=1:length(dates)
    cd(strcat(originalDirectory,'/Converted_Ephys_',dates{ii}))
    matrixNames = dir('*.mat');
    numFiles = length(matrixNames);
    load(matrixNames(1).name,'Fs')
    earlyCutOff = earlyCutOff*round(Fs);
    lateCutOff = lateCutOff*round(Fs);
    clear Fs;
    for jj=1:numFiles
        load(matrixNames(jj).name)
        numCombos = size(squareData,1);
        tData = squeeze(squareData(:,earlyCutOff:end-lateCutOff,2));
        sData = squeeze(spikeData(:,earlyCutOff:end-lateCutOff,2));
        timeSteps = size(tData,2);
        
        if mod(timeSteps,2) == 1
            tData = tData(:,1:end-1);
            sData = sData(:,1:end-1);
            timeSteps = timeSteps-1;
        end
        
        % IMPORTANT STEP FOR CREATION OF SPECTROGRAM
        % TIME AND FREQUENCY RESOLUTION OF THE RESULT
        T = 30; % data assumed stationary for T seconds, this should be an EVEN #
        N = round(T*Fs);
        if mod(N,2) == 1
            N = N+1;
        end
        R = 0.5; % desired spectral resolution (Hz)
        alpha = (T*R)/2; % must be greater than 1.25
        if alpha <= 1.25
            display(sprintf('alpha is equal to %2.4f',alpha))
            display('It needs to be greater than 1.25')
            display('Increase spectral resolution (R) or change stationarity time (T).')
            return;
        end
        K = N/2;
        times = N/2:K:timeSteps-N/2;
        realTimes = times./Fs;
        finalTime = realTimes(end);
        frequencySpectrogram = zeros(numCombos,length(times),N/2+1);
        IEIspectrogram = zeros(numCombos,length(times),N);
        ieis = (1:N)./Fs;
        x = linspace(0,finalTime/3600,length(realTimes));
        
        figure();
        plotcount = 1;
        for kk = 1:numCombos 
            % MAKE THE MULTI-TAPER SPECTROGRAM && IEI Spectrogram
            count = 1;
            for tt=times
                burstsAt = find(squeeze(sData(kk,(tt-(N/2-1)):(tt+(N/2)))));
                for zz=1:length(burstsAt)
                    if zz > 1
                        IEIspectrogram(kk,count,burstsAt(zz)-burstsAt(zz-1)) = ...
                            IEIspectrogram(kk,count,burstsAt(zz)-burstsAt(zz-1))+1;
                    end
                end
                [pxx,f] = pmtm(tData(kk,(tt-(N/2-1)):(tt+(N/2))),alpha,N,Fs);
                frequencySpectrogram(kk,count,:) = (squeeze(frequencySpectrogram(kk,count,:)))'+10*log10(pxx');
                count = count+1;
            end
            
            subplot(numCombos,2,plotcount)
            plotcount = plotcount+1;
            fspectro = squeeze(frequencySpectrogram(kk,:,:));
            imagesc(x,f,fspectro');
            hh = colorbar;
            ylabel(hh,'Power (dB/Hz)');
            title( ... 
                sprintf('Multi-taper Spectrogram for Electrode Combo #%i',kk));
                xlabel('Time (hours)');ylabel('Frequency (Hz)') 
            h = gca;
            h.YDir = 'normal';
            
            subplot(numCombos,2,plotcount)
            ieispectro = squeeze(IEIspectrogram(kk,:,:));
            imagesc(x,ieis,ieispectro')
            hh = colorbar;
            ylabel(hh,'Log Count')
            title( ... 
                sprintf('Inter-Event Interval Distribution for Electrode Combo #%i',kk));
                xlabel('Time (hours)');ylabel('Inter-Event Interval (seconds)') 
            h = gca;
            h.YDir = 'normal';
            plotcount = plotcount+1;
        end
        
        clear fspectro ieispectro;
        % MULTIVARIATE GAUSSIAN MAXIMUM LIKELIHOOD ESTIMATION - FREQUENCY
        fspectro = squeeze(mean(frequencySpectrogram,1));
        x_hat = (mean(fspectro,1))';
        sigma_hat = cov(fspectro);
        figure();
        subplot(1,2,1)
        plot(f,x_hat);title('Mean Power as a Function of Frequency');
        xlabel('Frequency (Hz)');ylabel('Power (dB/Hz)')
        subplot(1,2,2)
        imagesc(f,f,sigma_hat);title('Frequency Covariance Matrix \Sigma_f');
        xlabel('Frequency (Hz)');ylabel('Frequency (Hz)');colorbar
        
        % MULTIVARIATE GAUSSIAN - TIME
        x_hat = mean(fspectro,2);
        sigma_hat = cov(fspectro');
        figure();
        subplot(1,2,1)
        plot(x,x_hat);title('Mean Power as a Function of Time');
        xlabel('Time (hours)');ylabel('Power (dB/Hz)')
        subplot(1,2,2)
        imagesc(x,x,sigma_hat);title('Time Covariance Matrix \Sigma_t');
        xlabel('Time (hours)');ylabel('Time (hours)');colorbar
        
        % MULTIVARIATE GAUSSIAN MAXIMUM LIKELIHOOD ESTIMATION - FREQUENCY
        ieispectro = squeeze(mean(IEIspectrogram,1));
        x_hat = (mean(ieispectro,1))';
        sigma_hat = cov(ieispectro);
        figure();
        subplot(1,2,1)
        plot(ieis,x_hat);title('Mean Log Count as a Function of Inter-Event Interval');
        xlabel('Inter-Event Interval (seconds)');ylabel('Log Count')
        subplot(1,2,2)
        imagesc(ieis,ieis,sigma_hat);title('Inter-Event Interval Covariance Matrix \Sigma_IEI');
        xlabel('IEI (seconds)');ylabel('IEI (seconds)');colorbar
        
        % MULTIVARIATE GAUSSIAN - TIME
        x_hat = mean(ieispectro,2);
        sigma_hat = cov(ieispectro');
        figure();
        subplot(1,2,1)
        plot(x,x_hat);title('Mean Log Count as a Function of Time');
        xlabel('Time (hours)');ylabel('Power (dB/Hz)')
        subplot(1,2,2)
        imagesc(x,x,sigma_hat);title('Time Covariance Matrix \Sigma_t');
        xlabel('Time (hours)');ylabel('Time (hours)');colorbar
        clear fspectro x_hat sigma_hat Data squareData spikeData ieispectro;
    end
    clear fileNames numFiles numCombos times finalTime realTimes;
    earlyCutOff = earlyCutOff/round(Fs);
    lateCutOff = lateCutOff/round(Fs);
end
cd(originalDirectory)
end

