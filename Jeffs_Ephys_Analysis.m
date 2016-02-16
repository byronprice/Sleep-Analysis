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
%        earlyCutOff - amount of time (in hours) to exclude from the 
%         beginning of the night, if the bird was not asleep, for example
%        lateCutOff - amount of time (in hours) to exclude at the end of the
%         night, if the bird woke up, for example
% OUTPUT: pretty pictures

% Created: 2016/01/19 at 24 Cummington, Boston
%   Byron Price
% Updated: 2016/02/15
% By: Byron Price

if nargin < 2
    earlyCutOff = 1; % 1 hour chopped from the beginning and end of the night
    lateCutOff = 1;
end

originalDirectory = pwd;
for ii=1:length(dates)
    cd(strcat(originalDirectory,'/Converted_Ephys_',dates{ii}))
    fileNames = dir('*.mat');
    numFiles = length(fileNames);
    for jj=1:numFiles
        fileNames = dir('*.mat');
        load(fileNames(jj).name)
        numCombos = size(squareData,1);
        earlyCutOff = (earlyCutOff*60*60)*Fs;
        lateCutOff = (lateCutOff*60*60)*Fs;
        Data = squeeze(squareData(:,earlyCutOff:end-lateCutOff,2));
        
        timeSteps = size(Data,2);
        if mod(timeSteps,2) == 1
            Data = Data(:,1:end-1);
            timeSteps = timeSteps-1;
        end
        
        % IMPORTANT STEP FOR CREATION OF SPECTROGRAM
        % TIME AND FREQUENCY RESOLUTION OF THE RESULT
        T = 60; % data assumed stationary for T seconds, this should be an EVEN #
        N = T*Fs;
        R = 2; % desired spectral resolution (Hz)
        alpha = (N*R)/(2*Fs); % must be greater than 1.25
        if alpha <= 1.25
            display(sprintf('alpha is equal to %2.4f',alpha))
            display('It needs to be greater than 1.25')
            display('Increase spectral resolution (R) or change stationarity time (T).')
            return;
        end
        tau = T/2;
        K = (tau/T)*N;
        
        times = N/2:K:timeSteps-N/2;
        realTimes = times./Fs;
        finalSpectrogram = zeros(numCombos,length(times),N/2+1);
        figure();
        numRows = ceil(numCombos/2);
        for kk = 1:numCombos 
            % MAKE THE MULTI-TAPER SPECTROGRAM
            count = 1;
            for tt=times
                [pxx,f] = pmtm(Data(kk,(tt-(N/2-1)):(tt+(N/2))),alpha,N,Fs);
                finalSpectrogram(kk,count,:) = finalSpectrogram(kk,count,:)+10*log10(pxx');
                count = count+1;
            end
            x = linspace(0,realTimes(end)./3600,length(realTimes));
            
            subplot(numRows,2,kk)
            spectro = squeeze(finalSpectrogram(kk,:,:));
            imagesc(x,f,spectro');colorbar;title( ... 
                sprintf('Multi-taper Spectrogram for Electrode Combo #%i',kk));
                xlabel('Time (hours)');ylabel('Frequency (Hz)') 
            h = gca;
            h.YDir = 'normal';
        end
        
        clear spectro;
        % MULTIVARIATE GAUSSIAN MAXIMUM LIKELIHOOD ESTIMATION - FREQUENCY
%         timeSteps = size(finalSpectrogram,1);
%         numFreqs = size(finalSpectrogram,2);
        spectro = finalSpectrogram(1,:,:);
        x_hat = (mean(spectro,1))';
        sigma_hat = cov(spectro);
        figure();
        subplot(1,2,1)
        plot(f,x_hat);title('Mean Power as a Function of Frequency');
        xlabel('Frequency (Hz)');ylabel('Power (dB/Hz)')
        subplot(1,2,2)
        imagesc(f,f,sigma_hat);title('Frequency Covariance Matrix \Sigma_f');
        xlabel('Frequency (Hz)');ylabel('Frequency (Hz)');colorbar
        
        % MULTIVARIATE GAUSSIAN - TIME
        x_hat = mean(spectro,2);
        sigma_hat = cov(spectro');
        figure();
        subplot(1,2,1)
        plot(x,x_hat);title('Mean Power as a Function of Time');
        xlabel('Time (hours)');ylabel('Power (dB/Hz)')
        subplot(1,2,2)
        imagesc(x,x,sigma_hat);title('Time Covariance Matrix \Sigma_t');
        xlabel('Time (hours)');ylabel('Time (hours)');colorbar
        clearvars -except originalDirectory numFiles jj ii ...
            earlyCutOff lateCutOff dates;
    end
end
cd(originalDirectory)
end
