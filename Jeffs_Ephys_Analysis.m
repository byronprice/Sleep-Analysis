function [] = Jeffs_Ephys_Analysis(dates)
%Jeff_Ephys_Analysis.m
%   Once you have used the Jeff_Ephys_Conversion.m function to get the data
%    into a usable format, you use this function.  It will perform some
%    fairly simple analyses of the data, including creation of spectrograms of
%    the time-series over the night.  
% INPUT: dates - This is the data for the folder output from the conversion
%         code.
%         Input the dates as a cell array, if you want to analyze
%         data from multiple conversions.  The converted data is
%         automatically save to a folder named 'Converted_Ephys_YYYY-MM-DD'
%         so the code will grab those folders with the dates provided.
% OUTPUT: pretty pictures

% Created: 2016/01/19 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/01/21
% By: Byron Price

originalDirectory = pwd;
for i=1:length(dates)
    cd(strcat(originalDirectory,'/Converted_Ephys_',dates{i}))
    fileNames = dir('*.mat');
    numFiles = length(fileNames);
    for j=1:numFiles
        fileNames = dir('*.mat');
        load(fileNames(j).name)
        numCombos = size(squareData,1);
        Data = squeeze(squareData(:,10000:end-10000,2));
        
        timeSteps = size(Data,2);
        if mod(timeSteps,2) == 1
            Data = Data(:,1:end-1);
            timeSteps = timeSteps-1;
        end
        
        % IMPORTANT STEP FOR CREATION OF SPECTROGRAM
        % TIME AND FREQUENCY RESOLUTION OF THE RESULT
        T = 30; % data assumed stationary for T seconds, this should be an EVEN #
        N = T*Fs;
        R = 1; % desired spectral resolution (Hz)
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
        
        figure();
        numRows = ceil(numCombos/2);
        for k = 1:numCombos 
            finalSpectrogram = zeros(length(times),N/2+1);
            count = 1;
            for tt=times
                [pxx,f] = pmtm(Data(k,(tt-(N/2-1)):(tt+(N/2))),alpha,N,Fs);
                finalSpectrogram(count,:) = finalSpectrogram(count,:)+10*log10(pxx');
                count = count+1;
            end
            
            x = linspace(0,realTimes(end)/3600,length(realTimes));
            
            subplot(numRows,2,k)
            imagesc(x,f,finalSpectrogram');colorbar;title( ... 
                sprintf('Multi-taper Spectrogram for Electrode Combo #%d',k));
                xlabel('Time (hours from 10pm)');ylabel('Frequency (Hz)') 
            h = gca;
            h.YDir = 'normal';
        end
        
        
        % MULTIVARIATE GAUSSIAN MAXIMUM LIKELIHOOD ESTIMATION - FREQUENCY
        timeSteps = size(finalSpectrogram,1);
        numFreqs = size(finalSpectrogram,2);
        x_hat = zeros(numFreqs,1);
        sigma_hat = zeros(numFreqs,numFreqs);
        
        for z=1:numFreqs
            x_hat(z) = mean(finalSpectrogram(:,z));
        end
        for z=1:timeSteps
            x_t = finalSpectrogram(z,:);
            x_t = x_t';
            sigma_hat = sigma_hat + (1/timeSteps).*((x_t-x_hat)*(x_t-x_hat)');
        end
        figure();
        subplot(1,2,1)
        plot(f,x_hat);title('Mean Power as a Function of Frequency');
        xlabel('Frequency (Hz)');ylabel('Power (dB/Hz)')
        subplot(1,2,2)
        imagesc(f,f,sigma_hat);title('Frequency Covariance Matrix \Sigma_f');
        xlabel('Frequency (Hz)');ylabel('Frequency (Hz)');colorbar
        
        % MULTIVARIATE GAUSSIAN - TIME
        x_hat = zeros(timeSteps,1);
        sigma_hat = zeros(timeSteps,timeSteps);
        
        for z=1:timeSteps
            x_hat(z) = mean(finalSpectrogram(z,:));
        end
        for z=1:numFreqs
            x_f = finalSpectrogram(:,z);
            sigma_hat = sigma_hat + (1/numFreqs).*((x_f-x_hat)*(x_f-x_hat)');
        end
        x = linspace(0,realTimes(end)/3600,length(realTimes));
        figure();
        subplot(1,2,1)
        plot(x,x_hat);title('Mean Power as a Function of Time');
        xlabel('Time (hours)');ylabel('Power (dB/Hz)')
        subplot(1,2,2)
        imagesc(x,x,sigma_hat);title('Time Covariance Matrix \Sigma_t');
        xlabel('Time (hours)');ylabel('Time (hours)');colorbar
        clearvars -except originalDirectory numFiles squareData
    end
end
cd(originalDirectory)
end
