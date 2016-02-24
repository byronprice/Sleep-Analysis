function [] = Jeffs_Ephys_Conversion(cellArrayBirdNames,Fs)
%Jeff_Ephys_Conversion.m
%   Takes Jeff's data, combine the ~400 recordings in different files into a
%    single matrix, downsample to Fs Hz, rectify by squaring, subtract one
%    from another unless that step has already been done.
% 
% INPUT: cellArrayBirdNames - if you download the data as Jeff has it
%     saved, e.g. /lpi79/2011-03-04/sleep/sleepdata1_lpi79_110304_220008.mat
%     then, you can send in input as {'lpi79','lpi16'} or something
%     similar, whatever the bird names are, as long as all files are stored
%     in the same format
%        Fs - sampling frequency TO WHICH you want to downsample, due to
%        the way the function downsamples, the original sampling frequency
%        divided by the this new sampling frequency must equal a whole
%        number, if it does not, then the value of Fs will be rounded
% OUTPUT: no direct output, but all of the matrices will be stored in a
%  subdirectory of the original directory as 
%   Converted_Ephys_currentDate/combinedData_birdname_recordingDate.mat

% Then, you can send those files into the analysis function.

% Created: 2016/01/19 at 24 Cummington, Boston
%   Byron Price
% Updated: 2016/02/22
% By: Byron Price


% get original directory
originalDirectory = pwd;
currentDate = datetime('now','Format','yyyy-MM-dd');
resultDirectory = strcat('Converted_Ephys_',char(currentDate));
mkdir(resultDirectory)

% LOOP FOR THE DIFFERENT BIRDS
for bird = 1:length(cellArrayBirdNames)
    % FOR EACH BIRD, LOAD A SINGLE EXAMPLE FILE TO GET NUMBERS OF DIFFERENT ITEMS
    birdname = cellArrayBirdNames{bird};
    currentBirdDirectory = strcat(originalDirectory,'/',birdname);
    cd(currentBirdDirectory)
    subFolders = dir('*-*');
    example = subFolders(1).name;
    cd(strcat(pwd,'/',example,'/sleep'))
    files = dir('*.mat');
    load(files(1).name)
    numElectrodes = size(adc.names,2);
   
    FsOld = adc.fs;
    downSampleRate = round(FsOld/Fs); 
    Fs = FsOld/downSampleRate;
    
    clear FsOld adc audio parameters example ;
    
    % LOOP THROUGH EACH OF THE FOLDERS FOR A GIVEN BIRD, i.e. DIFFERENT
    %  DATES
    for ii=1:length(subFolders)
        cd(strcat(currentBirdDirectory,'/',subFolders(ii).name,'/sleep'))
        sleepFiles = dir('*.mat'); 
        numFiles = length(sleepFiles);
        fileNames = cell(1,numFiles);
       
        % CREATE SINGLE MATRIX WITH ALL THE DATA
        %  not all files are the same size, so we have to loop through each
        %   to figure out the total size and then add them into a final
        %   matrix one-by-one
        
        totalLength = 0;
        for jj=1:numFiles
            fileNames{jj} = sleepFiles(jj).name; 
            load(sleepFiles(jj).name,'adc')
            totalLength = totalLength+size(adc.data,1);
        end
        
        numCombos = nchoosek(numElectrodes,2);
        timeSteps = floor(totalLength/downSampleRate)+1;
        
        if ContainSubString(adc.names{1},'Gnd') == 1
            % SUBTRACT ELECTRODES FROM EACH OTHER, RECTIFY,
            %  SMOOTH, LOWPASS FILTER, DOWNSAMPLE
            TimeVec = zeros(numCombos,timeSteps);
            originalData = zeros(numCombos,timeSteps);
            squareData = zeros(numCombos,timeSteps);
            leftovers = zeros(numCombos,numFiles+1);
            leftovers(:,1) = downSampleRate*ones(numCombos,1);
            
            count = 1;
            for kk=1:numElectrodes
                for ll=kk+1:numElectrodes
                    index = 1;
                    for jj=1:numFiles
                        load(sleepFiles(jj).name,'adc','parameters')
                        stdev = std(adc.data(:,kk));
                        adc.data(adc.data(:,kk)>stdev*15,kk) = 0;
                        adc.data(adc.data(:,kk)<-stdev*15,kk) = 0;
                        smoothData = smooth((adc.data(:,kk)-adc.data(:,ll)).^2,downSampleRate);
                        unsmoothData = adc.data(:,kk)-adc.data(:,ll);
                        
                        currentLength = length(smoothData((1+downSampleRate-leftovers(kk,jj)):end));
                        newLength = ceil(currentLength/downSampleRate);
                        leftovers(kk,jj+1) = mod(currentLength-1,downSampleRate);
                        
                        % DOWNSAMPLE
                        TimeVec(count,index:index+newLength-1) = adc.t((1+downSampleRate-leftovers(kk,jj)):downSampleRate:end);
                        unsmoothData = unsmoothData((1+downSampleRate-leftovers(kk,jj)):downSampleRate:end);
                        newData = smoothData((1+downSampleRate-leftovers(kk,jj)):downSampleRate:end); 
                        
                        % FILTER
                        n = 30; 
                        lowpass = 0.99; % fraction of Nyquist frequency
                        blo = fir1(n,lowpass,'low',hamming(n+1));
                        newData = filter(blo,1,newData);
                        squareData(count,index:index+newLength-1) = newData;
                        
                        unsmoothData = filter(blo,1,unsmoothData);
                        originalData(count,index:index+newLength-1) = unsmoothData;
                        index = index+newLength;
                    end
                    count = count+1;
                end
            end  
        else
            TimeVec = zeros(numElectrodes,timeSteps);
            originalData = zeros(numElectrodes,timeSteps);
            squareData = zeros(numElectrodes,timeSteps);
            for kk=1:numElectrodes
                Data = zeros(totalLength,2);
                index = 1;
                for jj=1:numFiles
                    load(sleepFiles(jj).name,'adc','parameters')
                    stdev = std(adc.data(:,kk));
                    adc.data(adc.data(:,kk)>stdev*15,kk) = 0;
                    adc.data(adc.data(:,kk)<-stdev*15,kk) = 0;
                    currentLength = size(adc.data,1);
                    Data(index:index+currentLength-1,1) = adc.t;
                    Data(index:index+currentLength-1,2) = adc.data(:,kk);
                    index = index+currentLength+1;
                end
                
                %  RECTIFY, SMOOTH, LOWPASS FILTER
                %  DOWNSAMPLE (ELECTRODES ALREADY SUBTRACTED)
                
                TimeVec(kk,:) = squeeze(Data(1:downSampleRate:end,1));
                temp = squeeze(Data(:,2));
                temp = smooth(temp.^2,downSampleRate);
                newTemp = temp(1:downSampleRate:end);
                n = 30;
                lowpass = 0.99; % fraction of Nyquist frequency (FsOld/2)
                blo = fir1(n,lowpass,'low',hamming(n+1));
                newTemp = filter(blo,1,newTemp);
                squareData(kk,:) = newTemp;
                
                temp = squeeze(Data(:,2));
                newTemp = temp(1:downSampleRate:end);
                newTemp = filter(blo,1,newTemp);
                originalData(kk,:) = newTemp;
                clear Data;
            end
        end
     numCombos = size(squareData,1);
     startTime = parameters.rec_start_datenum;
     dataPoints = size(squareData,2);
     
     spikeData = zeros(numCombos,dataPoints);
     for zz =1:numCombos
         test = squeeze(squareData(zz,:));
         test = test-mean(test);
         thresh = std(test);
         test(test<thresh) = 0;
         test(test>0) = 1;
         spikeData(zz,:) = test;
     end
     filename = strcat('combinedData_',birdname,'_',subFolders(ii).name,'.mat');
     cd(strcat(originalDirectory,'/',resultDirectory))
     save(filename,'originalData','squareData','spikeData','TimeVec','Fs','fileNames','startTime')
        
    end
clearvars -except originalDirectory resultDirectory cellArrayBirdNames bird
cd(originalDirectory)

end
end

function [ trueorfalse ] = ContainSubString(MainString, SubString)
% ContainsSubstring.m
% Created: 06/22/2015
% Updated: 08/24/2015
% Byron Price

% Checks a string to determine whether or not a given substring forms part

% INPUT: MainString - a string
%        SubString - a substring
%
% OUTPUT: trueorfalse - Boolean value of 0 (false, MainString does not 
%       contain SubString) or 1 (true, MainString does contain SubString)

MainLength = length(MainString);
SubLength = length(SubString);

if SubLength > MainLength
    % SubString cannot be longer than MainString
    trueorfalse = false;
    return;
else
    % Call the recursive check between MainString and SubString
    % if RecursiveCheckBoolean returns true, SubString forms part of
    % MainString, and if false, it must keep searching
    for z = 1:(MainLength-(SubLength-1))
        trueorfalse = RecursiveCheckBoolean(MainString,SubString,z,1,SubLength);
        if trueorfalse == true
            return;
        else
            continue;
        end
    end

    return;

end

end

function [ aretrue ] = RecursiveCheckBoolean(Main,Sub,MainIndex,SubIndex,endcondition)
% Recursively compares each letter in Sub to letters in
% Main
    booleancheck = CheckLetter(Main(MainIndex),Sub(SubIndex));
    if booleancheck == true && SubIndex ~= endcondition
        aretrue = RecursiveCheckBoolean(Main,Sub,MainIndex+1,SubIndex+1,endcondition);
    elseif booleancheck == true && SubIndex == endcondition
        aretrue = true;
        return;
    else
        aretrue = false;
        return;
    end
end

function [ doescontain ] = CheckLetter(MainLetter, SubLetter)
    if MainLetter == SubLetter
        doescontain = true;
    else
        doescontain = false;
    end
end
