function [calib] = calibEar_swept

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to figure out correct attenuations for chin for SFOAE
% S. Hauser 3/13/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subj = input('Enter subject ID: ', 's');
earname = input('Enter L or R ear: ', 's');
calibNum = input('What NEL calib #? ', 's');

% Make directory to save results
paraDir = 'C:\Users\Heinz Lab - NEL2\Desktop\OAEs\Calib';
addpath(genpath(paraDir));
if(~exist(strcat(paraDir,'\',subj),'dir'))
    mkdir(strcat(paraDir,'\',subj));
end
respDir = strcat(paraDir,'\',subj,'\');

date = datetime('now');
date.Format = 'yyyy_MM_dd';
datetag = string(date);
calib.date = datetag;

% Get CalibFile from NEL DATA storage
ExpDataDir = 'C:\NEL\ExpData\';
cd(ExpDataDir);
findDir = dir(sprintf('*%s*%s*', datetag, subj));
if isempty(findDir)
    error('You need to run a NEL calibration first!\n');
elseif length(findDir)~=1
    error('Multiple Directories. I am confused. \n')
else
    datadir = [ExpDataDir findDir.name];
    cd(datadir);
    load(sprintf('coef_000%s_calib.mat', calibNum), 'b');
    calib.b = b; 
end

cd(paraDir);

% variables you can change
attenCh1 = 10;
attenCh2 = 10;
micSens = 50; % mV/Pa from standard of box 46.5mV is what comes out of calibrator (=94 dB)
freq = 1000; % play 1000 Hz tone to set level
fs = 48828.125;

t=0:1/fs:8; % 8 seconds long

% Initializing TDT
fig_num=99;
GB_ch=1;
FS_tag = 3;
Fs = 48828.125;
[f1RZ,RZ,~]=load_play_circuit_Nel2(FS_tag,fig_num,GB_ch);

%Set tones to play
tone = scaleSound(rampsound(cos(2*pi*freq*t), fs, 0.005));

for i = 1:2 %play tone in each channel
    
    buffdata = zeros(2, numel(t));
    buffdata(i, :) = tone;
    
    % Load the 2ch variable data into the RZ6:
    invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', buffdata(1, :));
    invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', buffdata(2, :)); %empty
    
    % Set the delay of the sound
    invoke(RZ, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
    playrecTrigger = 1;
    
    % Set attenuations
    rc = PAset([0, 0, attenCh1, attenCh2]);
    
    % Set total length of sample
    RZ6ADdelay = 97; % Samples
    resplength = size(buffdata,2) + RZ6ADdelay; % How many samples to read from OAE buffer
    invoke(RZ, 'SetTagVal', 'nsamps', resplength);
    
    %Start playing from the buffer:
    fprintf(1, 'Playing!\n');
    invoke(RZ, 'SoftTrg', playrecTrigger);
    currindex = invoke(RZ, 'GetTagVal', 'indexin');
    
    while(currindex < resplength)
        currindex=invoke(RZ, 'GetTagVal', 'indexin');
    end
    
    vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
        'F32','F64',1);
    
    % save the response
    resp(i,:) = vin((RZ6ADdelay + 1):end);
    
    % input RMS from oscilloscope
    RMS(i,1) = input('What is the RMS in mV on the oscilloscope? ');
    
    % Get ready for next trial
    invoke(RZ, 'SoftTrg', 8); % Stop and clear "OAE" buffer
    %Reset the play index to zero:
    invoke(RZ, 'SoftTrg', 5); %Reset Trigger
end

close_play_circuit(f1RZ, RZ);

%% Calculate RMS using input from oscilloscope
% calc dB w/ zero atten
P = RMS/micSens; % voltage to Pascal
dBSPL = 20*log10(P./20e-6); % Pascal to dBSPL

calib.oscil.Pa = P;
calib.oscil.dBSPL = dBSPL;

%% Calculate RMS using resp
% calc dB w/ zero atten
P2 = rms(resp,2)/(micSens/10e3); % voltage to Pascal
dBSPL2 = 20*log10(P2./20e-6); % Pascal to dBSPL

calib.resp.Pa = P2;
calib.resp.dBSPL = dBSPL2;
calib.attnCh1 = attenCh1;
calib.attnCh2 = attenCh2;
calib.subj = subj;
calib.ear = earname;
% % want levels at probeSPL and suppSPL (inputs to function)
% dropProbe = attenCh1 + (dB(1)-probeSPL)
% dropSupp = attenCh2 + (dB(2)-suppSPL)

fname = strcat(respDir,'Calib_',calibNum, '_' ,subj,'_', earname,'_',datetag, '.mat');

save(fname,'calib');
fprintf(1, 'Saved!\n');


end
