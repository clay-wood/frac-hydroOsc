(* clear all
close all
clc
warning off

% Script to anaylize acoustic data (verasonics device), once the sync has
% been performed.

% Load mechanical data (filename looks like 'pXXXX_data_bin')

runname=   'p5369d';
% addpath(genpath('/gpfs/group/cjm38/default/p5369/'));

%% read binary file (output of r_file)
[data,outname1] = ReadBinBiax(runname);


NormDisp        = data(:,4); 
LPDisp          = data(:,4); 
NormStress      = data(:,9); 
Sync            = data(:,13);
Time            = data(:,14);

clear data

% Filename:  p5369WGFracNSPPOsc 
% Number of records: 2551605
% Number of columns: 13
% Swp : 0
% dtime : 0
% 
% 
% ----------------------------------------------------------------------------
% |Column|             Name|             Unit|          Records|
% ----------------------------------------------------------------------------
% |     1|         Ver_disp|              bit|          2551605|
% |     2|         Ver_load|              bit|          2551605|
% |     3|         hor_disp|               um|          2551605|
% |     4|          Hstress|             MPai|          2551605|
% |     5|               Pc|              MPa|          2551605|
% |     6|          Pc_load|              bit|          2551605|
% |     7|         Ppa_disp|               um|          2551605|
% |     8|         Ppa_load|              bit|          2551605|
% |     9|              Ppb|              MPa|          2551605|
% |    10|         Ppb_load|              bit|          2551605|
% |    11|         Int_disp|               um|          2551605|
% |    12|             Sync|                V|          2551605|
% |    13|             Time|              sec|          2551605|
% ---------------------------------------------------------------------------- *)
(* %% define settings *)
(* acousticrun = 'run23'; % select the acoustic run to analyze after adjusting indexes (from figure 1) and the path

switch acousticrun 
    
    case 'run1' 
               
        AcSettingsfile = '/storage/home/c/cew52/research/p5369/matfiles/p5369_run1.mat'; % acoustic settings file (located in the current directory)
        AcSyncFile = '/storage/home/c/cew52/research/p5369/matfiles/p5369a_sync_run1.mat'; % sync file 
        WF_path = '/gpfs/group/cjm38/default/p5369/acData/run1/WF_'; % where the WFs are                                                 
        
(*         % Portion of the WF to analyze *)
        idx2analyze = zeros(2,3,3); % 1st dim => beg and end; 2nd dim => receivers; 3rd dim => transmitters 
        
        idx2analyze(:,:,1) = [249 374; 249 374; 249 374]'; %T1
        idx2analyze(:,:,2) = [249 374; 249 374; 249 374]'; %T2
        idx2analyze(:,:,3) = [249 374; 249 374; 249 374]'; %T3
        
(*         % Display waveforms sent by transmitter WhichTrans *)
        WhichTrans =1;           
(*         % Time Range *)
        TimeRange = []; % in seconds. Analyze acoustic data over that time range only. If empty, the whole run is analyzed
        AtWhichTimes = [4067 5300 6494];
        NtoStack = 3;
        ref.type = 'absref'; %'absref', 'relref' or 'mixref';
        
     

Filter.yes = 0; % 1 to filter, 0 not to filter
Filter.frq = [0.25 2]; % pass band filter 0.25MHz 2MHz
Filter.order = 256;
Filter.view = 1;

% number of waveforms to stack (either to display or when analyzing)

displayoptions = 1; % choose 0 to display all waveforms or 1 to display one set of waveforms over 100

ref.NrefWF = 50;

% offset waveforms by Offset when multiple channels are used
Offset1 = 10000;
Offset2 = 10000;

% used for 'relref' or 'mixref'
threshold = 0.95;

FreqQ = 400e3; % frequency at which the amplitude is monitored
NZPad = 2048; %2^18; %2048 number of samples for zero padding before fft 

% geometry
% [pos_L,pos_R,dist_s,dist_b] = GeoVesselAcBlocks(26.15,[1:3 5 7:8 10 12],[1:5 7:9],0,0);

%% plot mechanical data of interest within the considered run
load(AcSyncFile);
% Find sample number corresponding to the beginning and end of the acoustic run
FirstIdxAc = find(Time > acTime(1),1,'first'); 
LastIdxAc = find(Time > acTime(end),1,'first');
idxAc = FirstIdxAc:LastIdxAc;

FigRaw = figure;
ax1 = subplot(311);plot(Time(idxAc),NormDisp(idxAc)/1000);ylabel('Hori Disp (mm)');
ax2 = subplot(312);plot(Time(idxAc),LPDisp(idxAc)/1000);ylabel('LP Disp (mm)');hold on
ax3 = subplot(313);plot(Time(idxAc),NormStress(idxAc));ylabel('Normal Stress (MPa)');
xlabel('Time (s)');
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);

linkaxes([ax1,ax2,ax3],'x'); *)

(* %% show WFs at different times
ShowMeWFs_MC(WF_path,AcSettingsfile,AcSyncFile,AtWhichTimes,NtoStack,Offset1,WhichTrans,idx2analyze);
% return *)

(* % Uncomment when choosing the X-corr w  indows
%% process acoustic data (Time Shift, RmsAmp and Max Intercorrelation) *)
[MaxInter,TimeShift,RmsAmp,Amp,RmsAmpRef,AmpRef,fullWFref,LocalAcTime,freqQAmp,maxAmp,maxFreq] = ...
    ProcessAc_Tomo_MC_wAmp(WF_path,AcSettingsfile,AcSyncFile,...
    idx2analyze,ref,NtoStack,threshold,Offset2,displayoptions,Filter,TimeRange,NZPad,FreqQ);

%% save data
% filename of the resulting mat file with explicit name based on chosen paramters
if strcmp(ref.type,'relref') || strcmp(ref.type,'mixref')
    refname = [ref.type '_Th' num2str(threshold)]; % add threshold value to the name when using relative or mixed reference
else
    refname = ref.type;
end

if ~isempty(TimeRange)
   filenamedata = ['Results_' runname '_' acousticrun '_' num2str(LocalAcTime(1,1)) 's-' num2str(LocalAcTime(end,end)) 's_Stack' num2str(NtoStack) 'WFs_' refname '_wAmp.mat'];
else
   filenamedata = ['Results_' runname '_' acousticrun '_fullrun_Stack' num2str(NtoStack) 'WFs_' refname '_wAmp.mat'];    
end
       
save(filenamedata,...
    'LocalAcTime','RmsAmp','Amp','AmpRef','MaxInter','TimeShift',...        
    'RmsAmpRef','fullWFref','idx2analyze','NtoStack','ref','threshold','Filter','TimeRange', 'freqQAmp', 'maxAmp', 'maxFreq','FreqQ','NZPad');
       
return

