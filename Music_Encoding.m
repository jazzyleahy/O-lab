%Script to run acoustic portion of music encoding study

clear all
rng('shuffle','twister');

%%


lj = labJack('verbose',true);   % open labJack verbosely

p.sub_no = input('Subject ID (e.g. 7): '); %O-lab subject number
p.mussub_no = input('Subject ID (e.g. 7): '); %subject number within music encoding experiment
p.run_no = input('Enter run number : ');
p.experiment_ID = 'music_encoding';
p.timestamp = datestr(now, 'yyyymmddTHHMMSS');


idx_ex(p.mussub_no,:,:) %1-20
idx_task() % 1-4 1=fam 2=liking etc

if p.sub_no < 10
    p.sub_ID = ['SE0' num2str(p.sub_no) '_run' num2str(p.run_no) '_' p.experiment_ID '_' p.timestamp];
else
    p.sub_ID = ['SE' num2str(p.sub_no) '_run' num2str(p.run_no) '_' p.experiment_ID '_' p.timestamp];
end
diary([p.sub_ID '.txt'])

p.conds = {};
p.n_conditions = 1;
p.num_ex = 10;
p.n_trials_per_cond = 10; 

%create a queue for example presentation

%create matrix with randomized exemplars with index
p.ex_order = [];
for iSet = 1:(p.n_trials_per_cond/p.num_ex)
    thisSet = randperm(p.num_ex);
    if iSet > 1 %prevent examples from being presented back-to-back
        while thisSet(1,1) == p.ex_order(1,end)
            thisSet = randperm(p.num_ex);
        end
    end
    p.ex_order = [p.ex_order thisSet];
end

%Get stimulus names 
stimDir = dir([pwd '/Stimuli/Fritz/*_**.wav']);
stimFiles = {stimDir.name};
%stimFiles = stimFiles([4,5,1:3]); % reorder to keep numerical ordering

p.stim_dur = 30;
%ask SG what we want to do for onset
p.ISI = 1;
p.trial_onset = 2+[0:p.stim_dur+p.ISI:(p.stim_dur+p.ISI)*(1+p.n_trials_per_cond)];

p.trigger_codes = p.ex_order;

InitializePsychSound
srate = 44100;

PsychPortAudio('Verbosity', 10);
Dev = PsychPortAudio('GetDevices');  % get driver device ID
if strcmp(computer,'PCWIN64')
    pahandle = PsychPortAudio('Open', Dev(4).DeviceIndex, [], 2, srate, 2,[],[]);    % pahandle = PsychPortAudio('Open' [, deviceid][, mode][, reqlatencyclass][, freq][, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);
elseif strcmp(computer,'MACI64')
    for k = 1:length(Dev)
        if strcmp('Fire', Dev(k).DeviceName(1,1:4))
            pahandle = PsychPortAudio('Open', Dev(k).DeviceIndex, [], 2, srate, 2,[],[]);    % pahandle = PsychPortAudio('Open' [, deviceid][, mode][, reqlatencyclass][, freq][, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);
        end
    end
    if ~exist('pahandle','var')
        for k = 1:length(Dev)
            if strcmp('Built-in Output', Dev(k).DeviceName)
                pahandle = PsychPortAudio('Open', Dev(k).DeviceIndex, [], 2, srate, 2,[],[]);    % pahandle = PsychPortAudio('Open' [, deviceid][, mode][, reqlatencyclass][, freq][, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);
            end
        end
    end
end

%% load stimuli

disp('Loading sound files ...')
p.soundfilename = cell(length(p.conds),p.n_trials_per_cond);
cur_rep = ones(1,p.n_conds);


y = cell(1,p.n_trials);

 for k = 1:p.n_trials_per_cond
       c = p.cond_order(k);
      cur_rep(c) = cur_rep(c) + 1;

  p.soundfilename{k} = sprintf('%i_%i.wav', icond, iex);
 %   Lang{iLang}, Cond{iCond}, iEx); how should I alter this to reflect the
 %   conditions of our experiment?
  ls(p.soundfilename{k})
 end

 % load all audio files on memory:
bufferdata = zeros(1, p.n_trials);
for k = 1:p.n_trials
  try
    [y{k}, srate] = audioread(p.soundfilename{k});
  catch
    [y{k}, srate] = wavread(p.soundfilename{k});
    %why is there an error with wavread?
  end
  % matching # of output channels required:
  if size(y{k},2) ~= Dev(audioDeviceIndex).NrOutputChannels
    y{k} = repmat(y{k}, [1 Dev(audioDeviceIndex).NrOutputChannels]);
  end
  % assign buffer #:
  bufferdata(k) = PsychPortAudio('CreateBuffer', pahandle, y{k}'); 
end

% read boundary time and speaker change time:
p.speaker_change_sec = cell(1, p.n_trials);
for k = 1:p.n_trials
  [~,f1,~]=fileparts(p.soundfilename{k});
  load([f1,'.mat'],'SQ')
  p.speaker_change_sec{k} = SQ.speaker_change_sec;
  p.bndry_sec{k} = SQ.bndry_sec;
end

% load & play to redcue latency in initial playbacks
PsychPortAudio('Fillbuffer',pahandle, y{k}');  % load sound into buffer
PsychPortAudio('Volume', pahandle, 0); % muted
for i=1:3
  PsychPortAudio('Start', pahandle, 1, 0, 0);
  WaitSecs(1);
  PsychPortAudio('Stop', pahandle);
end
PsychPortAudio('Volume', pahandle, 1); % unmuted

%% INITIALIZE: timing control functions
% run these now once to avoid delays of first usage later on
priorityLevel = MaxPriority('GetSecs','KbWait','WaitSecs');
GetSecs;
WaitSecs(.01);
% KbWait([],1); % all keys to be released

%% prepare screen; the following are Psychtoolbox functions to display e.g. text on the screen

if iRun == 1 %setup the screen for the first run only
    Screen('Preference', 'SkipSyncTests', 1);
    ScreenNum = 0; 
    [wPtr, rect] = Screen('OpenWindow', ScreenNum, 0);%[0 0 800 500]);%0);%,[0 0 800 500]);%[0 0 1920 1100]); % [0 0 800 500] [0 0 1920 1100]
%     [wPtr, rect] = Screen('OpenWindow', ScreenNum, 0, [0 0 800 500]);%0);%,[0 0 800 500]);%[0 0 1920 1100]); % [0 0 800 500] [0 0 1920 1100]
    black = BlackIndex(wPtr);
    Screen('FillRect', wPtr, black); 
    breakfontsize = Screen('TextSize', wPtr, 25);
    HideCursor(ScreenNum);
    priorityLevel = MaxPriority(ScreenNum);
    Priority(priorityLevel);
end

WaitSecs(1);
GetSecs;
totalBreakTime = 0;
waitForSoundStart = 1;  % wait for sound onset confirmation

DrawFormattedText(wPtr,'Press any key to start the run', 'center', 'center', [255 255 255]);
Screen(wPtr,'flip');
KbWait;
Screen(wPtr,'flip');
WaitSecs(2);
p.Start = GetSecs;
p.starTime(iRun) = p.Start;

for trial = 1:p.n_trials_per_cond
    
    PsychPortAudio('Fillbuffer',pahandle,bufferdata(trial));  % load sound into buffer
    
    if trial < (p.n_trials_per_cond/4)+1
        WaitSecs('UntilTime',p.Start+p.trial_onset(1,trial));
        %do we want to add a break?
    elseif trial == (p.n_trials_per_cond/4)+1
        p.Pause1Start = GetSecs;
        DrawFormattedText(wPtr,'Take a Break - press any key to continue', 'center', 'center', [255 255 255]);
        Screen(wPtr,'Flip');
        KbWait;
        Screen(wPtr,'flip')
        WaitSecs(2);
        p.StartQ2 = GetSecs;
        totalBreakTime = p.StartQ2 - p.Pause1Start;
        WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial));
    elseif trial > (p.n_trials_per_cond/4)+1 && trial < (p.n_trials_per_cond/2)+1
        WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial));
    elseif trial == (p.n_trials_per_cond/2)+1
        p.Pause2Start = GetSecs;
        DrawFormattedText(wPtr,'Take a Break - press any key to continue', 'center', 'center', [255 255 255]);
        Screen(wPtr,'Flip');
        KbWait;
        Screen(wPtr,'flip')
        WaitSecs(2);
        p.StartQ3 = GetSecs;
        totalBreakTime = totalBreakTime + p.StartQ3 - p.Pause2Start;
        WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial));
    elseif trial > (p.n_trials_per_cond/2)+1 && trial < (3*p.n_trials_per_cond/4)+1
        WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial));
    elseif trial == (3*p.n_trials_per_cond/4)+1
        p.Pause3Start = GetSecs;
        DrawFormattedText(wPtr,'Take a Break - press any key to continue', 'center', 'center', [255 255 255]);
        Screen(wPtr,'Flip');
        KbWait;
        Screen(wPtr,'flip')
        WaitSecs(2);
        p.StartQ4 = GetSecs;
        totalBreakTime = totalBreakTime + p.StartQ4 - p.Pause3Start;
        WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial));
    elseif trial > (3*p.n_trials_per_cond/4)+1
        WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial));
    end

    DrawFormattedText(wPtr, ['Trial ' num2str(trial)], 'center', 'center', [255 255 255]);
    Screen(wPtr,'Flip');
%     p.sound_on(1,trial) = GetSecs;    % ONLY for testing without sounds
    p.sound_on(iRun,trial) = PsychPortAudio('Start', pahandle, 1, 0, waitForSoundStart);     % waitForSoundStart = 1

    % send sound condition trigger
    lj.setDIO([p.trigger_codes(iRun,trial),0,0]);
    WaitSecs(.01);
    lj.setDIO([0,0,0]);
    
%     disp(p.soundfilename{trial})

    [secs, keyCode, deltaSecs] = KbWait([],2,p.sound_on(iRun,trial)+p.stim_dur);
    
    save(p.sub_ID,'p')

    
    DrawFormattedText(wPtr,feedbackTxt, 'center', 'center', [255 255 255]);
    WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial)+p.stim_dur);
    Screen(wPtr,'Flip');
    
    WaitSecs('UntilTime',p.Start+totalBreakTime+p.trial_onset(1,trial)+p.stim_dur+0.5);
    Screen(wPtr,'Flip');
            
end
end

Screen('CloseAll');
PsychPortAudio('Close');
Priority(0);
close(lj);




