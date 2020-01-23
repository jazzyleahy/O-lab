%Script to run acoustic portion of music encoding study

clear all
rng('shuffle','twister');

%%

%do I need to open labjack?
% lj = labJack('verbose',true);   % open labJack verbosely

%what information do we want to ask here?
p.sub_no = input('Subject ID (e.g. 7): ');
p.run_no = input('Enter run number : ');
p.experiment_ID = 'music_encoding';
p.timestamp = datestr(now, 'yyyymmddTHHMMSS');

if p.sub_no < 10
    p.sub_ID = ['SE0' num2str(p.sub_no) '_run' num2str(p.run_no) '_' p.experiment_ID '_' p.timestamp];
else
    p.sub_ID = ['SE' num2str(p.sub_no) '_run' num2str(p.run_no) '_' p.experiment_ID '_' p.timestamp];
end
diary([p.sub_ID '.txt'])

%p.conds = {};
%p.n_conditions = 1;
%p.num_ex = 5;
%p.n_trials_per_cond = 300; 

%create a queue for example presentation
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
stimDir = dir([pwd '/Stimuli/ENTERSTIMNAME.wav']);
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


