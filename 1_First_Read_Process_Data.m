

clear all
close all
clc

% Choose the right subject directory for datadir variable

% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_4\SMD001'
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_4\SMD002'
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_4\SMD003';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_4\V035';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_4\SMD004\Day1';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_4\SMD004\Day2';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_4\SMD005';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD006';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD007';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD008';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD009\Day1';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD009\Day2';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD010\Day1';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD010\Day2';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD011\Day1';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD011\Day2';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD012';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD013';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD014';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\000\Fast';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5_b\SMD015\faster';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5_b\SMD016\faster';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD020';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD021';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD022';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD023';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD024';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD025';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD026';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD027';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\SMD028';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD029';
datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\SMD030';

% path(path,'C:\GoogleDrive\Synced Folder\Matlab\Data_Processing\NoForce4Targets')%Matlab function directory
path(path,'C:\GoogleDrive\Synced Folder\Matlab\Data_Processing\SMD_2Targ')

ExpData = zip_load('dir',datadir);  %creating the structure

%% extra steps

% subID = '000';
% load([subID,'_Processed'])

%% save unzipped data and save them

% % how many protocol blocks (number of zip files)        !!CHANGE IF NEEDED!!
prot_block = 1;

fs = 1000; % sampling rate Hz


tot_num_trials = 0;
for r=1:prot_block
    tot_num_trials= tot_num_trials + length(ExpData(r).c3d);
end

% Creating RawData cell
RawData.time_FS = cell(1,tot_num_trials);
RawData.time_trial = cell(1,tot_num_trials);
% RawData.time_fixed=cell(1,tot_num_trials);
% RawData.time_global=cell(1,tot_num_trials);
RawData.X = cell(1,tot_num_trials);
RawData.Y = cell(1,tot_num_trials);
RawData.VX = cell(1,tot_num_trials);
RawData.VY = cell(1,tot_num_trials);
RawData.Speed = cell(1,tot_num_trials);
RawData.Fx = cell(1,tot_num_trials);
RawData.Fy = cell(1,tot_num_trials);
% RawData.Theta1 = cell(1,tot_num_trials);
% RawData.Theta2 = cell(1,tot_num_trials);
RawData.reach_type = zeros(1,tot_num_trials);
RawData.target_num = zeros(1,tot_num_trials);
RawData.target_position = zeros(2,tot_num_trials); % 1=Xposition 2=Yposition
RawData.center_position = zeros(2,tot_num_trials); % 1=Xposition 2=Yposition
RawData.cursor_feedback = cell(1,tot_num_trials);
RawData.logical_radius = zeros(1,tot_num_trials);
RawData.score = zeros(1,tot_num_trials);     % Add score
RawData.event = zeros(1,tot_num_trials);
RawData.fs = fs;

% reading the data of the subject
fid = fopen('pat.dat','r');
Sub_Data = textscan(fid, '%s');
take_ID= strsplit(Sub_Data{1}{11},'\s*id=\s*','DelimiterType','RegularExpression');
subID=take_ID{end};


score_build = zeros(1,tot_num_trials); % to build to score vector

% Defining RawData values in the right trial order
for j=1:prot_block;
    
    for counter = 1:length(ExpData(j).c3d)
        
        if j == 1
            trial = ExpData(1).c3d(counter).TRIAL.TRIAL_NUM;
        end
        if j==2
            trial = ExpData(2).c3d(counter).TRIAL.TRIAL_NUM + length(ExpData(1).c3d);
        end
        if j==3
            trial = ExpData(3).c3d(counter).TRIAL.TRIAL_NUM + length(ExpData(1).c3d)+ length(ExpData(2).c3d);
        end
        if j==4
            trial = ExpData(4).c3d(counter).TRIAL.TRIAL_NUM + length(ExpData(1).c3d)+...
                length(ExpData(2).c3d)+ length(ExpData(3).c3d);
        end
        if j==5
            trial = ExpData(5).c3d(counter).TRIAL.TRIAL_NUM + length(ExpData(1).c3d)+ ...
                length(ExpData(2).c3d)+ length(ExpData(3).c3d)+ length(ExpData(4).c3d);
        end
        if j==6
            trial = ExpData(6).c3d(counter).TRIAL.TRIAL_NUM + length(ExpData(1).c3d)+ ...
                length(ExpData(2).c3d)+ length(ExpData(3).c3d)+ length(ExpData(4).c3d)+ length(ExpData(5).c3d);
        end
        
        
        RawData.X{trial} = ExpData(j).c3d(counter).Right_HandX;
        RawData.Y{trial} = ExpData(j).c3d(counter).Right_HandY - 0.15;
        RawData.Fx{trial} = ExpData(j).c3d(counter).Right_FS_ForceX;
        RawData.Fy{trial} = ExpData(j).c3d(counter).Right_FS_ForceY;
        %         RawData.Theta1{trial} = ExpData(j).c3d(counter).Theta1;
        %         RawData.Theta2{trial} = ExpData(j).c3d(counter).Theta2;
        
        % absolute time since experiment started
        time_FS = ExpData(j).c3d(counter).Right_FS_TimeStamp;   %Time from the Robot
        num_data_pts = length(RawData.X{trial});
        time = ((0:1:num_data_pts-1)/fs)';                    %Right time definition
        RawData.time_trial{trial} = time;
        RawData.time_FS{trial} = time_FS;
        
        
        Type = ExpData(j).c3d(counter).TRIAL.TP;
        RawData.target_num(trial) = mod(Type,4); %!!! CHANGE FOR NUM TARGETS !!!%this line defins the target number bassed on the tp row
        if mod(Type,4)==0
            RawData.target_num(trial) = 4;
        end
        
        %         Type = ExpData(j).c3d(counter).TRIAL.TP;
        %         RawData.target_num(trial) = mod(Type,2); %!!! CHANGE FOR NUM TARGETS !!!%this line defins the target number bassed on the tp row
        %         if mod(Type,4)==0
        %             RawData.target_num(trial) = 2;
        %         end
        
        %Cursor feedback info
        TP_Row=ExpData(j).c3d(counter).TRIAL.TP;
        RawData.cursor_feedback{trial} = ...
            ExpData(j).c3d(counter).TP_TABLE.Vision_or_Novision(ExpData(j).c3d(counter).TP_TABLE.Load(TP_Row));
        
        %         %target position info
        %         if RawData.target_num(trial) == 2 %for the second reach
        %             RawData.target_position(1,trial) =...
        %                 ExpData(j).c3d(counter).TARGET_TABLE.X(1)/100;
        %             RawData.target_position(2,trial) =...
        %                 ExpData(j).c3d(counter).TARGET_TABLE.Y(1)/100;
        %             RawData.logical_radius(1,trial) =...
        %                 ExpData(j).c3d(counter).TARGET_TABLE.Logical_Radius(1)/100;
        %         else
        %             RawData.target_position(1,trial) =...
        %                 ExpData(j).c3d(counter).TARGET_TABLE.X(2)/100;
        %             RawData.target_position(2,trial) =...
        %                 ExpData(j).c3d(counter).TARGET_TABLE.Y(2)/100;
        %             RawData.logical_radius(1,trial) =...
        %                 ExpData(j).c3d(counter).TARGET_TABLE.Logical_Radius(2)/100;
        %         end
        
        %target position info
        if RawData.target_num(trial) == 3
            RawData.target_position(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.X(8)/100;
            RawData.target_position(2,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Y(8)/100;
            RawData.logical_radius(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Logical_Radius(8)/100;
        elseif RawData.target_num(trial) == 4
            RawData.target_position(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.X(1)/100;
            RawData.target_position(2,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Y(1)/100;
            RawData.logical_radius(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Logical_Radius(1)/100;
        else
            RawData.target_position(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.X(RawData.target_num(trial)+1)/100;
            RawData.target_position(2,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Y(RawData.target_num(trial)+1)/100;
            RawData.logical_radius(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Logical_Radius(RawData.target_num(trial)+1)/100;
        end
        
        % Score info
        %         if isempty(find(strncmp(ExpData(j).c3d(counter).EVENTS.LABELS, 'Score event', 11)))
        if sum(~cellfun(@isempty, strfind(ExpData(j).c3d(counter).EVENTS.LABELS, 'Good')))
            score_build(trial)=1;
        else
            score_build(trial)=0;
        end
        
        % Event "Stopped moving"
        EventLabel = strfind(ExpData(j).c3d(counter).EVENTS.LABELS,'Stopped Moving');
        idx = find(~cellfun(@isempty,EventLabel));
        if length(idx)>1
            idx = idx(length(idx)-1);
        end
        Event_time = ExpData(j).c3d(counter).EVENTS.TIMES(idx);
        if isempty(idx)
            RawData.event(1,trial) = 0; %this trial is going to be thrown away. hand nerver came to a stop
            trial
            display('timed-out!')
        else
            RawData.event(1,trial) = Event_time;
        end
        
        
        %reach type
        RawData.reach_type(trial) = RawData.target_num(trial);
        
    end
    
end

% define score/ center position
RawData.score(1) = score_build(1);
for trial=2:tot_num_trials
    RawData.score(trial) = RawData.score(trial-1) + score_build(trial);
    RawData.center_position(:,trial) = RawData.target_position(:,trial-1);
end
%%
savefile = [subID,'_Raw'];
save(fullfile(datadir,savefile),'RawData')        %saving data


%% Process Data (New Version)

%Threshholds
vel_thresh_out=0.05;
vel_thresh_back=0.05;

% compute hand velocity & filter hand velocity
tot_num_trials = length(RawData.X);

figSpeed= figure;
velid_stop= zeros(1,tot_num_trials);
bad_index = [];


disp('Instructions:');
disp('-"d" to start');
disp('-Enter to exit');
disp('-"m" to maximize');
disp('-"a" to move back one trial');
disp('-"d" to move forward one trial');
disp('-"w" to mark as good');
disp('-"s" to mark as bad');

% trial = 0; %-1;
trial = 88;

while trial < tot_num_trials
    % for trial = [1:49, 51:tot_num_trials]
    
    ch = getkey();
    if ch == 109 % 'm'
        drawnow;
        set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    elseif ch == 13 % enter
        trial = tot_num_trials; % force end of while loop
    elseif ch == 119 % 'w'
        % if it exists in bad index find the index and remove
        % if it doesn't exist, do nothing
        disp(trial);
        if any(bad_index == trial)
            toRemove = find(bad_index == trial); % find index of trial to be removed
            bad_index(toRemove) = [];
        end
    elseif ch == 115 % 114 for 'r' or 115 for 's'
        % if it exists in bad index do nothing
        % if it doesn't exist, add it
        if ~any(bad_index == trial)
            bad_index = [bad_index trial]; % append bad trial to bad_index
        end
    elseif ch == 97  & trial > 1 % 'a'
        trial = trial - 1;
        time = RawData.time_trial{trial};
        if trial == 1
            time_global = time;
        else
            time_global = time + RawData.time_trial{trial-1}(end) + 1/fs;
        end
    elseif ch == 100 & trial < tot_num_trials % (115 for 's', 100 for 'd')
        trial = trial + 1;
        time = RawData.time_trial{trial};
        if trial == 1
            time_global = time;
        else
            time_global = time + RawData.time_trial{trial-1}(end) + 1/fs;
        end
    end
    
    % Defining the velocity
    
    V_X = diff(RawData.X{trial})*fs;
    V_X = [V_X; V_X(end)];
    V_Y = diff(RawData.Y{trial})*fs;
    V_Y = [V_Y; V_Y(end)];
    
    % Filtering the velocity
    
    omega = 20; % cuttoff Hz
    Wn = omega/fs/2;
    [c_b, c_a] = butter (5, Wn);
    
    RawData.VX{trial} = filtfilt (c_b, c_a, V_X);
    RawData.VY{trial} = filtfilt (c_b, c_a, V_Y);
    
    % Defining the Speed (magnitude of velocity)
    RawData.Speed{trial} = sqrt(RawData.VX{trial}.^2 + RawData.VY{trial}.^2);
    RawData.Magn{trial} = sqrt(RawData.Fx{trial}.^2 + RawData.Fy{trial}.^2);
    
    % check for velocity threshold crossing
    velid_start = find(RawData.Speed{trial}>vel_thresh_out,1);
    
    %debugging for bad trials
    if isempty(velid_start)
        velid_start = 1;
    end
    
    %     %find the second point that satisfices the threshhold conditions
    tempid = find(RawData.Speed{trial}((velid_start+1):end)<vel_thresh_back,1,'first')+velid_start;
    %     first_stopid(i) = tempid;
    %     tempid = find(RawData.Speed{trial}((tempid+5):end)>vel_thresh_back,1,'first')+tempid;
    %     tempid = find(RawData.Speed{trial}((tempid+5):end)<vel_thresh_back,1,'first')+tempid;
    
    velid_stop(trial)=tempid;
    
    temptime = RawData.time_trial{trial} - RawData.time_trial{trial}(velid_start);
    
    RawData.time_fixed{trial} = temptime;
    RawData.time_global{trial} = time_global;
    
    % Plotting the velocity threshold to collect the bad trials
    
    figure(figSpeed)
    clf
    %     subplot(3,1,1)
    subplot(2,1,1)
    hold on;
    plot(temptime,RawData.Speed{trial},'b')
    plot(temptime(velid_start),RawData.Speed{trial}(velid_start),'rx',...
        temptime(tempid),RawData.Speed{trial}(tempid),'ko')
    
    line([temptime(1) temptime(end)],[vel_thresh_out vel_thresh_out],'Color','r')
    line([temptime(1) temptime(end)],[vel_thresh_back vel_thresh_back],'Color','k')
    title([' Trial : ',num2str(trial)])
    grid
    
    subplot(2,1,2)
    hold on
    plot(RawData.X{trial},RawData.Y{trial},'r-')
    plot(RawData.X{trial}(velid_start:tempid),RawData.Y{trial}(velid_start:tempid),'k')
    plot(RawData.target_position(1,trial),RawData.target_position(2,trial),'kd')
    grid
    axis equal
    
    
    title([' Trial : ',num2str(trial)])
    disp(bad_index);
end
RawData.velid_stop = velid_stop;
RawData.bad_index = bad_index;
% RawData.first_stopid = first_stopid;

%%
savefile = [subID,'_Raw'];
save(fullfile(datadir,savefile),'RawData')


%% NOW FIX BAD TRIALS ALLIGNMENT


bad_index = RawData.bad_index;
% bad_index = [12];

figSpeed = figure;

% vel_thresh_out=0.05;
% vel_thresh_back=0.02;

for trial = bad_index
    
    % get this trial's data
    
    % check for velocity threshold crossing
    velid_start = find(RawData.Speed{trial}>vel_thresh_out,1);
    
    %     %debugging
    %     if isempty(velid_start)
    %         velid_start = 50;
    %     end
    
    temptime = RawData.time_trial{trial} - RawData.time_trial{trial}(velid_start);
    
    % make sure everything is synced well
    clf
    figure(figSpeed)
    subplot(3,1,2)
    hold on;
    plot(temptime,RawData.Y{trial},'b')
    plot(temptime(velid_start),RawData.Y{trial}(velid_start),'rx')
    plot(temptime(RawData.velid_stop(trial)),RawData.Y{trial}(RawData.velid_stop(trial)),'rx')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial) RawData.target_position(2,trial)],'Color','b')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial)-RawData.logical_radius(trial) RawData.target_position(2,trial)-RawData.logical_radius(trial)],'Color','r','LineStyle','--')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial)+RawData.logical_radius(trial) RawData.target_position(2,trial)+RawData.logical_radius(trial)],'Color','r','LineStyle','--')
    
    subplot(3,1,3)
    hold on
    plot(temptime,RawData.X{trial},'k')
    plot(temptime(velid_start),RawData.X{trial}(velid_start),'rx')
    plot(temptime(RawData.velid_stop(trial)),RawData.X{trial}(RawData.velid_stop(trial)),'rx')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial) RawData.target_position(1,trial)],'Color','k')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial)-RawData.logical_radius(trial) RawData.target_position(1,trial)-RawData.logical_radius(trial)],'Color','r','LineStyle','--')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial)+RawData.logical_radius(trial) RawData.target_position(1,trial)+RawData.logical_radius(trial)],'Color','r','LineStyle','--')
    
    subplot(3,1,1)
    hold on;
    plot(temptime,RawData.Speed{trial},'b')
    plot(temptime(velid_start),RawData.Speed{trial}(velid_start),'rx')
    line([temptime(1) temptime(end)],[vel_thresh_out vel_thresh_out])
    vel_stop_thresh = vel_thresh_back;
    line([temptime(1) temptime(end)],[vel_stop_thresh vel_stop_thresh],'Color','r')
    title([' Trial : ',num2str(trial)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %     subplot(1,2,2)
    %     hold on
    %     tempid = RawData.velid_stop(trial);
    %     plot(RawData.X{trial},RawData.Y{trial},'r-')
    %     plot(RawData.X{trial}(velid_start:tempid),RawData.Y{trial}(velid_start:tempid),'k')
    %     plot(RawData.target_position(1,trial),RawData.target_position(2,trial),'kd')
    %     grid
    %     axis equal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get new window for the start point
    [twindow, temp] = ginput(2);
    idx1 = find(temptime>twindow(1),1);
    idx2 = find(temptime>twindow(2),1);
    subplot(3,1,1)
    plot(temptime(idx2),RawData.Speed{trial}(idx2),'kx')
    plot(temptime(idx1),RawData.Speed{trial}(idx1),'kx')
    velid_start = find(RawData.Speed{trial}(idx1:idx2)>vel_thresh_out,1);
    velid_start = velid_start + idx1 - 1;
    plot(temptime(velid_start),RawData.Speed{trial}(velid_start),'ro')
    
    % get new window for the stop point
    [twindow, temp] = ginput(2);
    idx3 = find(temptime>twindow(1),1);
    idx4 = find(temptime>twindow(2),1);
    subplot(3,1,1)
    plot(temptime(idx4),RawData.Speed{trial}(idx4),'bx')
    plot(temptime(idx3),RawData.Speed{trial}(idx3),'bx')
    
    %save new stop index
    velid_stop_index = find(RawData.Speed{trial}(idx3:idx4)<vel_stop_thresh,1,'first')+idx3-1;
    velid_stop(trial)= velid_stop_index;
    
    plot(temptime(velid_stop(trial)),RawData.Speed{trial}(velid_stop(trial)),'ro')
    subplot(3,1,2)
    plot(temptime(velid_stop(trial)),RawData.Y{trial}(velid_stop(trial)),'ro')
    subplot(3,1,3)
    plot(temptime(velid_stop(trial)),RawData.X{trial}(velid_stop(trial)),'ro')
    
    %define again temptime and save it as time_fixed
    temptime = RawData.time_trial{trial} - RawData.time_trial{trial}(velid_start);
    RawData.time_fixed{trial}= temptime;
    
    pause
    
    
end

RawData.velid_stop = velid_stop;
%%
savefile = [subID,'_Raw'];
save(fullfile(datadir,savefile),'RawData')     %saving data

%% COLLATE, STORE GOOD DATA

% good_index = [1:43, 45:100]; %V035
% good_index = [1:11, 13:100]; % SMD004
% good_index = [1:100]; % SMD005
% good_index = [1:120]; % SMD010
% good_index = [1:50, 52:120]; % SMD011
% good_index = [1:57, 59:120]; %SMD007
good_index = [1:87, 89:120]; %SMD030
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_good_Trials = length(good_index);

Trial.Good_Indices = good_index;
Trial.XPosition = cell(1,num_good_Trials);
Trial.YPosition = cell(1,num_good_Trials);
Trial.Speed = cell(1,num_good_Trials);
Trial.XVelocity = cell(1,num_good_Trials);
Trial.YVelocity = cell(1,num_good_Trials);
Trial.XForce = cell(1,num_good_Trials);
Trial.YForce = cell(1,num_good_Trials);
% Trial.Theta1 = cell(1,num_good_Trials);
% Trial.Theta2 = cell(1,num_good_Trials);
Trial.Time_Trial = cell(1,num_good_Trials);
Trial.Time_Fixed = cell(1,num_good_Trials);
Trial.Time_Global = cell(1,num_good_Trials);
Trial.Start_Out_Index = zeros(1,num_good_Trials);
Trial.End_Out_Index = zeros(1,num_good_Trials);
% Trial.Subject_data = cell(1,3);

Trial.Reach_Type = cell(2,num_good_Trials);
Trial.Target_num = zeros(1,num_good_Trials);
Trial.Target_position = zeros(2,num_good_Trials);
Trial.Center_position = zeros(2,num_good_Trials);
Trial.Cursor_feedback = zeros(1,num_good_Trials);
Trial.score = zeros(1,num_good_Trials);
Trial.Reach_Type = zeros(1,num_good_Trials);
Trial.Subject_data = cell(1,5);

%saving the subject's sex/ID/DOB
Trial.Subject_data=Sub_Data{1}(10:end);

%Event "Stopped Moving"
Trial.Event_Idx = zeros(1,num_good_Trials);

% saving the subject's Height/Weight/sex/ID/Birthday
Trial.Subject_data{1}=Sub_Data{1}(3);
Trial.Subject_data{2}=Sub_Data{1}(5);
Trial.Subject_data{3}=Sub_Data{1}(10);
Trial.Subject_data{4}=subID;
Trial.Subject_data{5}=Sub_Data{1}(12);

% counter = 1;
%
% for i = 1:length(good_index)
%
%     Trial.XPosition{counter} = RawData.X{good_index(i)};
%     Trial.YPosition{counter} = RawData.Y{good_index(i)};
%     Trial.Speed{counter} = RawData.Speed{good_index(i)};
%     Trial.XVelocity{counter} = RawData.VX{good_index(i)};
%     Trial.YVelocity{counter} = RawData.VY{good_index(i)};
%     Trial.XForce{counter} =RawData.Fx{good_index(i)};
%     Trial.YForce{counter} =RawData.Fy{good_index(i)};
% %     Trial.Theta1{counter} =RawData.Theta1{good_index(i)};
% %     Trial.Theta2{counter} =RawData.Theta2{good_index(i)};
%     Trial.Time_Trial{counter} = RawData.time_trial{good_index(i)};
%     Trial.Time_Fixed{counter} = RawData.time_fixed{good_index(i)};
%     Trial.Time_Global{counter} = RawData.time_global{good_index(i)};
%     Trial.Start_Out_Index(counter) = find(Trial.Time_Fixed{counter}>0,1,'first')-50;
%     Trial.End_Out_Index(counter) = RawData.velid_stop(good_index(i))+50;
%     Trial.Reach_Type(counter) = RawData.reach_type(good_index(i));
%
%     if  RawData.cursor_feedback{good_index(i)} == 1
%         Trial.Cursor_feedback(counter) = false;
%     elseif RawData.cursor_feedback{good_index(i)} == 0
%         Trial.Cursor_feedback(counter) = true;
%     end
%
%     Trial.Target_num(counter) = RawData.target_num(good_index(i));
%     Trial.Target_position(1,counter) = RawData.target_position(1,good_index(i));
%     Trial.Target_position(2,counter) = RawData.target_position(2,good_index(i));
%     Trial.Center_position(1,counter) = RawData.center_position(1,good_index(i));
%     Trial.Center_position(2,counter) = RawData.center_position(2,good_index(i));
%     Trial.score(counter) = RawData.score(good_index(i));
%     Trial.Event_Idx(counter) = find(RawData.time_trial{good_index(i)}>=RawData.event(good_index(i)),1,'first');
%
%     counter = counter +1;
% end

% counter = 1;

for i = 1:length(good_index)
    
    Trial.XPosition{i} = RawData.X{good_index(i)};
    Trial.YPosition{i} = RawData.Y{good_index(i)};
    Trial.Speed{i} = RawData.Speed{good_index(i)};
    Trial.XVelocity{i} = RawData.VX{good_index(i)};
    Trial.YVelocity{i} = RawData.VY{good_index(i)};
    Trial.XForce{i} =RawData.Fx{good_index(i)};
    Trial.YForce{i} =RawData.Fy{good_index(i)};
    %     Trial.Theta1{counter} =RawData.Theta1{good_index(i)};
    %     Trial.Theta2{counter} =RawData.Theta2{good_index(i)};
    Trial.Time_Trial{i} = RawData.time_trial{good_index(i)};
    Trial.Time_Fixed{i} = RawData.time_fixed{good_index(i)};
    Trial.Time_Global{i} = RawData.time_global{good_index(i)};
    Trial.Start_Out_Index(i) = find(Trial.Time_Fixed{i}>0,1,'first')-50;
    Trial.End_Out_Index(i) = RawData.velid_stop(good_index(i))+50;
    Trial.Reach_Type(i) = RawData.reach_type(good_index(i));
    
    if  RawData.cursor_feedback{good_index(i)} == 1
        Trial.Cursor_feedback(i) = false;
    elseif RawData.cursor_feedback{good_index(i)} == 0
        Trial.Cursor_feedback(i) = true;
    end
    
    Trial.Target_num(i) = RawData.target_num(good_index(i));
    Trial.Target_position(1,i) = RawData.target_position(1,good_index(i));
    Trial.Target_position(2,i) = RawData.target_position(2,good_index(i));
    Trial.Center_position(1,i) = RawData.center_position(1,good_index(i));
    Trial.Center_position(2,i) = RawData.center_position(2,good_index(i));
    Trial.score(i) = RawData.score(good_index(i));
    Trial.Event_Idx(i) = find(RawData.time_trial{good_index(i)}>=RawData.event(good_index(i)),1,'first');
    
    %     counter = counter +1;
end

%%
savefile = [subID,'_Processed'];    %saving data
save(savefile,'Trial')                %saving data

%% compute perpendicular error, angular error and normalized path length
%

num_good_Trials= length(Trial.Good_Indices);

PerpError = zeros(1,num_good_Trials);
NormLength = zeros(1,num_good_Trials);
AbsPerpError = zeros(1,num_good_Trials);
EArea = zeros(1,num_good_Trials);
AngError = zeros(1,num_good_Trials);
peakF = zeros(1,num_good_Trials);
PeakSpeed = zeros (1, num_good_Trials);

for trial = 1:num_good_Trials % check for corecction of the last trial.
    
    num_data_pts = length(Trial.XPosition{trial});
    CenterPos = [Trial.Center_position(1,trial) Trial.Center_position(2,trial)];% x-y
    TargetPos = [Trial.Target_position(1,trial) Trial.Target_position(2,trial)];% x-y
    
    TargetDis = TargetPos-CenterPos;
    
    % rotate trajectory along axis of target direction
    reach_indices=Trial.Start_Out_Index(trial):Trial.End_Out_Index(trial);
    
    StartPos=[Trial.XPosition{trial}(reach_indices(1)),Trial.YPosition{trial}(reach_indices(1))];
    StopPos=[Trial.XPosition{trial}(reach_indices(end)),Trial.YPosition{trial}(reach_indices(end))];
    TargetLength = sqrt(sum((StartPos - StopPos).^2));
    
    %     targ_Ang = atan2(TargetDis(2),TargetDis(1));
    %
    %     HandDis = StopPos - CenterPos;
    %     hand_Ang = atan2(HandDis(2), HandDis(1));
    %     R = [cos(hand_Ang) -sin(hand_Ang); sin(hand_Ang) cos(hand_Ang)];
    %
    %     Xr = [Trial.XPosition{trial}-CenterPos(1)*ones(num_data_pts,1), ...
    %         Trial.YPosition{trial}-CenterPos(2)*ones(num_data_pts,1)]*R;
    
    %get peak speed
    PeakSpeed(trial) = max(Trial.Speed{trial}(reach_indices));
    
    % Perpendicular error version 2
    EndPos = [Trial.XPosition{trial}(Trial.Event_Idx(trial)), ...
        Trial.YPosition{trial}(Trial.Event_Idx(trial))];% x-y
    Trial.StopMove_Position(:,trial) = EndPos;
    
    % rotate trajectory along axis of target direction
    targ_Ang = atan2(TargetDis(2),TargetDis(1));
    R = [cos(targ_Ang) -sin(targ_Ang); sin(targ_Ang) cos(targ_Ang)];
    Xr2 = [Trial.XPosition{trial}-CenterPos(1)*ones(num_data_pts,1), ...
        Trial.YPosition{trial}-CenterPos(2)*ones(num_data_pts,1)]*R;
    
    % compute PERPENDICULAR ERROR version 2
    [AbsPerpError(trial), idx_Abs] = max(abs(Xr2(reach_indices,2)));
    idx_Abs = idx_Abs + reach_indices(1)-1;
    PerpError(trial) = Xr2(idx_Abs,2);
    
    [~, idx_Pos] = max(Xr2(reach_indices,2));
    idx_Pos = idx_Pos + reach_indices(1)-1;
    PerpError_Pos(trial) = Xr2(idx_Pos,2); %posotive values
    
    [~, idx_Neg] = min(Xr2(reach_indices,2));
    idx_Neg = idx_Neg + reach_indices(1)-1;
    PerpError_Neg(trial) = Xr2(idx_Neg,2); %negative values
    
    % compute normalized path length
    temp = 0;
    for i = reach_indices
        dx = Trial.XPosition{trial}(i+1) - Trial.XPosition{trial}(i);
        dy = Trial.YPosition{trial}(i+1) - Trial.YPosition{trial}(i);
        temp = temp + sqrt(dx^2 + dy^2);
    end
    NormLength(trial) = temp/TargetLength;
    
    % Enclosed Area
    for i = reach_indices
        EArea(trial) = EArea(trial) + abs((Xr2(i+1, 1) - Xr2(i,1)) * Xr2(i,2));
    end
    
    % compute angular error
    [temp, indx] = max(Trial.Speed{trial}(reach_indices));
    xm = Trial.XPosition{trial}(indx+reach_indices(1)-1);
    ym = Trial.YPosition{trial}(indx+reach_indices(1)-1);
    
    ang = atan2(ym-CenterPos(2), xm-CenterPos(1));
    temperr = ang-targ_Ang;
    % wrap ang errors
    if temperr > pi
        temperr = temperr - 2*pi;
    elseif temperr < -pi
        temperr = temperr + 2*pi;
    end
    
    AngError(trial) = temperr*180/pi; %change to degree
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get peak force & speed (avoid the dimple forces)
    
    dis2target= sqrt( (Trial.XPosition{trial}-TargetPos(1)*ones(num_data_pts,1)).^2 + ...
        (Trial.YPosition{trial}-TargetPos(2)*ones(num_data_pts,1)).^2 );
    dis2center= sqrt( (Trial.XPosition{trial}-CenterPos(1)*ones(num_data_pts,1)).^2 + ...
        (Trial.YPosition{trial}-CenterPos(2)*ones(num_data_pts,1)).^2 );
    
    idx = find(dis2center>=0.015,1,'first');    %after the dimple
    idxx = find(dis2target<=0.015,1,'first');   %after the dimple
    
    
    if isempty(idx)
        idx = Trial.Start_Out_Index(trial);
    elseif isempty(idxx)
        idxx = Trial.End_Out_Index(trial);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Fx = Trial.XForce{trial}(idx:idxx);
    Fy = Trial.YForce{trial}(idx:idxx);
    F = sqrt(Fx.^2 + Fy.^2);
    peakF(trial) = max(F);
    
    %covered area
    %why not use the rotated trajectory?
    Area = zeros(1,num_data_pts-1);
    for i = 1:num_data_pts-1
        ddy(i) = abs(Trial.YPosition{trial}(i+1)-Trial.YPosition{trial}(i));
        x = abs(Trial.XPosition{trial}(i));
        Area(i) = ddy(i)*x;
    end
    Trial.Area(trial) = sum(Area);
    
end

Trial.PeakForce = peakF;
Trial.AbsPerpError = AbsPerpError;
Trial.PerpError = PerpError;
Trial.PerpError_Pos = PerpError_Pos;
Trial.PerpError_Neg = PerpError_Neg;
Trial.AngError = AngError;
Trial.NormLength = NormLength;
Trial.EArea = EArea;
Trial.PeakSpeed = PeakSpeed;

%%
savefile = [subID,'_Processed'];    %saving data
save(savefile,'Trial')                %saving data

%% extention

clear all
close all
clc

subID = 'SMD028';
load([subID,'_Processed'])


%% Agerage Blocks

close all
%Define the Blocks

% Blocks = {[1:2:100], [2:2:100]};
Blocks = {[1:4:120], [2:4:120], [3:4:120], [4:4:120]};

for i = 1:length(Blocks)
    [a, Blocks{i}] = ismember(Blocks{i}, Trial.Good_Indices);
    Blocks{i} = Blocks{i}(Blocks{i}>=1);
    
end

%% Blocks

BlockAvg.Block_trial_num = cell(1,length(Blocks));

for i=1:length(Blocks)
    
    id = Blocks{i};
    trials= Blocks{i};
    
    success = [0,diff(Trial.score)]; %only successfull trials
    
    BlockAvg.Block_trial_num{i} = [];
    
    for j = trials
        if success(j)
            
            %store the block trial numbers
            BlockAvg.Block_trial_num{i} = [BlockAvg.Block_trial_num{i}, j];%id
            
            % compute and store averages, std's
            
            %Perpendicular error for stopped moving event
            BlockAvg.AbsPerpMean(i) = mean(Trial.AbsPerpError(id));
            BlockAvg.AbsPerpStd(i) = std(Trial.AbsPerpError(id));
            %Perpendicular error for stopped moving event
            BlockAvg.PerpMean(i) = mean(Trial.PerpError(id));
            BlockAvg.PerpStd(i) = std(Trial.PerpError(id));
            %Perpendicular error - negetive values
            BlockAvg.PerpNegMean(i) = mean(Trial.PerpError_Neg(id));
            BlockAvg.PerpNegStd(i) = std(Trial.PerpError_Neg(id));
            %Perpendicular error - positive values
            BlockAvg.PerpPosMean(i) = mean(Trial.PerpError_Pos(id));
            BlockAvg.PerpPosStd(i) = std(Trial.PerpError_Pos(id));
            % angular error
            BlockAvg.AngMean(i) = mean(Trial.AngError(id));
            BlockAvg.AngStd(i) = std(Trial.AngError(id));
            % normalized length
            BlockAvg.NormLMean(i) = mean(Trial.NormLength(id));
            BlockAvg.NormLStd(i) = std(Trial.NormLength(id));
            %Enclosed Area
            BlockAvg.EAreaMean(i) = mean(Trial.EArea(id));
            BlockAvg.EAreaStd(i) = std(Trial.EArea(id));
            % peek speed
            BlockAvg.PeakSpeedMean(i) = mean(Trial.PeakSpeed(id));
            BlockAvg.PeakSpeedStd(i) = std(Trial.PeakSpeed(id));
            % peek force
            BlockAvg.PeakForceMean(i) = mean(Trial.PeakForce(id));
            BlockAvg.PeakForceStd(i) = std(Trial.PeakForce(id));
            
            
            BlockAvg.Cnt(i) = length(id);
            
        end
    end
    
end


%%
savefile = [subID,'_BlockData_scc'];    %saving data
save(savefile,'BlockAvg')                %saving data

%% make some plots
close all

fig_trials = figure;
figure(fig_trials);


for blk = [1,2,3,4]
    
    trials = Blocks{blk};
    success = [0,diff(Trial.score)]; %only successfull trials
    for i = trials
        if success(i)
            id = Trial.Start_Out_Index(i):Trial.End_Out_Index(i);
            %
            %             subplot(3,2,blk)
            hold on
            plot(Trial.XPosition{i}(id), Trial.YPosition{i}(id))
            %
            %             subplot(3,2,blk+2)
            %             hold on
            %             plot(Trial.XVelocity{i}(id))
            %
            %             subplot(3,2,blk+4)
            %             hold on
            %             plot(Trial.YVelocity{i}(id))
        end
    end
end

% subplot(3,2,1)
axis equal
% subplot(3,2,2)
% axis equal
% subplot(3,2,3)
% title('x vel')
% subplot(3,2,5)
% title('y vel')
%
fig_vel = figure;
figure(fig_vel);

for blk = [1,2,3,4]
    
    trials = Blocks{blk};
    success = [0,diff(Trial.score)]; %only successfull trials
    for i = trials
        if success(i)
            id = Trial.Start_Out_Index(i):Trial.End_Out_Index(i);
            %
            subplot(2,2,1)
            hold on
            plot(Trial.XVelocity{i}(id))
            %
            subplot(2,2,2)
            hold on
            plot(Trial.YVelocity{i}(id))
            %
        end
    end
end
% axis equal
subplot(2,2,1)
title('x vel')
subplot(2,2,2)
title('y vel')

%% save fig
saveas(fig_trials, fullfile(datadir,['fig_trials_']),'fig')
saveas(fig_vel, fullfile(datadir,['fig_vel_']),'fig')

%% Add theta1 and theta2 to the processed data


l1 = ExpData(1).c3d(1).TP_TABLE.upper_arm(1)*0.01; %m
l2 = ExpData(1).c3d(1).TP_TABLE.lower_arm(1)*0.01; %m

shol_y = ExpData(1).c3d(1).TP_TABLE.sholder_distance(1)*0.01;

N = length(Trial.Good_Indices);

Trial.Theta1 = cell(1,N);
Trial.Theta2 = cell(1,N);

for i = 1:N
    
    id = Trial.Start_Out_Index(i):Trial.End_Out_Index(i);
    
    for j = 1:length(id)
        
        pos = [Trial.XPosition{i}(id(j));Trial.YPosition{i}(id(j)) + shol_y];
        THETA = inv_Position_2(pos,l1,l2);
        
        Trial.Theta1{i}(j) = THETA(1);
        Trial.Theta2{i}(j) = THETA(2);
    end
    
    
end

fig_theta = figure;
figure(fig_theta)
hold on

for i = 1:N
    
    plot(Trial.Theta1{i}, Trial.Theta2{i})
    
end

%%
saveas(fig_theta, fullfile(datadir,['fig_theta_']),'fig')

savefile = [subID,'_Processed'];    %saving data
save(savefile,'Trial')                %saving data

%% end code
