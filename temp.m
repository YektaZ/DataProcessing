clear
close all
clc

% remember to choose the right subject directory for datadir variable
% Do not change the path directory

datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets\NF04';
path(path,'C:\GoogleDrive\Synced Folder\Matlab\Data_Processing\NoForce4Targets')%Matlab function directory
ExpData = zip_load('dir',datadir);                            %creating the structure

%% save unzipped data and save them

fs = 1000; % sampling rate Hz
prot_block=2;   % how many protocol blocks         

tot_num_trials = 0;
for r=1:prot_block
    tot_num_trials= tot_num_trials + length(ExpData(r).c3d);
end

% Creating RawData cell
RawData.time_FS = cell(1,tot_num_trials);       %Time from the robot
RawData.time_trial = cell(1,tot_num_trials);    %Time elapsed on the trial
RawData.time_fixed=cell(1,tot_num_trials);      %Time synched to the Vel thresholds
RawData.time_global=cell(1,tot_num_trials);     %Total elapsed time of experiment
RawData.X = cell(1,tot_num_trials);
RawData.Y = cell(1,tot_num_trials);
RawData.VX = cell(1,tot_num_trials);
RawData.VY = cell(1,tot_num_trials);
RawData.Speed = cell(1,tot_num_trials);
RawData.Fx = cell(1,tot_num_trials);
RawData.Fy = cell(1,tot_num_trials);
RawData.Magn = cell(1,tot_num_trials);
RawData.reach_type = zeros(1,tot_num_trials);   %Load condition 1-6
RawData.target_num = zeros(1,tot_num_trials);
RawData.target_position = zeros(2,tot_num_trials); % 1=Xposition 2=Yposition 
RawData.center_position = zeros(2,tot_num_trials); % 1=Xposition 2=Yposition
RawData.cursor_feedback = cell(1,tot_num_trials);   %1: Visible, 0: Invisible
RawData.logical_radius = zeros(1,tot_num_trials);
RawData.score = zeros(1,tot_num_trials);            % Add score
RawData.event = zeros(1,tot_num_trials);

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
            trial = ExpData(1).c3d(counter).TRIAL.TRIAL_NUM;        %%%order?
        end
        if j==2
            trial = ExpData(2).c3d(counter).TRIAL.TRIAL_NUM + length(ExpData(1).c3d);
        end      
        
        RawData.X{trial} = ExpData(j).c3d(counter).Right_HandX;
        RawData.Y{trial} = ExpData(j).c3d(counter).Right_HandY - 0.15;
        RawData.Fx{trial} = ExpData(j).c3d(counter).Right_FS_ForceX;
        RawData.Fy{trial} = ExpData(j).c3d(counter).Right_FS_ForceY;
      
        % absolute time since experiment started
        time_FS = ExpData(j).c3d(counter).Right_FS_TimeStamp;   %Time from the Robot
        num_data_pts = length(RawData.X{trial});
        time = ((0:1:num_data_pts-1)/fs)';                    %Right time definition
        RawData.time_trial{trial} = time;
        RawData.time_FS{trial} = time_FS;
        
        % Target and reach information
        Type = ExpData(j).c3d(counter).TRIAL.TP-1;
        if Type <= 7 
            RawData.target_num(trial) = Type+1;
        else
            RawData.target_num(trial) = mod(Type,8)+1;
        end
        
        %Cursor feedback info
%         TP_Row=ExpData(j).c3d(counter).TRIAL.TP;
%         TP_Load = ExpData(j).c3d(counter).TP_TABLE.x3_Load(TP_Row);
%         RawData.cursor_feedback{trial} = ...
%             ExpData(j).c3d(counter).LOAD_TABLE.Cursor_Feedback(TP_Load);
        
        %target & center position info
        RawData.target_position(1,trial) =...
            ExpData(j).c3d(counter).TARGET_TABLE.X(RawData.target_num(trial)+1)/100;
        RawData.target_position(2,trial) =...
            ExpData(j).c3d(counter).TARGET_TABLE.Y(RawData.target_num(trial)+1)/100;
        RawData.logical_radius(1,trial) =...
            ExpData(j).c3d(counter).TARGET_TABLE.Logical_Radius(RawData.target_num(trial)+1)/100;
%         RawData.reach_type(trial) = ExpData(j).c3d(counter).TP_TABLE.x3_Load(TP_Row);
        RawData.center_position(1,trial) = ExpData(j).c3d(counter).TARGET_TABLE.X(1);
        RawData.center_position(2,trial) = ExpData(j).c3d(counter).TARGET_TABLE.Y(1);
        
        % Event "Stopped moving"
        EventLabel = strfind(ExpData(j).c3d(counter).EVENTS.LABELS,'Stoped Moving');
        idx = find(~cellfun(@isempty,EventLabel));
        if length(idx)>1
        idx = idx(length(idx)-1);
        end
        Event_time = ExpData(j).c3d(counter).EVENTS.TIMES(idx);
        RawData.event(1,trial) = Event_time;
        
        % Score info
        if isempty(find(strncmp(ExpData(j).c3d(counter).EVENTS.LABELS, 'Score event', 11),1))
            score_build(trial)=0;
        else
            score_build(trial)=1;
        end
        
    end
    
end



% define score
RawData.score(1) = score_build(1);
for trial=2:tot_num_trials
   RawData.score(trial) = RawData.score(trial-1) + score_build(trial);
end

%% Raw Velocity before Processing
% compute hand velocity & filter hand velocity

velid_stop= zeros(1,tot_num_trials);
bad_index = [];

trial = 0;

while trial < tot_num_trials
    
    ch = 100; % assume pressing d for all, if all trials are good
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
            bad_index = [bad_index trial] % append bad trial to bad_index
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
    
    vel_thresh_out=0.05;
    vel_thresh_back=0.01;
    
    % check for velocity threshold crossing
    velid_start = find(RawData.Speed{trial}>vel_thresh_out,1);
    
    if length(velid_start) == 0   %in case the trial never pass the vel thresh!
    bad_index = [bad_index trial]
    continue
    end
    
    tempid = find(RawData.Speed{trial}((velid_start+1):end)<vel_thresh_back,1,'first')+velid_start;
    velid_stop(trial)=tempid;
    temptime = RawData.time_trial{trial} - RawData.time_trial{trial}(velid_start);
    
    RawData.time_fixed{trial} = temptime;
    RawData.time_global{trial} = time_global;
  
  
end
RawData.velid_stop = velid_stop;

%% process data
  
% compute hand velocity & filter hand velocity
  
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
 
trial = 0;

while trial < tot_num_trials
        
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
            bad_index = [bad_index trial] % append bad trial to bad_index
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
    
    vel_thresh_out=0.05;
    vel_thresh_back=0.01;
    
    % check for velocity threshold crossing
    velid_start = find(RawData.Speed{trial}>vel_thresh_out,1);
    
    if length(velid_start) == 0   %in case the trial never pass the vel thresh!
    bad_index = [bad_index trial]
    continue
    end
    
    tempid = find(RawData.Speed{trial}((velid_start+1):end)<vel_thresh_back,1,'first')+velid_start;
    velid_stop(trial)=tempid;
    temptime = RawData.time_trial{trial} - RawData.time_trial{trial}(velid_start);
    
    RawData.time_fixed{trial} = temptime;
    RawData.time_global{trial} = time_global;
    
    % Plotting the velocity threshold to collect the bad trials
    
    figure(figSpeed)
    clf
    subplot(3,2,1)
    hold on;
    plot(temptime,RawData.Speed{trial},'b')
    plot(temptime(velid_start),RawData.Speed{trial}(velid_start),'rx',...
        temptime(tempid),RawData.Speed{trial}(tempid),'ko')
    
    line([temptime(1) temptime(end)],[vel_thresh_out vel_thresh_out],'Color','r')
    line([temptime(1) temptime(end)],[vel_thresh_back vel_thresh_back],'Color','k')
    
    title([' Trial : ',num2str(trial)]) 

    grid
    
    %
    subplot(3,2,3)
    hold on;
    plot(temptime,RawData.Y{trial},'b')
    plot(temptime(velid_start),RawData.Y{trial}(velid_start),'rx',...
        temptime(tempid),RawData.Y{trial}(tempid),'ko')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial) RawData.target_position(2,trial)],'Color','b')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial)-RawData.logical_radius(trial) RawData.target_position(2,trial)-RawData.logical_radius(trial)],'Color','b','LineStyle','--')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial)+RawData.logical_radius(trial) RawData.target_position(2,trial)+RawData.logical_radius(trial)],'Color','b','LineStyle','--')
    
    if ((RawData.Y{trial}(tempid))>=(RawData.target_position(2,trial)-RawData.logical_radius(trial))) && ((RawData.Y{trial}(tempid))<=(RawData.target_position(2,trial)+RawData.logical_radius(trial)))
    title([' Trial : ',num2str(trial)],'color','g')
    else
    title([' Trial : ',num2str(trial)],'color','r')
    end
    grid
    
    subplot(3,2,5)
    hold on
    plot(temptime,RawData.X{trial},'k')
    plot(temptime(velid_start),RawData.X{trial}(velid_start),'rx',...
        temptime(tempid),RawData.X{trial}(tempid),'ko')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial) RawData.target_position(1,trial)],'Color','k')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial)-RawData.logical_radius(trial) RawData.target_position(1,trial)-RawData.logical_radius(trial)],'Color','k','LineStyle','--')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial)+RawData.logical_radius(trial) RawData.target_position(1,trial)+RawData.logical_radius(trial)],'Color','k','LineStyle','--') 

    if ((RawData.X{trial}(tempid))>=(RawData.target_position(1,trial)-RawData.logical_radius(trial))) && ((RawData.X{trial}(tempid))<=(RawData.target_position(1,trial)+RawData.logical_radius(trial)))
    title([' Trial : ',num2str(trial)],'color','g')
    else
    title([' Trial : ',num2str(trial)],'color','r')
    end
    grid
    
    subplot(3,2,[2,4,6])
    hold on
    plot(RawData.X{trial},RawData.Y{trial})
    plot(RawData.X{trial}(velid_start),RawData.Y{trial}(velid_start),'rx')
    plot(RawData.X{trial}(tempid),RawData.Y{trial}(tempid),'ko')
    axis equal
    grid
    
    
    %disp(bad_index);
end
RawData.velid_stop = velid_stop;
%% NOW FIX BAD TRIALS ALLIGNMENT

vel_thresh_out = 0.05;
vel_stop_thresh = 0.02;
figSpeed = figure;
velid_stop = RawData.velid_stop;
remove_index=[];
remove_counter=1;

disp('Instructions:');
disp('- Press "s" to remove');
disp('- Press any key to proceed');
  temp_bad_index=bad_index(45:end);
for trial = temp_bad_index %bad_index%
%     if trial == 336
%         remove_index(remove_counter) = [trial]
%         remove_counter=remove_counter+1;
%     else
    % get this trial's data
   
    
    % check for velocity threshold crossing
    velid_start = find(RawData.Speed{trial}>vel_thresh_out,1);
    temptime = RawData.time_trial{trial} - RawData.time_trial{trial}(velid_start);
     
    
        
    % make sure everything is synced well
    clf
    figure(figSpeed)
    subplot(3,1,2)
    hold on;
    plot(temptime,RawData.Y{trial},'b')
    plot(temptime(velid_start),RawData.Y{trial}(velid_start),'rx')
    plot(temptime(velid_stop(trial)),RawData.Y{trial}(velid_stop(trial)),'rx')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial) RawData.target_position(2,trial)],'Color','b')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial)-RawData.logical_radius(trial) RawData.target_position(2,trial)-RawData.logical_radius(trial)],'Color','r','LineStyle','--')
    line([temptime(1) temptime(end)],[RawData.target_position(2,trial)+RawData.logical_radius(trial) RawData.target_position(2,trial)+RawData.logical_radius(trial)],'Color','r','LineStyle','--')

    subplot(3,1,3)
    hold on
    plot(temptime,RawData.X{trial},'k')
    plot(temptime(velid_start),RawData.X{trial}(velid_start),'rx')
    plot(temptime(velid_stop(trial)),RawData.X{trial}(velid_stop(trial)),'rx')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial) RawData.target_position(1,trial)],'Color','k')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial)-RawData.logical_radius(trial) RawData.target_position(1,trial)-RawData.logical_radius(trial)],'Color','r','LineStyle','--')
    line([temptime(1) temptime(end)],[RawData.target_position(1,trial)+RawData.logical_radius(trial) RawData.target_position(1,trial)+RawData.logical_radius(trial)],'Color','r','LineStyle','--')
     
    subplot(3,1,1)
    hold on;
    plot(temptime,RawData.Speed{trial},'b') 
    plot(temptime(velid_start),RawData.Speed{trial}(velid_start),'rx')
    line([temptime(1) temptime(end)],[vel_thresh_out vel_thresh_out])
    line([temptime(1) temptime(end)],[vel_stop_thresh vel_stop_thresh],'Color','r')
    title([' Trial : ',num2str(trial)])
    
        
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
    RawData.velid_stop(trial)= velid_stop_index;
    
    plot(temptime(RawData.velid_stop(trial)),RawData.Speed{trial}(RawData.velid_stop(trial)),'ro')
    subplot(3,1,2)
    plot(temptime(RawData.velid_stop(trial)),RawData.Y{trial}(RawData.velid_stop(trial)),'ro')
    subplot(3,1,3)
    plot(temptime(RawData.velid_stop(trial)),RawData.X{trial}(RawData.velid_stop(trial)),'ro')
    
    %define again temptime and save it as time_fixed
    temptime = RawData.time_trial{trial} - RawData.time_trial{trial}(velid_start);
    RawData.time_fixed{trial}= temptime;
    
    ch2 = getkey(); % collect the trials to remove after fixing
    if ch2 == 115
        remove_index(remove_counter) = [trial]
        remove_counter=remove_counter+1;

    end
    
    pause
    

%     end
end

%
%  SAVE DATA
%
savefile = [subID,'_Raw'];
save(fullfile(datadir,savefile),'RawData')        %saving data

%% COLLATE, STORE GOOD DATA

temp_all_index= [1:tot_num_trials];
good_index=temp_all_index;
for temp_remove=remove_index
    temp_removal=find(good_index==temp_remove);
    good_index(temp_removal)=[];
end

num_good_Trials = length(good_index);
 
Trial.Good_Indices = good_index;
Trial.Bad_Indices = bad_index;
Trial.trial_num = cell(1,num_good_Trials);
Trial.XPosition = cell(1,num_good_Trials);
Trial.YPosition = cell(1,num_good_Trials);
Trial.Speed = cell(1,num_good_Trials);
Trial.XVelocity = cell(1,num_good_Trials);
Trial.YVelocity = cell(1,num_good_Trials);
Trial.XForce = cell(1,num_good_Trials);
Trial.YForce = cell(1,num_good_Trials);
Trial.Magn = cell(1,num_good_Trials);
Trial.Time_Trial = cell(1,num_good_Trials);
Trial.Time_Fixed = cell(1,num_good_Trials);
Trial.Time_Global = cell(1,num_good_Trials);
Trial.Start_Out_Index = zeros(1,num_good_Trials);
Trial.End_Out_Index = zeros(1,num_good_Trials);
Trial.Subject_data = cell(1,5);

Trial.Reach_Type = cell(2,num_good_Trials);
Trial.Target_num = zeros(1,num_good_Trials);
Trial.Target_position = zeros(2,num_good_Trials);
Trial.Center_position = zeros(2,num_good_Trials);
Trial.Cursor_feedback = zeros(1,num_good_Trials);
Trial.score = zeros(1,num_good_Trials);
Trial.remove_index=remove_index;             %%%To save the removed trials
%Event "Stopped Moving"
Trial.Event_Idx = zeros(1,num_good_Trials);

%saving the subject's Height/Weight/sex/ID/Birthday
Trial.Subject_data{1}=Sub_Data{1}(3);
Trial.Subject_data{2}=Sub_Data{1}(5);
Trial.Subject_data{3}=Sub_Data{1}(10);
Trial.Subject_data{4}=subID;
Trial.Subject_data{5}=Sub_Data{1}(12);
   
count=1;
for i = good_index
    
    Trial.trial_num{count} = i;
    Trial.XPosition{count} = RawData.X{i};
    Trial.YPosition{count} = RawData.Y{i};
    Trial.Speed{count} = RawData.Speed{i};
    Trial.XVelocity{count} = RawData.VX{i};
    Trial.YVelocity{count} = RawData.VY{i};
    Trial.XForce{count} =RawData.Fx{i};
    Trial.YForce{count} =RawData.Fy{i};
    Trial.Magn{count} = RawData.Magn{i};
    Trial.Time_Trial{count} = RawData.time_trial{i};
    Trial.Time_Fixed{count} = RawData.time_fixed{i};
    Trial.Time_Global{count} = RawData.time_global{i};
%          Trial.Start_Out_Index(count) = find(Trial.Time_Fixed{count}>0,1,'first')-100; %%Reason?
%          Trial.End_Out_Index(count) = RawData.velid_stop(i)+ 100;
    
    Trial.Start_Out_Index(count) = find(Trial.Time_Fixed{count}>0,1,'first')-50; %%Reason?
    Trial.End_Out_Index(count) = RawData.velid_stop(i)+50;
    
    
    % Target and load information
    if RawData.reach_type(i) == 1 || RawData.reach_type(i) == 4
        Trial.Reach_Type{1,count} = RawData.reach_type(i);
        Trial.Reach_Type{2,count} = 'NULL FIELD';
    elseif RawData.reach_type(i) == 2 || RawData.reach_type(i) == 5
        Trial.Reach_Type{1,count} = RawData.reach_type(i);
        Trial.Reach_Type{2,count} = 'CURL FIELD';
    elseif RawData.reach_type(i) == 3 || RawData.reach_type(i) == 6
        Trial.Reach_Type{1,count} = RawData.reach_type(i);
        Trial.Reach_Type{2,count} = 'CHANNEL FIELD';
    end
    
    if  RawData.cursor_feedback{i} == 1
        Trial.Cursor_feedback(count) = false;   % That means VISIBLE
    elseif RawData.cursor_feedback{i} == 0
        Trial.Cursor_feedback(count) = true;    % That means INVISIBLE
    end
    
    Trial.Target_num(count) = RawData.target_num(i);
    Trial.Target_position(1,count) = RawData.target_position(1,i);
    Trial.Target_position(2,count) = RawData.target_position(2,i);
    Trial.Center_position(1,count) = RawData.center_position(1,i);
    Trial.Center_position(2,count) = RawData.center_position(2,i);
    Trial.score(count) = RawData.score(i);
    
    Trial.Event_Idx(count) = find(RawData.time_trial{i}>=RawData.event(i),1,'first');
    
    count=count+1;
    
end

%
% SAVE DATA
%

savefile = [Trial.Subject_data{4},'_Processed'];
save(savefile,'Trial')

%%  compute perpendicular error, angular error and normalized path length

tot_num_trials = length(Trial.Good_Indices);        %change the number

Reach_Type = nan(1,tot_num_trials);
PerpError_Pos = nan(1,tot_num_trials);
PerpError_Neg = nan(1,tot_num_trials);
PerpError_Abs = nan(1,tot_num_trials);
NormLength = nan(1,tot_num_trials);
AngError = nan(1,tot_num_trials);
peakF = nan(1,tot_num_trials);
Max_Speed = nan(1,tot_num_trials);
Trial.StopMove_Position = nan(2,tot_num_trials);

% for debugging
% errorfig = figure;
 
for trial = 1:1:tot_num_trials
    
    num_data_pts = length(Trial.XPosition{trial});
    CenterPos = [Trial.Center_position(1,trial) Trial.Center_position(2,trial)];% x-y
    
    % !!!!!
%     %Old Target Position
%     TargetPos = [Trial.XPosition{trial}(Trial.End_Out_Index(trial)), ...
%         Trial.YPosition{trial}(Trial.End_Out_Index(trial))];% x-y
%     
%     TargetLength = sqrt(sum((CenterPos - TargetPos).^2));
%     TargetDis = TargetPos-CenterPos;
%     
%     % rotate trajectory along axis of target direction
%     targ_Ang = atan2(TargetDis(2),TargetDis(1));
%     R = [cos(targ_Ang) -sin(targ_Ang); sin(targ_Ang) cos(targ_Ang)];
%     Xr = [Trial.XPosition{trial}-CenterPos(1)*ones(num_data_pts,1), ...
%         Trial.YPosition{trial}-CenterPos(2)*ones(num_data_pts,1)]*R;
%     
%     reach_indices = Trial.Start_Out_Index(trial):Trial.End_Out_Index(trial);
    
    
    
    % !!!!!
    % New Target Position
    TargetPos = [Trial.XPosition{trial}(Trial.Event_Idx(trial)), ...
        Trial.YPosition{trial}(Trial.Event_Idx(trial))];% x-y
    Trial.StopMove_Position(:,trial) = TargetPos;
    
    TargetLength = sqrt(sum((CenterPos - TargetPos).^2));
    TargetDis = TargetPos-CenterPos;
    
    % rotate trajectory along axis of target direction
    targ_Ang = atan2(TargetDis(2),TargetDis(1));
    R = [cos(targ_Ang) -sin(targ_Ang); sin(targ_Ang) cos(targ_Ang)];
    Xr = [Trial.XPosition{trial}-CenterPos(1)*ones(num_data_pts,1), ...
        Trial.YPosition{trial}-CenterPos(2)*ones(num_data_pts,1)]*R;
    
    reach_indices = Trial.Start_Out_Index(trial):Trial.End_Out_Index(trial);            %Refer to the reason of definition?
    
    % compute PERPENDICULAR ERROR
    [~, idx] = max(abs(Xr(reach_indices,2)));
    idx = idx + reach_indices(1)-1;
    PerpError(trial) = Xr(idx,2);

%     [~, idx_Abs] = max(abs(Xr(reach_indices,2)));
%     idx_Abs = idx_Abs + reach_indices(1)-1;
%     PerpError_Abs(trial) = Xr(idx_Abs,2);
        
    [~, idx_Pos] = max(Xr(reach_indices,2));
    idx_Pos = idx_Pos + reach_indices(1)-1;
    PerpError_Pos(trial) = Xr(idx_Pos,2);
    
    [~, idx_Neg] = min(Xr(reach_indices,2));
    idx_Neg = idx_Neg + reach_indices(1)-1;
    PerpError_Neg(trial) = Xr(idx_Neg,2);
    
    
    % compute NORMALIZED PATH LENGTH
    temp = 0;
    for i = reach_indices                   %%%no data points
        dx = Trial.XPosition{trial}(i+1) - Trial.XPosition{trial}(i);
        dy = Trial.YPosition{trial}(i+1) - Trial.YPosition{trial}(i);
        temp = temp + sqrt(dx^2 + dy^2);
    end
    NormLength(trial) = temp/TargetLength;
    
    % compute ANGULAR ERROR
    [temp, indx] = max(Trial.Speed{trial}(reach_indices));
    xm = Trial.XPosition{trial}(reach_indices(1)+indx-1);
    ym = Trial.YPosition{trial}(reach_indices(1)+indx-1);
    Max_Speed(trial) = temp;
    
    ang = atan2(ym-CenterPos(2), xm-CenterPos(1));
    temperr = ang-targ_Ang;
    % wrap ang errors
    if temperr > pi
        temperr = temperr - 2*pi;
    elseif temperr < -pi
        temperr = temperr + 2*pi;
    end
    AngError(trial) = temperr;
    
    % get peak force
    Fx = Trial.XForce{trial}(reach_indices);
    Fy = Trial.YForce{trial}(reach_indices);
    F = sqrt(Fx.^2 + Fy.^2);
    peakF(trial) = max(F);
    

% debugging


%         figure(errorfig)
%         clf
%         plot(Trial.XPosition{trial},Trial.YPosition{trial})
%         hold on
%         plot(Xr(:,1),Xr(:,2),'r',Xr(reach_indices,1),Xr(reach_indices,2),'k.')
%         plot(xm,ym,'ro')
%         plot(Xr(idx,1),Xr(idx,2),'md')
%         plot(Xr(idx_Pos,1),Xr(idx_Pos,2),'bd')
%         plot(Xr(idx_Neg,1),Xr(idx_Neg,2),'rd')
%         plot(Xr(idx_Abs,1),Xr(idx_Abs,2),'gd')
%         plot(TargetPos(1),TargetPos(2),'b*')
%         grid
%         axis equal
%         title([' perp : ',num2str(PerpError_Abs(trial)),', ang : ',num2str(temperr*180/pi),...
%             ', norm : ',num2str(NormLength(trial)), ' Trial : ' num2str(Trial.Good_Indices(trial))])
%         pause
 
end

Trial.PeakForce = peakF;
Trial.PerpError = PerpError;
Trial.PerpError_Abs = PerpError_Abs;
Trial.PerpError_Pos = PerpError_Pos;
Trial.PerpError_Neg = PerpError_Neg;
Trial.AngError = AngError;
Trial.NormLength = NormLength;
Trial.Max_Speed = Max_Speed;


%% SAVE DATA
%
savefile = [Trial.Subject_data{4},'_Processed'];    %saving data
save(savefile,'Trial')                %saving data
  

%% end code
