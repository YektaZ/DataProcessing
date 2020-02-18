% load all individual subject data
% average across velocities

% average across subjects and get the group average.

%%


clear all
close all
clc

SubID = 'SMD020';

% path(path,'C:\GoogleDrive\Synced Folder\Matlab\Data_Processing\SMD_2Targ')
% path(path, 'C:\Users\Yeki\Google Drive\Synced Folder\Matlab\Data_Processing\SMD_2Targ');
path(path,'/Users/yektazahed/Google Drive/Synced Folder/Matlab/Data_Processing/SMD_2Targ')

% load(['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\',SubID,'_Processed.mat'])
% load(['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\',SubID,'_BlockData.mat'])

% load(['C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\',SubID,'_Processed.mat'])
% load(['C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\',SubID,'_BlockData.mat'])

load(['/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2/',SubID,'/',SubID,'_Processed.mat'])
load(['/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2/',SubID,'/',SubID,'_BlockData.mat'])

fs=1000; %Hz\
theta = 0:pi/12:2*pi;
x = cos(theta);
y = sin(theta);


%%%%% which trials to average over
block_trials = BlockAvg.Block_trial_num;

%% Compute Subject Averages For Blocks

% fig_vel = figure;
fig_speed = figure;

success = [0,diff(Trial.score)]; %only successfull trials
depth = zeros(1,4);

success(success == 2) = 1;

for i = 1:length(block_trials)
    
    trials = block_trials{i};
    trials_successful = nonzeros(trials.*success(trials));
    
    [Xvel_avg, Yvel_avg, Xvel_SEM, Yvel_SEM, Avg_Time] = AVGES_4(trials_successful',Trial,fs);
    
    Speed_avg = sqrt(Xvel_avg.^2 + Yvel_avg.^2);
    dt = Avg_Time(2)- Avg_Time(1);
    d_Speed = diff(Speed_avg)./dt;
    dd_Speed = diff(d_Speed)./dt;
   
    %%%%%%%%%%%%%% find the depth of the bump in speed profile
    
    % make sure the max and min point have zero slop and +/- acceleration
    
    numpt = length(Speed_avg);
    id_speed = floor(numpt*0.21):floor(numpt*0.65);
    
    thresh_Max = 0.001;
    thresh_Min = 100;
    
    for k = id_speed(3:end)
        if Speed_avg(k)>thresh_Max & d_Speed(k-1)<0.05 & dd_Speed(k-2)<0
            thresh_Max = Speed_avg(k);
            idx_max = k;
        end
        if Speed_avg(k)<thresh_Min & d_Speed(k-1)<0.05 & dd_Speed(k-2)>0
            thresh_Min = Speed_avg(k);
            id_min = k;
        end
        
    end
   
    
    depth(i) = abs(Speed_avg(idx_max) - Speed_avg(id_min));
    
    figure(fig_speed)
    subplot(2,2,i)
    hold on
    plot(Avg_Time, Speed_avg,'LineWidth',1.5,'Color',[0,0,0])
%     plot(Avg_Time(idx_max), Speed_avg(idx_max),'ro')
%     plot(Avg_Time(id_min),Speed_avg(id_min),'ko')
    title([SubID, ' Speed'])
    axis([-0.5,2,0,0.5])
    grid on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    AVG_Vel.X{i} = Xvel_avg;
    AVG_Vel.X_SEM{i} = Xvel_SEM;
    AVG_Vel.Y{i} = Yvel_avg;
    AVG_Vel.Y_SEM{i} = Yvel_SEM;
    AVG_Vel.Time{i} = Avg_Time;
    AVG_Vel.Speed{i} = Speed_avg;
    
%     figure(fig_vel)
%     subplot(2,2,i)
%     hold on
%     plot(Xvel_avg,'LineWidth',1.5,'Color','b')
%     plot(Yvel_avg,'LineWidth',1.5,'Color','k')
%     plot(Xvel_avg + Xvel_SEM,'LineWidth',1.5,'Color',[0.7,0.8,0.8])
%     plot(Xvel_avg - Xvel_SEM,'LineWidth',1.5,'Color',[0.7,0.8,0.8])
%     plot(Yvel_avg + Yvel_SEM,'LineWidth',1.5,'Color',[0.7,0.8,0.8])
%     plot(Yvel_avg - Yvel_SEM,'LineWidth',1.5,'Color',[0.7,0.8,0.8])
%     title('average velocities')
%     legend('Xvel','Yvel')
    
end

AVG_Vel.depth = depth;

%% saving for each subject

% datadir = ['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\',SubID,'\'];
% datadir = ['C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\'];
datadir = ['/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2/', SubID,'/'];
saveas(fig_speed, fullfile(datadir,['fig_speed']),'fig')
save(fullfile(datadir,[SubID, '_AVG_Vel']),'AVG_Vel')


%%
clear all
% close all
clc

% Define the Subjects

% Group = 'LimbFeedback';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5';
% Subjects = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014', 'SMD028'};

Group = 'EnvFeedback';
datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
Subjects = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

Subfile = Subjects;

num_subs = length(Subjects);


%% Load data into big matrix
% compute work

num_blocks = 4;

for sub = 1:num_subs
    
    S = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_AVG_Vel']);
    
    for i = 1:num_blocks
        BlockXvelAvg{i,sub} = S.AVG_Vel.X{i};
        BlockYvelAvg{i,sub} = S.AVG_Vel.Y{i};
        BlockXvelSEM{i,sub} = S.AVG_Vel.X_SEM{i};
        BlockYvelSEM{i,sub} = S.AVG_Vel.Y_SEM{i};
        BlockAvgTime{i,sub} = S.AVG_Vel.Time{i};
    end
end

%% Compute the averages (April 2019)
% in this section I have a larger averaging time.
% I keep the count of the subjects in each data point

% I am syncing the data from time = 0 upto the smalest maximum time

fs = 1000;

fig_vel = figure;
fig_speed = figure;

theta = 0:pi/12:2*pi;

for block_cnt = 1:num_blocks
    
    min_time = 0 -50/fs;;
    max_time = 100;
    
    
    for i = 1:num_subs
        temp = BlockAvgTime{block_cnt,i}(end);
        if temp <= max_time
            max_time = temp; %get the smalest max_time
        end
    end
    
    Avg_Time = min_time:1/1000:max_time;
    numpt = length(Avg_Time);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tempX = zeros(numpt,1);
    tempY = zeros(numpt,1);
    tempSp = zeros(numpt,1);
    tempX_SS = zeros(numpt,1);
    tempY_SS = zeros(numpt,1);
    tempSp_SS = zeros(numpt,1);
    tempCnt = zeros(numpt,1);
    
    count = 0;
    for j = 1:num_subs
        
        count = count +1;
        
        id_start = 1;
        id_end = length(BlockAvgTime{block_cnt,j}); %this is specialized for each subject
        Start_idx = find(BlockAvgTime{block_cnt,j}>=min_time,1,'first');
        
        idx = Start_idx:numpt;
        
        if length(idx)>1
            
            IDX(count,:) = [idx(1), idx(end)];
            
            l = length(idx);
            
            tempX(1:l) = tempX(1:l) + BlockXvelAvg{block_cnt,j}(idx);
            tempY(1:l) = tempY(1:l) + BlockYvelAvg{block_cnt,j}(idx);
            tempCnt(1:l) = tempCnt(1:l) + 1;
        end
    end
    
    Xvel_avg{block_cnt} = tempX./tempCnt;
    Yvel_avg{block_cnt} = tempY./tempCnt;
    Time_output{block_cnt} = Avg_Time;
    X_Cnt{block_cnt} = tempCnt;
    Speed_avg{block_cnt} = sqrt(Xvel_avg{block_cnt}.^2 + Yvel_avg{block_cnt}.^2);
    
    for j = 1:num_subs
        
        idx = IDX(j,1):IDX(j,2);
        l = length(idx);
        if l>1
            tempX_SS(1:l) = tempX_SS(1:l) + (BlockXvelAvg{block_cnt,j}(idx) - Xvel_avg{block_cnt}(1:l)).^2;
            tempY_SS(1:l) = tempY_SS(1:l) + (BlockYvelAvg{block_cnt,j}(idx) - Yvel_avg{block_cnt}(1:l)).^2;
        end
    end
    
    X_SEM{block_cnt} = sqrt(tempX_SS./X_Cnt{block_cnt});
    Y_SEM{block_cnt} = sqrt(tempY_SS./X_Cnt{block_cnt});
    
    
    id = X_Cnt{block_cnt}>=(num_subs)/2;
    
    figure(fig_vel)
    subplot(2,2,block_cnt)
    hold on
    
    plot(Avg_Time(id), Xvel_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    plot(Avg_Time(id),Xvel_avg{block_cnt}(id)+X_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    plot(Avg_Time(id),Xvel_avg{block_cnt}(id)-X_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    
    plot(Avg_Time(id), Yvel_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    plot(Avg_Time(id),Yvel_avg{block_cnt}(id)+Y_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    plot(Avg_Time(id),Yvel_avg{block_cnt}(id)-Y_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    title([Group, ' average velocities'])
    axis([-0.5,1.5,-0.3,0.3])
    grid on
    
    %     id_time(block_cnt) = find(Avg_Time>0,1);
    
    figure(fig_speed)
    subplot(2,2,block_cnt)
    hold on
    
    plot(Avg_Time(id), Speed_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    title([Group, ' average speed'])
    axis([-0.5,1.5,0,0.3])
    grid on
end

BlockData.X_avg = Xvel_avg;
BlockData.Y_avg = Yvel_avg;
BlockData.X_SEM = X_SEM;
BlockData.Y_SEM = Y_SEM;
BlockData.Time = Time_output;
BlockData.X_Cnt = X_Cnt;
BlockData.Speed = Speed_avg;

%% saving

datadir_temp = [datadir,'\Averages'];
saveas(fig_vel, fullfile(datadir_temp, 'fig_successfull_velocity'),'fig')
saveas(fig_speed, fullfile(datadir_temp, 'fig_successfull_speed'),'fig')
save(fullfile(datadir_temp,[Group, '_Vel_Successful_Averages']), 'BlockData')

%% make a plot of the trajectories on top of each other


clear all
% close all
clc

LimbFB_avg = load('C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\Averages\LimbFeedback_Averages_succ.mat');
EnvFB_avg = load('C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\Averages\EnvFeedback_Averages_succ.mat');

for block_cnt = 1:4
    
    figure
    hold on
    plot(LimbFB_avg.BlockData.Time{block_cnt}, LimbFB_avg.BlockData.X_avg{block_cnt},'k')
    plot(EnvFB_avg.BlockData.Time{block_cnt}, EnvFB_avg.BlockData.X_avg{block_cnt},'b')
    plot(LimbFB_avg.BlockData.Time{block_cnt}, LimbFB_avg.BlockData.X_avg{block_cnt} + LimbFB_avg.BlockData.X_SEM{block_cnt},'--','Color',[0.7,0.8,0.8])
    plot(LimbFB_avg.BlockData.Time{block_cnt}, LimbFB_avg.BlockData.X_avg{block_cnt} - LimbFB_avg.BlockData.X_SEM{block_cnt},'--','Color',[0.7,0.8,0.8])
    plot(EnvFB_avg.BlockData.Time{block_cnt}, EnvFB_avg.BlockData.X_avg{block_cnt} + EnvFB_avg.BlockData.X_SEM{block_cnt},'--','Color',[0.,0.8,0.8])
    plot(EnvFB_avg.BlockData.Time{block_cnt}, EnvFB_avg.BlockData.X_avg{block_cnt} - EnvFB_avg.BlockData.X_SEM{block_cnt},'--','Color',[0.,0.8,0.8])
    axis equal
    legend('limb fb')
    title('average trajectories')
    
    figure
    hold on
    plot(LimbFB_avg.BlockData.Time{block_cnt}, LimbFB_avg.BlockData.Y_avg{block_cnt},'k')
    plot(EnvFB_avg.BlockData.Time{block_cnt}, EnvFB_avg.BlockData.Y_avg{block_cnt},'b')
    plot(LimbFB_avg.BlockData.Time{block_cnt}, LimbFB_avg.BlockData.Y_avg{block_cnt} + LimbFB_avg.BlockData.Y_SEM{block_cnt},'--','Color',[0.7,0.8,0.8])
    plot(LimbFB_avg.BlockData.Time{block_cnt}, LimbFB_avg.BlockData.Y_avg{block_cnt} - LimbFB_avg.BlockData.Y_SEM{block_cnt},'--','Color',[0.7,0.8,0.8])
    plot(EnvFB_avg.BlockData.Time{block_cnt}, EnvFB_avg.BlockData.Y_avg{block_cnt} + EnvFB_avg.BlockData.Y_SEM{block_cnt},'--','Color',[0.,0.8,0.8])
    plot(EnvFB_avg.BlockData.Time{block_cnt}, EnvFB_avg.BlockData.Y_avg{block_cnt} - EnvFB_avg.BlockData.Y_SEM{block_cnt},'--','Color',[0.,0.8,0.8])
    axis equal
    legend('limb fb')
    title('average trajectories')
    
end


%%

close all
clear all
clc

% Define the Subjects

Group = 'LimbFeedback';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5';
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5';
datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_5';
Subjects = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014', 'SMD028'};

% Group = 'EnvFeedback';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2'
% Subjects = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

Subfile = Subjects;

num_subs = length(Subjects);

%%

fig_traj = figure;
figure(fig_traj)
hold on

num_blocks = 4;

for sub = 1:num_subs
    
    S = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_AVG_Data_']);
    
    for i = 1:num_blocks
        plot(S.AVG_Data.X{i},S.AVG_Data.Y{i})
    end
end

axis equal
title(Group)

%% saving

datadir_temp = [datadir,'\Averages'];
saveas(fig_traj, fullfile(datadir_temp, 'fig_sub_trajectories'),'fig')
% save(fullfile(datadir_temp,[Group, '_Vel_Successful_Averages']), 'BlockData')


%% load subject's avg_vel --> get depth

clear all
close all
clc

% Define the Subjects

Group_LFB = 'LimbFeedback';
% datadir_LFB = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5';
% datadir_LFB = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5';
datadir_LFB = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_5';
Subjects_LFB = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014', 'SMD028'};

Group_EFB = 'EnvFeedback';
% datadir_EFB = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% datadir_EFB = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
datadir_EFB = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2';
Subjects_EFB = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

num_subs_LFB = length(Subjects_LFB);
num_subs_EFB = length(Subjects_EFB);

%% 

Depths_LFB = [];

for sub = 1:num_subs_LFB
    
    V_LFB = load([datadir_LFB,'/',Subjects_LFB{sub},'/',Subjects_LFB{sub},'_AVG_Vel']);
    
    Depths_LFB = [Depths_LFB; V_LFB.AVG_Vel.depth];
end

Depths_EFB = [];

for sub = 1:num_subs_EFB
    
    V_EFB = load([datadir_EFB,'/',Subjects_EFB{sub},'/',Subjects_EFB{sub},'_AVG_Vel']);
    
    Depths_EFB = [Depths_EFB; V_EFB.AVG_Vel.depth];
end


[h, p] = ttest2(nonzeros(Depths_LFB(:)), nonzeros(Depths_EFB(:)))
[h, p] = ttest(nonzeros(Depths_LFB(:)), 0)
[h, p] = ttest(nonzeros(Depths_EFB(:)), 0)

%% load the data for each subject
% count the number of minimums for each subject (is it 1 or 2? )
% compare the two groups
% this code doesn't work, I ended up counting them with hand, much easier

%%%%%% notes: Sub = 5 LimbFB group has no sucsessful reach in block 1

clear all
close all
clc

% Define the Subjects

% Group = 'LimbFeedback';
% % datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5';
% datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_5';
% % /Users/yektazahed/Google Drive/Synced Folder/Matlab/Data_Processing/SMD_2Targ
% Subjects = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014', 'SMD028'};

Group = 'EnvFeedback';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2';
Subjects = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

Subfile = Subjects;

num_subs = length(Subjects);

min_count = zeros(num_subs,4);
max_count = zeros(num_subs,4);

thresh_Max0 = 0.001;
thresh_Min0 = 100;

for sub = 1%:num_subs
    
    clf; fig_speed = figure;
    
    Data = load([datadir,'/',Subjects{sub},'/',Subjects{sub},'_AVG_Vel.mat']);
    
    for i = 1:4
        
        if Data.AVG_Vel.X{i}
            
            id_max = [];
            id_min = [];
            
            Speed_avg = sqrt(Data.AVG_Vel.X{i}.^2 + Data.AVG_Vel.Y{i}.^2);
            dt = Data.AVG_Vel.Time{i}(2)- Data.AVG_Vel.Time{i}(1);
            d_Speed = diff(Speed_avg)./dt;
            dd_Speed = diff(d_Speed)./dt;
            
            %%%%%%%%%%%%%% find the depth of the bump in speed profile
            %%%%%%%%%%%%%% count the number of max and mins
            
            % make sure the max and min point have zero slop and +/- acceleration
            
            numpt = length(Speed_avg);
            id_speed = floor(numpt*0.2):floor(numpt*0.80);
            
            thresh_Max = thresh_Max0;
            thresh_Min = thresh_Min0;
            
            for k = id_speed(3:end)
                if Speed_avg(k)>thresh_Max & 0<d_Speed(k-1) & d_Speed(k-1)<0.01 & dd_Speed(k-2)<0
                    Max = Speed_avg(k);
                    thresh_Max = Max;
                    thresh_Min = thresh_Min0;
                    id_max = [id_max, k];
                    max_count(sub,i) = max_count(sub, i) + 1;
                end
                if Speed_avg(k)<thresh_Min & 0<d_Speed(k-1) & d_Speed(k-1)<0.01 & dd_Speed(k-2)>0
                    Min = Speed_avg(k);
                    thresh_Min = Min;
                    thresh_Max =thresh_Max0;
                    id_min = [id_min, k];
                    min_count(sub,i) = min_count(sub, i) + 1;
                end
                
            end
            
            
            figure(fig_speed)
            subplot(2,2,i)
            hold on
            plot(Data.AVG_Vel.Time{i}, Speed_avg,'LineWidth',1.5,'Color',[0,0,0])
            plot(Data.AVG_Vel.Time{i}(id_max), Speed_avg(id_max),'ro')
            plot(Data.AVG_Vel.Time{i}(id_min),Speed_avg(id_min),'ko')
            title([Subjects{sub}, ' Speed'])
            grid on
            pause(0.5)
            
        end
    end
end

display(min_count)
display(max_count)

%% THE END

