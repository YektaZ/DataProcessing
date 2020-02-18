

%% INITIALIZATION

clear all
close all
clc

SubID = 'SMD020';

% path(path,'C:\GoogleDrive\Synced Folder\Matlab\Data_Processing\SMD_2Targ')
path(path,'/Users/yektazahed/Google Drive/Synced Folder/Matlab/Data_Processing/SMD_2Targ')

% load(['C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Target_mix\',SubID,'\',SubID,'_Processed.mat']);
% load(['C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Target_mix\',SubID,'\',SubID,'_BlockData.mat']);
% load(['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\',SubID,'\Day1\',SubID,'_Processed.mat'])
% load(['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\',SubID,'\Day1\',SubID,'_BlockData.mat'])
% load(['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\',SubID,'_Processed.mat'])
% load(['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\',SubID,'_BlockData.mat'])
load(['/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2/',SubID,'/',SubID,'_Processed.mat'])
load(['/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2/',SubID,'/',SubID,'_BlockData.mat'])

% path(path,'C:\Users\Yeki\Google Drive\Synced Folder\Matlab\Data_Processing\NoForce4Targets')

fs=1000; %Hz\
theta = 0:pi/12:2*pi;
x = cos(theta);
y = sin(theta);


%%%%% which trials to average over
block_trials = BlockAvg.Block_trial_num;
% Blocks = {[1:2:100], [2:2:100]};
% block_trials = {[81:4:120], [81+1:4:120], [81+2:4:120], [81+3:4:120]};
% num_trials = '40';

%% Compute Subject Averages For Blocks
% close all

fig_speed = figure;
% fig_force = figure;
fig_avg = figure;

for i = 1:length(block_trials)
    
%     [X_avg, X_sem, Y_avg, Y_sem, Speed_avg, Speed_se, XForce_avg, XForce_se, YForce_avg, YForce_se, Time] =...
%         AVGES_2(block_trials{i},Trial,fs);
    [X_avg, Y_avg, Speed_avg, X_SEM, Y_SEM, Speed_SEM, Time] =...
        AVGES_3(block_trials{i},Trial,fs);
    
    AVG_Data.X{i} = X_avg;
    AVG_Data.X_SEM{i} = X_SEM;
    AVG_Data.Y{i} = Y_avg;
    AVG_Data.Y_SEM{i} = Y_SEM;
    AVG_Data.Speed{i} = Speed_avg;
    AVG_Data.Speed_SEM{i} = Speed_SEM;
    AVG_Data.Time{i} = Time;
%     AVG_Data.Cnt{i} = X_cnt;
    
    figure(fig_speed)
    subplot(2,2,i)
    hold on
    plot(Time, Speed_avg + Speed_SEM,'--','Color','g')
    plot(Time, Speed_avg - Speed_SEM,'--','Color','g')
    plot(Time, Speed_avg, 'LineWidth',1.5,'Color','k')
    title('average speed')
    
% id = X_cnt>8;
    figure(fig_avg)
%     subplot(2,2,i)
    hold on
    for j = 1:10:length(X_avg)
        patch(X_avg(j) + X_SEM(j)*cos(theta), Y_avg(j) + Y_SEM(j)*sin(theta), [0.7,0.8,0.8],'EdgeColor','none')
    end
    plot(X_avg,Y_avg,'LineWidth',1.5,'Color',[0,0,0])
    title('average trajectories')
    axis equal
    
end

% figure(fig_avg)

% subplot(2,2,1)
% axis([-0.15, 0.15, -0.1, 0])
% plot([-0.1, 0.1], [-0.05, -0.05],'o','LineWidth',2)
% 
% subplot(2,2,2)
% axis([-0.15, 0.15, -0.1, 0])
% plot([-0.1, 0.1], [-0.05, -0.05],'o','LineWidth',2)
% 
% subplot(2,2,3)
% axis([-0.15, 0.15, -0.1, 0])
% plot([-0.1, 0.1], [-0.05, -0.05],'o','LineWidth',2)
% 
% subplot(2,2,4)
% axis([-0.15, 0.15, -0.1, 0])
% plot([-0.1, 0.1], [-0.05, -0.05],'o','LineWidth',2)

%% saving

datadir = ['C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\',SubID,'\'];
saveas(fig_avg, fullfile(datadir,['fig_avg_']),'fig')
saveas(fig_speed, fullfile(datadir, ['fig_speed_']),'fig')
save([SubID, '_AVG_Data_'],'AVG_Data')


%% THE END





