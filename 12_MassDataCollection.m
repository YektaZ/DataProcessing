
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

%%

prot_block = 1;

fs = 1000; % sampling rate Hz


tot_num_trials = 0;
for r=1:prot_block
    tot_num_trials= tot_num_trials + length(ExpData(r).c3d);
end

% Creating RawData cell
RawData.MassX = cell(1,tot_num_trials);
RawData.MassY = cell(1,tot_num_trials);
RawData.MassVx = cell(1,tot_num_trials);
RawData.MassVy = cell(1,tot_num_trials);

fid = fopen('pat.dat','r');
Sub_Data = textscan(fid, '%s');
take_ID= strsplit(Sub_Data{1}{11},'\s*id=\s*','DelimiterType','RegularExpression');
subID=take_ID{end};


% Defining RawData values in the right trial order
for j=1:prot_block
    
    for counter = 1:length(ExpData(j).c3d)
        
        if j == 1
            trial = ExpData(1).c3d(counter).TRIAL.TRIAL_NUM;
        end
        RawData.MassX{trial} = ExpData(j).c3d(counter).MassX;
        RawData.MassY{trial} = ExpData(j).c3d(counter).MassY;
        RawData.MassVx{trial} = ExpData(j).c3d(counter).MassVX;
        RawData.MassVy{trial} = ExpData(j).c3d(counter).MassVY;
   
    end
    
end

%%%%%% load the indicies and get the processed data

Processed_Data = load([datadir,'\',subID,'_Processed']);

counter = 1;

for trial = Processed_Data.Trial.Good_Indices
    
    MassProcessed.XPosition{counter} = RawData.MassX{trial}; 
    MassProcessed.XVelocity{counter} = RawData.MassVx{trial}; 
    MassProcessed.YPosition{counter} = RawData.MassY{trial};
    MassProcessed.YVelocity{counter} = RawData.MassVy{trial};
    MassProcessed.Speed{counter} = sqrt(RawData.MassVx{trial}.^2 + RawData.MassVy{trial}.^2);
    
    
    counter = counter + 1;
    
end
MassProcessed.Time_Fixed = Processed_Data.Trial.Time_Fixed;

%%%%% average across the trials of each block

Block_Data = load([datadir,'\',subID,'_BlockData']);

Blocks = Block_Data.BlockAvg.Block_trial_num;

fs = 1/(MassProcessed.Time_Fixed{1}(2) - MassProcessed.Time_Fixed{1}(1));
theta = 0:pi/12:2*pi;
x = cos(theta);
y = sin(theta);

fig_speed = figure;
fig_avg = figure;


for i = 1:length(Blocks)

    [X_avg, Y_avg, Speed_avg, X_SEM, Y_SEM, Speed_SEM, Time] =...
        AVGES_3(Blocks{i},MassProcessed,fs);
 
    Mass_AVG_Data.X{i} = X_avg;
    Mass_AVG_Data.X_SEM{i} = X_SEM;
    Mass_AVG_Data.Y{i} = Y_avg;
    Mass_AVG_Data.Y_SEM{i} = Y_SEM;
    Mass_AVG_Data.Speed{i} = Speed_avg;
    Mass_AVG_Data.Speed_SEM{i} = Speed_SEM;
    Mass_AVG_Data.Time{i} = Time;
    
    figure(fig_speed)
    subplot(2,2,i)
    hold on
    plot(Time, Speed_avg + Speed_SEM,'--','Color','g')
    plot(Time, Speed_avg - Speed_SEM,'--','Color','g')
    plot(Time, Speed_avg, 'LineWidth',1.5,'Color','k')
    title('average speed')
    

    figure(fig_avg)
    hold on
    for j = 1:10:length(X_avg)
        patch(X_avg(j) + X_SEM(j)*cos(theta), Y_avg(j) + Y_SEM(j)*sin(theta), [0.7,0.8,0.8],'EdgeColor','none')
    end
    plot(X_avg,Y_avg,'LineWidth',1.5,'Color',[0,0,0])
    title('average trajectories')
    axis equal
end

%% saving

saveas(fig_avg, fullfile(datadir,'\',['fig_mass_avg']),'fig')
saveas(fig_speed, fullfile(datadir,'\', ['fig_mass_speed']),'fig')
save([subID, '_Mass_AVG_Data_'],'Mass_AVG_Data')

%% now average across subjects

clear all
close all
clc

% Define the Subjects

Group = 'LimbFeedback'; 
datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5'; 
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5';
Subjects = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014', 'SMD028'}; 

% Group = 'EnvFeedback';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2'; 
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% Subjects = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

Subfile = Subjects;

num_subs = length(Subjects);
num_blocks = 4;


for sub = 1:num_subs
    
%     M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData_scc']);
    S = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_Mass_AVG_Data_']);
    
for i = 1:num_blocks
        BlockXrAvg{i,sub} = S.Mass_AVG_Data.X{i};
        BlockYrAvg{i,sub} = S.Mass_AVG_Data.Y{i} - 0.15;
        BlockSpeedAvg{i,sub} = S.Mass_AVG_Data.Speed{i};
        BlockXrSEM{i,sub} = S.Mass_AVG_Data.X_SEM{i};
        BlockYrSEM{i,sub} = S.Mass_AVG_Data.Y_SEM{i};
        BlockSpeedSEM{i,sub} = S.Mass_AVG_Data.Speed_SEM{i};
        BlockAvgTime{i,sub} = S.Mass_AVG_Data.Time{i};
    end
end

fs = 1000;
fig_avg = figure;
fig_vel = figure;

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
        
%         idx = Start_idx:Start_idx + numpt-1;
        idx = Start_idx:numpt;
        
       IDX(count,:) = [idx(1), idx(end)];

        l = length(idx);
        
        tempX(1:l) = tempX(1:l) + BlockXrAvg{block_cnt,j}(idx);
        tempY(1:l) = tempY(1:l) + BlockYrAvg{block_cnt,j}(idx);
        tempSp(1:l) = tempSp(1:l) + BlockSpeedAvg{block_cnt,j}(idx);
        tempCnt(1:l) = tempCnt(1:l) + 1;
    end
    
    X_avg{block_cnt} = tempX./tempCnt;
    Y_avg{block_cnt} = tempY./tempCnt;
    Speed_avg{block_cnt} = tempSp./tempCnt;
    Time_output{block_cnt} = Avg_Time; %min_time:1/fs:max_time;
    X_Cnt{block_cnt} = tempCnt;
    
    for j = 1:num_subs
        
        idx = IDX(j,1):IDX(j,2);
        l = length(idx);
        tempX_SS(1:l) = tempX_SS(1:l) + (BlockXrAvg{block_cnt,j}(idx) - X_avg{block_cnt}(1:l)).^2;
        tempY_SS(1:l) = tempY_SS(1:l) + (BlockYrAvg{block_cnt,j}(idx) - Y_avg{block_cnt}(1:l)).^2;
        tempSp_SS(1:l) = tempSp_SS(1:l) + (BlockSpeedAvg{block_cnt,j}(idx) - Speed_avg{block_cnt}(1:l)).^2;

    end
    
    X_SEM{block_cnt} = sqrt(tempX_SS./X_Cnt{block_cnt});
    Y_SEM{block_cnt} = sqrt(tempY_SS./X_Cnt{block_cnt});
    Speed_SEM{block_cnt} = sqrt(tempSp_SS./X_Cnt{block_cnt});
    
    id = X_Cnt{block_cnt}>=5; %(num_subs)/2;
    
    figure(fig_avg)
    hold on
    for j = 1:10:length(X_avg{block_cnt}(id))
        if id(j)
            patch(X_avg{block_cnt}(j) + X_SEM{block_cnt}(j)*cos(theta), Y_avg{block_cnt}(j) + Y_SEM{block_cnt}(j)*sin(theta), [0.7,0.8,0.8],'EdgeColor','none')
        end
    end
    plot(X_avg{block_cnt}(id),Y_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    title('mass average trajectories')
    axis equal
    
    figure(fig_vel)
    subplot(2,2,block_cnt)
    hold on

    plot(Avg_Time(id), Speed_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    plot(Avg_Time(id),Speed_avg{block_cnt}(id)+Speed_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    plot(Avg_Time(id),Speed_avg{block_cnt}(id)-Speed_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    title('mass average speed')
    axis([-0.5,1.5,0,0.4])
    
    id_time(block_cnt) = find(Avg_Time>0,1);
end

figure(fig_avg)
axis([-0.15, 0.15, -0.15, 0.15])
plot([-0.1, 0.1, -0.1, 0.1],[-0.1, -0.1, 0.1, 0.1],'o','LineWidth',2)

MassBlockData.X_avg = X_avg;
MassBlockData.Y_avg = Y_avg;
MassBlockData.Speed_avg = Speed_avg;
MassBlockData.X_SEM = X_SEM;
MassBlockData.Y_SEM = Y_SEM;
MassBlockData.Speed_SEM = Speed_SEM;
MassBlockData.Time = Time_output;
MassBlockData.X_Cnt = X_Cnt;

%% Saving

datadir_temp = [datadir,'\Averages'];
saveas(fig_avg, fullfile(datadir_temp, 'fig_mass_avg'),'fig')
saveas(fig_vel, fullfile(datadir_temp, 'fig_mass_speed'),'fig')
save(fullfile(datadir_temp,[Group, '_Mass_Averages']), 'MassBlockData')


%% THE END

