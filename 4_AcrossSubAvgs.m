
%%
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


%% Load data into big matrix
% compute work 

num_blocks = 4;

BlockAbsPerp = zeros(num_blocks,num_subs);
BlockPerp = zeros(num_blocks,num_subs);
BlockPerpNeg = zeros(num_blocks,num_subs);
BlockPerpPos = zeros(num_blocks,num_subs);
BlockAng = zeros(num_blocks,num_subs);
BlockNorm = zeros(num_blocks,num_subs);
BlockEArea = zeros(num_blocks,num_subs);
BlockPeakSpeed = zeros(num_blocks,num_subs);
BlockPeakForce = zeros(num_blocks,num_subs);
BlockXrAvg = cell(num_blocks,num_subs);
BlockYrAvg = cell(num_blocks,num_subs);
BlockVelAvg = cell(num_blocks,num_subs);
BlockNetWork = cell(num_blocks,num_subs);

for sub = 1:num_subs
    
    M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData_scc']);
    S = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_AVG_Data_']);
    
    BlockAbsPerp(:,sub) = M.BlockAvg.AbsPerpMean;
    BlockPerp(:,sub) = M.BlockAvg.PerpMean;
    BlockPerpNeg(:,sub) = M.BlockAvg.PerpNegMean;
    BlockPerpPos(:,sub) = M.BlockAvg.PerpPosMean;
    BlockAng(:,sub) = M.BlockAvg.AngMean;
    BlockNorm(:,sub) = M.BlockAvg.NormLMean;
    BlockEArea(:,sub) = M.BlockAvg.EAreaMean;
    BlockPeakSpeed(:,sub) = M.BlockAvg.PeakSpeedMean;
    BlockPeakForce(:,sub) = M.BlockAvg.PeakForceMean;
    
for i = 1:num_blocks
        BlockXrAvg{i,sub} = S.AVG_Data.X{i};
        BlockYrAvg{i,sub} = S.AVG_Data.Y{i};
        BlockSpeedAvg{i,sub} = S.AVG_Data.Speed{i};
        BlockXrSEM{i,sub} = S.AVG_Data.X_SEM{i};
        BlockYrSEM{i,sub} = S.AVG_Data.Y_SEM{i};
        BlockSpeedSEM{i,sub} = S.AVG_Data.Speed_SEM{i};
        BlockAvgTime{i,sub} = S.AVG_Data.Time{i};
    end
end

%% Compute the averages new version (April 2019)
% in this section I have a larger averaging time.
% I keep the count of the subjects in each data point

% I am syncing the data from time = 0 upto the smalest maximum time

% close all

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
    
    id = X_Cnt{block_cnt}>=(num_subs)/2;
    
    figure(fig_avg)
    hold on
    for j = 1:10:length(X_avg{block_cnt}(id))
        if id(j)
            patch(X_avg{block_cnt}(j) + X_SEM{block_cnt}(j)*cos(theta), Y_avg{block_cnt}(j) + Y_SEM{block_cnt}(j)*sin(theta), [0.7,0.8,0.8],'EdgeColor','none')
        end
    end
    plot(X_avg{block_cnt}(id),Y_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    title('average trajectories')
    axis equal
    
    figure(fig_vel)
    subplot(2,2,block_cnt)
    hold on

    plot(Avg_Time(id), Speed_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    plot(Avg_Time(id),Speed_avg{block_cnt}(id)+Speed_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    plot(Avg_Time(id),Speed_avg{block_cnt}(id)-Speed_SEM{block_cnt}(id),'Color',[0.7,0.8,0.8])
    title('average speed')
    axis([-0.5,1.5,0,0.4])
    
    id_time(block_cnt) = find(Avg_Time>0,1);
end

figure(fig_avg)
axis([-0.15, 0.15, -0.15, 0.15])
plot([-0.1, 0.1, -0.1, 0.1],[-0.1, -0.1, 0.1, 0.1],'o','LineWidth',2)

BlockData.X_avg = X_avg;
BlockData.Y_avg = Y_avg;
BlockData.Speed_avg = Speed_avg;
BlockData.X_SEM = X_SEM;
BlockData.Y_SEM = Y_SEM;
BlockData.Speed_SEM = Speed_SEM;
BlockData.Time = Time_output;
BlockData.X_Cnt = X_Cnt;

%%%%%%%%%%%%%% average the metrices

BlockData.AbsPerpMean = mean(BlockAbsPerp,2);
BlockData.AbsPerpSEM = std(BlockAbsPerp,0,2)/sqrt(num_subs);
BlockData.BlockAbsPerp= BlockAbsPerp;

BlockData.PerpMean = mean(BlockPerp,2);
BlockData.PerpSEM = std(BlockPerp,0,2)/sqrt(num_subs);
BlockData.BlockPerp= BlockPerp;

BlockData.PerpNegMean = mean(BlockPerpNeg,2);
BlockData.PerpNegSEM = std(BlockPerpNeg,0,2)/sqrt(num_subs);
BlockData.BlockPerpNeg= BlockPerpNeg;

BlockData.PerpPosMean = mean(BlockPerpPos,2);
BlockData.PerpPosSEM = std(BlockPerpPos,0,2)/sqrt(num_subs);
BlockData.BlockPerpPos= BlockPerpPos;

BlockData.AngMean = mean(BlockAng,2);
BlockData.AngSEM = std(BlockAng,0,2)/sqrt(num_subs);
BlockData.BlockAng = BlockAng;

BlockData.NormMean = mean(BlockNorm,2);
BlockData.NormSEM = std(BlockNorm,0,2)/sqrt(num_subs);
BlockData.BlockNorm= BlockNorm;

BlockData.EAreaMean = mean(BlockEArea,2);
BlockData.EAreaSEM = std(BlockEArea,0,2)/sqrt(num_subs);
BlockData.BlockEArea= BlockEArea;

BlockData.PeakSpeedMean = mean(BlockPeakSpeed,2);
BlockData.PeakSpeedSEM = std(BlockPeakSpeed,0,2)/sqrt(num_subs);
BlockData.BlockPeakSpeed= BlockPeakSpeed;

BlockData.PeakForceMean = mean(BlockPeakForce,2);
BlockData.PeakForceSEM = std(BlockPeakForce,0,2)/sqrt(num_subs);
BlockData.BlockPeakForce= BlockPeakForce;

%% saving

datadir_temp = [datadir,'\Averages'];
saveas(fig_avg, fullfile(datadir_temp, 'fig_avg_40'),'fig')
saveas(fig_vel, fullfile(datadir_temp, 'fig_speed'),'fig')
save(fullfile(datadir_temp,[Group, '_Averages_succ']), 'BlockData')


%% THE END

