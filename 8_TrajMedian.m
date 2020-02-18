
%%
clear all
% close all
clc

% Define the Subjects

% Group = 'LimbFeedback'; 
% % datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets'; 
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets';
% Subjects = {'NF04', 'NF05', 'NF06', 'NF07', 'NF08', 'NF15', 'NF16', 'NF30', 'NF31', 'NF32'}; 

% Group = 'CursorFeedback';
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control';
% % datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets_control';
% Subjects = {'NF20', 'NF21', 'NF22', 'NF23', 'NF24', 'NF25', 'NF26', 'NF27', 'NF33', 'NF34'}; 

Group = 'MixFeedback';
datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Target_mix';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Target_mix';
Subjects = {'NF42', 'NF43', 'NF44', 'NF45', 'NF46', 'NF47', 'NF48', 'NF49', 'NF50', 'NF51'}; 


Subfile = Subjects;

num_subs = length(Subjects);


%% Load data into big matrix
%

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

for sub = 1:num_subs
    
%     M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData3']);
%     M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData2']);
    M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData']);
    S = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_AVG_Data_40']);
    
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

%% Compute median
% in this section I have a larger averaging time.
% I keep the count of the subjects in each data point

% I am syncing the data from time = 0 upto the smalest maximum time

% close all

fs = 1000;
fig_median = figure;

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
    tempX_SS = zeros(numpt,1);
    tempY_SS = zeros(numpt,1);
%     tempSp_SS = zeros(numpt,1);
    tempCnt = zeros(numpt,1);
    
    count = 0;
    for j = 1:num_subs
        
        count = count +1;
        
%         id_start = 1;
%         id_end = length(BlockAvgTime{block_cnt,j}); %this is specialized for each subject
        Start_idx = find(BlockAvgTime{block_cnt,j}>=min_time,1,'first');
        
        idx = Start_idx:numpt;
        
       IDX(count,:) = [idx(1), idx(end)];

        l = length(idx);
        
        tempX(1:l) = tempX(1:l) + BlockXrAvg{block_cnt,j}(idx);
        tempY(1:l) = tempY(1:l) + BlockYrAvg{block_cnt,j}(idx);
        tempCnt(1:l) = tempCnt(1:l) + 1;
    end
    
    X_avg{block_cnt} = tempX./tempCnt;
    Y_avg{block_cnt} = tempY./tempCnt;
    Time_output{block_cnt} = Avg_Time; %min_time:1/fs:max_time;
    X_Cnt{block_cnt} = tempCnt;
        
    % for when all the subjects have data, compute the median
    id1 = X_Cnt{block_cnt}>=9;
    for pnt = 1:sum(id1)
        temp1 = [];
        temp2 = [];
        for j = 1:num_subs
            
            id_start = IDX(j,1);
            id_end = IDX(j,2);
            
%             if length(BlockXrAvg{block_cnt,j})>=(pnt)
            if (id_end-id_start)>=(pnt)
                temp1 = [temp1, BlockXrAvg{block_cnt,j}(id_start + pnt-1)];
                temp2 = [temp2, BlockYrAvg{block_cnt,j}(id_start + pnt-1)];
            end
        end
        X_median{block_cnt}(pnt) = median(temp1);
        Y_median{block_cnt}(pnt) = median(temp2);
    end
    
    
    
    for j = 1:num_subs
        
        idx = IDX(j,1):IDX(j,2);
        l = length(idx);
        tempX_SS(1:l) = tempX_SS(1:l) + (BlockXrAvg{block_cnt,j}(idx) - X_avg{block_cnt}(1:l)).^2;
        tempY_SS(1:l) = tempY_SS(1:l) + (BlockYrAvg{block_cnt,j}(idx) - Y_avg{block_cnt}(1:l)).^2;

    end
    
    X_SEM{block_cnt} = sqrt(tempX_SS./X_Cnt{block_cnt});
    Y_SEM{block_cnt} = sqrt(tempY_SS./X_Cnt{block_cnt});
    
    figure(fig_median)
    hold on
    
    id = X_Cnt{block_cnt}>=1;
    
    for j = 1:10:length(X_avg{block_cnt})
        if id(j)
%             patch(X_median{block_cnt}(j) + X_SEM{block_cnt}(j)*cos(theta), Y_median{block_cnt}(j) + Y_SEM{block_cnt}(j)*sin(theta), [0.7,0.8,0.8],'EdgeColor','none')
%             patch(X_avg{block_cnt}(j) + X_SEM{block_cnt}(j)*cos(theta), Y_avg{block_cnt}(j) + Y_SEM{block_cnt}(j)*sin(theta), [0.7,0.8,0.8],'EdgeColor','none')
        end
    end
    plot(X_median{block_cnt},Y_median{block_cnt},'LineWidth',1.5,'Color',[0,0,0])
%     plot(X_avg{block_cnt}(id),Y_avg{block_cnt}(id),'LineWidth',1.5,'Color',[0,0,0])
    title('average trajectories')
    axis equal
    
end

BlockData.X_avg = X_avg;
BlockData.Y_avg = Y_avg;
BlockData.X_SEM = X_SEM;
BlockData.Y_SEM = Y_SEM;
BlockData.Time = Time_output;
BlockData.X_Cnt = X_Cnt;
BlockData.X_median = X_median;
BlockData.Y_median = Y_median;


%% saving

datadir_temp = [datadir,'\Averages'];
saveas(fig_median, fullfile(datadir_temp, 'fig_median_40'),'fig')
save(fullfile(datadir_temp,[Group, '_Median_40']), 'BlockData')

%% make a plot of all the averages

num_blocks = 4;

fig_reaches = figure;
hold on

for sub = 1:num_subs
    
    S = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_AVG_Data_40']);
    
    for i = 1:num_blocks
        plot(S.AVG_Data.X{i}, S.AVG_Data.Y{i})
    end
end

title(Group)
axis equal

%% saving
datadir_temp = [datadir,'\Averages'];
saveas(fig_reaches, fullfile(datadir_temp, 'fig_reaches'),'fig')

%% THE END

