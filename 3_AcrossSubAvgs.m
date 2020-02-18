%%
% script to average numbers across subjects
%

%%
clear all
close all
clc

% Define the Subjects

Group = 'LimbFeedback'; 
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets'; 
datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets';
Subjects = {'NF04', 'NF05', 'NF06', 'NF07', 'NF08', 'NF15', 'NF16', 'NF30', 'NF31', 'NF32'}; 

Subfile = Subjects;

num_subs = length(Subjects);

%% Load data into big matrix
%

num_blocks = 4;

BlockPerp = zeros(num_blocks,num_subs);
BlockPerpNeg = zeros(num_blocks,num_subs);
BlockPerpPos = zeros(num_blocks,num_subs);
BlockAng = zeros(num_blocks,num_subs);
BlockNorm = zeros(num_blocks,num_subs);
BlockPeakSpeed = zeros(num_blocks,num_subs);
BlockPeakForce = zeros(num_blocks,num_subs);
BlockXrAvg = cell(num_blocks,num_subs);
BlockYrAvg = cell(num_blocks,num_subs);
BlockVelAvg = cell(num_blocks,num_subs);

for sub = 1:num_subs
    
    M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData']);
    S = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_AVG_Data']);
    

    BlockPerp(:,sub) = M.BlockAvg.PerpMean;
    BlockPerpNeg(:,sub) = M.BlockAvg.PerpNegMean;
    BlockPerpPos(:,sub) = M.BlockAvg.PerpPosMean;
    BlockAng(:,sub) = M.BlockAvg.AngMean;
    BlockNorm(:,sub) = M.BlockAvg.NormLMean;
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

%% Compute the averages new version (March 2019)

close all

fs = 1000;
fig_avg = figure;

theta = 0:pi/12:2*pi;

for block_cnt = 1:num_blocks
    
    min_time = -100;
    max_time = 100;
    
    
    for i = 1:num_subs
        Time = BlockAvgTime{block_cnt,i};
        temp = min(Time);
        if temp >= min_time
            min_time = temp; %get the largest min_time
        end
        
        temp = max(Time);
        if temp <= max_time
            max_time = temp; %get the smallest max_time
        end
    end
    
    % get rid of small roundoff errors?
    max_time = max_time - mod(max_time,1/fs);
    min_time = min_time - mod(min_time,1/fs);
    numpt = ceil((max_time-min_time)*fs)+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tempX = [];
    tempY = [];
    tempSp = [];
    
    for j = 1:num_subs
        Start_idx = find(BlockAvgTime{block_cnt,j}>=(min_time - 0.5/fs),1);
        Stop_idx  = Start_idx + numpt-5; %-1; %thorw some final points out
        tempX = [tempX; BlockXrAvg{block_cnt,j}(Start_idx : Stop_idx)'];
        tempY = [tempY; BlockYrAvg{block_cnt,j}(Start_idx : Stop_idx)'];
        tempSp = [tempSp; BlockSpeedAvg{block_cnt,j}(Start_idx : Stop_idx)'];
    end
    
    X_avg{block_cnt} = mean(tempX);
    X_sem{block_cnt} = std(tempX)/sqrt(num_subs);
    Y_avg{block_cnt} = mean(tempY);
    Y_sem{block_cnt} = std(tempY)/sqrt(num_subs);
    Speed_avg{block_cnt} = mean(tempSp);
    Speed_sem{block_cnt} = std(tempSp)/sqrt(num_subs);
    Time_output{block_cnt} = min_time:1/fs:max_time;
    Time_output{block_cnt} = Time_output{block_cnt}(1:length(tempX)); % for the last points
    
    
    figure(fig_avg)
    hold on
    for j = 1:10:length(X_avg{block_cnt})
        patch(X_avg{block_cnt}(j) + X_sem{block_cnt}(j)*cos(theta), Y_avg{block_cnt}(j) + Y_sem{block_cnt}(j)*sin(theta), [0.7,0.8,0.8],'EdgeColor','none')
    end
    plot(X_avg{block_cnt},Y_avg{block_cnt},'LineWidth',1.5,'Color',[0,0,0])
    title('average trajectories')
    axis equal
    
end

BlockData.X_avg = X_avg;
BlockData.Y_avg = Y_avg;
BlockData.Speed_avg = Speed_avg;
BlockData.X_SEM = X_sem;
BlockData.Y_SEM = Y_sem;
BlockData.Speed_SEM = Speed_sem;
BlockData.Time = Time_output;

%%%%%%%%%%%%%% average the metrices

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

BlockData.PeakSpeedMean = mean(BlockPeakSpeed,2);
BlockData.PeakSpeedSEM = std(BlockPeakSpeed,0,2)/sqrt(num_subs);
BlockData.BlockPeakSpeed= BlockPeakSpeed;

BlockData.PeakForceMean = mean(BlockPeakForce,2);
BlockData.PeakForceSEM = std(BlockPeakForce,0,2)/sqrt(num_subs);
BlockData.BlockPeakForce= BlockPeakForce;

%% saving

datadir_temp = [datadir,'\Averages'];
saveas(fig_avg, fullfile(datadir, 'fig_avg'),'fig')
save(fullfile(datadir_temp,[Group, '_Averages']), 'BlockData')

%% Compute averages and SEM's and plot the AVG Traj
close all

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

BlockData.PeakSpeedMean = mean(BlockPeakSpeed,2);
BlockData.PeakSpeedSEM = std(BlockPeakSpeed,0,2)/sqrt(num_subs);
BlockData.BlockPeakSpeed= BlockPeakSpeed;

BlockData.PeakForceMean = mean(BlockPeakForce,2);
BlockData.PeakForceSEM = std(BlockPeakForce,0,2)/sqrt(num_subs);
BlockData.BlockPeakForce= BlockPeakForce;

%

% BinData.PerpMean = mean(BinPerp,2);
% BinData.PerpSEM = std(BinPerp,0,2)/sqrt(num_subs);
% BinData.BinPerp= BinPerp;
% 
% BinData.PerpNegMean = mean(BinPerpNeg,2);
% BinData.PerpNegSEM = std(BinPerpNeg,0,2)/sqrt(num_subs);
% BinData.BinPerpNeg= BinPerpNeg;
% 
% BinData.PerpPosMean = mean(BinPerpPos,2);
% BinData.PerpPosSEM = std(BinPerpPos,0,2)/sqrt(num_subs);
% BinData.BinPerpPos= BinPerpPos;
% 
% BinData.AngMean = mean(BinAng,2);
% BinData.AngSEM = std(BinAng,0,2)/sqrt(num_subs);
% BinData.BinAng = BinAng;
% 
% BinData.NormMean = mean(BinNorm,2);
% BinData.NormSEM = std(BinNorm,0,2)/sqrt(num_subs);
% BinData.BinNorm= BinNorm;
% 
% BinData.PeakSpeedMean = mean(BinPeakSpeed,2);
% BinData.PeakSpeedSEM = std(BinPeakSpeed,0,2)/sqrt(num_subs);
% BinData.BinPeakSpeed= BinPeakSpeed;
% 
% BinData.PeakForceMean = mean(BinPeakForce,2);
% BinData.PeakForceSEM = std(BinPeakForce,0,2)/sqrt(num_subs);
% BinData.BinPeakForce= BinPeakForce;

% Trajectory Average

ang = 0:.1:2*pi;
circle_x = cos(ang);
circle_y = sin(ang);

traj_avg = figure;
% vel_avg = figure;

for block_cnt = 1:num_blocks
    
    min_t = 10^4;
    max_t = -10^4;
    
    for i = 1:num_subs
        temp1 = BlockAvgTime{block_cnt,i}(1);
        if temp1 <= min_t
            min_t = temp1;
        end
        temp2 = BlockAvgTime{block_cnt,i}(end);
        if temp2 >= max_t
            max_t = temp2;
        end
    end
    
    Avg_Time = min_t:1/1000:max_t;
    N = length(Avg_Time);
    
    
    Xr_Sum = zeros(N,1);
    Yr_Sum = zeros(N,1);
    Vel_Sum = zeros(N,1);
    Xr_Cnt = zeros(N,1);
    Xr_SS = zeros(N,1);
    Yr_SS = zeros(N,1);
    Vel_SS = zeros(N,1);
    
    count = 0;
    for j = 1:num_subs
        
        count = count +1;
        
        id_start = 1;
        id_end = length(BlockAvgTime{block_cnt,j});
        
        t0 = BlockAvgTime{block_cnt,j}(1);
        t0 = str2num(sprintf('%.2f',t0)); %decimal precision error
        id1 = find(Avg_Time>=t0,1,'first'); %sync the data
        idx = id1:(id1+id_end-id_start);
%         lgth = length (idx);
        
        IDX(count,:) = [idx(1),idx(end)]; %save the idx for the next loop
        
        Xr_Sum(idx) = Xr_Sum(idx) + BlockXrAvg{block_cnt,j}(idx);
        Yr_Sum(idx) = Yr_Sum(idx) + BlockYrAvg{block_cnt,j}(idx);
%         Vel_Sum(idx) = Vel_Sum(idx) + BlockVelAvg{block_cnt,j}(id_start:id_end);
        Xr_Cnt(idx) = Xr_Cnt(idx) + 1;
        
    end
    
    BlockData.XrAvg{block_cnt} = Xr_Sum./Xr_Cnt;
    BlockData.YrAvg{block_cnt} = Yr_Sum./Xr_Cnt;
%     BlockData.VelAvg{block_cnt} = Vel_Sum./Xr_Cnt;
    BlockData.Time{block_cnt} = Avg_Time;
    
    count = 0;
    for j=1:num_subs
        
        count = count +1;
        
        id_start= 1;
        id_end = length(BlockAvgTime{block_cnt,j});
        
        idx = IDX(count,1):IDX(count,2);
        
        Xr_SS(idx) = Xr_SS(idx) + ( BlockXrAvg{block_cnt,j}(idx)-BlockData.XrAvg{block_cnt}(idx) ).^2;
        Yr_SS(idx) = Yr_SS(idx) + ( BlockYrAvg{block_cnt,j}(idx)-BlockData.YrAvg{block_cnt}(idx) ).^2;
%         Vel_SS(idx) = Vel_SS(idx) + ( BlockVelAvg{block_cnt,j}(id_start:id_end)-BlockData.VelAvg{block_cnt}(idx) ).^2;
        
    end
    
   XrStd{block_cnt} = sqrt( Xr_SS./Xr_Cnt );
   YrStd{block_cnt} = sqrt( Yr_SS./Xr_Cnt );
%    VelStd = sqrt(Vel_SS./Xr_Cnt);

BlockData.XrCnt{block_cnt} = Xr_Cnt';    
    
    %ploting the rotated trajectories
    figure(traj_avg)
%     subplot(3,4,block_cnt)
subplot(2,3,block_cnt)
    plot(0, 0.1,'o')
    hold on
    
    id = Xr_Cnt>=(num_subs);
    
    for dd = 1:5:length(BlockData.XrAvg{block_cnt})
        if id(dd)
        patch(BlockData.XrAvg{block_cnt}(dd)+circle_x*XrStd{block_cnt}(dd), ...
            BlockData.YrAvg{block_cnt}(dd)+circle_y*YrStd{block_cnt}(dd),[.8 .8 .8],'EdgeColor','none')
        end
    end
    
    if block_cnt==11 || block_cnt==12
    plot(BlockData.XrAvg{block_cnt}(id),BlockData.YrAvg{block_cnt}(id),'r')
    else
        plot(BlockData.XrAvg{block_cnt}(id),BlockData.YrAvg{block_cnt}(id),'b')
    end
    
%     figure(vel_avg)
%     subplot(3,4,block_cnt)
%     plot(BlockData.Time{block_cnt}(id),BlockData.VelAvg{block_cnt}(id))
%     hold on
%     plot(BlockData.Time{block_cnt}(id),BlockData.VelAvg{block_cnt}(id)+VelStd(id),'color',[0.5 0.5 0.5])
%     plot(BlockData.Time{block_cnt}(id),BlockData.VelAvg{block_cnt}(id)-VelStd(id),'color',[0.5 0.5 0.5])
%     
end

BlockData.XrStd = XrStd;
BlockData.YrStd = YrStd;

figure(traj_avg)
% subplot(3,4,1)
% subplot(2,3,1)
% axis([-0.04 0.04 0 .12])
title([Group, ' Traj Avgs'])
grid
% % subplot(3,4,2)
% subplot(2,3,2)
% axis([-0.04 0.04 0 .12])
% title('early adp D1')
% grid
% % subplot(3,4,3)
% subplot(2,3,3)
% axis([-0.04 0.04 0 .12])
% title('late adp D1')
% grid
% % subplot(3,4,4)
% subplot(2,3,4)
% axis([-0.04 0.04 0 .12])
% title('early wash D1')
% grid
% % subplot(3,4,5)
% subplot(2,3,5)
% axis([-0.04 0.04 0 .12])
% title('late wash D1')
% grid
% % subplot(3,4,6)
% subplot(2,3,6)
% axis([-0.04 0.04 0 .12])
% title('baseline D2')
% grid

%%
directory = [datadir, '/Averages'];
filename = [Group,'_traj_avg'];
% filename2 = [Group,'_vel_avg'];
saveas(traj_avg,fullfile(directory,filename),'fig')
saveas(traj_avg,fullfile(directory,filename),'eps')
% saveas(vel_avg,filename2,'fig')

file = [Group,'_Block&Bin_AVGS'];
save(fullfile(directory, file),'BinData', 'BlockData')

%% Plots
close all

error_perp_blocks = figure ;
subplot(3,1,1)
bar(1:num_blocks,BlockData.PerpMean)
hold on
errorbar(1:num_blocks,BlockData.PerpMean,BlockData.PerpSEM,'.')
title('PerpError Avg on Blocks')
grid
subplot(3,1,2)
bar(1:num_blocks,BlockData.PerpNegMean)
hold on
errorbar(1:num_blocks,BlockData.PerpNegMean,BlockData.PerpNegSEM,'.')
title('PerpError Negetive Avg on Blocks')
grid
subplot(3,1,3)
bar(1:num_blocks,BlockData.PerpPosMean)
hold on
errorbar(1:num_blocks,BlockData.PerpPosMean,BlockData.PerpPosSEM,'.')
title('PerpError Positive Avg on Blocks')
grid


error_blocks = figure;
subplot(4,1,1)
bar(BlockData.AngMean)
hold on
errorbar(1:num_blocks,BlockData.AngMean,BlockData.AngSEM,'.')
title('AngError Avg on Blocks')
grid
subplot(4,1,2)
bar(BlockData.NormMean)
hold on
errorbar(1:num_blocks,BlockData.NormMean,BlockData.NormSEM,'.')
title('NormLength Avg on Blocks')
grid
subplot(4,1,3)
bar(BlockData.PeakSpeedMean)
hold on
errorbar(1:num_blocks,BlockData.PeakSpeedMean,BlockData.PeakSpeedSEM,'.')
title('Peak Speed Avg on Blocks')
grid
subplot(4,1,4)
bar(BlockData.PeakForceMean)
hold on
errorbar(1:num_blocks,BlockData.PeakForceMean,BlockData.PeakForceSEM,'.')
title('Peak Force Avg on Blocks')
grid

%
error_perp_bins = figure;
subplot(3,1,1)
errorbar(BinData.PerpMean, BinData.PerpSEM,'.')
title('PerpError Avg on Bins')
% axis([0 160 -0.04 0.04])
grid
subplot(3,1,2)
errorbar(BinData.PerpNegMean, BinData.PerpNegSEM,'.')
title('PerpError Negetive Avg on Bins')
% axis([0 160 -0.04 0])
grid
subplot(3,1,3)
errorbar(BinData.PerpPosMean, BinData.PerpPosSEM,'.')
title('PerpError Positive Avg on Bins')
% axis([0 160 0 0.04])
grid
%
error_bins = figure;
subplot(4,1,1)
errorbar(BinData.AngMean, BinData.AngSEM,'.')
title('AngError Avg on Bins')
grid
subplot(4,1,2)
errorbar(BinData.NormMean, BinData.NormSEM,'.')
title('NormLength Avg on Bins')
grid
subplot(4,1,3)
errorbar(BinData.PeakSpeedMean, BinData.PeakSpeedSEM,'.')
title('Peak Speed Avg on Bins')
grid
subplot(4,1,4)
errorbar(BinData.PeakForceMean, BinData.PeakForceSEM,'.')
title('Peak Force Avg on Bins')
grid

% %
% figure
% title(Group)
% subplot(2,3,1)
% plot(BlockData.XrAvg{1},BlockData.YrAvg{1})
% grid
% axis equal
% subplot(2,3,2)
% plot(BlockData.XrAvg{2},BlockData.YrAvg{2})
% grid
% axis equal
% subplot(2,3,3)
% plot(BlockData.XrAvg{3},BlockData.YrAvg{3})
% grid
% axis equal
% subplot(2,3,4)
% plot(BlockData.XrAvg{4},BlockData.YrAvg{4})
% grid
% axis equal
% subplot(2,3,5)
% plot(BlockData.XrAvg{5},BlockData.YrAvg{5})
% grid
% axis equal


%% save across-subject averages

% path(path,...
%     '/Users/Yekta/Documents/PhD/Lab/KinArm Data/Bar Feedback/four target/Agv across subjects');

directory = [datadir, '/Averages'];
file = [Group,'_Block&Bin_AVGS'];
save(fullfile(directory, file),'BinData', 'BlockData')

filename = [Group,'_blocks_1'];
saveas(error_perp_blocks,fullfile(directory,filename),'fig')

filename = [Group,'_blocks_2'];
saveas(error_blocks,fullfile(directory,filename),'fig')

filename = [Group,'_bins_1'];
saveas(error_perp_bins,fullfile(directory,filename),'fig')

filename = [Group,'_bins_2'];
saveas(error_bins,fullfile(directory,filename),'fig')

%% error figures for the paper
%
close all
%
error_perp_bins = figure;
subplot(2,1,2)
hold on
plot(BinData.AngMean(1:16),'b')
plot(17:60, BinData.AngMean(17:60),'r')
plot(61:75, BinData.AngMean(61:75),'b')
plot(76:80, BinData.AngMean(76:80),'r')
plot(81:124, BinData.AngMean(81:124),'b')
plot(125:154, BinData.AngMean(125:154),'r')
%
plot(BinData.AngMean(1:16)+BinData.AngSEM(1:16),'c')
plot(BinData.AngMean(1:16)-BinData.AngSEM(1:16),'c')
plot(17:60,BinData.AngMean(17:60)+BinData.AngSEM(17:60),'c')
plot(17:60,BinData.AngMean(17:60)-BinData.AngSEM(17:60),'c')
plot(61:75,BinData.AngMean(61:75)+BinData.AngSEM(61:75),'c')
plot(61:75,BinData.AngMean(61:75)-BinData.AngSEM(61:75),'c')
plot(76:80,BinData.AngMean(76:80)+BinData.AngSEM(76:80),'c')
plot(76:80,BinData.AngMean(76:80)-BinData.AngSEM(76:80),'c')
plot(81:124,BinData.AngMean(81:124)+BinData.AngSEM(81:124),'c')
plot(81:124,BinData.AngMean(81:124)-BinData.AngSEM(81:124),'c')
plot(125:154,BinData.AngMean(125:154)+BinData.AngSEM(125:154),'c')
plot(125:154,BinData.AngMean(125:154)-BinData.AngSEM(125:154),'c')
title('AngError Avg on Bins')
axis([0 155 -30 30])
line([0, 155], [0, 0])


subplot(2,1,1)
hold on
plot(BinData.PerpNegMean(1:16),'b')
plot(17:60, BinData.PerpNegMean(17:60),'r')
plot(61:75, BinData.PerpNegMean(61:75),'b')
plot(76:80, BinData.PerpNegMean(76:80),'r')
plot(81:124, BinData.PerpNegMean(81:124),'b')
plot(125:154, BinData.PerpNegMean(125:154),'r')
%
plot(BinData.PerpNegMean(1:16)+BinData.PerpNegSEM(1:16),'c')
plot(BinData.PerpNegMean(1:16)-BinData.PerpNegSEM(1:16),'c')
plot(17:60,BinData.PerpNegMean(17:60)+BinData.PerpNegSEM(17:60),'c')
plot(17:60,BinData.PerpNegMean(17:60)-BinData.PerpNegSEM(17:60),'c')
plot(61:75,BinData.PerpNegMean(61:75)+BinData.PerpNegSEM(61:75),'c')
plot(61:75,BinData.PerpNegMean(61:75)-BinData.PerpNegSEM(61:75),'c')
plot(76:80,BinData.PerpNegMean(76:80)+BinData.PerpNegSEM(76:80),'c')
plot(76:80,BinData.PerpNegMean(76:80)-BinData.PerpNegSEM(76:80),'c')
plot(81:124,BinData.PerpNegMean(81:124)+BinData.PerpNegSEM(81:124),'c')
plot(81:124,BinData.PerpNegMean(81:124)-BinData.PerpNegSEM(81:124),'c')
plot(125:154,BinData.PerpNegMean(125:154)+BinData.PerpNegSEM(125:154),'c')
plot(125:154,BinData.PerpNegMean(125:154)-BinData.PerpNegSEM(125:154),'c')
% title('PerpError Negetive Avg on Bins')
% axis([0 155 -0.05 0])
% subplot(3,1,1)
% hold on
line([0, 155], [0, 0])

plot(BinData.PerpPosMean(1:16),'b')
plot(17:60, BinData.PerpPosMean(17:60),'r')
plot(61:75, BinData.PerpPosMean(61:75),'b')
plot(76:80, BinData.PerpPosMean(76:80),'r')
plot(81:124, BinData.PerpPosMean(81:124),'b')
plot(125:154, BinData.PerpPosMean(125:154),'r')
%
plot(BinData.PerpPosMean(1:16)+BinData.PerpPosSEM(1:16),'c')
plot(BinData.PerpPosMean(1:16)-BinData.PerpPosSEM(1:16),'c')
plot(17:60,BinData.PerpPosMean(17:60)+BinData.PerpPosSEM(17:60),'c')
plot(17:60,BinData.PerpPosMean(17:60)-BinData.PerpPosSEM(17:60),'c')
plot(61:75,BinData.PerpPosMean(61:75)+BinData.PerpPosSEM(61:75),'c')
plot(61:75,BinData.PerpPosMean(61:75)-BinData.PerpPosSEM(61:75),'c')
plot(76:80,BinData.PerpPosMean(76:80)+BinData.PerpPosSEM(76:80),'c')
plot(76:80,BinData.PerpPosMean(76:80)-BinData.PerpPosSEM(76:80),'c')
plot(81:124,BinData.PerpPosMean(81:124)+BinData.PerpPosSEM(81:124),'c')
plot(81:124,BinData.PerpPosMean(81:124)-BinData.PerpPosSEM(81:124),'c')
plot(125:154,BinData.PerpPosMean(125:154)+BinData.PerpPosSEM(125:154),'c')
plot(125:154,BinData.PerpPosMean(125:154)-BinData.PerpPosSEM(125:154),'c')
title('PerpError Positive and Negetive Avg on Bins')
% axis([0 155 0 0.05])
axis([0 155 -0.05 0.05])
% %
error_perp_blocks = figure;
subplot(2,1,2)
hold on
% bar([16,17,60,61,75,80, 81, 124,125,154],[20,20,20,20,20,20,20,20,20,20],'w')
% bar([16,17,60,61,75,80, 81, 124,125,154],[-20,-20,-20,-20,-20,-20,-20,-20,-20,-20],'w')
errorbar([16,17,60,61,75,80, 81, 124,125,154],BlockData.AngMean, BlockData.AngSEM,'.')
title('AngError Avg on Blocks')
axis([0 155 -30 30])


subplot(2,1,1)
hold on
% bar([16,17,60,61,75,80, 81, 124,125,154],[0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04],'w')
% bar([16,17,60,61,75,80, 81, 124,125,154],[-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04],'w')
errorbar([16,17,60,61,75,80, 81, 124,125,154],BlockData.PerpNegMean, BlockData.PerpNegSEM,'.')
% title('PerpError Negetive Avg on Blocks')
% axis([0 155 -0.05 0.05])
% subplot(3,1,1)
% hold on
line([0, 155], [0, 0])
% bar([16,17,60,61,75,80, 81, 124,125,154],[0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04],'w')
% bar([16,17,60,61,75,80, 81, 124,125,154],[-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04],'w')
errorbar([16,17,60,61,75,80, 81, 124,125,154],BlockData.PerpPosMean, BlockData.PerpPosSEM,'.')
title('PerpError Positive and Negetive Avg on Blocks')
axis([0 155 -0.05 0.05])
%
%% Saving for the paper

directory = [datadir, '/Averages'];
file = [Group,'_errors'];
figure(error_perp_bins);
saveas(gca, fullfile(directory, file),'fig')
saveas(gca, fullfile(directory, file),'eps')
%
file = [Group,'_block_errors'];
figure(error_perp_blocks);
saveas(gca, fullfile(directory, file),'fig')
saveas(gca, fullfile(directory, file),'eps')

%% END CODE
