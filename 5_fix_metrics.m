% fix the target_pos/ enclosed area / Abs perp error

% load raw data
% get the target positions/ center position
% load processed data
% replace the target positions/ center position
% get all the metrics again
% save processed data

%% 
clear all
close all
clc

% Choose the right subject directory for datadir variable

% datadir = '/Users/Yekta/Documents/PhD/Lab/KinArm Data/SpringMassDamper_11/Timothy'; %zip file directory
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF04'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF05'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF06'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF07'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF08'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF15'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF16'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF30'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF31'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets\NF32'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF20'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF21'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF22'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF23'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF24'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF25'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF26'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF27'
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF33'
datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control\NF34'

% path(path,'C:\GoogleDrive\Synced Folder\Matlab\Data_Processing\NoForce4Targets')%Matlab function directory
path(path,'C:\Users\Yeki\Google Drive\Synced Folder\Matlab\Data_Processing\NoForce4Targets')

ExpData = zip_load('dir',datadir);  %creating the structure

% reading the data of the subject
fid = fopen('pat.dat','r');
Sub_Data = textscan(fid, '%s');
take_ID= strsplit(Sub_Data{1}{11},'\s*id=\s*','DelimiterType','RegularExpression');
subID=take_ID{end};

load([subID,'_Raw'])

%%
prot_block = 1;

fs = 1000; % sampling rate Hz

tot_num_trials = 0;
for r=1:prot_block
    tot_num_trials= tot_num_trials + length(ExpData(r).c3d);
end


RawData.reach_type = zeros(1,tot_num_trials);
RawData.target_num = zeros(1,tot_num_trials);
RawData.target_position = zeros(2,tot_num_trials); % 1=Xposition 2=Yposition
RawData.center_position = zeros(2,tot_num_trials); % 1=Xposition 2=Yposition


% Defining RawData values in the right trial order
for j=1:prot_block
    
    for counter = 1:length(ExpData(j).c3d)
        
        if j == 1
            trial = ExpData(1).c3d(counter).TRIAL.TRIAL_NUM;
        end
        
        Type = ExpData(j).c3d(counter).TRIAL.TP;
        RawData.target_num(trial) = mod(Type,4); %!!! CHANGE FOR NUM TARGETS !!!%this line defins the target number bassed on the tp row
        if mod(Type,4)==0
            RawData.target_num(trial) = 4;
        end
        
 %target position info
        if RawData.target_num(trial) == 3     
            RawData.target_position(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.X(8)/100;
            RawData.target_position(2,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Y(8)/100;
        elseif RawData.target_num(trial) == 4
            RawData.target_position(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.X(1)/100;
            RawData.target_position(2,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Y(1)/100;
        else
            RawData.target_position(1,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.X(RawData.target_num(trial)+1)/100;
            RawData.target_position(2,trial) =...
                ExpData(j).c3d(counter).TARGET_TABLE.Y(RawData.target_num(trial)+1)/100;
        end

     
    end
    
end

% define score/ center position
for trial=2:tot_num_trials
    RawData.center_position(:,trial) = RawData.target_position(:,trial-1);
end

%% load Processed Data
load([subID,'_Processed2'])

% get all the metrics

num_good_Trials= length(Trial.Good_Indices);
good_index = Trial.Good_Indices;

Trial.Target_num = zeros(1,num_good_Trials);
Trial.Target_position = zeros(2,num_good_Trials);
Trial.Center_position = zeros(2,num_good_Trials);
Trial.Reach_Type = zeros(1,num_good_Trials);

counter = 1;

for i = 1:length(good_index)

   Trial.Target_num(counter) = RawData.target_num(good_index(i));
    Trial.Target_position(1,counter) = RawData.target_position(1,good_index(i));
    Trial.Target_position(2,counter) = RawData.target_position(2,good_index(i));
    Trial.Center_position(1,counter) = RawData.center_position(1,good_index(i));
    Trial.Center_position(2,counter) = RawData.center_position(2,good_index(i));

    counter = counter +1;
end


% compute perpendicular error, angular error and normalized path length
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
    for i = 1:num_data_pts-1% reach_indices
        EArea_temp(i) = abs((Xr2(i+1, 1) - Xr2(i,1)) * Xr2(i,2));
    end
    
    EArea(trial) = sum(EArea_temp);
    
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
    Trial.Area(trial) = sum(Area); % it can't be right
    
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

%% Agerage Blocks

close all
%Define the Blocks

Blocks = {[1:4:120], [2:4:120], [3:4:120], [4:4:120]};

for i = 1:length(Blocks)
    [a, Blocks{i}] = ismember(Blocks{i}, Trial.Good_Indices);
    Blocks{i} = Blocks{i}(Blocks{i}>=1);
    
end

% Blocks

BlockAvg.Block_trial_num = cell(1,length(Blocks));

for i=1:length(Blocks)
        
    id= Blocks{i};
    
    %store the block trial numbers
    BlockAvg.Block_trial_num{i} = id;
    
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

%% save 
savefile = [subID,'_Raw2'];
save(fullfile(datadir,savefile),'RawData')     %saving data

savefile = [subID,'_Processed2'];    %saving data
save(savefile,'Trial')                %saving data

savefile = [subID,'_BlockData2'];    %saving data
save(savefile,'BlockAvg')                %saving data

display('...done...')

%% THE END

