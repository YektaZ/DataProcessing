
% compute the work done by each subject and compare across groups

clear all
close all
clc

% Define the Subjects

Group = 'LimbFeedback';
datadir_LFB = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5';
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5';
% datadir_LFB = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_5';
Subjects_LFB = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014', 'SMD028'};

% Group = 'EnvFeedback';
datadir_EFB = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% datadir_EFB = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2';
Subjects_EFB = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

num_subs_LFB = length(Subjects_LFB);
num_subs_EFB = length(Subjects_EFB);

%% compute the work done by each subject in both groups

num_blocks = 4;

Work_x_LFB = zeros(num_subs_LFB,4);
Work_y_LFB = zeros(num_subs_LFB,4);

% fig_LFB = figure;
% hold on

for sub = 1:num_subs_LFB
    
    S = load([datadir_LFB,'/',Subjects_LFB{sub},'/',Subjects_LFB{sub},'_Processed']);
    B = load([datadir_LFB,'/',Subjects_LFB{sub},'/',Subjects_LFB{sub},'_BlockData']);
    
    Blocks = B.BlockAvg.Block_trial_num;
    
    for i = 1:num_blocks
        
        dt = S.Trial.Time_Fixed{i}(2)-S.Trial.Time_Fixed{i}(1);
        
        trials = Blocks{i};
        
        wx = 0;
        wy = 0;
        
        for t = trials
            
            %compute the work done for that subject and that block
            id1 = S.Trial.Start_Out_Index(t);
            id2 = S.Trial.End_Out_Index(t);
            idx = id1:id2;
            
            wx = wx + sum(S.Trial.XForce{t}(idx).*S.Trial.XVelocity{t}(idx))*dt; % absolout?
            wy = wy + sum(S.Trial.YForce{t}(idx).*S.Trial.YVelocity{t}(idx))*dt;
        end
        
        Work_x_LFB(sub, i) = wx/length(trials); %compute the average
        Work_y_LFB(sub, i) = wy/length(trials);
        
%         plot(S.Trial.XForce{t}(idx).*S.Trial.XVelocity{t}(idx))
%         pause(0.05)
    end
end


Work_LFB = Work_x_LFB + Work_y_LFB;

%%%%%% 

Work_x_EFB = zeros(num_subs_EFB,4);
Work_y_EFB = zeros(num_subs_EFB,4);

% figure
% hold on

for sub = 1:num_subs_EFB
    
    S = load([datadir_EFB,'/',Subjects_EFB{sub},'/',Subjects_EFB{sub},'_Processed']);
    B = load([datadir_EFB,'/',Subjects_EFB{sub},'/',Subjects_EFB{sub},'_BlockData']);
    
    Blocks = B.BlockAvg.Block_trial_num;
    
    for i = 1:num_blocks
        
        dt = S.Trial.Time_Fixed{i}(2)-S.Trial.Time_Fixed{i}(1);
        
        trials = Blocks{i};
        
        wx = 0;
        wy = 0;
        
        for t = trials
            
            %compute the work done for that subject and that block
            id1 = S.Trial.Start_Out_Index(t);
            id2 = S.Trial.End_Out_Index(t);
            idx = id1:id2;
            
            wx = wx + sum(S.Trial.XForce{t}(idx).*S.Trial.XVelocity{t}(idx))*dt; 
            wy = wy + sum(S.Trial.YForce{t}(idx).*S.Trial.YVelocity{t}(idx))*dt;
        end
        
        Work_x_EFB(sub, i) = wx/length(trials); %get average
        Work_y_EFB(sub, i) = wy/length(trials);
        
%         plot(S.Trial.XForce{t}(idx).*S.Trial.XVelocity{t}(idx))
%         pause(0.05)
    end
end

Work_EFB = Work_x_EFB + Work_y_EFB;

display(Work_LFB)
display(Work_EFB)

%% significance tests
clc

[h,p] = ttest2(Work_LFB(:), Work_EFB(:))

g1 = [ones(num_subs_LFB*4,1); zeros(num_subs_EFB*4,1)]; % zero is environment feedback
g2 = [ones(num_subs_LFB, 1); 2*ones(num_subs_LFB, 1); 3*ones(num_subs_LFB, 1); 4* ones(num_subs_LFB, 1); ...
    ones(num_subs_EFB, 1); 2*ones(num_subs_EFB, 1); 3*ones(num_subs_EFB, 1); 4* ones(num_subs_EFB, 1)]; %repmat([1; 2; 3; 4], num_subs_LFB + num_subs_EFB, 1);

[p, tbl, stats] = anovan([Work_LFB(:); Work_EFB(:)], {g1, g2}, 'varnames', {'FB','reach type'});

multcompare(stats, 'Dimension', [1])

%% plot bar plot

n = sqrt(length(Work_LFB(:)));

figure
hold on
bar([1,2], [mean(Work_EFB(:)),mean(Work_LFB(:))])
errorbar([1,2], [mean(Work_EFB(:)),mean(Work_LFB(:))],...
    [std(Work_EFB(:))/n,std(Work_LFB(:))/n]);

axis([0,3,0,0.26])
xticks([1 2])
xticklabels({'EFB','LFB'})
title('average work')

%% THE END

