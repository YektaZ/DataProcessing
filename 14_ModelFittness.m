
clear all
close all
clc


% Group = 'LimbFeedback';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5';
% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5';
datadir_LFB = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_5';
Subjects_LFB = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014', 'SMD028'};

% Group = 'EnvFeedback';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
datadir_EFB = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2';
Subjects_EFB = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

num_subs_LFB = length(Subjects_LFB);
num_subs_EFB = length(Subjects_EFB);


sho_x = 0;
sho_y = -0.4;

%% load the model solutions

addpath('/Users/yektazahed/Google Drive/Synced Folder/Matlab/Spring_Mass_Damper')

Environment_model = load('Simulations_Spline_N10.mat');
Coupled_model = load('Simulations_Coupled_dT.mat');


Coupled_model.Sol.Parameters = Coupled_model.Parameters;

Coupled_model.Sol.sho_x = 0;
Coupled_model.Sol.sho_y = -0.4;
Environment_model.sho_x = 0;
Environment_model.sho_y = 0;


%% model comparisons for Environment model and Coupled model

%% Timeless RMS

% each group is compared seperately

figure
hold on

Error_Env = [];
Error_Coupled = [];

for sub = 1:num_subs_LFB
    
    Data = load([datadir_LFB,'/',Subjects_LFB{sub},'/',Subjects_LFB{sub},'_AVG_Data_']);
    
    S = load([datadir_LFB,'/',Subjects_LFB{sub},'/',Subjects_LFB{sub},'_Processed']);
    B = load([datadir_LFB,'/',Subjects_LFB{sub},'/',Subjects_LFB{sub},'_BlockData']);
    
    Blocks = B.BlockAvg.Block_trial_num;
    
    for i = 1:length(Blocks)
        plot(Environment_model.X{i}(:,1)-sho_x,Environment_model.X{i}(:,2)-sho_y,'b')
        plot(Coupled_model.Sol.X{i}(:,1),Coupled_model.Sol.X{i}(:,2),'b')
        plot(Data.AVG_Data.X{i}-sho_x, Data.AVG_Data.Y{i}-sho_y,'k--')
    end
    legend('Env_dF','Coupled dT','Data');
    axis equal
    
    Error_Env = [Error_Env; TimeLessError_indv_coupled_fun(Environment_model, Data)];
    Error_Coupled = [Error_Coupled; TimeLessError_indv_coupled_fun(Coupled_model.Sol, Data)];
    
    
end
display(Error_Env)
display(Error_Coupled)

[h,p] = ttest2(Error_Env(:),Error_Coupled(:))

%%%%%%%%% EFB
figure
hold on

Error_Env = [];
Error_Coupled = [];

for sub = 1:num_subs_EFB
    
    Data = load([datadir_EFB,'/',Subjects_EFB{sub},'/',Subjects_EFB{sub},'_AVG_Data_']);
    
    S = load([datadir_EFB,'/',Subjects_EFB{sub},'/',Subjects_EFB{sub},'_Processed']);
    B = load([datadir_EFB,'/',Subjects_EFB{sub},'/',Subjects_EFB{sub},'_BlockData']);
    
    Blocks = B.BlockAvg.Block_trial_num;
    
    for i = 1:length(Blocks)
        plot(Environment_model.X{i}(:,1)-sho_x,Environment_model.X{i}(:,2)-sho_y,'b')
        plot(Coupled_model.Sol.X{i}(:,1),Coupled_model.Sol.X{i}(:,2),'b')
        plot(Data.AVG_Data.X{i}-sho_x, Data.AVG_Data.Y{i}-sho_y,'k--')
    end
    legend('Env_dF','Coupled dT','Data');
    axis equal
    
    Error_Env = [Error_Env; TimeLessError_indv_coupled_fun(Environment_model, Data)];
    Error_Coupled = [Error_Coupled; TimeLessError_indv_coupled_fun(Coupled_model.Sol, Data)];
    
    
end
display(Error_Env)
display(Error_Coupled)

[h,p] = ttest2(Error_Env(:),Error_Coupled(:))




%% THE END


