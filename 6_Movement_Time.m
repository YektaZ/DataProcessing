% compute the movement time for each subject on each trial and average
% across them


%% LimbFeedback 
% close all
clear all
clc


% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets';
% Subjects = {'NF04', 'NF05', 'NF06', 'NF07', 'NF08', 'NF15', 'NF16', 'NF30', 'NF31', 'NF32'}; 

datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_5';
Subjects = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014','SMD028'};

Subfile = Subjects;
num_subs = length(Subjects);

% block_trials = {[81:4:120], [81+1:4:120], [81+2:4:120], [81+3:4:120]};

Blocks = {[1:4:120], [2:4:120], [3:4:120], [4:4:120]};

num_blocks = length(Blocks);

Time = zeros(num_subs,num_blocks);

for sub = 1:num_subs
   
    T = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_Processed']);
    
    for i = 1:length(Blocks)
        [a, Blocks{i}] = ismember(Blocks{i}, T.Trial.Good_Indices);
        Blocks{i} = Blocks{i}(Blocks{i}>=1);
        
    end
    
    for block = 1:num_blocks
        
        counter = 1;
        
        for trial = Blocks{block}
            
            % id0 = T.Trial.Start_Out_Index(trial)+50;
            idf = T.Trial.End_Out_Index(trial);
            
            time(counter) = T.Trial.Time_Fixed{trial}(idf);
            
            counter = counter + 1;
        end
        
        Time(sub,block) = mean(time);
        
    end
    
end

figure
hold on
bar( mean(Time))
errorbar(mean(Time),std(Time))
title('movement time Limb feedback group')
xlabel('reach type')

display(mean(mean(Time)))

%% Cursor Feedback group

% clear all

% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control';
% Subjects = {'NF20', 'NF21', 'NF22', 'NF23', 'NF24', 'NF25', 'NF26', 'NF27', 'NF33', 'NF34'}; 

% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2';
Subjects = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

Subfile = Subjects;
num_subs = length(Subjects);

Blocks = {[81:4:120], [81+1:4:120], [81+2:4:120], [81+3:4:120]};
num_blocks = length(Blocks);

Time_Env = zeros(num_subs,num_blocks);

for sub = 1:num_subs
   
    T = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_Processed']);
    
    for i = 1:length(Blocks)
        [a, Blocks{i}] = ismember(Blocks{i}, T.Trial.Good_Indices);
        Blocks{i} = Blocks{i}(Blocks{i}>=1);
        
    end
    
    for block = 1:num_blocks
        
        counter = 1;
        
        for trial = Blocks{block}
            
            % id0 = T.Trial.Start_Out_Index(trial)+50;
            idf = T.Trial.End_Out_Index(trial);
            
            time(counter) = T.Trial.Time_Fixed{trial}(idf);
            
            counter = counter + 1;
        end
        
        Time_Env(sub,block) = mean(time);
        
    end
    
end

figure
hold on
bar( mean(Time_Env))
errorbar(mean(Time_Env),std(Time_Env))
title('movement time Cursor feedback group')
xlabel('reach type')

display(mean(mean(Time_Env)))

%% t-test
close all

[h, p] = ttest2(Time_Env(:), Time(:), 'Alpha',1e-9)

figure
hold on
bar([mean(Time_Env(:)), mean(Time(:))])
errorbar([mean(Time_Env(:)), mean(Time(:))], [std(Time_Env(:)), std(Time(:))])

%% THE END


