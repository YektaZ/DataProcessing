% compute the work done by each individual subject
% average across them and creat bar plots

% W = integral(Force*Velocity)
% W = sum(XForce*XVel)  or  W = sum(YForce*YVel)

%% computet the work for block trials of all subs

close all
clear all
clc

Group = 'LimbFeedback';
datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5';
Subjects = {'SMD006', 'SMD007', 'SMD008', 'SMD009', 'SMD010', 'SMD011', 'SMD012', 'SMD013', 'SMD014','SMD028'};

Subfile = Subjects;

num_subs = length(Subjects);

dt = 1/1000;

for sub = 1:num_subs
    
    M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData']);
    T = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_Processed']);
    
    for block = 1:4
        trials = M.BlockAvg.Block_trial_num{block};
        Work = zeros(1,length(trials));
        
        counter = 1;
        
        for t = trials
            
            id = T.Trial.Start_Out_Index(t):T.Trial.End_Out_Index(t);

                Work(counter) = sum(T.Trial.XForce{t}(id).*T.Trial.XVelocity{t}(id))*dt + ...
                    sum(T.Trial.YForce{t}(id).*T.Trial.YVelocity{t}(id))*dt; %shouldn't be negative!

            counter = counter + 1;
        end
        
        NetWork_block_LimbFB(sub,block) = mean(Work);
        
    end
end

display(NetWork_block_LimbFB)

%% 


Group = 'EnvFeedback';
datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2';
Subjects = {'SMD020', 'SMD021', 'SMD022', 'SMD023', 'SMD024', 'SMD025', 'SMD026', 'SMD027', 'SMD029', 'SMD030'};

Subfile = Subjects;

num_subs = length(Subjects);

dt = 1/1000;

for sub = 1:num_subs
    
    M = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_BlockData']);
    T = load([datadir,'/',Subjects{sub},'/',Subfile{sub},'_Processed']);
    
    for block = 1:4
        trials = M.BlockAvg.Block_trial_num{block};
        Work = zeros(1,length(trials));
        
        counter = 1;
        
        for t = trials
            
            id = T.Trial.Start_Out_Index(t):T.Trial.End_Out_Index(t);

                Work(counter) = sum(T.Trial.XForce{t}(id).*T.Trial.XVelocity{t}(id))*dt + ...
                    sum(T.Trial.YForce{t}(id).*T.Trial.YVelocity{t}(id))*dt; %shouldn't be negative!

            counter = counter + 1;
        end
        
        NetWork_block_EnvFB(sub,block) = mean(Work);
        
    end
end

display(NetWork_block_EnvFB)

%%
close all

figure
hold on
bar(1:2:8,mean(NetWork_block_LimbFB),0.3,'r')
errorbar(1:2:8,mean(NetWork_block_LimbFB), std(NetWork_block_LimbFB))
bar(2:2:8,mean(NetWork_block_EnvFB),0.3,'k')
errorbar(2:2:8,mean(NetWork_block_EnvFB), std(NetWork_block_EnvFB))
legend('LimbFB')
title('Net Work Avg Across Subs for different Blocks')

%% THE END


