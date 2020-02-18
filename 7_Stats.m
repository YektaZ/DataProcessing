

%%
clear all
close all
clc

% load averages

% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_5\Averages\'; 
datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_5\Averages\';
% datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_5/Averages/';
LimbF_Data = load([datadir,'LimbFeedback_Averages_succ']);

% datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\Averages\';
datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\SMD_2Targ_Env_2\Averages\';
% datadir = '/Users/yektazahed/Google Drive/Synced Folder/Experimental Data/SMD_2Targ_Env_2/Averages/';
EnvF_Data = load([datadir,'EnvFeedback_Averages_succ']);

%% 

AbsMean = [mean(EnvF_Data.BlockData.BlockAbsPerp(:)*1e2),
    mean(LimbF_Data.BlockData.BlockAbsPerp(:)*1e2)];

NormLMean = [mean(EnvF_Data.BlockData.BlockNorm(:)),
    mean(LimbF_Data.BlockData.BlockNorm(:))];

EAreaMean = [mean(EnvF_Data.BlockData.BlockEArea(:)*1e4),
    mean(LimbF_Data.BlockData.BlockEArea(:)*1e4)];

PForceMean = [mean(EnvF_Data.BlockData.BlockPeakForce(:)),
    mean(LimbF_Data.BlockData.BlockPeakForce(:))];

PSpeedMean = [mean(EnvF_Data.BlockData.BlockPeakSpeed(:)),
    mean(LimbF_Data.BlockData.BlockPeakSpeed(:))];

%%%%%

n1 = length(EnvF_Data.BlockData.BlockAbsPerp(:));
n2 = length(LimbF_Data.BlockData.BlockAbsPerp(:));

AbsSTE = [std(EnvF_Data.BlockData.BlockAbsPerp(:)*1e2/sqrt(n1)),
    std(LimbF_Data.BlockData.BlockAbsPerp(:)*1e2)/sqrt(n2)];

NormLSTE = [std(EnvF_Data.BlockData.BlockNorm(:)/sqrt(n1)),
    std(LimbF_Data.BlockData.BlockNorm(:))/sqrt(n2)];

EAreaSTE = [std(EnvF_Data.BlockData.BlockEArea(:)*1e4/sqrt(n1)),
    std(LimbF_Data.BlockData.BlockEArea(:)*1e4)/sqrt(n2)];

PForceSTE = [std(EnvF_Data.BlockData.BlockPeakForce(:)/sqrt(n1)),
    std(LimbF_Data.BlockData.BlockPeakForce(:))/sqrt(n2)];

PSpeedSTE = [std(EnvF_Data.BlockData.BlockPeakSpeed(:)/sqrt(n1)),
    std(LimbF_Data.BlockData.BlockPeakSpeed(:))/sqrt(n2)];


figure
hold on
bar(AbsMean)
errorbar(AbsMean, AbsSTE)
xticks([1 2])
xticklabels({'EFB','LFB'})
title('maximum abs perp error')
axis([0,3,0,4])
ylabel('cm')

figure
hold on
bar(NormLMean)
errorbar(NormLMean, NormLSTE)
xticks([1 2])
xticklabels({'EFB','LFB'})
title('Normalized path length')
axis([0,3,1,1.2])

figure
hold on
bar(EAreaMean)
errorbar(EAreaMean, EAreaSTE)
xticks([1 2])
xticklabels({'EFB','LFB'})
title('Enclosed Area')
axis([0,3,0,50])
ylabel('cm2')

figure
hold on
bar(PForceMean)
errorbar(PForceMean, PForceSTE)
xticks([1 2])
xticklabels({'EFB','LFB'})
title('Peak Force')
axis([0,3,0,6])
ylabel('N')

figure
hold on
bar(PSpeedMean)
errorbar(PSpeedMean, PSpeedSTE)
xticks([1 2])
xticklabels({'EFB','LFB'})
title('Peak Speed')
axis([0,3,0,0.3])
ylabel('m/s2')

%% ttests: two sample t-tests
% change to one way anova?

[h, p] = ttest2(EnvF_Data.BlockData.BlockAbsPerp(:), LimbF_Data.BlockData.BlockAbsPerp(:),'Alpha',0.05/5) % p = 2.4e-4

[h, p] = ttest2(EnvF_Data.BlockData.BlockNorm(:), LimbF_Data.BlockData.BlockNorm(:),'Alpha',0.05/5) % p = 0.093

[h, p] = ttest2(EnvF_Data.BlockData.BlockEArea(:), LimbF_Data.BlockData.BlockEArea(:),'Alpha',0.05/5) % p = 1.0e-4

[h, p] = ttest2(EnvF_Data.BlockData.BlockPeakForce(:), LimbF_Data.BlockData.BlockPeakForce(:),'Alpha',0.05/5) % p = 3.8e-3

[h, p] = ttest2(EnvF_Data.BlockData.BlockPeakSpeed(:), LimbF_Data.BlockData.BlockPeakSpeed(:),'Alpha',0.05/5) %p = 0.0276

%% one-way ANOVA: on blocks of limb feedback

[p,tbl,stats] = anova1(LimbF_Data.BlockData.BlockAbsPerp')
multcompare(stats)
[p,tbl,stats] = anova1(LimbF_Data.BlockData.BlockNorm')
multcompare(stats)
[p,tbl,stats] = anova1(LimbF_Data.BlockData.BlockEArea')
multcompare(stats)

%% compare the final scores in two groups

Limbfb_scores = [59, 52, 32, 45, 34, 54, 50, 77, 48, 70];
Envfb_scores = [45, 69, 70, 81, 62, 74, 93, 85, 54, 87];

figure
hold on
bar([1,2], [mean(Limbfb_scores), mean(Envfb_scores)])
errorbar([1,2], [mean(Limbfb_scores), mean(Envfb_scores)],...
    [std(Limbfb_scores)/length(Limbfb_scores), std(Envfb_scores)/length(Envfb_scores)])
xticklabels('Limb FB')
title('scores')
grid on

[h, p] = ttest2(Limbfb_scores, Envfb_scores)

%% THE END

