
% figure out if a subject is an outlier or not

%%
clear all
% close all
clc

% Define the Subjects

% Group = 'LimbFeedback'; 
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets'; 
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets';
% Subjects = {'NF04', 'NF05', 'NF06', 'NF07', 'NF08', 'NF15', 'NF16', 'NF30', 'NF31', 'NF32'}; 
% Avg_Data = load([datadir,'\Averages\','LimbFeedback_Averages_40b']);

Group = 'CursorFeedback';
datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Targets_control';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Targets_control'; 
Subjects = {'NF20', 'NF21', 'NF22', 'NF23', 'NF24', 'NF25', 'NF26', 'NF27', 'NF33', 'NF34'}; 
Avg_Data = load([datadir,'\Averages\','CursorFeedback_Averages_40b_10subs']);

% Group = 'MixFeedback';
% % datadir = 'C:\Users\Yeki\Google Drive\Synced Folder\Experimental Data\NoForce4Target_mix';
% datadir = 'C:\GoogleDrive\Synced Folder\Experimental Data\NoForce4Target_mix';
% Subjects = {'NF42', 'NF43', 'NF44', 'NF45', 'NF46', 'NF47', 'NF48', 'NF49', 'NF50', 'NF51'}; 
% Avg_Data = load([datadir,'\Averages\','MixFeedback_Averages_40b']);


Subfile = Subjects;

num_subs = length(Subjects);


%%

clc
B = zeros(1,num_subs);

for i = 1:4
    
    [~,b1] = hampel(Avg_Data.BlockData.BlockAbsPerp(i,:));
    [~,b2] = hampel(Avg_Data.BlockData.BlockEArea(i,:));
    [~,b3] = hampel(Avg_Data.BlockData.BlockNorm(i,:));
    
    B = B + b1+b2+b3;
end

B
outliers = (B==12);
display('Outliers: ')
display('...', Subjects{outliers});


% clc
% B = zeros(1,num_subs);
% 
% for i = 1:4
%     
%     b1 = isoutlier(Avg_Data.BlockData.BlockAbsPerp(i,:),'mean');
%     b2 = isoutlier(Avg_Data.BlockData.BlockEArea(i,:),'mean');
%     b3 = isoutlier(Avg_Data.BlockData.BlockNorm(i,:),'mean');
%     
%     B = B + b1+b2+b3;
% end
% 
% B
% outliers = (B==12);
% display('Outliers: ')
% display('...', Subjects{outliers});

%% THE END

