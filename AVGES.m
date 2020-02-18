function [X_avg, X_sem, Y_avg, Y_sem, Speed_avg, Speed_sem, XForce_avg, XForce_sem, YForce_avg, YForce_sem, Time_output] = AVGES(indices,Data,fs)

min_time = -100;        %%change?
max_time = 100;


for i = indices
    Time = Data.Time_Fixed{i};
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
tempXF = []; 
tempYF = [];

for j = indices
    Start_idx = find(Data.Time_Fixed{j}>=(min_time - 0.5/fs),1);
    Stop_idx  = Start_idx + numpt-5; %-1; %thorw some points out
    tempX = [tempX, Data.XPosition{j}(Start_idx : Stop_idx)];
    tempY = [tempY, Data.YPosition{j}(Start_idx : Stop_idx)];
    tempSp = [tempSp, Data.Speed{j}(Start_idx : Stop_idx)];
    tempXF = [tempXF, Data.XForce{j}(Start_idx : Stop_idx)];
    tempYF = [tempYF, Data.YForce{j}(Start_idx : Stop_idx)];
end

X_avg = mean(tempX');
X_sem = std(tempX'); %/sqrt(length(indices));
Y_avg = mean(tempY');
Y_sem = std(tempY'); %/sqrt(length(indices));
Speed_avg = mean(tempSp');
Speed_sem = std(tempSp')/sqrt(length(indices));
XForce_avg = mean(tempXF');
XForce_sem = std(tempXF')/sqrt(length(indices));
YForce_avg = mean(tempYF');
YForce_sem = std(tempYF')/sqrt(length(indices));
Time_output = min_time:1/fs:max_time;
Time_output = Time_output(1:length(tempX)); % for the last points

end




