function [X_avg, Y_avg, X_SEM, Y_SEM, Speed_avg, Speed_SEM, Xr_Cnt] = AVGES_2(indices,Data,fs)

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
% numpt = ceil((max_time-min_time)*fs)+1;

Avg_Time = min_time:1/1000:max_time;
N = length(Avg_Time);
    

Xr_Sum = zeros(N,1);
Yr_Sum = zeros(N,1);
Speed_Sum = zeros(N,1);
Xr_Cnt = zeros(N,1);
Xr_SS = zeros(N,1);
Yr_SS = zeros(N,1);
Speed_SS = zeros(N,1);

count = 0;
for j = indices
    
    t0 = Data.Time_Fixed{j}(1);
    id_time = find(Avg_Time>=t0,1,'first');
    id_speed = find(Data.Speed{j}(id_time:end)>=0.05,1,'first')+id_time;
    
    npt = min(N, length(Data.XPosition{j}(id_speed:end))); % get the most data points you can 
    
    Xr_Sum(1:npt) = Xr_Sum(1:npt) + Data.XPosition{j}(id_speed:id_speed+npt-1);
    Yr_Sum(1:npt) = Yr_Sum(1:npt) + Data.YPosition{j}(id_speed:id_speed+npt-1);
    Speed_Sum(1:npt) = Speed_Sum(1:npt) + Data.Speed{j}(id_speed:id_speed+npt-1);
     
    count = count +1;
    
    IDX(count,:) = [id_speed, id_speed+npt-1];
    
    Xr_Cnt(1:npt) = Xr_Cnt(1:npt) + 1;
    
end

X_avg = Xr_Sum./Xr_Cnt;
Y_avg = Yr_Sum./Xr_Cnt;
Speed_avg = Speed_Sum./Xr_Cnt;

count = 0;

for j = indices
    count = count + 1;
    
   idx = IDX(count,1):IDX(count,2);
   npt = length(idx);
        
   Xr_SS(1:npt) = Xr_SS(1:npt) + ( Data.XPosition{j}(idx) - X_avg(1:npt) ).^2;
   Yr_SS(1:npt) = Yr_SS(1:npt) + ( Data.YPosition{j}(idx) - Y_avg(1:npt) ).^2;
   Speed_SS(1:npt) = Speed_SS(1:npt) + ( Data.Speed{j}(idx) - Speed_avg(1:npt) ).^2;
   
end

X_SEM = sqrt( Xr_SS./Xr_Cnt );
Y_SEM = sqrt( Yr_SS./Xr_Cnt );
Speed_SEM = sqrt( Speed_SS./Xr_Cnt );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_avg = mean(tempX');
% X_sem = std(tempX'); %/sqrt(length(indices));
% Y_avg = mean(tempY');
% Y_sem = std(tempY'); %/sqrt(length(indices));
% Speed_avg = mean(tempSp');
% Speed_sem = std(tempSp')/sqrt(length(indices));
% XForce_avg = mean(tempXF');
% XForce_sem = std(tempXF')/sqrt(length(indices));
% YForce_avg = mean(tempYF');
% YForce_sem = std(tempYF')/sqrt(length(indices));
% Time_output = min_time:1/fs:max_time;
% Time_output = Time_output(1:length(tempX)); % for the last points

end




