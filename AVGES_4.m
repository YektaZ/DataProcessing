
% find the velocity average

function [Vx_avg, Vy_avg, Xvel_SEM, Yvel_SEM, Avg_Time] = AVGES_4(indices,Data,fs)

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

Avg_Time = min_time:1/1000:max_time;
N = length(Avg_Time);
    

Xr_Sum = zeros(N,1);
Yr_Sum = zeros(N,1);
% Speed_Sum = zeros(N,1);
% Xr_Cnt = zeros(N,1);
Xvel_SS = zeros(N,1);
Yvel_SS = zeros(N,1);
% Speed_SS = zeros(N,1);

count = 0;

for j = indices
    
    count = count +1;

    id1 = find(Data.Time_Fixed{j}>=(min_time - 0.5/fs),1);
    id2 = id1+N-1;
    idx = id1:id2;
    
    IDX(count,:) = [id1, id2];
    
    Xr_Sum = Xr_Sum + Data.XVelocity {j}(idx);
    Yr_Sum = Yr_Sum + Data.YVelocity{j}(idx);
    
end

Vx_avg = Xr_Sum./length(indices);
Vy_avg = Yr_Sum./length(indices);

count = 0;

for j = indices
    
    count = count +1;
    
    idx = IDX(count,1):IDX(count,2);
    
    Xvel_SS = Xvel_SS + (Data.XVelocity{j}(idx) - Vx_avg).^2;
    Yvel_SS = Yvel_SS + (Data.YVelocity{j}(idx) - Vy_avg).^2;
    
end

Xvel_SEM = sqrt(Xvel_SS./length(indices));
Yvel_SEM = sqrt(Yvel_SS./length(indices));

end




