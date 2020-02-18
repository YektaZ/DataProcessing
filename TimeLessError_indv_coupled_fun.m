
% This metric computes the error between two paths

% The error for each data point is computed as the shortest distance
% between that point and the other path

function Error = TimeLessError_indv_coupled_fun(Simulations, Data)

% global sho_x sho_y
sho_x = Simulations.sho_x;
sho_y = Simulations.sho_y;

Error = zeros(1,4);

for i = 1:4
    
    N = length(Data.AVG_Data.Time{i}); 
    %%%%%%%%%%%%%%%%%     interpolate
    
    Speed_sim = sqrt(Simulations.u{i}(1,:).^2 + Simulations.u{i}(2,:).^2);
    
    if isfield(Simulations.Parameters, 'TT')
        T = Simulations.Parameters.TT;
    else
        T = Simulations.Parameters.t;
    end
    
    TT = Data.AVG_Data.Time{i}(51:N);
    dt = TT(2)-TT(1);
    
    count = 1;
    for j = 0.001:dt:T(end)
        Speed_temp(count) = interp1(T,Speed_sim,j);
        X_temp(count) = interp1(T,Simulations.X{i}(:,1),j);
        Y_temp(count) = interp1(T,Simulations.X{i}(:,2),j);
        count = count + 1;
    end
    
    id0_data = find(Data.AVG_Data.Time{i}>0, 1);
    id0_sim = find(Speed_temp>0.05, 1);
    
    numpt = length(Speed_temp); % numpt for simulations
    
    id_data = floor(linspace(id0_data, N, 100)); % only consider 100 data points
    id_sim = floor(linspace(id0_sim,numpt, 100));
    
    X_data = Data.AVG_Data.X{i}(id_data) - sho_x;
    Y_data = Data.AVG_Data.Y{i}(id_data) - sho_y;
    
    
    for pt = id_sim
        dis = sqrt(sum(([X_data,Y_data] - [X_temp(pt),Y_temp(pt)]).^2,2)); % do I need to take the sqrt?
        Error(i) = Error(i) + min(dis);
    end
    
    Error(i) = Error(i)/length(id_sim); % to get the mean
    
end

end