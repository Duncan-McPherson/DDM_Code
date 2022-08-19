
function [d_Data,Data]=Get_Sim_Data(ODE,state0,u,tspan)
%% Get the size of the state and control
[N1,M1]=size(state0);
[N2,M2]=size(u);
Noise=0.05;
Shuffle=0;
%% Get simulation data by simulating the system using ODE113

% Determine the left hand side derivative

    for i=2:length(u)
        [t_1,y_1] = ode113(@(t_1,y_1)ODE(t_1,y_1',u(i-1,:)),tspan(1,i-1:i),state0');
        y_list(i,:)=y_1(end,:);
        d_y_list(i,:)=ODE(0,y_list(i,:),u(i,:));
        state0=y_list(i,:)';
%         y_list(i,:)
    end


%% Add some noise to the system
for i=1:N1
    Data(:,i)=y_list(:,i)+Noise*randn(size(y_list(:,i)));
end
%
for i=1:N1
    d_Data(:,i)=d_y_list(:,i)+Noise*randn(size(d_y_list(:,i)));
end

%% Shuffle the data
if Shuffle==1
    Sequence=randperm(size(Data,1));
    Data=Data(Sequence,:);
    d_Data=d_Data(Sequence,:);
end
