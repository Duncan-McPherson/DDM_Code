clear; clc; close all;

%%True Values
%Mass and Gain
J = 7.65e-3;
Ka =  6.4898;
Kt = 0.4769;
B = 0.0321;
Ku = (15:21);

%Kernel Hyperparameters
omega = 500;
gamma = 1;
theta1 = 1;
theta2 = 1;
noise = 0.1;

%Time
N0 = 110;
Tsamp = 0.001;
Tspan = (0:(6*N0)-1)*Tsamp;
Tend = 2*Tspan(end)+Tsamp;

%Create input for Sim
tmat=reshape(Tspan(1:end),6,N0);
cnst=0.01;
for j = 1:size(tmat,1)
    if rem(j,2)==0
        cnst=cnst+0.1;
    end
    if j~=7
        u_hlf((j-1)*size(tmat,2)+1:(j)*size(tmat,2))=(-1)^(j-1)*(cnst);
    else
        u_hlf((j-1)*size(tmat,2)+1:(j)*size(tmat,2)/2)=(-1)^(j-1)*(cnst);
    end
end

%Cross Validation paramters
signal = [u_hlf,-flip(u_hlf)]';
tspan0 = (0:length(signal)-1)*Tsamp;
scale = Ku(1);
test_signal = signal*scale;

sim("Machine_Learning.slx");  
Sim_Output = ans.Velocity(:,2);
Sim_Input = ans.Control(:,2);

K = Kernel_Function(omega, theta1, theta2);
gpr = GPR(Sim_Input, Sim_Output, K, noise);
[avg, gpr] = gpr.Prediction(0.5);

