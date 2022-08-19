clear; clc; close all;

%%True Values
%Mass and Gain
J = 7.65e-3;
Ka = 6.4898;
Kt = 0.4769;
%B = 0.0321;
rg = 1.5915;

%Recording Values
Ku = (5:2.5:15);

%Friction Values
Tp_stat = 3.647;
Tn_stat = -2.702;
Tp_coul = 2.205;
Tn_coul = -1.656;
B = 0.0321;
Om1p = 2;
Om2p = 2;
Om1n = -2;
Om2n = -2;
Tp_visc = 0;
Tn_visc = 0;

%Hyperparameters
gamma = 0.1;
theta1 = 1;
theta2 = 0.25;
omega = 500;
%Frequency of Tool
omega = 1;

%ODE(Keivan Code)
state0 = 0.0050;

%Machine Learning Variables
Sim_Output = 0;
Outputs = zeros(length(Ku),1);
Sim_Input = 0;

%Time
N0 = 110;
Tsamp = 0.001;
Tspan = (0:(6*N0)-1)*Tsamp;
Tend = 2*Tspan(end)+Tsamp;

%Plots
% Keivan_J = zeros(1,5);
% Keivan_B = zeros(1,5);
% Keivan_UJ = zeros(1,5);
% Keivan_UB = zeros(1,5);
% Duncan_J = zeros(1,5);
% Duncan_B = zeros(1,5);
% Duncan_UJ = zeros(1,5);
% Duncan_UB = zeros(1,5);
Duncan_KJ = zeros(1,5);
Duncan_KB = zeros(1,5);

%% Simulation

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

%% Different Gains
for i = (1:length(Ku))
    %Scale Input
    signal = [u_hlf,-flip(u_hlf)]'*Ku(i);
    T0 = signal*Kt*Ka;
    tspan0 = (0:length(signal)-1)*Tsamp;
%     figure(4);
%     plot(tspan0, signal);
%     title('Control Signal')
%     xlabel('Voltage [V]')
%     ylabel('Time [sec]')
       
    %Execute the simulations and get the data
    sim("Machine_Learning.slx");  
    Sim_Output = ans.Velocity(:,2);
    Sim_Input = ans.Control(:,2);
%     %T0 = ans.Torque(:,2);
%     [dData0,Data0]=Get_Sim_Data(@(t,omga,T0)Drive_ODE(t,omga,T0,J,B,Tp_stat,Om1p,Tp_coul,Om2p,Tp_visc,Tn_stat,Om1n,Tn_coul,Om2n,Tn_visc),state0,T0,tspan0); 

    %Graph differences
    output = Sim_Output;
    %Outputs(i,1) = output;
    input = Sim_Input;
%     figure(1);
%     subplot(2,3,i)
%     plot(tspan0, Sim_Output, '--r',tspan0, Data0);
%     sgtitle('Simulink and ODE113 Simulations') 
%     legend('Simulink Velocity', 'ODE Velocity')
%     xlabel('Time [sec]')
%     ylabel('Velocity [rad/sec]')
%     figure(2);
%     subplot(2,3,i)
%     plot(tspan0, Data0-Sim_Output);
%     sgtitle('Error between Simulink and ODE113 Velocities')  
%     xlabel('Time [sec]')
%     ylabel('Velocity [rad/sec]')
    
%     %% Keivan's Code
%     %LS fitting
%     X=[Data0(1:end-1),signal(1:end-1)];
%     y=Data0(2:end);
% 
%     PHI=X(:,1:2);
%     LS_theta= (PHI'*PHI)\(PHI'*y);
% 
%     Keivan_J(i)=(LS_theta(1)-1)*Kt*Ka*Tsamp/LS_theta(2)/log(LS_theta(1));
%     Keivan_B(i)=(1-LS_theta(1))*Kt*Ka/LS_theta(2);
% 
%     %unbiased LS fitting
%     U_PHI=[X(:,1:2),0.5*(sign(X(:,1))+1),-0.5*(sign(X(:,1))-1)];
%     ULS_theta= (U_PHI'*U_PHI)\(U_PHI'*y);
% 
%     Keivan_UJ(i)=(ULS_theta(1)-1)*Kt*Ka*Tsamp/ULS_theta(2)/log(ULS_theta(1));
%     Keivan_UB(i)=(1-ULS_theta(1))*Kt*Ka/ULS_theta(2);
%     
    %% Duncan's Code
    %Simplified LS
    %Regressor Matrix, Paramter Vector, and Output Vector
    Y = output(2:end);
    phi = [output(1:end-1), input(1:end-1)];
    ls = (phi'*phi)\(phi'*Y);

    %Find J and B values
    Duncan_J(i) = ((ls(1)-1)*Kt*Ka*Tsamp)/ls(2)/log(ls(1));
    Duncan_B(i) = (1-ls(1))*Kt*Ka/ls(2);

    %Unbiased LS
    %Regressor Matrix, Paramter Vector, and Output Vector
    u_phi = [phi, -(0.5*sign(output(1:end-1)).*(1+sign(output(1:end-1)))), -(-0.5*sign(output(1:end-1)).*(1-sign(output(1:end-1))))];
    uls_inv=pinv(u_phi)*Y;
    uls=(u_phi'*u_phi)\(u_phi'*Y);
    
    %Find J and B values
    Duncan_UJ(i) = ((uls(1)-1)*Kt*Ka*Tsamp)/uls(2)/log(uls(1));
    Duncan_UB(i) = (1-uls(1))*Kt*Ka/uls(2);

    %% Reproducing Kernel Hilbert Spaces (RKHS)  
    %Building a Kernel function
    K = zeros(length(output)-1,length(input)-1);
    for k=(1:length(output)-1)
        for j=(1:length(input)-1)
            K(k,j) = exp((-1*abs(k-j)^2)/2*theta1)*exp(-1*(sin((pi()*abs(k-j))/omega))^2/(2*theta2));
        end
    end
    u_phi = [output(1:end-1), input(1:end-1), -(0.5*sign(output(1:end-1)).*(1+sign(output(1:end-1)))), -(-0.5*sign(output(1:end-1)).*(1-sign(output(1:end-1))))];
    Niner = [K+(gamma\eye(length(output)-1)), ones(length(output)-1, 1), u_phi;
             ones(1,length(output)-1),        zeros(1),                  zeros(1,4);
             transpose(u_phi),                zeros(4,1),                zeros(4)];
    Threer = [output(2:end);0;zeros(4,1)];
    ls_svr = (Niner'*Niner)\(Niner'*Threer);
    J_hat = ls_svr(1321);
    B_hat = ls_svr(1322);
    
    Duncan_KJ(i) = ((J_hat-1)*Kt*Ka*Tsamp)/B_hat/log(J_hat);
    Duncan_KB(i) = (1-J_hat)*Kt*Ka/B_hat;
end    
    
% %% Plot output
% figure(3);
% subplot(2,2,1)
% sgtitle('Inertia and Damping of Values of Simulink and ODE113 Simulations') 
% plot(Ku, Duncan_KJ)%Ku, Duncan_J, Ku, Duncan_UJ,
% legend('SIM-LS','SIM-UBLS', 'LS-SVR')
% title('Simulink J')
% xlabel('Gain [V/V]')
% ylabel('Inertia [kg*m^2]')
% subplot(2,2,2)
% plot(Ku, Duncan_KB) %Ku, Duncan_B, Ku, Duncan_UB, 
% legend('SIM-LS','SIM-UBLS', 'LS-SVR')
% title('Simulink B')
% xlabel('Gain [V/V]')
% ylabel('Damping [kg*m^2/sec]')
% subplot(2,2,3)
% plot(Ku, Keivan_J, Ku, Keivan_UJ)
% legend('ODE-LS','ODE-UBLS')
% title('ODE J')
% xlabel('Gain [V/V]')
% ylabel('Inertia [kg*m^2]')
% subplot(2,2,4)
% plot(Ku, Keivan_B, Ku, Keivan_UB)
% legend('ODE-LS','ODE-UBLS')
% title('ODE B')
% xlabel('Gain [V/V]')
% ylabel('Damping [kg*m^2/sec]')

    
    