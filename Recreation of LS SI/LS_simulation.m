clear; clc; close all;

%True Values
%Mass and Gain
J = 0.00765;
Ka = 6.4898;
Kt = 0.4769;
rg = 1.5915;
%B = 0.0321;
Ts = 0.001;

%Recording Values
x = [1, 2.5, 5, 7.5, 10, 12.5, 15];
LSUJ = x;
LSUB = x;
LSSJ = x;
LSSB = x;

for i = (1:length(x))
    %Set the control signal for the simulation
    Gain = x(i);    
   
    %Execute the simulink sim and get the data
    sim("Machine_Learning.slx");
    output = ans.Velocity(:,2);
    input = ans.Control(:,2);
    N=length(output);

    %Simplified LS
    %Regressor Matrix, Paramter Vector, and Output Vector
    phi = [output(1:N-1), input(1:N-1)];
    Y = output(2:N);
    LSE=((phi'*phi)\eye(2))*phi'*Y;
    p=LSE(1);
    k=LSE(2);

    %Find J and B values
    LSSJ(i) =((p-1)*Kt*Ka*Ts)/(k*log(p));
    LSSB(i) =(1-p)*Kt*Ka/k;

    %Unbiased LS
    %Regressor Matrix, Paramter Vector, and Output Vector
    phi = [output(1:N-1), input(1:N-1), -1*pv(output(1:N-1), 0.5), -1*nv(output(1:N-1), 0.5)];
    LS=((phi'*phi)\eye(4))*phi'*Y;
    P=LS(1);
    K=LS(2);
   
    %Find J and B values
    LSUJ(i) =((P-1)*Kt*Ka*Ts)/(K*log(P));
    LSUB(i) =(1-P)*Kt*Ka/K;
end    
   
%Plot outputs
subplot(1,2,1)
plot(x, LSSJ, x, LSUJ)
legend('LS','UBLS')
subplot(1,2,2)
plot(x, LSSB, x, LSUB)
legend('LS','UBLS')

%Functions
function y = pv(w, Ohm)
    y = w;
    for i = (1:length(w))
        if abs(w(i))<=Ohm
            sig = 0;
        elseif w(i)>Ohm
            sig = 1;
        else
            sig = -1;
        end
        y(i) = 1/2*(sig)*(1+sig);
    end
end
function y = nv(w, Ohm)
    y = w;
    for i = (1:length(w))
        if abs(w(i))<=Ohm
            sig = 0;
        elseif w(i)>Ohm
            sig = 1;
        else
            sig = -1;
        end
        y(i) = -1/2*(sig)*(1-sig);
    end
end