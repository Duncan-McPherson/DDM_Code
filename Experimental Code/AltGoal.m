%  The RSPKLS - Recursive Semi-Parametric Kernel Least Squares
%
%  Alternate goal of predicting force given all hyperparameters including
%  pyhsics parameters
%
%  [Fc, Y] = DK(Xi, Yi, M, beta, hp, lambda, gam, del)
%    IN
%    Xi, Yi:      X and Y identification data
%    M:           Max size of the kernel model
%    beta:        Pyhsical Parameters of the System
%    lambda:      Forgetting factor                                             
%    gam:         Regularization factor
%    del:         Error factor/jitter
%
%    OUT
%    Fc:          Force Predictions
%    Y:           Velocity Outputs
%
%   For further details, send email to jdmcpher@uvic.ca
%   2023 - Duncan McPherson

function [Fc, Y] = AltGoal(Xi, Yi, M, beta, hp, lambda, sign, del)
%% Constants
omega = hp(1); %Kernel HP #1
sig1 = hp(2);  %Kernel HP #2
sig2 = hp(3);  %Kernel HP #3

n = size(Xi,1); %Size of data
Y = zeros(n,1); %Output prediction
Fc = zeros(n,1);%Force prediction

kss = 1 + del;                              %Covariance of first data
mu = (Yi(1)-Xi(1,1:5)*beta)*kss/(sign+kss); %mean alpha values
Sigma = kss - kss^2/(sign+kss);             %alpha varaince
Q = 1/kss;                                  %Inverted Kernel Matrix

Hyp = {sig1, sig2, omega};
Xb = Xi(1,6); %Dictionary of Time Steps
m = 1;          %Size of model
%% Loop through data and determine output
for t = 1:n 
    X = Xi(t,:);
    Xl = X(1:5);
    Xnl = X(6);
    y = Yi(t,:);

    % Forget a bit
    Sigma = lambda*Sigma + (1-lambda)*(DK("LP",Xb,Xb,Hyp) + del*eye(m));
    mu = sqrt(lambda)*mu;
   
    % Predict distrubance at sample t
    kbs = DK("LP",Xb,Xnl,Hyp); kss = 1 + del;
    q = Q*kbs; % O(n^2)s
    Fc(t) = q'*mu;
    Y(t) = Fc(t) + Xl*beta;
    h = Sigma*q;
    gammat = kss - q'*kbs; gammat(gammat<0)=0;
    sf2 = gammat + q'*h; sf2(sf2<0)=0;
    sy2 = sign + sf2;

    % Include sample t and add a new basis
    Qold = Q;
    p = [q;-1];
    Q = [Q zeros(m,1);zeros(1,m) 0] + 1/gammat*(p*p'); % O(n^2)
  
    p = [h;sf2];
    mu = [mu;Fc(t)] + ((y - Xl*beta - Fc(t))/sy2)*p;
    Sigma = [Sigma h; h' sf2] - 1/sy2*(p*p');
    Xb = [Xb; Xnl]; m = m + 1; 
        
    %----- Delete a basis if necessary
    if m > M || gammat < del
        if gammat < del % To avoid roundoff error
            if gammat < del/10
                %warning('Numerical roundoff error too high, you should increase jitter noise')
            end
            criterium = [ones(1,m-1) 0];
        else
            errors = (Q*mu)./diag(Q);
            criterium = abs(errors);
        end
        % Remove element r, which incurs in the minimum error 
        [~, r] = min(criterium);            
        small = 1:m; small(r) = [];

        if  r == m   % If we must remove the element we just added, perform reduced update instead
            Q = Qold;
        else
            Qs = Q(small, r); qs = Q(r,r); Q = Q(small, small);
            Q = Q - (Qs*Qs')/qs; 
        end
        %Reduce other components
        mu = mu(small);
        Sigma = Sigma(small, small);
        Xb = Xb(small); m = m - 1;
    end
end
end