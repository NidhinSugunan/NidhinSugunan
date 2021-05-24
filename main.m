%% Hypnosis control system to administer intravenous anesthetic agent, a comparative study using LQR, PI and MPC control strategies
% Nidhin Sugunan (5026040) 

clc; clear; close all;

%% System Model
BW = 60; % body weight
age = 40;
V1 = 1.72*(BW^0.71)*(age^(-0.39));
V2 = 3.32*(BW^0.61);
V3 = 266;
V4 = V1/100;
K1 = 0.0595*(BW^0.75);
K2 = 0.0969*(BW^0.62);
K3 = 0.0889*(BW^0.55);
K4 = 0.12;

% System matrix
A = [-(K1+K2+K3+K4)/V1     K2/V1   K3/V1   K4/V1;
       K2/V2               -K2/V2  0       0;
       K3/V3               0       -K3/V3  0;
       K4/V4               0       0       -K4/V4];

B = [1/V1;  0;  0;  0];

C = [0  0   0   1];

Mcon = ss(A,B,C,0); % continuous model

% Discrete model
Minfo = stepinfo(Mcon);
Mrise = Minfo.RiseTime;
Ts = Mrise/15; % rule of thumb is to take 1/10 of rise time
M = c2d(Mcon,Ts,'zoh'); % Discrete model

figure
pzmap(M) 

% LTI system definition
LTI.A = M.A;
LTI.B = M.B;       
LTI.C = M.C;
LTI.D = M.D;
LTI.Bd = [1;  1;  0;  0];
LTI.Cd = 1;
LTI.x0 = [0; 0;  0;  0];
LTI.d = 0.5;
LTI.yref = 3.6; 

% Definition of system dimension
dim.nx = 4;     % state dimension
dim.nu = 1;     % input dimension
dim.ny = 1;     % output dimension
dim.nd = 1;     % disturbance dimension
dim.N = 20;     % prediction horizon

% LTI extended system definition
LTIe.A = [LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd)];
LTIe.B = [LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C = [LTI.C LTI.Cd];
LTIe.Bd = [1;  1;  0;  0]; %[0; zeros(dim.nx-1,dim.nd)]; 
LTIe.Cd = ones(dim.ny,dim.nd);
LTIe.x0 = [LTI.x0; LTI.d];
LTIe.yref = LTI.yref; 

% Definition of extended system dimension
dime.nx = 5;     % state dimension
dime.nu = 1;     % input dimension
dime.ny = 1;     % output dimension
dime.nd = 1;     % disturbance dimension
dime.N = dim.N;     % prediction horizon

% Definition of quadratic cost function
weight.Q = diag([2, 2, 5, 150]);                                  % weight on output
weight.R = 1;                                                   % weight on input
weight.P = idare(LTI.A , LTI.B , weight.Q , weight.R , [], []);   % terminal cost
weighte.Q = blkdiag(weight.Q,zeros(dim.nd));            % weight on output
weighte.R = weight.R;                                   % weight on input
weighte.P = blkdiag(weight.P,zeros(dim.nd));            % terminal cost


T = 100;     % simulation horizon

%% Offset-free MPC from output
predmode = predmodgen(LTIe,dime);         % Generation of prediction model 
[He,he] = costgen(predmode,weighte,dime);    % Generate cost

% Kalman weight matrices
Qkal = diag([0.5, 0.5, 0.5, 0.5, 0.5]); 
Rkal = 50; 

% Discrete Riccati eqn to find Kalmain gain
[Pkal, Kkal, poles, info] = idare(LTIe.A',LTIe.C',Qkal,Rkal,[], []);
K = ss(LTIe.A-Kkal'*LTIe.C, [LTIe.B Kkal'],  LTIe.C,  [], Ts);

% Initial State and Signals
xhat = zeros(dime.nx, T);
yhat = zeros(dime.ny, T);
xe = zeros(dime.nx, T);
ue = zeros(dime.nu, T);
ye = zeros(dime.ny, T);
dhat = zeros(dim.nd, T);

xe(:,1) = LTIe.x0;
ye(:,1) = LTIe.C*LTIe.x0;
yrefe = 3.6*ones(T+1,dim.ny);
yref = yrefe;
% Stage Cost and Final Cost weights
Q = diag([2, 2, 20, 20]);                             
R = 10;   
P = idare(LTI.A, LTI.B, Q, R, [], []);

for k =1:T
    
    x0hat = xhat(:,k);
    dhat =  xhat(end,k);
    x0 = xe(:,k);
    [A,b] = constraintgen(dime, predmode, x0hat);
    
    % Compute optimal ss (online, at every iteration)
    eqconstraints = eqconstraintsgen(LTI,dim,dhat);
    [xr,urKe] = optimalss(LTI,dim,weight,[],eqconstraints); 
    xrKe = [xr; dhat];
    
    uostar = sdpvar(dime.nu*dime.N,1);                             % define optimization variable
    options = sdpsettings('verbose',0);
    Constraint= A*uostar<=b;                                       % define constraints
    Objective = 0.5*uostar'*He*uostar+(he*[x0hat; xrKe; urKe])'*uostar;% define cost function
    optimize(Constraint,Objective,options);                      % solve the problem
    uostar = value(uostar);        
    
    % Select the first input only
    ue(:,k)= uostar(1:dim.nu);
    
    % Compute the state/output evolution
    xe(:,k+1) = LTIe.A*xe(:,k) + LTIe.B*ue(:,k);
    ye(:,k+1) = LTIe.C*xe(:,k+1);
    l(k) = 0.5*((xe(:,k+1)-xrKe)'*weighte.Q*(xe(:,k+1)-xrKe)+((ue(:,k)-urKe)'*weighte.R*(ue(:,k)-urKe)));
    vff(k) = Objective;
    vf(k) = 0.5*((xe(:,k+1)-xrKe)'*weighte.P*(xe(:,k+1)-xrKe));
    clear u_uncon
    
    
        if (k > 50) % disturbance
        ye(:,k+1) = ye(:,k+1) + 0.5;
        end
    
    % Update extended-state estimation
    xhat(:,k+1) = LTIe.A*xhat(:,k)+LTIe.B*ue(:,k)-Kkal'*(yhat(:,k)-ye(:,k));                
    yhat(:,k+1) = LTIe.C*xhat(:,k);
    dhat(:,k) = xhat(end,k);
    
end

%% LQR controller

% Augmented System
Aaug = [M.A , zeros(size(M.A,1),1);...
        zeros(1,size(M.A,1)) , 1    ]; 
Baug = [M.B ; 0];
Caug = [M.C , 1];
SYSdaug = ss(Aaug,Baug,Caug,[],Ts);

bolus = 2*BW/2; % 2 mgs divided equally into 2 compartiments
x0 = [ bolus/V1 ; bolus/V2; 0 ; 0];

% Definition of system dimension
dim_aug.nx=5;     % state dimension
dim_aug.ny=1;     % output dimension
dim_aug.nu=1;     % input dimension
dim_aug.N=10;     % horizon

% define Kalman weight matrices
R_aug = 35 ;
Q_aug = diag([0.2 , 0.1 , 0.1 , 25, 1]);

% Solve Discrete Riccati Equation
[Pkalman_aug , K_kalman_aug , poles_aug, info] = idare(SYSdaug.A',SYSdaug.C',Q_aug,R_aug,[], []);
KALMAN_aug = ss(SYSdaug.A-K_kalman_aug'*SYSdaug.C  ,[SYSdaug.B K_kalman_aug'] ,  SYSdaug.C  ,  [SYSdaug.D 0]  , Ts);



% LQR controller 
Q = diag([1 1 1 100]); 
R = 0.01;
[Klqr,~,Eigen7_1] = dlqr(M.A, M.B, Q, R, []);

SYSdlqr = ss(M.A-M.B*Klqr,M.B,M.C,0,Ts);
[y,~]=step(SYSdlqr);
FF = 1/y(end);

%% Plot

simOut = sim('PIDnLQR'); % LQR and PID controller designed using SIMULINK
iplqr = simOut.ulqr;
ippid = simOut.upid;
oplqr = simOut.ylqr;
oppid = simOut.ypid;

%% Plot PID vs MPC comparision
t11 = Ts*(0:1:T-1);
t12 = Ts*(0:1:T);
figure
    tiledlayout(2,1)
    nexttile
    plot(iplqr)
    hold on
    plot(ippid)
    hold on
    plot(t11,ue)
    xlim([0 4000])
    title('\textbf{Performance comparison of PID and MPC}', 'interpreter','latex')
    xlabel('$t$[sec]','interpreter','latex')
    ylabel('Input [$mg/kg/h$]','interpreter','latex')
    legend({'$u_{LQR}$','$u_{PID}$','$u_{MPC}$'},'interpreter','latex');
    nexttile
    xlim([0 4000])
    plot(oplqr,'m')
    hold on
    plot(oppid,'g')
    hold on
    plot(t12,ye,'b')
    hold on
    plot(t12,yref,'--r')
    xlim([0 5259])
    legend({'$y_{LQR}$','$y_{PID}$','$y_{MPC}$','$y_{ref}$'},'interpreter','latex');
    xlabel('$t$[sec]','interpreter','latex')
    ylabel('Output [$\mu g/mL$]','interpreter','latex')
    
