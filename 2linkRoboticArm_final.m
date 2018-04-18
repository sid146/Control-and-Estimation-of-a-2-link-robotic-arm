%% Test
clc,clear all, close all

% CLTI System Matrices
A = [0 0 1 0
    0 0 0 1
    2 -1 0 0
    -2 2 0 0];

B = [0 0;
    0 0;
    1 0;
    0 1];

C = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

D = [0 0;
    0 0;
    0 0;
    0 0];

% Setting up a system
sysc = ss(A,B,C,D);

% setting up time
dt = .1; % choosing a time step size
tfin = 5;
tspan = 0:dt:tfin; % time span

% discretizing the system
sysd = c2d(sysc,dt); 

Ad = sysd.a;
Bd = sysd.b;
Cd = sysd.c;
Dd = sysd.d;

% Ad = expm(A*(tfin/dt));

% Checking for cotrollability and obervability
co = ctrb(Ad,Bd);
ob = obsv(Ad,Cd);

% Initial and final coditions
x0 = [pi/4 pi/3 3 -1]';
xf = [0 0 0 0]';

%% obtaining input to go from x0 to xf
co = [];
ob = [];

% costructing the cotrollability matrix
for k = 1:tfin/dt
    co = [co (Ad^((tfin/dt)-k))*Bd];
    ob = [ob;Cd*(Ad^(k-1))];
end

% Extracting the inputs
xbar = xf - (Ad^(tfin/dt))*x0;
[U1,S1,V1] = svd(co);
uco = V1*pinv(S1)*U1'*xbar;

% Separating the inputs
u1 = uco(1:2:length(uco));
u2 = uco(2:2:length(uco));

% Simulation
xsc = x0;
for k = 1:length(u1)
    xsc(:,k+1) = Ad*xsc(:,k)+ Bd*[u1(k);u2(k)];
end

% Plotting
figure
plot(tspan,xsc,'LineWidth',2)

title('Open Loop Control steering system to stability')
xlabel('time (s)')
ylabel('States')
legend('Theta_1','Theta_2','Ang Velocity_1','Ang Velocity_2')


% Now the system has been driven from X_0 to X-F using open loop control
% using controllability matrix

%initial condition

temp = zeros(size(C,1),1);
ybar = [];
for i = 1:tfin/dt
    for j = 1:i-1   % summation
        temp = temp + C*Ad^(i-1-j)*Bd* [u1(j);u2(j)];
    end
    cx = C*xsc(:,i);    % y(t) is c*x(t)
    ybar = [ybar; cx-temp];
    temp = zeros(size(C,1),1);
end

[U2,S2,V2] = svd(ob);
u_ini = V2*pinv(S2)*U2'*ybar;

%closed loop
lambda=[-.7 .6 .9 -.8];
h=place(Ad,Bd,lambda);
xscc = x0;
for k = 1:length(u1)
    xscc(:,k+1) = (Ad-Bd*h)*xscc(:,k);
end

figure
plot(tspan,xscc,'LineWidth',2)

title('Closed Loop Control stabilizing using Feedback Law')
xlabel('time (s)')
ylabel('States')
legend('Theta_1','Theta_2','Ang Velocity_1','Ang Velocity_2')

% Now the system has been stabilised using closed loop control
% using the feedback law 

%observer
l2=[-0.5; 0.6; 0.3; -0.6];
L=place(Ad,C,l2)


yba=[];
for k = 1:length(u1)
    xscc(:,k+1) = (Ad-Bd*h)*xscc(:,k);
    yba(:,k+1)=C*xscc(:,k);



end

e=xscc-yba;
e_dot=(Ad-C*L)*e;

lambda=[-1 .06 .03 -3];
z=place(Ad,Bd,lambda);
xscn = [];
xscn(:,1)=[3.14,3.14,3,3];
xnl=[];
xnld=[];
xnl(:,1)=[3.14,3.14,3,3];
xnld(:,1)=[3.14,3.14,3,3];

for k = 1:length(u1)
    %linear closed loop estimator/cpntroller
    xscn(:,k+1) = (Ad-Bd*z)*yba(:,k);
    
    
    %cnl
    xnl(1,k+1)=yba(1,k);
    xnl(2,k+1)=yba(2,k);
    xnl(3,k+1)=yba(3,k);
    xnl(4,k+1)=yba(4,k);
    
   %cnl with parameter variations
    xnld(1,k+1)=yba(1,k)+0.5*sin(yba(2,k));
    xnld(2,k+1)=yba(2,k)+ sin(yba(1,k));
    xnld(3,k+1)=yba(3,k)+ 2*(yba(1,k));
    xnld(4,k+1)=yba(4,k) + 3*(yba(2,k));
    
end

%{figure
 plot(tspan,xnl(1,:),'k')
 
  
 hold on;
plot(tspan,xnld(1,:),'r')
title('Change in behaviour due to perturbance ')
xlabel('time (s)')
ylabel('theta_1')
legend('CNL','CNL with disturbance')

 hold on;
plot(tspan,xnld(1,:),'r')
legend('cnl','cnl with disturbance')

figure
plot(tspan,e_dot,'LineWidth',1)
title('Closed Loop Estimator/Controller')
xlabel('time (s)')
ylabel('States')
legend('Theta_1','Theta_2','Ang Velocity_1','Ang Velocity_2')
 
 figure
  plot(tspan,xscc(3,:),'r')
  hold on;
  plot(tspan,e(3,:),'k')
 legend('actual','estimated')
 title('Estimation of  Angular Velocity-1')
xlabel('time (s)')
ylabel('Angular Velocity_1')
 
 figure
  plot(tspan,xscc(2,:),'r')
  hold on;
  plot(tspan,e(2,:),'k')
 legend('actual','estimated')
  title('Estimation of Angle Theta_2')
xlabel('time (s)')
ylabel('Theta_2')
 
 figure
  plot(tspan,xscc(1,:),'r')
  hold on;
  plot(tspan,e(1,:),'k')
 legend('actual','estimated')
  title('Estimation of Angle Theta_1')
xlabel('time (s)')
ylabel('Theta_1')
 
 
 figure
  plot(tspan,xscc(4,:),'r')
  hold on;
  plot(tspan,e(4,:),'k')
 legend('actual','estimated')
 title('Estimation of Angular Velocity-2')
xlabel('time (s)')
ylabel('Angular Velocity_2')

  figure
 plot(tspan,xsc(1,:),'LineWidth',2)
 hold on;
 plot(tspan,xscc(1,:),'r','LineWidth',2)
 xlabel('time')
 ylabel('Theta_1')
 legend('open loop','closed loop')
 title('Open Loop and Closed loop control of Angle-1')
 
 figure
 plot(tspan,xsc(2,:),'LineWidth',2)
 hold on;
 plot(tspan,xscc(2,:),'r','LineWidth',2)
 xlabel('time')
 ylabel('Theta_2')
 legend('open loop','closed loop')
 title('Open Loop and Closed loop control of Angle-2')
 
 figure
  plot(tspan,xscc(1,:),'r')
  hold on;
  plot(tspan,e(1,:),'k')
 legend('actual','estimated')
 title('Estimation of the State - Theta_1')
xlabel('time (s)')
ylabel('Theta_1')
 
 
 figure
 plot(tspan,xsc(3,:),'LineWidth',2)
 hold on;
 plot(tspan,xscc(3,:),'r','LineWidth',2)
 xlabel('time')
 ylabel('Angular Velocity_1')
 legend('open loop','closed loop')
 title('Open Loop and Closed loop control of Anglular velocity-1')
 
 figure
 plot(tspan,xsc(4,:),'LineWidth',2)
 hold on;
 plot(tspan,xscc(4,:),'r','LineWidth',2)
 xlabel('time')
 ylabel('Angular Velocity_2')
 legend('open loop','closed loop')
 title('Open Loop and Closed loop control of Anglular velocity-2')
 
 figure
 plot(tspan,xnl,'LineWidth',1)
 xlabel('time')
 ylabel('States')
 legend('Theta_1','Theta_2','Ang Velocity_1','Ang Velocity_2')
 title('Estimator/Controller applied to a Non-Linear System')
%}

figure
 plot(tspan,xnl(1,:),'k')
 
  
 hold on;
plot(tspan,xnld(1,:),'r')
title('Change in behaviour due to perturbance for Theta-1 ')
xlabel('time (s)')
ylabel('theta_1')
legend('CNL','CNL with disturbance')

figure
 plot(tspan,xnl(2,:),'k')
 
  
 hold on;
plot(tspan,xnld(2,:),'r')
title('Change in behaviour due to perturbance for Theta-2 ')
xlabel('time (s)')
ylabel('theta_2')
legend('CNL','CNL with disturbance')

figure
 plot(tspan,xnl(3,:),'k')
 
  
 hold on;
plot(tspan,xnld(3,:),'r')
title('Change in behaviour due to perturbance for Angular velocity-1 ')
xlabel('time (s)')
ylabel('Ang velocity-1')
legend('CNL','CNL with disturbance')


figure
 plot(tspan,xnl(4,:),'k')
 
  
 hold on;
plot(tspan,xnld(4,:),'r')
title('Change in behaviour due to perturbance for Angular velocity-2 ')
xlabel('time (s)')
ylabel('Ang velocity-2')
legend('CNL','CNL with disturbance')





