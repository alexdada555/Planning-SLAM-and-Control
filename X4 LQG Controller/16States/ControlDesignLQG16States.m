close all; % close all figures
clear;     % clear workspace variables
clc;       % clear command window
format short;

%% Octave Packeges
pkg load control;

%% Mass of the Multirotor in Kilograms as taken from the CAD

M = 0.857945; 
g = 9.81;

%% Dimensions of Multirotor

L = 0.16319; % along X-axis and Y-axis Distance from left and right motor pair to center of mass

%%  Mass Moment of Inertia as Taken from the CAD

Ixx = 1.061E+05/(1000*10000);
Iyy = 1.061E+05/(1000*10000);
Izz = 2.011E+05/(1000*10000);

%% Motor Thrust and Torque Constants (To be determined experimentally)

Ktau =  7.708e-10 * 2;
Kthrust =  1.812e-07;
Kthrust2 = 0.0007326;
Mtau = 1/44.22;
Ku = 515.5;

%% Equilibrium Input

W_e = ((-4*Kthrust2) + sqrt((4*Kthrust2)^2 - (4*(-M*g)*(4*Kthrust))))/(2*(4*Kthrust))*ones(4,1);
U_e = (W_e/(Ku*Mtau));


%% Define Discrete-Time BeagleBone Dynamics

T = 0.01; % Sample period (s)- 100Hz
ADC = 3.3/((2^12)-1); % 12-bit ADC Quantization
DAC =  3.3/((2^12)-1); % 12-bit DAC Quantization

%% Define Linear Continuous-Time Multirotor Dynamics: x_dot = Ax + Bu, y = Cx + Du         

% A = 16x16 matrix
A = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, -g, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, g, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*Kthrust*W_e(1)/M, -2*Kthrust*W_e(2)/M, -2*Kthrust*W_e(3)/M, -2*Kthrust*W_e(4)/M;
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, L*2*Kthrust*W_e(1)/Ixx, L*2*Kthrust*W_e(2)/Ixx, -L*2*Kthrust*W_e(3)/Ixx, -L*2*Kthrust*W_e(4)/Ixx;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, L*2*Kthrust*W_e(1)/Iyy, -L*2*Kthrust*W_e(2)/Iyy, L*2*Kthrust*W_e(3)/Iyy, -L*2*Kthrust*W_e(4)/Iyy;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*Ktau*W_e(1)/Izz, -2*Ktau*W_e(2)/Izz, -2*Ktau*W_e(3)/Izz, 2*Ktau*W_e(4)/Izz;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/Mtau];

% B = 16x4 matrix
B = [0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     Ku, 0, 0, 0;
     0, Ku, 0, 0;
     0, 0, Ku, 0;
     0, 0, 0, Ku];

% C = 4x16 matrix
C = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];

% D = 4x4 matrix
D = [0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0];

%% Discrete-Time System

sysdt = c2d(ss(A,B,C,D),T,'zoh');  % Generate Discrete-Time System

Adt = sysdt.a; 
Bdt = sysdt.b; 
Cdt = sysdt.c; 
Ddt = sysdt.d;

%% System Characteristics

poles = eig(Adt);
% System Unstable

figure(1);
plot(poles,'*');
grid on;
title('Discrete System Eigenvalues');

cntr = rank(ctrb(Adt,Bdt));
% Fully Reachable

obs = rank(obsv(Adt,Cdt));
% Fully Observable

%% Discrete-Time Full Integral Augmaneted System 

Cr = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];   

u = size(B,2);                        % number of controlled inputs
r = size(Cr,1);                       % number of reference inputs
n = size(A,2);                        % number of states
q = size(Cr,1);                       % number of controlled outputs

Dr = zeros(q,u);

Adtaug = [Adt zeros(n,r); 
          -Cr*Adt eye(q,r)];

Bdtaug = [Bdt; 
        -Cr*Bdt];

Cdtaug = [C zeros(r,r)];

%% Discrete-Time Full State-Feedback Control
% State feedback control design with integral control via state augmentation
% X Y Z Psi are controlled outputs

Q = diag([5000,0,5000,0,5000,0,1000,0,1000,0,5000,0,0,0,0,0,5,5,5,0.4]); % State penalty
%Q = diag([5,0,5,0,5,0,10,0,10,0,5,0,0,0,0,0,1,1,1,0.1]); % State penalty
R = eye(4,4)*(10^-4);  % Control penalty

Kdtaug = dlqr(Adtaug,Bdtaug,Q,R); % DT State-Feedback Controller Gains
Kdt = Kdtaug(:,1:n);        % LQR Gains
Kidt = Kdtaug(:,n+1:end);  % Integral Gains

%% Discrete-Time Kalman Filter Design x_dot = A*x + B*u + G*w, y = C*x + D*u + H*w + v

Gdt = 1e-1*eye(n);

Rw = diag([0.1,1,0.1,1,0.1,1,1,0.1,1,0.1,1,0.1,1000,1000,1000,1000]);   % Process noise covariance matrix
Rv = diag([10^-5,10^-5,10^-5,10^-5]);     % Measurement noise covariance matrix Note: use low gausian noice for Rv

Ldt = dlqe(Adt,Gdt,Cdt,Rw,Rv);

%Hdt = zeros(size(Cy,1),size(Gdt,2)); % No process noise on measurements
%sys4kf = ss(Adt,[Bdt Gdt],Cy,[Dy Hdt],T);
%sys4kf = ss(Adt,Bdt,Cy,Dy,T);
%[kalm,Ldt] = kalman(sys4kf,Rw,Rv);  

%% Dynamic Simulation

Time = 1;
kT = round(Time/T);

X = zeros(16,kT);
Xreal = zeros(16,kT);
Xest = zeros(16,kT);
Xe = zeros(4,kT);
Y = zeros(4,kT);
e = zeros(4,kT);

U = ones(4,kT);
U(:,1) = U_e;
Ref = [0;0;0;0];

t_span = [0,T];

profile on;
PROT = profile("info");

for k = 2:kT-1

    %Reference Setting
    if k == 10/T
        Ref(1) = 0.2;
    end
    if k == 15/T
        Ref(1) = 0;
    end
    if k == 20/T
        Ref(2) = 0.2;
    end
    if k == 25/T
        Ref(2) = 0;
    end
    if k == 30/T
        Ref(3) = -0.2;
    end
    if k == 35/T
        Ref(3) = -0.2;
    end
    if k == 40/T
        Ref(4) = 90*pi/180;
    end
    if k == 50/T
        Ref(4) = 0;
    end

    %Estimation
%    Xest(:,k) = Adt*Xest(:,k-1) + Bdt*(U(:,k-1)-U_e); % No KF Linear Prediction   

%    Xest(:,k) = Xreal(:,k);  % No KF Non Linear Prediction

%    Y(:,k) = Xreal([1,3,5,11],k);
%    xkf = [Xest(:,k-1)];
%    xode = ode45(@(t,X) Quad_Dynamics(t,X,U(:,k-1)),t_span,xkf); % Nonlinear Prediction
%    Xest(:,k) = xode.y(:,end);
%    e(:,k) = [Y(:,k) - Xest([1,3,5,11],k)];
%    Xest(:,k) = Xest(:,k) + Ldt*e(:,k);

    Y(:,k) = Xreal([1,3,5,11],k);
    Xest(:,k) = Adt*Xest(:,k-1) + Bdt*(U(:,k-1)-U_e);   % Linear Prediction
    e(:,k) = [Y(:,k) - Xest([1,3,5,11],k)];
    Xest(:,k) = Xest(:,k) + Ldt*e(:,k);

    %Control
    Xe(:,k) = Xe(:,k-1) + (Ref - Xest([1,3,5,11],k));   % Integrator 
    U(:,k) = min(800, max(0, U_e - [Kdt,Kidt]*[Xest(:,k) ;Xe(:,k)])); % Constraint Saturation 

    %Simulation    
    t_span = [0,T];
    xode = ode45(@(t,X) Quad_Dynamics(t,X,U(:,k)),t_span,Xreal(:,k)); % Runge-Kutta Integration Nonlinear Dynamics
    Xreal(:,k+1) = xode.y(:,end);

%    t=t_span;
%    [dX] = Quad_Dynamics(t,Xreal(:,k),U(:,k)); % Forward Euler Integration Nonlinear Dynamics
%    Xreal(:,k+1) = Xreal(:,k) + T*dX;

%    Xreal(:,k+1) = Adt*Xreal(:,k) + Bdt*(U(:,k)-U_e);  % Fully Linear Dynamics
end

PROT = profile("info");
profile off;

Rad2Deg = [180/pi,180/pi,180/pi]';

%Plots
t = (0:kT-1)*T;
figure(2);
subplot(2,1,1);
plot(t,Xreal([1,3,5],:).*[1,1,-1]');
legend('X','Y','Z');
title('Real Position');
xlabel('Time(s)');
ylabel('Meters(m)');

subplot(2,1,2);
plot(t,Xreal([7,9,11],:).*Rad2Deg);
legend('\phi','\theta','\psi');
title('Real Attitude');
xlabel('Time(s)');
ylabel('Degrees(d)');

figure(3);
subplot(2,1,1);
plot(t,Xest([1,3,5],:).*[1,1,-1]');
legend('X_e','Y_e','Z_e');
title('Estimated Position');
xlabel('Time(s)');
ylabel('Meters(m)');

subplot(2,1,2);
plot(t,Xest([7,9,11],:).*Rad2Deg);
legend('\phi_e','\theta_e','\psi_e');
title('Estimated Attitude');
xlabel('Time(s)');
ylabel('Degrees(d)');

figure(4);
subplot(2,1,1);
plot(t,e([1,2,3],:));
legend('e_x','e_y','e_z');
title('Position prediction error');
xlabel('Time(s)');
ylabel('Error meters(m)');

subplot(2,1,2);
plot(t,e(4,:)*Rad2Deg(1));
legend('e_\psi');
title('Hedding prediction error');
xlabel('Time(s)');
ylabel('Error degrees(d)');

figure(5);
plot(t,U);
legend('U1','U2','U3','U4');
title('Inputs PWM  Signal');
xlabel('Time(s)');
ylabel('Micro Seconds(ms)');

Cpoles = eig(Adtaug - (Bdtaug*[Kdt,Kidt]));
figure(6);
plot(Cpoles,'*');
grid on;
title('Closed-Loop Eigenvalues');

Obspoles = eig(Adt - (Ldt*Cdt));
figure(7);
plot(Obspoles,'*');
grid on;
title('Closed-Loop Estimator Eigenvalues');

%% PRINT TO CONFIGURATION FILES

dlmwrite("Adt.txt", Adt,',', 0, 0);
dlmwrite("Bdt.txt", Bdt,',', 0, 0);
dlmwrite("Cdt.txt", Cdt,',', 0, 0);
dlmwrite("Ddt.txt", Ddt,',', 0, 0);
dlmwrite("Kdt.txt", Kdt,',', 0, 0);
dlmwrite("Kidt.txt", Kidt,',', 0, 0);
dlmwrite("Ldt.txt", Ldt,',', 0, 0);
dlmwrite("U_e.txt", U_e,',', 0, 0);