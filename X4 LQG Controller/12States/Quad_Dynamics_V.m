function [dX] = Quad_Dynamics_V(t,X,U)
%% Mass of the Multirotor in Kilograms as taken from the CAD
M = 1.157685; 
g = 9.81;

%% Dimensions of Multirotor
L = 0.16319; % along X-axis and Y-axis Distance from left and right motor pair to center of mass

%%  Mass Moment of Inertia as Taken from the CAD
% Inertia Matrix and Diagolalisation CAD model coordinate system rotated 90 degrees about X
Ixx = 1.129E+05/(1000*10000);
Iyy = 1.129E+05/(1000*10000);
Izz = 2.033E+05/(1000*10000);

%% Motor Thrust and Torque Constants (To be determined experimentally)
Ktau =  7.708e-10 * 2;
Kthrust =  1.812e-07;
Kthrust2 = 0.0007326;
Mtau = (1/44.22);
Ku = 515.5;

%% Air resistance damping coeeficients
Dxx = 0.01212;
Dyy = 0.01212;
Dzz = 0.0648;                          

%% X = [x xdot y ydot z zdot phi p theta q psi r w1 w2 w3 w4]
X = num2cell(X);
[x xdot y ydot z zdot phi p theta q psi r w1 w2 w3 w4] = X{:};

%% Initialise Output and state vectors
dX = zeros(16,1);

%% Motor Forces and Torques
F = zeros(4,1);
T = zeros(4,1);

F = Kthrust*[w1^2;w2^2;w3^2;w4^2] + Kthrust2*[w1;w2;w3;w4];
T = [Ktau,-Ktau,-Ktau,Ktau]*[w1^2;w2^2;w3^2;w4^2];

Fn = sum(F);
Tn = sum(T);

%% First Order Direvatives dX = [xdot ydot zdot phidot thetadot psidot]
dX([1,3,5]) = [xdot;ydot;zdot];

dX([7,9,11]) = [1,(sin(phi)*tan(theta)),cos(phi)*tan(theta);0,cos(phi),-sin(phi);0,sin(phi)/cos(theta),cos(phi)/cos(theta)] * [p;q;r];               

%% Second Order Direvatives: dX = [xddot yddot zddot pdot qdot rdot]
dX([2,4,6]) = (Fn/M)*[cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);-cos(phi)*cos(theta)] + [0;0;g] - [Dxx/M,0,0;0,Dyy/M,0;0,0,Dzz/M]*[xdot;ydot;zdot]; 

dX([8,10,12]) = [(L/Ixx)*((F(1)+F(2)) - (F(3)+F(4)));(L/Iyy)*((F(1)+F(3)) - (F(2)+F(4)));Tn/Izz] - [(((Izz-Iyy)/Ixx)*(r*q));(((Izz-Ixx)/Iyy)*(p*r));(((Iyy-Ixx)/Izz)*(p*q))];

%% Motor Dynamics: dX = [w1dot w2dot w3dot w4dot], U = PWM Pulse Width between (0-1000)us
dX(13:end) = -(1/Mtau)*[w1;w2;w3;w4] + Ku*U;
end