function [out,uin,t,B] = rdid_nonlnr_mdl_msd(fs,tval,C,refin)
%  Date-Driven Flatness Output Searching and ADRC Control of 
%  Underactuated Nonlinear System 
%
%  Date : 07 - 12 - 2019
%  Shangjie Ma
%  --------------------------------------------------------------------
%  -Nonlinear System Case Model - 4D Minimum-phase case
%   - Parameter Def.
%   - Open-loop Matrix Def. 
%   - Stability & Controller
%   - ODE Solution

% Parameter Def.
%========================================================================
dt = 1/fs;
t = (0:dt:tval-dt)'; % time vector
N = size(t,1);

% Open-loop Matrix Def.
%========================================================================
% A = [0  0  0  1  0  0;
%      0  0  0  0  1  0;
%      0  0  0  0  0  1;
%   -49.4521 -34.4652 14.3340 -0.6580 -0.2186 -0.1607;
%   -30.4147 -47.2673 27.2193 -0.3723 -0.1025 -0.1753;
%   -42.5996 -30.3623 7.1132  -0.8928 -0.4002 -0.5570];
% 

% 
% D = [0 0;0 0;0 0];
% 
% rank(ctrb(A,B))
% eig(A)
% Stability & Controller
%========================================================================
% LM_LQR_Q =diag([1 1 1 1]);
% LM_LQR_R = 1;
% LM_LQR_K =lqr(A,B,LM_LQR_Q,LM_LQR_R);
% Ke = LM_LQR_K;

% % ODE Solution
% %========================================================================
% x1(1) = 0.01*rand
% x2(1) = 0.01*rand
% x3(1) = 0.01*rand
% x4(1) = 0.01*rand
% x5(1) = 0.01*rand
% x6(1) = 0.01*rand
x1(1) = 0.0066
x2(1) = 4.6711e-05
x3(1) = 0.0093
x4(1) = 0.0057
x5(1) = 0.0064
x6(1) = 0.0045

u(:,1)= refin(:,1);

% c1 = 0.7363
% c2 = 0.3592
% c3 = 0.9588
% c1 = 0
% c2 = 0
% c3 = 0
m1 = 100;
m2 = 200;
m3 = 300;
k1 = 100;
k2 = 150;
k3 = 200;
c1 = 20;
c2 = 20;
c3 = 20;
b11 = 50;
b12 = 15;
b21 = 22;
b22 = 13;
b31 = 16;
b32 = 38; 
B = [0 0;0 0;0 0;b11/m1 b12/m1;b21/m2 b22/m2;b31/m3 b32/m3]
S = -inv([B(5,1) B(6,1);B(5,2) B(6,2)])*[B(4,1);B(4,2)]
for i = 1:N-1
    x = [x1(i);x2(i);x3(i);x4(i);x5(i);x6(i)];
    x1(i+1) = x1(i) + x4(i)*dt;
    x2(i+1) = x2(i) + x5(i)*dt;
    x3(i+1) = x3(i) + x6(i)*dt;
    x4(i+1) = x4(i) + ((-k1-k2)/m1 * x1(i) + k2/m1 * x2(i) - c2/m1 * x4(i) + c2/m1 * x5(i) - c1/m1*abs(x4(i))*x4(i) + b11/m1*u(1,i) + b12/m1 * u(2,i) )*dt;
    x5(i+1) = x5(i) + (k2/m2 * x1(i) - (k2+k3)/m2 * x2(i) + k3/m2*x3(i) + c2/m2 * x4(i) - c2/m2 * x5(i) -c3/m2 * x5(i) + c3/m2*x6(i) + b21/m2*u(1,i) + b22/m2 * u(2,i) )*dt;
    x6(i+1) = x6(i) + (k3/m3 * x2(i) - k3/m3* x3(i) + abs(x5(i))*c3/m3*x5(i)-abs(x5(i))*c3/m3*x6(i)+b31/m3*u(1,i) + b32/m3 * u(2,i))*dt;
    u(:,i+1) = refin(:,i+1);
end

plot(t,x1,t,x2,t,x3,t,x4,t,x5,t,x6)
%title('Response of System')
xlabel('Time(/s)')
ylabel('Response')
legend('$x_1$','$x_2$','$x_3$','$\dot{x}_1$','$\dot{x}_2$','$\dot{x}_3$','Interpreter','Latex');
x = [x1;x2;x3;x4;x5;x6];
out = C * x;
uin = refin;




    