function [out,uin,t,B] = rdid_nonlnr_mdl_vessel(fs,tval,C,refin)
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
% B = [0 0;0 0;0 0;0.5604 0.1862;0.3717 0.0233;0.2747 0.2576];
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
x1(1) = 0.0367
x2(1) = 0.0477
x3(1) = 0.0166
x4(1) = 0.0615
x5(1) = 0.0018
x6(1) = 0.0738
u(:,1)= refin(:,1);
m11=200;
m22=250;
m33=80;
d11=70;
d22=100;
d33=50;

% c1 = 0.7363
% c2 = 0.3592
% c3 = 0.9588
% B
for i = 1:N-1
    x = [x1(i);x2(i);x3(i);x4(i);x5(i);x6(i)];
    x1(i+1) = x1(i) + x2(i)*dt;
    x2(i+1) = x2(i) + (cos(x5(i))/m11 * ((m22-m11)*x6(i)-d11*x2(i)+u(1,i))-sin(x5(i))/m11*m11/m22*((m22-m11)*x2(i)*x6(i)-d22*x4(i)))*dt;
    x3(i+1) = x3(i) + x4(i)*dt;
    x4(i+1) = x4(i) + (sin(x5(i))/m11 * ((m22-m11)*x6(i)-d11*x2(i))+cos(x5(i))/m11*m11/m22*((m22-m11)*x2(i)*x6(i)-d22*x4(i)))*dt;
    x5(i+1) = x5(i) + x6(i)*dt;
    x6(i+1) = x6(i) + ((-(m22-m11)*x2(i)*x4(i)-d33*x6(i))/m33 + u(2,i)/m33 )*dt;
    u(:,i+1) = refin(:,i+1);
end

%plot(t,x1,t,x2,t,x3,t,x4,t,x5,t,x6)
plot(t,x1,t,x3,t,x5)
title('Response of System')
xlabel('Time(/s)')
ylabel('Response')
legend('$x_1$','$x_2$','$x_3$','$\dot{x}_1$','$\dot{x}_2$','$\dot{x}_3$','Interpreter','Latex');
x = [x1;x2;x3;x4;x5;x6];
out = C * x;
uin = refin;




    