clc
clear all

%% 给定系统参数,参数来自 《Mode-Independent H2-Control of a DC Motor Modeled as a Markov Jump Linear System》
%% 控制系统参数
modes = 3;
A(:,:,1) = [-0.4799 5.1546 0;-3.8162 14.4723 0;0.1399 0 -0.9925]; %被控系统状态方程状态转移矩阵
A(:,:,2) = [-1.6026 9.1632 0;-0.5918 3.0317 0;0.0740 0 -0.4338];
A(:,:,3) = [0.6439 0.9178 0;-0.5056 2.4811 0;0.3865 0 0.0982];

B(:,:,1) = [5.8705 15.5010 0]'; %被控系统状态方程控制输入增益矩阵
B(:,:,2) = [10.2851 2.2282 0]';
B(:,:,3) = [0.7874 1.5302 0]';

C(:,:,1) = [1.0230 2.1100 0.9500]; %被控系统输出方程输出矩阵
C(:,:,2) = [0.9800 2.0500 1.1000];
C(:,:,3) = [1.0000 2.0000 1.0500];

D = [1.000 0.5000 -0.5000]; %被控系统输出方程控制输入增益矩阵

Pr(:,:) = [0.95 0.05 0;0.36 0.6 0.04;0.1 0.1 0.8]; % 被控系统模态转移概率

%% 跟踪系统参数
alpha = 0.15;
a = 0.4;
b = 0.5;
c = sqrt(1-a^2-b^2); %由以上参数构造三维正弦波
hat_A = [cos(alpha)+(1-cos(alpha))*a^2   (1-cos(alpha))*a*b-sin(alpha)*c  (1-cos(alpha))*a*c+sin(alpha)*b;
    (1-cos(alpha))*a*b+sin(alpha)*c  cos(alpha)+(1-cos(alpha))*b^2    (1-cos(alpha))*b*c-sin(alpha)*a;
    (1-cos(alpha))*a*c-sin(alpha)*b  (1-cos(alpha))*b*c+sin(alpha)*a  cos(alpha)+(1-cos(alpha))*c^2;];

hat_C = [1 1 1]; %被跟踪系统状态方程噪声增益矩阵

%% 给定性能指标参数正定矩阵Q、R、衰减因子gamma
[A_row,~] = size(A(:,:,1));  %size函数默认返回矩阵维度,获取被控系统状态的维数
[~,B_col] = size(B(:,:,1));  %获取被控系统控制输入的维数
[C_row,~] = size(C(:,:,1));  %获取被控系统输出的维数
[hat_A_row,~] = size(hat_A); %获取被跟踪系统状态的维数

Q(:,:,1) = 10*eye(C_row); %系统跟踪误差权重矩阵
Q(:,:,2) = 10*eye(C_row);
Q(:,:,3) = 10*eye(C_row);

R = 0.2*[1 1 1];  %系统控制输入权重矩阵

gamma = 0.99;  %衰减因子gamma

%% 根据[x r]增广系统参数生成新的系统参数
for mode=1:modes
    tilde_A(:,:,mode) = blkdiag(A(:,:,mode),hat_A);  %blkdiag函数根据已有矩阵生成准对角矩阵
    tilde_B(:,:,mode) = [B(:,:,mode)' zeros(B_col,hat_A_row)]';
    tilde_C(:,:,mode) = [C(:,:,mode) -hat_C];
end

%% 直接迭代参数
P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % 给定解的初始值
sigmma_P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %定义按概率加权求和矩阵
S_u(:,:,:) = zeros(B_col,A_row+hat_A_row,modes);  % 给定控制增益K_u的初始值
norm = 0;
episodes = 2000; %给定求解迭代次数

for episode = 1:episodes
    for mode = 1:modes
        sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3); %由当前P构造sigmma_P
        P(:,:,mode) = tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) + gamma*tilde_A(:,:,mode)'*sigmma_P(:,:,mode)*tilde_A(:,:,mode) - (tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode))*inv(D(mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_B(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode) + R(mode))*(D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P(:,:,mode)*tilde_A(:,:,mode));
    end
    norm(episode) = trace(P(:,:,1)'*P(:,:,1) + P(:,:,2)'*P(:,:,2) + P(:,:,3)'*P(:,:,3));
end

for mode = 1:modes
    sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3); 
    S_u(:,:,mode) = -inv(gamma*tilde_B(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode) + R(mode))*(D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P(:,:,mode)*tilde_A(:,:,mode));
end

figure(1)
plot(0:1:(episodes-1),norm)