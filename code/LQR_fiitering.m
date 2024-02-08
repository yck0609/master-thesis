clc
clear all

%% 给定系统参数,参数来自 《Mode-Independent H2-Control of a DC Motor Modeled as a Markov Jump Linear System》
%% 控制系统参数
modes = 3; %系统模态数量
A(:,:,1) = [-0.4799 5.1546 0;-3.8162 14.4723 0;0.1399 0 -0.9925]; %被控系统状态方程状态转移矩阵
A(:,:,2) = [-1.6026 9.1632 0;-0.5918 3.0317 0;0.0740 0 -0.4338];
A(:,:,3) = [0.6436 0.9178 0;-0.5056 2.4811 0;0.3865 0 0.0982];

B(:,:,1) = [5.8705 15.5010 0]'; %被控系统状态方程控制输入增益矩阵
B(:,:,2) = [10.2851 2.2282 0]';
B(:,:,3) = [0.7874 1.5302 0]';

F(:,:,1) = [0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0]; %被控系统状态方程噪声增益矩阵
F(:,:,2) = [0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0];
F(:,:,3) = [0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0];

C(:,:,1) = [1.0230 2.1100 0.9500]; %被控系统输出方程输出矩阵
C(:,:,2) = [0.9800 2.0500 1.1000];
C(:,:,3) = [1.0000 2.0000 1.0500];

D = [1.000 0.5000 -0.5000]; %被控系统输出方程控制输入增益矩阵

H(:,:,1) = [0 0 0 0 0 0 0 0.10 0.15]; %被控系统输出方程噪声增益矩阵
H(:,:,2) = [0 0 0 0 0 0 0 0.10 0.15]; 
H(:,:,3) = [0 0 0 0 0 0 0 0.10 0.15]; 

E(:,:,1) = [2 0 0;0 1.05 0]; %被控系统量测方程量测矩阵
E(:,:,2) = [2 0 0;0 1.05 0];
E(:,:,3) = [2 0 0;0 1.05 0];

G(:,:,1) = [0 0 0 0.5 0 0 0 0 0;0 0 0 0 0.6 0 0 0 0]; %被控系统量测方程噪声增益矩阵
G(:,:,2) = [0 0 0 0.5 0 0 0 0 0;0 0 0 0 0.6 0 0 0 0]; 
G(:,:,3) = [0 0 0 0.5 0 0 0 0 0;0 0 0 0 0.6 0 0 0 0]; 

Pr(:,:) = [0.95 0.05 0;0.36 0.6 0.04;0.1 0.1 0.8]; % 被控系统模态转移概率

%% 跟踪系统参数
alpha = 0.15;
a = 0.4;
b = 0.5;
c = sqrt(1-a^2-b^2); %由以上参数构造三维正弦波
hat_A = [cos(alpha)+(1-cos(alpha))*a^2   (1-cos(alpha))*a*b-sin(alpha)*c  (1-cos(alpha))*a*c+sin(alpha)*b;
    (1-cos(alpha))*a*b+sin(alpha)*c  cos(alpha)+(1-cos(alpha))*b^2    (1-cos(alpha))*b*c-sin(alpha)*a;
    (1-cos(alpha))*a*c-sin(alpha)*b  (1-cos(alpha))*b*c+sin(alpha)*a  cos(alpha)+(1-cos(alpha))*c^2;]; %被跟踪系统状态方程状态转移矩阵

hat_F = [0.2 0 0 0 0 0 0 0 0;0 0.2 0 0 0 0 0 0 0;0 0 0.2 0 0 0 0 0 0]; %被跟踪系统状态方程噪声增益矩阵

hat_C = [1 1 1]; %被跟踪系统输出方程输出矩阵

hat_E = [1.5 0 0;0 2 0]; %被跟踪系统量测方程量测矩阵

hat_G = [0 0 0 0.25 0 0 0 0 0;0 0 0 0 0.20 0 0 0 0];%被跟踪系统量测方程噪声增益矩阵

%% 构造增广系统参数
[A_row,~] = size(A(:,:,1));  %size函数默认返回矩阵维度,获取被控系统状态的维数
[~,B_col] = size(B(:,:,1));  %获取被控系统控制输入的维数
[~,F_col] = size(F(:,:,1));  %获取被控系统噪声的维数
[C_row,~] = size(C(:,:,1));  %获取被控系统输出的维数
[E_row,~] = size(E(:,:,1));  %获取被控系统量测的维数
[hat_A_row,~] = size(hat_A); %获取被跟踪系统状态的维数
[~,hat_F_col] = size(hat_F); %获取被跟踪系统随机输入的维数
[hat_E_row,~] = size(hat_E); %获取被跟踪系统量测的维数

tilde_A = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %根据系统维数定义增广矩阵
tilde_B = zeros(A_row + hat_A_row,B_col,modes);
tilde_F = zeros(A_row + hat_A_row,F_col + hat_F_col,modes);
tilde_C = zeros(C_row,A_row + hat_A_row,modes);
tilde_H = zeros(C_row,F_col + hat_F_col,modes);
tilde_E = zeros(E_row + hat_E_row,A_row + hat_A_row,modes);
tilde_G = zeros(E_row + hat_E_row,F_col + hat_F_col,modes);

for mode=1:modes
    tilde_A(:,:,mode) = blkdiag(A(:,:,mode),hat_A);  %blkdiag函数根据已有矩阵生成准对角矩阵
    tilde_B(:,:,mode) = [B(:,:,mode)' zeros(B_col,hat_A_row)]';
    tilde_F(:,:,mode) = blkdiag(F(:,:,mode),hat_F);
    tilde_C(:,:,mode) = [C(:,:,mode) -hat_C];
    tilde_H(:,:,mode) = [H(:,:,mode) zeros(C_row,hat_F_col)];
    tilde_E(:,:,mode) = blkdiag(E(:,:,mode),hat_E);
    tilde_G(:,:,mode) = blkdiag(G(:,:,mode),hat_G);
end

%% 给定性能指标参数正定矩阵Q、R、衰减因子gamma
Q(:,:,1) = 10*eye(C_row); %系统跟踪误差权重矩阵
Q(:,:,2) = 10*eye(C_row);
Q(:,:,3) = 10*eye(C_row);

R = 0.5*[1 1 1];  %系统控制输入权重矩阵

gamma = 0.99;  %衰减因子gamma
theta = 22.35; %H无穷控制L2增益
theta_filtering = 3.60; %H无穷滤波L2增益

%% 求解H无穷滤波器滤波器
P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % 给定解的初始值
sigmma_P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %定义按概率加权求和矩阵
L_y(:,:,:) = zeros(A_row + hat_A_row,E_row + hat_E_row,modes);  % 给定控制增益K_u_LQR的初始值
norm_P = [];
episodes = 21; %给定LQR控制器求解迭代次数

for episode = 1:episodes
    for mode = 1:modes
        sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3); %由当前P构造sigmma_P
        P(:,:,mode) = tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_A(:,:,mode)' + tilde_F(:,:,mode)*tilde_F(:,:,mode)' - (tilde_F(:,:,mode)*tilde_G(:,:,mode)' ...
            + tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)')*inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)') ...
            *(tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)')';
        P(:,:,mode) = (P(:,:,mode) + P(:,:,mode)')/2;
    end
    norm_P(episode) = trace(P(:,:,1)'*P(:,:,1) + P(:,:,2)'*P(:,:,2) + P(:,:,3)'*P(:,:,3));
end

for mode = 1:modes
    sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3); %计算按概率加权求和的P
    L_y(:,:,mode) = -(tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)' - tilde_A(:,:,mode)*sigmma_P(:,:,mode)*inv(sigmma_P(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P(:,:,mode)*tilde_E(:,:,mode)') ...
        *inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)' - tilde_E(:,:,mode)*sigmma_P(:,:,mode)*inv(sigmma_P(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P(:,:,mode)*tilde_E(:,:,mode)');
end

figure(6)
plot(0:1:(episodes-1),norm_P)
