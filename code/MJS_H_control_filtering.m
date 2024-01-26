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

F(:,:,1) = [0.10 0 0;0 0.10 0;0 0 0.10]; %被控系统状态方程噪声增益矩阵
F(:,:,2) = [0.10 0 0;0 0.10 0;0 0 0.10];
F(:,:,3) = [0.10 0 0;0 0.10 0;0 0 0.10];

C(:,:,1) = [1.0230 2.1100 0.9500]; %被控系统输出方程输出矩阵
C(:,:,2) = [0.9800 2.0500 1.1000];
C(:,:,3) = [1.0000 2.0000 1.0500];

D = [1.000 0.5000 -0.5000]; %被控系统输出方程控制输入增益矩阵

H(:,:,1) = [0 0.10 0.15]; %被控系统输出方程噪声增益矩阵
H(:,:,2) = [0 0.10 0.15];
H(:,:,3) = [0 0.10 0.15];

E(:,:,1) = [1 0 0;0 1.05 0]; %被控系统量测方程量测矩阵
E(:,:,2) = [1 0 0;0 1.05 0];
E(:,:,3) = [1 0 0;0 1.05 0];

G(:,:,1) = [0.5 0 0;0 0.6 0]; %被控系统量测方程噪声增益矩阵
G(:,:,2) = [0.5 0 0;0 0.6 0];
G(:,:,3) = [0.5 0 0;0 0.6 0];

Pr(:,:) = [0.95 0.05 0;0.36 0.6 0.04;0.1 0.1 0.8]; % 被控系统模态转移概率

%% 跟踪系统参数
alpha = 0.15;
a = 0.4;
b = 0.5;
c = sqrt(1-a^2-b^2); %由以上参数构造三维正弦波
hat_A = [cos(alpha)+(1-cos(alpha))*a^2   (1-cos(alpha))*a*b-sin(alpha)*c  (1-cos(alpha))*a*c+sin(alpha)*b;
    (1-cos(alpha))*a*b+sin(alpha)*c  cos(alpha)+(1-cos(alpha))*b^2    (1-cos(alpha))*b*c-sin(alpha)*a;
    (1-cos(alpha))*a*c-sin(alpha)*b  (1-cos(alpha))*b*c+sin(alpha)*a  cos(alpha)+(1-cos(alpha))*c^2;]; %被跟踪系统状态方程状态转移矩阵

hat_F = [0.2 0 0;0 0.2 0;0 0 0.2]; %被跟踪系统状态方程噪声增益矩阵

hat_C = [1 1 1]; %被跟踪系统输出方程输出矩阵

hat_E = [1 0 0;0 1 0]; %被跟踪系统量测方程量测矩阵

hat_G = [0.25 0 0;0 0.20 0];%被跟踪系统量测方程噪声增益矩阵

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

R = 0.02*[1 1 1];  %系统控制输入权重矩阵

gamma = 0.99;  %衰减因子gamma
theta = 4.55; %H无穷控制L2增益
theta_filtering = 3.60; %H无穷滤波L2增益

%% 控制器求解
%% LQR控制器迭代参数
P_LQR = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % 给定解的初始值
sigmma_P_LQR = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %定义按概率加权求和矩阵
K_u_LQR(:,:,:) = zeros(B_col,A_row+hat_A_row,modes);  % 给定控制增益K_u_LQR的初始值
norm_LQR = 0;
episodes_LQR = 1000; %给定LQR控制器求解迭代次数

for episode_LQR = 1:episodes_LQR
    for mode = 1:modes
        sigmma_P_LQR(:,:,mode) = Pr(mode,1)*P_LQR(:,:,1) + Pr(mode,2)*P_LQR(:,:,2) + Pr(mode,3)*P_LQR(:,:,3); %由当前P构造sigmma_P
        P_LQR(:,:,mode) = tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_A(:,:,mode) - (tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_B(:,:,mode))*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_A(:,:,mode));
    end
    norm_LQR(episode_LQR) = trace(P_LQR(:,:,1)'*P_LQR(:,:,1) + P_LQR(:,:,2)'*P_LQR(:,:,2) + P_LQR(:,:,3)'*P_LQR(:,:,3));
end

for mode = 1:modes
    sigmma_P_LQR(:,:,mode) = Pr(mode,1)*P_LQR(:,:,1) + Pr(mode,2)*P_LQR(:,:,2) + Pr(mode,3)*P_LQR(:,:,3);  %计算按概率加权求和的P
    K_u_LQR(:,:,mode) = -inv(gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_B(:,:,mode) + R(mode))*(D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_A(:,:,mode));
end

figure(1)
plot(0:1:(episodes_LQR-1),norm_LQR)

%% H无穷控制器迭代参数
P_H = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % 给定解的初始值
sigmma_P_H = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %定义按概率加权求和矩阵
sigmma_V_H = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);

% K_u_H(:,:,1) = [0.1910  -0.9150   0.0320   0.0380   0.0360   0.0370];
% K_u_H(:,:,2) = [0.1700  -0.9900   0.0260   0.0690   0.0630   0.0680];
% K_u_H(:,:,3) = [0.0730  -1.4850   0.0050   0.2100   0.1900   0.2100];

K_u_H = K_u_LQR; % 将LQR控制器作为迭代初始值
K_w_H(:,:,:) = zeros(F_col + hat_F_col,A_row+hat_A_row,modes);  % 给定噪声增益K_w的初始值

K_H(:,:,1) = [eye(A_row+hat_A_row);K_u_H(:,:,1);K_w_H(:,:,1)];
K_H(:,:,2) = [eye(A_row+hat_A_row);K_u_H(:,:,2);K_w_H(:,:,2)];
K_H(:,:,3) = [eye(A_row+hat_A_row);K_u_H(:,:,3);K_w_H(:,:,3)];

P_H_episode = zeros(A_row+hat_A_row,A_row+hat_A_row,modes); % 记录迭代过程中每一幕P的值
K_H_episode = zeros(A_row+hat_A_row+B_col+F_col+hat_F_col,A_row+hat_A_row,modes);

norm_P_H_1_episode = [];  %记录解矩阵P的2范数每一幕的变化
norm_P_H_2_episode = [];
norm_P_H_3_episode = [];
norm_K_H_1_episode = [];
norm_K_H_2_episode = [];
norm_K_H_3_episode = [];

Gamma_H = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
Upsilon_H = zeros(A_row+hat_A_row+B_col+F_col+hat_F_col,A_row+hat_A_row+B_col+F_col+hat_F_col,modes);

episodes_H = 11; %给定迭代次数

norm_P_H_1_episode(1) = trace(P_H(:,:,1)'*P_H(:,:,1));
norm_P_H_2_episode(1) = trace(P_H(:,:,2)'*P_H(:,:,2));
norm_P_H_3_episode(1) = trace(P_H(:,:,3)'*P_H(:,:,3));
norm_K_H_1_episode(1) = trace(K_H(:,:,1)'*K_H(:,:,1));
norm_K_H_2_episode(1) = trace(K_H(:,:,2)'*K_H(:,:,2));
norm_K_H_3_episode(1) = trace(K_H(:,:,3)'*K_H(:,:,3));

for episode_H = 1:episodes_H
    episode_H
    for mode = 1:modes %基于K_u以及K_w求解P
        Gamma_H(:,:,mode) = sqrt(gamma)*(tilde_A(:,:,mode) + tilde_B(:,:,mode)*K_u_H(:,:,mode) + tilde_F(:,:,mode)*K_w_H(:,:,mode));  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
        Upsilon_H(:,:,mode) = [tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode); ...
            D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode) D(mode)'*Q(:,:,mode)*D(mode)+R(mode) D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode); ...
            tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) tilde_H(:,:,mode)'*Q(:,:,mode)*D(mode) tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)-theta^(2)*eye(F_col + hat_F_col)];
        K_H_episode(:,:,mode) = [eye(A_row+hat_A_row);K_u_H(:,:,mode);K_w_H(:,:,mode)];
    end
    
    n = 1;
    V_H = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
    while n < 200  %进行耦合Lyapunov方程的求解
        for mode = 1:modes   %对三个模态进行迭代求解
            sigmma_V_H(:,:,mode) = Pr(mode,1)*V_H(:,:,1)+Pr(mode,2)*V_H(:,:,2)+Pr(mode,3)*V_H(:,:,3); %计算按概率加权求和的P
            V_H(:,:,mode) = Gamma_H(:,:,mode)'*sigmma_V_H(:,:,mode)*Gamma_H(:,:,mode) + K_H_episode(:,:,mode)'*Upsilon_H(:,:,mode)*K_H_episode(:,:,mode);  % 迭代求解P
        end
        n = n + 1;
    end
    P_H = V_H;
    
    for mode = 1:modes %基于求解的P更新K_u以及K_w
        sigmma_P_H(:,:,mode) = Pr(mode,1)*P_H(:,:,1)+Pr(mode,2)*P_H(:,:,2)+Pr(mode,3)*P_H(:,:,3);  %计算按概率加权求和的P
        
        K_u_H(:,:,mode) = inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))') ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))');
        K_w_H(:,:,mode) = inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col+hat_F_col)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))) ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))');
        K_H(:,:,mode) = [eye(A_row+hat_A_row);K_u_H(:,:,mode);K_w_H(:,:,mode)];
    end
    norm_P_H_1_episode(episode_H + 1) = log(trace(P_H(:,:,1)'*P_H(:,:,1)));
    norm_P_H_2_episode(episode_H + 1) = log(trace(P_H(:,:,2)'*P_H(:,:,2)));
    norm_P_H_3_episode(episode_H + 1) = log(trace(P_H(:,:,3)'*P_H(:,:,3)));
    norm_K_H_1_episode(episode_H + 1) = log(trace(K_H(:,:,1)'*K_H(:,:,1)));
    norm_K_H_2_episode(episode_H + 1) = log(trace(K_H(:,:,2)'*K_H(:,:,2)));
    norm_K_H_3_episode(episode_H + 1) = log(trace(K_H(:,:,3)'*K_H(:,:,3)));
end

for mode = 1:modes %得到最终结果
    sigmma_P_H(:,:,mode) = Pr(mode,1)*P_H(:,:,1)+Pr(mode,2)*P_H(:,:,2)+Pr(mode,3)*P_H(:,:,3);  %计算按概率加权求和的P
    
    K_u_H(:,:,mode) = inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))') ...
        *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))');
    K_w_H(:,:,mode) = inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col+hat_F_col)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))) ...
        *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))');
    K_H(:,:,mode) = [eye(A_row+hat_A_row);K_u_H(:,:,mode);K_w_H(:,:,mode)];
    
    delta_H(:,:,mode) = P_H(:,:,mode) - tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode)-gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_A(:,:,mode)+[(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)) tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)]*inv([(gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode)+D(mode)'*Q(:,:,mode)*D(mode)) (D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode));(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))' (tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))])*[(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)) tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)]';
end

delta_H = trace(delta_H(:,:,1)'*delta_H(:,:,1)+delta_H(:,:,2)'*delta_H(:,:,2)+delta_H(:,:,3)'*delta_H(:,:,3))

figure(2)
plot(0:1:(episodes_H-1),norm_P_H_1_episode(1:episodes_H),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_P_H_2_episode(1:episodes_H),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_P_H_3_episode(1:episodes_H),'Color','g','LineWidth',1.5)
legend('$log(\left\|P_{1}\right\|_{2})$','$log(\left\|P_{2}\right\|_{2})$','$log(\left\|P_{3}\right\|_{2})$','Interpreter','latex'); %legend在坐标区上添加图例
axis([0 episodes_H-1 0 10.5]) %调整坐标轴范围 axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes_H-1)/10:episodes_H-1]) %设置 x 轴刻度值
ylabel('$log(\left\|P_{i}\right\|_{2})$','interpreter','latex')
yticks([0:1:10])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

figure(3)
plot(0:1:(episodes_H-1),norm_K_H_1_episode(1:episodes_H),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_K_H_2_episode(1:episodes_H),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_K_H_3_episode(1:episodes_H),'Color','g','LineWidth',1.5)
legend('$log(\left\|K_{1}\right\|_{2})$','$log(\left\|K_{2}\right\|_{2})$','$log(\left\|K_{3}\right\|_{2})$','Interpreter','latex');
axis([0 episodes_H-1 0 9]) %调整坐标轴范围axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes_H-1)/10:(episodes_H-1)])
ylabel('$log(\left\|K_{i}\right\|_{2})$','interpreter','latex')
yticks([0:1:9])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5
annotation(figure(3),'ellipse',[0.7015625 0.311320754716981 0.0427083333333333 0.0451215932914054]); % 创建 ellipse
annotation(figure(3),'arrow',[0.7 0.646354166666667],[0.362683438155136 0.419287211740042]); % 创建 arrow

%% TP未知下求解H无穷控制器迭代参数
sigmma_P_TD(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
sigmma_V_TD(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);

%给定初始镇定控制器
K_u_TD = K_u_LQR;
K_w_TD(:,:,:) = zeros(F_col + hat_F_col,A_row+hat_A_row,modes);  % 给定噪声增益K_w的初始值

%给定初始镇定噪声增益
K_TD(:,:,1) = [eye(A_row+hat_A_row);K_u_TD(:,:,1);K_w_TD(:,:,1)];
K_TD(:,:,2) = [eye(A_row+hat_A_row);K_u_TD(:,:,2);K_w_TD(:,:,2)];
K_TD(:,:,3) = [eye(A_row+hat_A_row);K_u_TD(:,:,3);K_w_TD(:,:,3)];

Gamma_TD = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
Upsilon_TD = zeros(A_row+hat_A_row+B_col+F_col+hat_F_col,A_row+hat_A_row+B_col+F_col+hat_F_col,modes);

delta_sigmma_P_TD = [0 0 0]; %用于储存sigmma_P距离最优的差值
delta_K_TD = [0 0 0]; %用于储存K_u距离最优的差值

episodes_TD = 16; %定义迭代的幕数
steps_TD = 100; %定义每一幕的步数
mu = 0.95; %定义回报权重

for mode = 1:modes
    delta_sigmma_P_TD(1,mode,1) = log(1+(trace((sigmma_P_H(:,:,mode)-sigmma_P_TD(:,:,mode))'*(sigmma_P_H(:,:,mode)-sigmma_P_TD(:,:,mode))))/(trace(sigmma_P_H(:,:,mode)'*sigmma_P_H(:,:,mode))));
    delta_K_TD(1,mode,1) = log(1+(trace((K_H(:,:,mode)-K_TD(:,:,mode))'*(K_H(:,:,mode)-K_TD(:,:,mode))))/(trace(K_H(:,:,mode)'*K_H(:,:,mode))));
end

for episode_TD = 1:episodes_TD
    episode_TD
    
    lambda = 1/(episode_TD); %更新步长
    sigmma_V_TD = sigmma_P_TD;
    
    for mode = 1:modes
        Gamma_TD(:,:,mode) = sqrt(gamma)*(tilde_A(:,:,mode)+tilde_B(:,:,mode)*K_u_TD(:,:,mode)+tilde_F(:,:,mode)*K_w_TD(:,:,mode));  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
        Upsilon_TD(:,:,mode) = [tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode); ...
            D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode) D(mode)'*Q(:,:,mode)*D(mode)+R(mode) D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode); ...
            tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) tilde_H(:,:,mode)'*Q(:,:,mode)*D(mode) tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)-theta^(2)*eye(F_col + hat_F_col)];
        K_TD(:,:,mode) = [eye(A_row+hat_A_row);K_u_TD(:,:,mode);K_w_TD(:,:,mode)];
    end
    
    for mode = 1:modes
        mode_now = mode; %保留每一幕的第一个模态
        TD_sum(:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row);
        for step_TD=1:steps_TD
            mode_old = mode_now; %储存当前模态
            mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %基于转移概率更新模态
            
            if step_TD == 1
                Lambda(:,:,mode) = eye(A_row+hat_A_row); %定义Lambda矩阵
            else
                Lambda(:,:,mode) = Gamma_TD(:,:,mode_old)*Lambda(:,:,mode); %更新Lambda矩阵
            end
            
            TD = Lambda(:,:,mode)'*(Gamma_TD(:,:,mode_now)'*sigmma_V_TD(:,:,mode_now)*Gamma_TD(:,:,mode_now)+K_TD(:,:,mode_now)'*Upsilon_TD(:,:,mode_now)*K_TD(:,:,mode_now)-sigmma_V_TD(:,:,mode_old))*Lambda(:,:,mode);
            
            TD_sum = TD_sum + mu^(step_TD-1)*TD; %进行第step步的sigmma_P更新
        end
        sigmma_P_TD(:,:,mode) = sigmma_V_TD(:,:,mode) + lambda*TD_sum; %进行第step步的sigmma_P更新
        K_u_TD(:,:,mode) = inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode)+R(mode)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))') ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode))');
        K_w_TD(:,:,mode) = inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col+hat_F_col)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))) ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))');
    end
    for mode = 1:modes
        delta_sigmma_P_TD(1,mode,episode_TD+1) = log(1 + (trace((sigmma_P_H(:,:,mode)-sigmma_P_TD(:,:,mode))'*(sigmma_P_H(:,:,mode)-sigmma_P_TD(:,:,mode))))/(trace(sigmma_P_H(:,:,mode)'*sigmma_P_H(:,:,mode))));
        delta_K_TD(1,mode,episode_TD+1) = log(1 + (trace((K_H(:,:,mode)-K_TD(:,:,mode))'*(K_H(:,:,mode)-K_TD(:,:,mode))))/(trace(K_H(:,:,mode)'*K_H(:,:,mode))));
    end
end

for i = 1:modes
    P_TD(:,:,i) = tilde_C(:,:,i)'*Q(:,:,i)*tilde_C(:,:,i) + gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_A(:,:,i) - [(tilde_C(:,:,i)'*Q(:,:,i)*D(i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_B(:,:,i)) tilde_C(:,:,i)'*Q(:,:,i)*tilde_H(:,:,i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i)]*inv([(gamma*tilde_B(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_B(:,:,i)+R(i)+D(i)'*Q(:,:,i)*D(i)) (D(i)'*Q(:,:,i)*tilde_H(:,:,i) + gamma*tilde_B(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i));(D(i)'*Q(:,:,i)*tilde_H(:,:,i) + gamma*tilde_B(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i))' (tilde_H(:,:,i)'*Q(:,:,i)*tilde_H(:,:,i)+gamma*tilde_F(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i)-theta^(2)*eye(F_col + hat_F_col))])*[(tilde_C(:,:,i)'*Q(:,:,i)*D(i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_B(:,:,i)) tilde_C(:,:,i)'*Q(:,:,i)*tilde_H(:,:,i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i)]';
end
for i = 1:modes
    delta_TD(:,:,i) = sigmma_P_TD(:,:,i) - Pr(i,1)*P_TD(:,:,1) - Pr(i,2)*P_TD(:,:,2) - Pr(i,3)*P_TD(:,:,3);
end
delta_TD_norm = trace(delta_TD(:,:,1)'*delta_TD(:,:,1)+delta_TD(:,:,2)'*delta_TD(:,:,2)+delta_TD(:,:,3)'*delta_TD(:,:,3))

delta_sigmma_P1 = delta_sigmma_P_TD(1,1,:);
delta_sigmma_P1 = delta_sigmma_P1(:);  %通过（：）将多维数组变为一维数组
delta_sigmma_P2 = delta_sigmma_P_TD(1,2,:);
delta_sigmma_P2 = delta_sigmma_P2(:);
delta_sigmma_P3 = delta_sigmma_P_TD(1,3,:);
delta_sigmma_P3 = delta_sigmma_P3(:);

delta_K_TD_1 = delta_K_TD(1,1,:);
delta_K_TD_1 = delta_K_TD_1(:);  %通过（：）将多维数组变为一维数组
delta_K_TD_2 = delta_K_TD(1,2,:);
delta_K_TD_2 = delta_K_TD_2(:);
delta_K_TD_3 = delta_K_TD(1,3,:);
delta_K_TD_3 = delta_K_TD_3(:);

figure(4)
plot(0:1:(episodes_TD-1),delta_sigmma_P1(1:episodes_TD),'--','Color','b',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_sigmma_P2(1:episodes_TD),'-.','Color','r',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_sigmma_P3(1:episodes_TD),'Color','g',"LineWidth",1.5)
legend('$\Delta^{(l)}_{1}$','$\Delta^{(l)}_{2}$','$\Delta^{(l)}_{3}$','Interpreter','latex','Position', ...
    [0.691145833333333 0.212788259958072 0.188194444444442 0.205170875242091]); %legend在坐标区上添加图例
xlabel('$Iteration$','interpreter','latex')
ylabel('error of value function','interpreter','latex')
axis([0 episodes_TD-1 -0.05 0.8])
xticks([0:1:episodes_TD-1])
% yticks([0:0.1:0.8])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

figure(5)
plot(0:1:(episodes_TD-1),delta_K_TD_1(1:episodes_TD),'--','Color','b',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_K_TD_2(1:episodes_TD),'-.','Color','r',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_K_TD_3(1:episodes_TD),'Color','g',"LineWidth",1.5)
legend('$\delta^{(l)}_{1}$','$\delta^{(l)}_{2}$','$\delta^{(l)}_{3}$','Interpreter','latex','Position', ...
    [0.691145833333333 0.212788259958072 0.188194444444442 0.205170875242091]); %legend在坐标区上添加图例
xlabel('$Iteration$','interpreter','latex')
ylabel('error of controller','interpreter','latex')
axis([0 episodes_TD-1 -0.0005 0.0040])
xticks([0:1:episodes_TD-1])
% yticks([0:0.00005:0.0004])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

%% 滤波器求解
%% 求解H无穷滤波器初始镇定滤波器
P_filtering_initial = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % 给定解的初始值
sigmma_P_filtering_initial = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %定义按概率加权求和矩阵
L_y_H_filtering_initial(:,:,:) = zeros(A_row + hat_A_row,E_row + hat_E_row,modes);  % 给定控制增益K_u_LQR的初始值
norm_P_filtering_initial = [];
episodes_P_filtering_initial = 51; %给定LQR控制器求解迭代次数

for episode_P_filtering_initial = 1:episodes_P_filtering_initial
    for mode = 1:modes
        sigmma_P_filtering_initial(:,:,mode) = Pr(mode,1)*P_filtering_initial(:,:,1) + Pr(mode,2)*P_filtering_initial(:,:,2) + Pr(mode,3)*P_filtering_initial(:,:,3); %由当前P构造sigmma_P
        P_filtering_initial(:,:,mode) = tilde_A(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*tilde_A(:,:,mode)' + tilde_F(:,:,mode)*tilde_F(:,:,mode)' - (tilde_F(:,:,mode)*tilde_G(:,:,mode)' ...
            + tilde_A(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*tilde_E(:,:,mode)')*inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*tilde_E(:,:,mode)') ...
            *(tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*tilde_E(:,:,mode)')';
        P_filtering_initial(:,:,mode) = (P_filtering_initial(:,:,mode) + P_filtering_initial(:,:,mode)')/2;
    end
    norm_P_filtering_initial(episode_P_filtering_initial) = trace(P_filtering_initial(:,:,1)'*P_filtering_initial(:,:,1) + P_filtering_initial(:,:,2)'*P_filtering_initial(:,:,2) + P_filtering_initial(:,:,3)'*P_filtering_initial(:,:,3));
end

for mode = 1:modes
    sigmma_P_filtering_initial(:,:,mode) = Pr(mode,1)*P_filtering_initial(:,:,1) + Pr(mode,2)*P_filtering_initial(:,:,2) + Pr(mode,3)*P_filtering_initial(:,:,3); %计算按概率加权求和的P
    L_y_H_filtering_initial(:,:,mode) = -(tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*tilde_E(:,:,mode)' - tilde_A(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*inv(sigmma_P_filtering_initial(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_filtering_initial(:,:,mode)*tilde_E(:,:,mode)') ...
        *inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*tilde_E(:,:,mode)' - tilde_E(:,:,mode)*sigmma_P_filtering_initial(:,:,mode)*inv(sigmma_P_filtering_initial(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_filtering_initial(:,:,mode)*tilde_E(:,:,mode)');
end

figure(6)
plot(0:1:(episodes_P_filtering_initial-1),norm_P_filtering_initial)

%% H无穷滤波器迭代参数
P_H_filtering = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % 给定解的初始值
sigmma_P_H_filtering = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %定义按概率加权求和矩阵
sigmma_V_H_filtering = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);

L_y_H_filtering = L_y_H_filtering_initial; % 将LQR控制器作为迭代初始值
L_w_H_filtering(:,:,:) = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);  % 给定噪声增益K_w的初始值

L_H_filtering(:,:,1) = [eye(A_row + hat_A_row) L_y_H_filtering(:,:,1) L_w_H_filtering(:,:,1)];
L_H_filtering(:,:,2) = [eye(A_row + hat_A_row) L_y_H_filtering(:,:,2) L_w_H_filtering(:,:,2)];
L_H_filtering(:,:,3) = [eye(A_row + hat_A_row) L_y_H_filtering(:,:,3) L_w_H_filtering(:,:,3)];

norm_P_H_filtering_1_episode = [];  %记录解矩阵P的2范数每一幕的变化
norm_P_H_filtering_2_episode = [];
norm_P_H_filtering_3_episode = [];

norm_L_H_filtering_1_episode = [];
norm_L_H_filtering_2_episode = [];
norm_L_H_filtering_3_episode = [];

Gamma_H_filtering = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
Upsilon_H_filtering = zeros(A_row + hat_A_row + A_row + hat_A_row + E_row + hat_E_row,A_row + hat_A_row + A_row + hat_A_row + E_row + hat_E_row,modes);

episodes_H_filtering = 101; %给定迭代次数

norm_P_H_filtering_1_episode(1) = log( 1 + trace(P_H_filtering(:,:,1)'*P_H_filtering(:,:,1)));
norm_P_H_filtering_2_episode(1) = log( 1 + trace(P_H_filtering(:,:,2)'*P_H_filtering(:,:,2)));
norm_P_H_filtering_3_episode(1) = log( 1 + trace(P_H_filtering(:,:,3)'*P_H_filtering(:,:,3)));
norm_L_H_filtering_1_episode(1) = log( 1 + trace(L_H_filtering(:,:,1)'*L_H_filtering(:,:,1)));
norm_L_H_filtering_2_episode(1) = log( 1 + trace(L_H_filtering(:,:,2)'*L_H_filtering(:,:,2)));
norm_L_H_filtering_3_episode(1) = log( 1 + trace(L_H_filtering(:,:,3)'*L_H_filtering(:,:,3)));

for episode_H_filtering = 1:episodes_H_filtering
    episode_H_filtering
    for mode = 1:modes %基于K_u以及K_w求解P
        Gamma_H_filtering(:,:,mode) = tilde_A(:,:,mode) + L_y_H_filtering(:,:,mode)*tilde_E(:,:,mode) + L_w_H_filtering(:,:,mode);  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
        Upsilon_H_filtering(:,:,mode) = [tilde_F(:,:,mode)*tilde_F(:,:,mode)' tilde_F(:,:,mode)*tilde_G(:,:,mode)' zeros(A_row + hat_A_row,A_row + hat_A_row); ...
            tilde_G(:,:,mode)*tilde_F(:,:,mode)' tilde_G(:,:,mode)*tilde_G(:,:,mode)' zeros(E_row + hat_E_row,A_row + hat_A_row); ...
            zeros(A_row + hat_A_row,A_row + hat_A_row) zeros(A_row + hat_A_row,E_row + hat_E_row) -theta_filtering^(2)*eye(A_row + hat_A_row)];
        L_H_filtering(:,:,mode) = [eye(A_row+hat_A_row) L_y_H_filtering(:,:,mode) L_w_H_filtering(:,:,mode)];
    end
    
    n = 1;
    V_H_filtering = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);
    while n < 200  %进行耦合Lyapunov方程的求解
        for mode = 1:modes   %对三个模态进行迭代求解
            sigmma_V_H_filtering(:,:,mode) = Pr(mode,1)*V_H_filtering(:,:,1) + Pr(mode,2)*V_H_filtering(:,:,2) + Pr(mode,3)*V_H_filtering(:,:,3); %计算按概率加权求和的P
            V_H_filtering(:,:,mode) = Gamma_H_filtering(:,:,mode)*sigmma_V_H_filtering(:,:,mode)*Gamma_H_filtering(:,:,mode)' + L_H_filtering(:,:,mode)*Upsilon_H_filtering(:,:,mode)*L_H_filtering(:,:,mode)';  % 迭代求解P
            V_H_filtering(:,:,mode) = (V_H_filtering(:,:,mode) + V_H_filtering(:,:,mode)')/2;
        end
        n = n + 1;
    end
    P_H_filtering = V_H_filtering;
    
    for mode = 1:modes %基于求解的P更新K_u以及K_w
        sigmma_P_H_filtering(:,:,mode) = Pr(mode,1)*P_H_filtering(:,:,1) + Pr(mode,2)*P_H_filtering(:,:,2) + Pr(mode,3)*P_H_filtering(:,:,3);  %计算按概率加权求和的P
        
        L_y_H_filtering(:,:,mode) = (tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*inv(sigmma_P_H_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_F(:,:,mode)*tilde_G(:,:,mode)') ...
            *inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*inv(sigmma_P_H_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)');
        
        L_w_H_filtering(:,:,mode) = ((tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)')*inv(tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' + tilde_G(:,:,mode)*tilde_G(:,:,mode)')*tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode) - tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)) ...
            *inv(sigmma_P_H_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row) - sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)'*inv(tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' + tilde_G(:,:,mode)*tilde_G(:,:,mode)')*tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode));
        
        L_H_filtering(:,:,mode) = [eye(A_row+hat_A_row) L_y_H_filtering(:,:,mode) L_w_H_filtering(:,:,mode)];
    end
    norm_P_H_filtering_1_episode(episode_H_filtering + 1) = log(1 + trace(P_H_filtering(:,:,1)'*P_H_filtering(:,:,1)));
    norm_P_H_filtering_2_episode(episode_H_filtering + 1) = log(1 + trace(P_H_filtering(:,:,2)'*P_H_filtering(:,:,2)));
    norm_P_H_filtering_3_episode(episode_H_filtering + 1) = log(1 + trace(P_H_filtering(:,:,3)'*P_H_filtering(:,:,3)));
    norm_L_H_filtering_1_episode(episode_H_filtering + 1) = log(1 + trace(L_H_filtering(:,:,1)'*L_H_filtering(:,:,1)));
    norm_L_H_filtering_2_episode(episode_H_filtering + 1) = log(1 + trace(L_H_filtering(:,:,2)'*L_H_filtering(:,:,2)));
    norm_L_H_filtering_3_episode(episode_H_filtering + 1) = log(1 + trace(L_H_filtering(:,:,3)'*L_H_filtering(:,:,3)));
end

for mode = 1:modes %得到最终结果
    sigmma_P_H_filtering(:,:,mode) = Pr(mode,1)*P_H_filtering(:,:,1) + Pr(mode,2)*P_H_filtering(:,:,2) + Pr(mode,3)*P_H_filtering(:,:,3);  %计算按概率加权求和的P
    L_y_H_filtering(:,:,mode) = (tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*inv(sigmma_P_H_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_F(:,:,mode)*tilde_G(:,:,mode)') ...
        *inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*inv(sigmma_P_H_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)');
    L_w_H_filtering(:,:,mode) = ((tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)')*inv(tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' + tilde_G(:,:,mode)*tilde_G(:,:,mode)')*tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode) - tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)) ...
        *inv(sigmma_P_H_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row) - sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)'*inv(tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' + tilde_G(:,:,mode)*tilde_G(:,:,mode)')*tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode));
    L_H_filtering(:,:,mode) = [eye(A_row+hat_A_row) L_y_H_filtering(:,:,mode) L_w_H_filtering(:,:,mode)];
    delta_H_filtering(:,:,mode) = P_H_filtering(:,:,mode) - (tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_A(:,:,mode)' + tilde_F(:,:,mode)*tilde_F(:,:,mode)' - ([tilde_F(:,:,mode)*tilde_G(:,:,mode)' ...
        + tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)])*inv([tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' tilde_E(:,:,mode)*sigmma_P_H_filtering(:,:,mode); ...
        sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' sigmma_P_H_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row)]) ...
        *([tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)*tilde_E(:,:,mode)' tilde_A(:,:,mode)*sigmma_P_H_filtering(:,:,mode)])');
end

delta_H_filtering = trace(delta_H_filtering(:,:,1)'*delta_H_filtering(:,:,1) + delta_H_filtering(:,:,2)'*delta_H_filtering(:,:,2) + delta_H_filtering(:,:,3)'*delta_H_filtering(:,:,3))

figure(7)
plot(0:1:(episodes_H_filtering-1),norm_P_H_filtering_1_episode(1:episodes_H_filtering),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes_H_filtering-1),norm_P_H_filtering_2_episode(1:episodes_H_filtering),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes_H_filtering-1),norm_P_H_filtering_3_episode(1:episodes_H_filtering),'Color','g','LineWidth',1.5)
legend('$log(\left\|P_{1}\right\|_{2})$','$log(\left\|P_{2}\right\|_{2})$','$log(\left\|P_{3}\right\|_{2})$','Interpreter','latex'); %legend在坐标区上添加图例
axis([0 episodes_H_filtering-1 0 10.5]) %调整坐标轴范围 axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes_H_filtering-1)/10:episodes_H_filtering-1]) %设置 x 轴刻度值
ylabel('$log(\left\|P_{i}\right\|_{2})$','interpreter','latex')
yticks([0:1:10])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

figure(8)
plot(0:1:(episodes_H_filtering-1),norm_L_H_filtering_1_episode(1:episodes_H_filtering),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes_H_filtering-1),norm_L_H_filtering_2_episode(1:episodes_H_filtering),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes_H_filtering-1),norm_L_H_filtering_3_episode(1:episodes_H_filtering),'Color','g','LineWidth',1.5)
legend('$log(\left\|K_{1}\right\|_{2})$','$log(\left\|K_{2}\right\|_{2})$','$log(\left\|K_{3}\right\|_{2})$','Interpreter','latex');
axis([0 episodes_H_filtering-1 0 9]) %调整坐标轴范围axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes_H_filtering-1)/10:(episodes_H_filtering-1)])
ylabel('$log(\left\|K_{i}\right\|_{2})$','interpreter','latex')
yticks([0:1:9])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5
annotation(figure(8),'ellipse',[0.7015625 0.311320754716981 0.0427083333333333 0.0451215932914054]); % 创建 ellipse
annotation(figure(8),'arrow',[0.7 0.646354166666667],[0.362683438155136 0.419287211740042]); % 创建 arrow

%% TP未知求解H无穷滤波器
sigmma_P_TD_filtering(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
sigmma_V_TD_filtering(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);

L_y_TD_filtering = L_y_H_filtering_initial; %给定初始镇定滤波器
L_w_TD_filtering(:,:,:) = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);  % 给定噪声增益K_w的初始值

L_TD_filtering(:,:,1) = [eye(A_row+hat_A_row) L_y_TD_filtering(:,:,1) L_w_TD_filtering(:,:,1)];
L_TD_filtering(:,:,2) = [eye(A_row+hat_A_row) L_y_TD_filtering(:,:,2) L_w_TD_filtering(:,:,2)];
L_TD_filtering(:,:,3) = [eye(A_row+hat_A_row) L_y_TD_filtering(:,:,3) L_w_TD_filtering(:,:,3)];

Gamma_TD_filtering = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);
Upsilon_TD_filtering = zeros(A_row + hat_A_row + A_row + hat_A_row + E_row + hat_E_row,A_row + hat_A_row + A_row + hat_A_row + E_row + hat_E_row,modes);

delta_sigmma_P_TD_filtering = zeros(1,modes); %用于储存sigmma_P距离最优的差值
delta_L_TD_filtering = zeros(1,modes); %用于储存K_u距离最优的差值

episodes_TD_filtering = 56; %定义迭代的幕数
steps_TD_filtering = 300; %定义每一幕的步数
mu_filtering = 0.05; %定义回报权重

for mode = 1:modes
    delta_sigmma_P_TD_filtering(1,mode,1) = log(1 + (trace((sigmma_P_H_filtering(:,:,mode) - sigmma_P_TD_filtering(:,:,mode))'*(sigmma_P_H_filtering(:,:,mode) - sigmma_P_TD_filtering(:,:,mode))))/(trace(sigmma_P_H_filtering(:,:,mode)'*sigmma_P_H_filtering(:,:,mode))));
    delta_L_TD_filtering(1,mode,1) = log(1 + (trace((L_H_filtering(:,:,mode) - L_TD_filtering(:,:,mode))'*(L_H_filtering(:,:,mode) - L_TD_filtering(:,:,mode))))/(trace(L_H_filtering(:,:,mode)'*L_H_filtering(:,:,mode))));
    Gamma_TD_filtering(:,:,mode) = tilde_A(:,:,mode) + L_y_TD_filtering(:,:,mode)*tilde_E(:,:,mode) + L_w_TD_filtering(:,:,mode);  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
    Upsilon_TD_filtering(:,:,mode) = [tilde_F(:,:,mode)*tilde_F(:,:,mode)' tilde_F(:,:,mode)*tilde_G(:,:,mode)' zeros(A_row + hat_A_row,A_row + hat_A_row); ...
        tilde_G(:,:,mode)*tilde_F(:,:,mode)' tilde_G(:,:,mode)*tilde_G(:,:,mode)' zeros(E_row + hat_E_row,A_row + hat_A_row); ...
        zeros(A_row + hat_A_row,A_row + hat_A_row) zeros(A_row + hat_A_row,E_row + hat_E_row) -theta_filtering^(2)*eye(A_row + hat_A_row)];
    L_TD_filtering(:,:,mode) = [eye(A_row + hat_A_row) L_y_TD_filtering(:,:,mode) L_w_TD_filtering(:,:,mode)];
end

for episode_TD_filtering = 1:episodes_TD_filtering
    episode_TD_filtering
    
    lambda_filtering = 1/(episode_TD_filtering); %更新步长
    sigmma_V_TD_filtering = sigmma_P_TD_filtering;
    
    for mode = 1:modes
        mode_now = mode; %保留每一幕的第一个模态
        TD_filtering_sum(:,:) = zeros(A_row + hat_A_row,A_row + hat_A_row);
        for step_TD_filtering = 1:steps_TD_filtering
            mode_old = mode_now; %储存当前模态
            mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %基于转移概率更新模态
            
            if step_TD_filtering == 1
                Lambda_filtering(:,:,mode) = eye(A_row + hat_A_row); %定义Lambda矩阵
            else
                Lambda_filtering(:,:,mode) = Gamma_TD_filtering(:,:,mode_old)'*Lambda_filtering(:,:,mode); %更新Lambda矩阵
            end
            
            TD_filtering = Lambda_filtering(:,:,mode)'*(Gamma_TD_filtering(:,:,mode_now)*sigmma_V_TD_filtering(:,:,mode_now)*Gamma_TD_filtering(:,:,mode_now)' + L_TD_filtering(:,:,mode_now)*Upsilon_TD_filtering(:,:,mode_now)*L_TD_filtering(:,:,mode_now)' - sigmma_V_TD_filtering(:,:,mode_old))*Lambda_filtering(:,:,mode);
            TD_filtering = (TD_filtering + TD_filtering' + 0.0001*eye(A_row+hat_A_row))/2;
            TD_filtering_sum = TD_filtering_sum + mu_filtering^(step_TD_filtering - 1)*TD_filtering; %进行第step步的sigmma_P更新
            TD_filtering_sum = (TD_filtering_sum + TD_filtering_sum')/2;
        end
        sigmma_P_TD_filtering(:,:,mode) = sigmma_V_TD_filtering(:,:,mode) + lambda_filtering*TD_filtering_sum; %进行第step步的sigmma_P更新
        L_y_TD_filtering(:,:,mode) = (tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*inv(sigmma_P_TD_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_F(:,:,mode)*tilde_G(:,:,mode)') ...
            *inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' - tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*inv(sigmma_P_TD_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)');
        L_w_TD_filtering(:,:,mode) = ((tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)')*inv(tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' + tilde_G(:,:,mode)*tilde_G(:,:,mode)')*tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode) - tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)) ...
            *inv(sigmma_P_TD_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row) - sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)'*inv(tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' + tilde_G(:,:,mode)*tilde_G(:,:,mode)')*tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode));
        Gamma_TD_filtering(:,:,mode) = tilde_A(:,:,mode) + L_y_TD_filtering(:,:,mode)*tilde_E(:,:,mode) + L_w_TD_filtering(:,:,mode);  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
        L_TD_filtering(:,:,mode) = [eye(A_row + hat_A_row) L_y_TD_filtering(:,:,mode) L_w_TD_filtering(:,:,mode)];
    end
    
    for mode = 1:modes
        delta_sigmma_P_TD_filtering(1,mode,episode_TD_filtering + 1) = log(1 + (trace((sigmma_P_H_filtering(:,:,mode) - sigmma_P_TD_filtering(:,:,mode))'*(sigmma_P_H_filtering(:,:,mode) - sigmma_P_TD_filtering(:,:,mode))))/(trace(sigmma_P_H_filtering(:,:,mode)'*sigmma_P_H_filtering(:,:,mode))));
        delta_L_TD_filtering(1,mode,episode_TD_filtering + 1) = log(1 + (trace((L_H_filtering(:,:,mode) - L_TD_filtering(:,:,mode))'*(L_H_filtering(:,:,mode) - L_TD_filtering(:,:,mode))))/(trace(L_H_filtering(:,:,mode)'*L_H_filtering(:,:,mode))));
    end
end

for i = 1:modes
    P_TD_filtering(:,:,i) = tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_A(:,:,mode)' + tilde_F(:,:,mode)*tilde_F(:,:,mode)' - ([tilde_F(:,:,mode)*tilde_G(:,:,mode)' ...
        + tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)])*inv([tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' tilde_E(:,:,mode)*sigmma_P_TD_filtering(:,:,mode); ...
        sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' sigmma_P_TD_filtering(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row)]) ...
        *([tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)*tilde_E(:,:,mode)' tilde_A(:,:,mode)*sigmma_P_TD_filtering(:,:,mode)])';
end
for i = 1:modes
    delta_TD_filtering(:,:,i) = sigmma_P_TD_filtering(:,:,i) - Pr(i,1)*P_TD_filtering(:,:,1) - Pr(i,2)*P_TD_filtering(:,:,2) - Pr(i,3)*P_TD_filtering(:,:,3);
end
delta_TD_filtering_norm = trace(delta_TD_filtering(:,:,1)'*delta_TD_filtering(:,:,1) + delta_TD_filtering(:,:,2)'*delta_TD_filtering(:,:,2) + delta_TD_filtering(:,:,3)'*delta_TD_filtering(:,:,3))

delta_sigmma_P_TD_filtering_1 = delta_sigmma_P_TD_filtering(1,1,:);
delta_sigmma_P_TD_filtering_1 = delta_sigmma_P_TD_filtering_1(:);  %通过（：）将多维数组变为一维数组
delta_sigmma_P_TD_filtering_2 = delta_sigmma_P_TD_filtering(1,2,:);
delta_sigmma_P_TD_filtering_2 = delta_sigmma_P_TD_filtering_2(:);
delta_sigmma_P_TD_filtering_3 = delta_sigmma_P_TD_filtering(1,3,:);
delta_sigmma_P_TD_filtering_3 = delta_sigmma_P_TD_filtering_3(:);

delta_L_TD_filtering_1 = delta_L_TD_filtering(1,1,:);
delta_L_TD_filtering_1 = delta_L_TD_filtering_1(:);  %通过（：）将多维数组变为一维数组
delta_L_TD_filtering_2 = delta_L_TD_filtering(1,2,:);
delta_L_TD_filtering_2 = delta_L_TD_filtering_2(:);
delta_L_TD_filtering_3 = delta_L_TD_filtering(1,3,:);
delta_L_TD_filtering_3 = delta_L_TD_filtering_3(:);

figure(9)
plot(0:1:(episodes_TD_filtering-1),delta_sigmma_P_TD_filtering_1(1:episodes_TD_filtering),'--','Color','b',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD_filtering-1),delta_sigmma_P_TD_filtering_2(1:episodes_TD_filtering),'-.','Color','r',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD_filtering-1),delta_sigmma_P_TD_filtering_3(1:episodes_TD_filtering),'Color','g',"LineWidth",1.5)
legend('$\Delta^{(l)}_{1}$','$\Delta^{(l)}_{2}$','$\Delta^{(l)}_{3}$','Interpreter','latex','Position', ...
    [0.691145833333333 0.212788259958072 0.188194444444442 0.205170875242091]); %legend在坐标区上添加图例
xlabel('$Iteration$','interpreter','latex')
ylabel('error of value function','interpreter','latex')
% axis([0 episodes_TD_filtering-1 -0.05 0.8])
xticks([0:1:episodes_TD_filtering-1])
% yticks([0:0.1:0.8])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

figure(10)
plot(0:1:(episodes_TD_filtering-1),delta_L_TD_filtering_1(1:episodes_TD_filtering),'--','Color','b',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD_filtering-1),delta_L_TD_filtering_2(1:episodes_TD_filtering),'-.','Color','r',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD_filtering-1),delta_L_TD_filtering_3(1:episodes_TD_filtering),'Color','g',"LineWidth",1.5)
legend('$\delta^{(l)}_{1}$','$\delta^{(l)}_{2}$','$\delta^{(l)}_{3}$','Interpreter','latex','Position', ...
    [0.691145833333333 0.212788259958072 0.188194444444442 0.205170875242091]); %legend在坐标区上添加图例
xlabel('$Iteration$','interpreter','latex')
ylabel('error of controller','interpreter','latex')
% axis([0 episodes_TD_filtering-1 -0.0005 0.0040])
xticks([0:1:episodes_TD_filtering-1])
% yticks([0:0.00005:0.0004])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

%% 进行控制器测试
%% 系统模拟仿真参数
episodes_test = 1; %给定迭代次数
steps_test = 300;
hat_x = [20 -20 0]'; % 被跟踪系统参数
hat_w = zeros(hat_F_col,1); % 被跟踪系统参数
hat_y = zeros(hat_E_row,1);
hat_z = zeros(C_row,1);

%% 系统模拟仿真参数处于LQR控制器
x_LQR = [0 0 0]'; % 被控系统参数
u_LQR = zeros(B_col,1);
w_LQR = zeros(F_col,1);
z_LQR = zeros(C_row,1);
y_LQR = zeros(E_row,1);
tilde_x_LQR = [x_LQR' hat_x']'; % 增广系统参数
estimate_tilde_x_LQR = 0*tilde_x_LQR; % 增广系统参数
tilde_w_LQR = zeros(F_col + hat_F_col,1);
tilde_z_LQR = zeros(C_row,1);
tilde_y_LQR = [y_LQR' hat_y']';
estimate_tilde_y_LQR = 0*tilde_y_LQR;
delta_tilde_x_LQR = 0*tilde_x_LQR;

%% 系统模拟仿真参数处于H无穷控制器
x_H = [0 0 0]'; % 被控系统参数
u_H = zeros(B_col,1);
w_H = zeros(F_col,1);
z_H = zeros(C_row,1);
y_H = zeros(E_row,1);
tilde_x_H = [x_H' hat_x']'; % 增广系统参数
estimate_tilde_x_H = 0*tilde_x_H;
tilde_w_H = zeros(F_col + hat_F_col,1);
tilde_z_H = zeros(C_row,1);
tilde_y_H = [y_H' hat_y']';
estimate_tilde_y_H = 0*tilde_y_H;
delta_tilde_x_H = 0*tilde_x_H;

%% 系统模拟仿真参数处于TD控制器
x_TD = [0 0 0]'; % 被控系统参数
u_TD = zeros(B_col,1);
w_TD = zeros(F_col,1);
z_TD = zeros(C_row,1);
y_TD = zeros(E_row,1);
tilde_x_TD = [x_TD' hat_x']'; % 增广系统参数
estimate_tilde_x_TD = 0*tilde_x_TD;
tilde_w_TD = zeros(F_col + hat_F_col,1);
tilde_z_TD = zeros(C_row,1);
tilde_y_TD = [y_TD' hat_y']';
estimate_tilde_y_TD = 0*tilde_y_TD;
delta_tilde_x_TD = 0*tilde_x_TD;

%% 仿真跟踪系统
for episode_test = 1:episodes_test
    mode_now = randsrc(1,1,[1 2 3;1/3 1/3 1/3]);
    for step_test = 1:steps_test
        mode_old = mode_now; %储存当前模态
        mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %基于转移概率更新模态

        delta_tilde_x_LQR(:,step_test) = tilde_x_LQR(:,step_test) - estimate_tilde_x_LQR(:,step_test);
        delta_tilde_x_H(:,step_test) = tilde_x_H(:,step_test) - estimate_tilde_x_H(:,step_test);
        delta_tilde_x_TD(:,step_test) = tilde_x_TD(:,step_test) - estimate_tilde_x_TD(:,step_test);

        tilde_w_LQR(:,step_test) = normrnd(0,0.005,[F_col + hat_F_col 1]);
        tilde_w_H(:,step_test) = normrnd(0,0.005,[F_col + hat_F_col 1]);
        tilde_w_TD(:,step_test) = normrnd(0,0.005,[F_col + hat_F_col 1]);
        
        %         tilde_w_LQR(:,step_test) = K_w_H(:,:, mode_now)*tilde_x_LQR(:,step_test);
        %         tilde_w_H(:,step_test) = K_w_H(:,:, mode_now)*tilde_x_H(:,step_test);
        %         tilde_w_TD(:,step_test) = K_w_H(:,:, mode_now)*tilde_x_TD(:,step_test);
        
        u_LQR(:,step_test) = K_u_LQR(:,:, mode_now)*estimate_tilde_x_LQR(:,step_test); %基于当前控制器给出反馈控制律
        u_H(:,step_test) = K_u_H(:,:, mode_now)*estimate_tilde_x_H(:,step_test);
        u_TD(:,step_test) = K_u_TD(:,:, mode_now)*estimate_tilde_x_TD(:,step_test);
        
        %         u_LQR(:,step_test) = K_u_LQR(:,:, mode_now)*tilde_x_LQR(:,step_test); %基于当前控制器给出反馈控制律
        %         u_H(:,step_test) = K_u_H(:,:, mode_now)*tilde_x_H(:,step_test);
        %         u_H(:,step_test) = K_u_TD(:,:, mode_now)*tilde_x_TD(:,step_test);
        
        tilde_y_LQR(:,step_test) = tilde_E(:,:,mode_now)*tilde_x_LQR(:,step_test) + tilde_G(:,:, mode_now)*tilde_w_LQR(:,step_test); %基于当前模态以及反馈控制器、状态更新跟踪误差
        tilde_y_H(:,step_test) = tilde_E(:,:,mode_now)*tilde_x_H(:,step_test) + tilde_G(:,:, mode_now)*tilde_w_H(:,step_test);
        tilde_y_TD(:,step_test) = tilde_E(:,:,mode_now)*tilde_x_TD(:,step_test) + tilde_G(:,:, mode_now)*tilde_w_TD(:,step_test);
        
        tilde_z_LQR(:,step_test) = tilde_C(:,:,mode_now)*tilde_x_LQR(:,step_test) + D(mode_now)*u_LQR(:,step_test) + tilde_H(:,:, mode_now)*tilde_w_LQR(:,step_test); %基于当前模态以及反馈控制器、状态更新跟踪误差
        tilde_z_H(:,step_test) = tilde_C(:,:,mode_now)*tilde_x_H(:,step_test) + D(mode_now)*u_H(:,step_test) + tilde_H(:,:, mode_now)*tilde_w_H(:,step_test);
        tilde_z_TD(:,step_test) = tilde_C(:,:,mode_now)*tilde_x_TD(:,step_test) + D(mode_now)*u_TD(:,step_test) + tilde_H(:,:, mode_now)*tilde_w_TD(:,step_test);

        tilde_x_LQR(:,step_test+1) = tilde_A(:,:,mode_now)*tilde_x_LQR(:,step_test) + tilde_B(:,:,mode_now)*u_LQR(:,step_test) + tilde_F(:,:,mode_now)*tilde_w_LQR(:,step_test); %基于当前模态以及反馈控制器更新状态
        tilde_x_H(:,step_test+1) = tilde_A(:,:,mode_now)*tilde_x_H(:,step_test) + tilde_B(:,:,mode_now)*u_H(:,step_test) + tilde_F(:,:,mode_now)*tilde_w_H(:,step_test);
        tilde_x_TD(:,step_test+1) = tilde_A(:,:,mode_now)*tilde_x_TD(:,step_test) + tilde_B(:,:,mode_now)*u_TD(:,step_test) + tilde_F(:,:,mode_now)*tilde_w_TD(:,step_test);
        
        estimate_tilde_y_LQR(:,step_test) = tilde_E(:,:,mode_now)*estimate_tilde_x_LQR(:,step_test); %基于当前模态以及反馈控制器、状态更新跟踪误差
        estimate_tilde_y_H(:,step_test) = tilde_E(:,:,mode_now)*estimate_tilde_x_H(:,step_test);
        estimate_tilde_y_TD(:,step_test) = tilde_E(:,:,mode_now)*estimate_tilde_x_TD(:,step_test);
        
        estimate_tilde_x_LQR(:,step_test+1) = tilde_A(:,:,mode_now)*estimate_tilde_x_LQR(:,step_test) + tilde_B(:,:,mode_now)*u_LQR(:,step_test) - L_y_H_filtering(:,:, mode_now)*(tilde_y_LQR(:,step_test) - estimate_tilde_y_LQR(:,step_test)); %基于当前模态以及反馈控制器更新状态
        estimate_tilde_x_H(:,step_test+1) = tilde_A(:,:,mode_now)*estimate_tilde_x_H(:,step_test) + tilde_B(:,:,mode_now)*u_H(:,step_test) - L_y_H_filtering(:,:, mode_now)*(tilde_y_H(:,step_test) - estimate_tilde_y_H(:,step_test));
        estimate_tilde_x_TD(:,step_test+1) = tilde_A(:,:,mode_now)*estimate_tilde_x_TD(:,step_test) + tilde_B(:,:,mode_now)*u_TD(:,step_test) - L_y_TD_filtering(:,:, mode_now)*(tilde_y_TD(:,step_test) - estimate_tilde_y_TD(:,step_test));

        z_LQR(:,step_test) = C(:,:, mode_now)*tilde_x_LQR(1:A_row,step_test) + D(mode_now)*u_LQR(:,step_test) + H(:,:, mode_now)*tilde_w_LQR(1:F_col,step_test); %基于当前模态以及反馈控制器、状态更新输出
        z_H(:,step_test) = C(:,:, mode_now)*tilde_x_H(1:A_row,step_test) + D(mode_now)*u_H(:,step_test) + H(:,:, mode_now)*tilde_w_H(1:F_col,step_test);
        z_TD(:,step_test) = C(:,:, mode_now)*tilde_x_TD(1:A_row,step_test) + D(mode_now)*u_TD(:,step_test) + H(:,:, mode_now)*tilde_w_TD(1:F_col,step_test);
        hat_z(:,step_test) = hat_C*tilde_x_LQR(A_row+1:A_row+hat_A_row,step_test);
    end
    episode_test
end

figure(11)
subplot(3,1,1)
plot(1:1:(steps_test),delta_tilde_x_H(1,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),delta_tilde_x_TD(1,1:steps_test),'-.','Color','r',"LineWidth",1)
legend('delta-$\tilde{x}_{H}$','delta-$\tilde{x}_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,2)
plot(1:1:(steps_test),delta_tilde_x_H(2,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),delta_tilde_x_TD(2,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
legend('delta-$\tilde{x}_{LQR}$','delta-$\tilde{x}_H$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,3)
plot(1:1:(steps_test),delta_tilde_x_H(3,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),delta_tilde_x_TD(3,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
legend('delta-$\tilde{x}_{LQR}$','delta-$\tilde{x}_H$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

figure(12)
subplot(3,1,1)
plot(1:1:(steps_test),delta_tilde_x_H(4,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),delta_tilde_x_TD(4,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
legend('delta-$\tilde{x}_{LQR}$','delta-$\tilde{x}_H$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,2)
plot(1:1:(steps_test),delta_tilde_x_H(5,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),delta_tilde_x_TD(5,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
legend('delta-$\tilde{x}_{LQR}$','delta-$\tilde{x}_H$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,3)
plot(1:1:(steps_test),delta_tilde_x_H(6,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),delta_tilde_x_TD(6,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
legend('delta-$\tilde{x}_{LQR}$','delta-$\tilde{x}_H$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

figure(13)
subplot(3,1,1)
plot(1:1:(steps_test),tilde_x_H(1,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_H(1,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_TD(1,1:steps_test),'Color','r',"LineWidth",1)
legend('$\tilde{x}$','$\hat_x$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,2)
plot(1:1:(steps_test),tilde_x_H(2,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_H(2,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_TD(2,1:steps_test),'Color','r',"LineWidth",1)
legend('$\tilde{x}$','$\hat_x$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,3)
plot(1:1:(steps_test),tilde_x_H(3,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_H(3,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_TD(3,1:steps_test),'Color','r',"LineWidth",1)
legend('$\tilde{x}$','$\hat_x$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

figure(14)
subplot(3,1,1)
plot(1:1:(steps_test),tilde_x_H(4,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_H(4,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_TD(4,1:steps_test),'Color','r',"LineWidth",1)
legend('$\tilde{x}$','$\hat_x$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,2)
plot(1:1:(steps_test),tilde_x_H(5,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_H(5,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_TD(5,1:steps_test),'Color','r',"LineWidth",1)
legend('$\tilde{x}$','$\hat_x$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

subplot(3,1,3)
plot(1:1:(steps_test),tilde_x_H(6,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_H(6,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),estimate_tilde_x_TD(6,1:steps_test),'Color','r',"LineWidth",1)
legend('$\tilde{x}$','$\hat_x$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

figure(15)
plot(1:1:(steps_test),tilde_z_LQR(1,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),tilde_z_H(1,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),tilde_z_TD(1,1:steps_test),'Color','g',"LineWidth",1)
legend('$\tilde{z}_{LQR}$','$\tilde{z}_H$','$\tilde{z}_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controller ','interpreter','latex')
% axis([1 steps_test -1 1])
xticks([0:10:steps_test])
% yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

figure(16)
plot(1:1:(steps_test),z_H(1,1:steps_test),'--','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),z_TD(1,1:steps_test),'-.','Color','g',"LineWidth",1)
hold on
plot(1:1:(steps_test),hat_z(1,1:steps_test),'Color','k',"LineWidth",1)
legend('$z_H$','$z_{TD}$','$\hat{z}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$z$ and $\hat{z}$','interpreter','latex')
axis([1 steps_test -10 6])
xticks([0:10:steps_test])
yticks([-10:2:6])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);