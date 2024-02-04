clc
clear all

%% 给定系统参数,参数来自 《Mode-Independent H2-Control of a DC Motor Modeled as a Markov Jump Linear System》
%% 控制系统参数
modes = 3;
A(:,:,1) = [-0.4799 5.1546 0;-3.8162 14.4723 0;0.1399 0 -0.9925]; %被控系统状态方程状态转移矩阵
A(:,:,2) = [-1.6026 9.1632 0;-0.5918 3.0317 0;0.0740 0 -0.4338];
A(:,:,3) = [0.6436 0.9178 0;-0.5056 2.4811 0;0.3865 0 0.0982];

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

R = 0.5*[1 1 1];  %系统控制输入权重矩阵

gamma = 0.99;  %衰减因子gamma

%% 根据[x r]增广系统参数生成新的系统参数
for mode=1:modes
    tilde_A(:,:,mode) = blkdiag(A(:,:,mode),hat_A);  %blkdiag函数根据已有矩阵生成准对角矩阵
    tilde_B(:,:,mode) = [B(:,:,mode)' zeros(B_col,hat_A_row)]';
    tilde_C(:,:,mode) = [C(:,:,mode) -hat_C];
end

%% Kleinman迭代参数
P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % 给定解的初始值
sigmma_P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %定义按概率加权求和矩阵
sigmma_V = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);

S_u_initial(:,:,1) = [0.250  -0.850   0.050   0.050   0.050   0.050];
S_u_initial(:,:,2) = [0.201  -1.502   0.100   0.100   0.000   0.100];
S_u_initial(:,:,3) = [0.100  -0.998   0.100   0.100   0.100   0.100];

S_u = S_u_initial; % 将LQR控制器作为迭代初始值

S(:,:,1) = [eye(A_row+hat_A_row);S_u(:,:,1)];
S(:,:,2) = [eye(A_row+hat_A_row);S_u(:,:,2)];
S(:,:,3) = [eye(A_row+hat_A_row);S_u(:,:,3)];

P_episode = zeros(A_row+hat_A_row,A_row+hat_A_row,modes); % 记录迭代过程中每一幕P的值
S_episode = zeros(A_row+hat_A_row+B_col,A_row+hat_A_row,modes);

delta = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);

norm_P_1_episode = [];  %记录解矩阵P的2范数每一幕的变化
norm_P_2_episode = [];
norm_P_3_episode = [];

norm_S_1_episode = [];
norm_S_2_episode = [];
norm_S_3_episode = [];

Gamma = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
Upsilon = zeros(A_row+hat_A_row+B_col,A_row+hat_A_row+B_col,modes);

episodes = 21; %给定迭代次数

%% 求解H无穷控制器
norm_P_1_episode(1) = trace(P(:,:,1)'*P(:,:,1));
norm_P_2_episode(1) = trace(P(:,:,2)'*P(:,:,2));
norm_P_3_episode(1) = trace(P(:,:,3)'*P(:,:,3));
norm_S_1_episode(1) = trace(S_u(:,:,1)'*S_u(:,:,1));
norm_S_2_episode(1) = trace(S_u(:,:,2)'*S_u(:,:,2));
norm_S_3_episode(1) = trace(S_u(:,:,3)'*S_u(:,:,3));

for episode = 1:episodes
    episode
    for mode = 1:modes %基于K_u以及K_w求解P
        Gamma(:,:,mode) = sqrt(gamma)*(tilde_A(:,:,mode) + tilde_B(:,:,mode)*S_u(:,:,mode));  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
        Upsilon(:,:,mode) = [tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode); D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode) D(mode)'*Q(:,:,mode)*D(mode)+R(mode)];
        S_episode(:,:,mode) = [eye(A_row+hat_A_row);S_u(:,:,mode)];
    end
    
    n = 1;
    V = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
    while n < 150  %进行耦合Lyapunov方程的求解
        for mode = 1:modes   %对三个模态进行迭代求解
            sigmma_V(:,:,mode) = Pr(mode,1)*V(:,:,1) + Pr(mode,2)*V(:,:,2) + Pr(mode,3)*V(:,:,3); %计算按概率加权求和的P
            V(:,:,mode) = Gamma(:,:,mode)'*sigmma_V(:,:,mode)*Gamma(:,:,mode) + S_episode(:,:,mode)'*Upsilon(:,:,mode)*S_episode(:,:,mode);  % 迭代求解P
            V(:,:,mode) = (V(:,:,mode)' + V(:,:,mode))/2;
        end
        n = n + 1;
    end
    P = V;
    
    for mode = 1:modes %基于求解的P更新K_u以及K_w
        sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3);  %计算按概率加权求和的P
        
        S_u(:,:,mode) = -inv(D(mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_B(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode) + R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode))';
        S(:,:,mode) = [eye(A_row+hat_A_row);S_u(:,:,mode)];
    end
    P_episode(:,:,:,episode+1) = P;  %更新求解的P
    S_episode(:,:,:,episode+1) = S;  %更新求解的P
    
    norm_P_1_episode(episode+1) = log(trace(P(:,:,1)'*P(:,:,1)));
    norm_P_2_episode(episode+1) = log(trace(P(:,:,2)'*P(:,:,2)));
    norm_P_3_episode(episode+1) = log(trace(P(:,:,3)'*P(:,:,3)));
    norm_S_1_episode(episode+1) = log(trace(S(:,:,1)'*S(:,:,1)));
    norm_S_2_episode(episode+1) = log(trace(S(:,:,2)'*S(:,:,2)));
    norm_S_3_episode(episode+1) = log(trace(S(:,:,3)'*S(:,:,3)));
end

for mode = 1:modes %得到最终结果
    sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3);  %计算按概率加权求和的P
    
    S_u(:,:,mode) = -inv(D(mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_B(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode) + R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode))';
    S(:,:,mode) = [eye(A_row+hat_A_row);S_u(:,:,mode)];
    
    delta(:,:,mode) = P(:,:,mode)-tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode)-gamma*tilde_A(:,:,mode)'*sigmma_P(:,:,mode)*tilde_A(:,:,mode) + (tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode))*inv(gamma*tilde_B(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode)+R(mode)+D(mode)'*Q(:,:,mode)*D(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P(:,:,mode)*tilde_B(:,:,mode))';
end

delta1 = trace(delta(:,:,1)'*delta(:,:,1) + delta(:,:,2)'*delta(:,:,2) + delta(:,:,3)'*delta(:,:,3))

sigmma_P_optimal = sigmma_P; %给定sigmma_P的最优值
S_u_optimal = S_u; %给定控制器的最优值
S_optimal = S;

%% 画图
figure(1)
plot(0:1:(episodes-1),norm_P_1_episode(1:episodes),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes-1),norm_P_2_episode(1:episodes),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes-1),norm_P_3_episode(1:episodes),'Color','g','LineWidth',1.5)
legend('$log(\left\|P_{1}\right\|_{2})$','$log(\left\|P_{2}\right\|_{2})$','$log(\left\|P_{3}\right\|_{2})$','Interpreter','latex'); %legend在坐标区上添加图例
axis([0 episodes-1 0 200]) %调整坐标轴范围 axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes-1)/10:episodes-1]) %设置 x 轴刻度值
ylabel('$log(\left\|P_{i}\right\|_{2})$','interpreter','latex')
yticks([0:20:200])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

figure(2)
plot(0:1:(episodes-1),norm_S_1_episode(1:episodes),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes-1),norm_S_2_episode(1:episodes),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes-1),norm_S_3_episode(1:episodes),'Color','g','LineWidth',1.5)
legend('$log(\left\|K_{1}\right\|_{2})$','$log(\left\|K_{2}\right\|_{2})$','$log(\left\|K_{3}\right\|_{2})$','Interpreter','latex');
axis([0 episodes-1 0 9]) %调整坐标轴范围axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes-1)/10:(episodes-1)])
ylabel('$log(\left\|K_{i}\right\|_{2})$','interpreter','latex')
yticks([0:10:100])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5
annotation(figure(2),'ellipse',[0.7015625 0.311320754716981 0.0427083333333333 0.0451215932914054]); % 创建 ellipse
annotation(figure(2),'arrow',[0.7 0.646354166666667],[0.362683438155136 0.419287211740042]); % 创建 arrow

%% TP未知下求解H无穷控制器迭代参数
sigmma_P_TD(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
sigmma_V_TD(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);

%给定初始镇定控制器
S_u_TD = S_u_initial;

%给定初始镇定噪声增益
S_TD(:,:,1) = [eye(A_row+hat_A_row);S_u_TD(:,:,1)];
S_TD(:,:,2) = [eye(A_row+hat_A_row);S_u_TD(:,:,2)];
S_TD(:,:,3) = [eye(A_row+hat_A_row);S_u_TD(:,:,3)];

Gamma_TD = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
Upsilon_TD = zeros(A_row+hat_A_row+B_col,A_row+hat_A_row+B_col,modes);

delta_sigmma_P_TD = [0 0 0]; %用于储存sigmma_P距离最优的差值
delta_S_TD = [0 0 0]; %用于储存K_u距离最优的差值

episodes_TD = 51; %定义迭代的幕数
steps_TD = 200; %定义每一幕的步数
mu = 0.10; %定义回报权重

%% TP未知下求解H无穷控制器
for mode = 1:modes
    Gamma_TD(:,:,mode) = sqrt(gamma)*(tilde_A(:,:,mode) + tilde_B(:,:,mode)*S_u_TD(:,:,mode));  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
    delta_sigmma_P_TD(1,mode) = log(1+(trace((sigmma_P_optimal(:,:,mode)-sigmma_P_TD(:,:,mode))'*(sigmma_P_optimal(:,:,mode)-sigmma_P_TD(:,:,mode))))/(trace(sigmma_P_optimal(:,:,mode)'*sigmma_P_optimal(:,:,mode))));
    delta_S_TD(1,mode) = log(1+(trace((S_optimal(:,:,mode)-S_TD(:,:,mode))'*(S_optimal(:,:,mode)-S_TD(:,:,mode))))/(trace(S_optimal(:,:,mode)'*S_optimal(:,:,mode))));
end

for episode_TD = 1:episodes_TD
    episode_TD
    
    lambda = 1/(episode_TD); %更新步长
    
    sigmma_V_TD = sigmma_P_TD;
    for mode = 1:modes
        mode_now = mode; %保留每一幕的第一个模态
        TD_sum(:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row);
        for step_test=1:steps_TD
            mode_old = mode_now; %储存当前模态
            mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %基于转移概率更新模态
            
            Gamma_TD(:,:,mode_now) = sqrt(gamma)*(tilde_A(:,:,mode_now)+tilde_B(:,:,mode_now)*S_u_TD(:,:,mode_now));  %记录这一幕各个模态的A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
            Upsilon_TD(:,:,mode_now) = [tilde_C(:,:,mode_now)'*Q(:,:,mode_now)*tilde_C(:,:,mode_now) tilde_C(:,:,mode_now)'*Q(:,:,mode_now)*D(mode_now);D(mode_now)'*Q(:,:,mode_now)*tilde_C(:,:,mode_now) D(mode_now)'*Q(:,:,mode_now)*D(mode_now)+R(mode_now)];
            S_TD(:,:,mode_now) = [eye(A_row+hat_A_row);S_u_TD(:,:,mode_now)];
            
            if step_test == 1
                Lambda(:,:,mode) = eye(A_row+hat_A_row); %定义Lambda矩阵
            else
                Lambda(:,:,mode) = Gamma_TD(:,:,mode_old)*Lambda(:,:,mode); %更新Lambda矩阵
            end
            
            TD = Lambda(:,:,mode)'*(Gamma_TD(:,:,mode_now)'*sigmma_V_TD(:,:,mode_now)*Gamma_TD(:,:,mode_now)+S_TD(:,:,mode_now)'*Upsilon_TD(:,:,mode_now)*S_TD(:,:,mode_now)-sigmma_V_TD(:,:,mode_old))*Lambda(:,:,mode);
            TD_sum = TD_sum + mu^(step_test-1)*TD; %进行第step步的sigmma_P更新
        end
        sigmma_P_TD(:,:,mode) = sigmma_V_TD(:,:,mode) + lambda*TD_sum; %进行第step步的sigmma_P更新
        S_u_TD(:,:,mode) = -inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_A(:,:,mode));
    end
    for mode = 1:modes
        delta_sigmma_P_TD(1,mode,episode_TD+1) = log(1 + (trace((sigmma_P_optimal(:,:,mode) - sigmma_P_TD(:,:,mode))'*(sigmma_P_optimal(:,:,mode)-sigmma_P_TD(:,:,mode))))/(trace(sigmma_P_optimal(:,:,mode)'*sigmma_P_optimal(:,:,mode))));
        delta_S_TD(1,mode,episode_TD + 1) = log(1 + (trace((S_optimal(:,:,mode) - S_TD(:,:,mode))'*(S_optimal(:,:,mode) - S_TD(:,:,mode))))/(trace(S_optimal(:,:,mode)'*S_optimal(:,:,mode))));
    end
end

delta_sigmma_P1 = delta_sigmma_P_TD(1,1,:);
delta_sigmma_P1 = delta_sigmma_P1(:);  %通过（：）将多维数组变为一维数组
delta_sigmma_P2 = delta_sigmma_P_TD(1,2,:);
delta_sigmma_P2 = delta_sigmma_P2(:);
delta_sigmma_P3 = delta_sigmma_P_TD(1,3,:);
delta_sigmma_P3 = delta_sigmma_P3(:);

delta_S_TD_1 = delta_S_TD(1,1,:);
delta_S_TD_1 = delta_S_TD_1(:);  %通过（：）将多维数组变为一维数组
delta_S_TD_2 = delta_S_TD(1,2,:);
delta_S_TD_2 = delta_S_TD_2(:);
delta_S_TD_3 = delta_S_TD(1,3,:);
delta_S_TD_3 = delta_S_TD_3(:);

figure(3)
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
xticks([0:10:episodes_TD-1])
% yticks([0:0.1:0.8])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

figure(4)
plot(0:1:(episodes_TD-1),delta_S_TD_1(1:episodes_TD),'--','Color','b',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_S_TD_2(1:episodes_TD),'-.','Color','r',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_S_TD_3(1:episodes_TD),'Color','g',"LineWidth",1.5)
legend('$\delta^{(l)}_{1}$','$\delta^{(l)}_{2}$','$\delta^{(l)}_{3}$','Interpreter','latex','Position', ...
    [0.691145833333333 0.212788259958072 0.188194444444442 0.205170875242091]); %legend在坐标区上添加图例
xlabel('$Iteration$','interpreter','latex')
ylabel('error of controller','interpreter','latex')
axis([0 episodes_TD-1 -0.0005 0.040])
xticks([0:10:episodes_TD-1])
% yticks([0:0.00005:0.0004])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %设置坐标轴字体为Times New Roman，大小为26，线宽0.5

%% 进行控制器测试
%% 系统模拟仿真参数
episodes_test = 1; %给定迭代次数
steps_test = 150; % 运行步长

% 被控系统参数
x = [0 0 0]';
u = zeros(1);
z = [];

% 被跟踪系统参数
hat_x = [10 -10 0]'; 
hat_z = [hat_C*hat_x];

% 增广系统参数
tilde_x = [x' hat_x']';
tilde_z = [0];

%% 系统模拟仿真参数处于TD控制器
% 被控系统参数
x_TD = [0 0 0]';
u_TD = zeros(1);
z_TD = [];

% 增广系统参数
tilde_x_TD = [x_TD' hat_x']';
tilde_z_TD = [0];

%% 仿真跟踪系统
mode_now = randsrc(1,1,[1 2 3;1/3 1/3 1/3]);
for step_test = 1:steps_test
    mode_old = mode_now; %储存当前模态
    mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %基于转移概率更新模态

    u(:,step_test) = S_u(:,:, mode_now)*tilde_x(:,step_test); %基于当前控制器给出反馈控制律
    u_TD(:,step_test) = S_u_TD(:,:, mode_now)*tilde_x_TD(:,step_test);

    tilde_x(:,step_test+1) = tilde_A(:,:,mode_now)*tilde_x(:,step_test) + tilde_B(:,:,mode_now)*u(:,step_test); %基于当前模态以及反馈控制器更新状态
    tilde_x_TD(:,step_test+1) = tilde_A(:,:,mode_now)*tilde_x_TD(:,step_test) + tilde_B(:,:,mode_now)*u_TD(:,step_test);

    tilde_z(:,step_test) = tilde_C(:,:,mode_now)*tilde_x(:,step_test) + D(mode_now)*u(:,step_test); %基于当前模态以及反馈控制器、状态更新跟踪误差
    tilde_z_TD(:,step_test) = tilde_C(:,:,mode_now)*tilde_x_TD(:,step_test) + D(mode_now)*u_TD(:,step_test) ;

    z(:,step_test) = C(:,:,mode_now)*tilde_x(1:A_row,step_test) + D(mode_now)*u(:,step_test); %基于当前模态以及反馈控制器、状态更新输出
    z_TD(:,step_test) = C(:,:,mode_now)*tilde_x_TD(1:A_row,step_test) + D(mode_now)*u_TD(:,step_test);
    hat_z(:,step_test) = hat_C*tilde_x(A_row+1:A_row+hat_A_row,step_test);
end

figure(5)
plot(1:1:(steps_test),z(1,1:steps_test),'--','Color','r',"LineWidth",1)
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