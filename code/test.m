%% 进行控制器测试
%% 系统模拟仿真参数
episodes_test = 1; %给定迭代次数
steps_test = 150;
hat_x = [10 -10 0]'; % 被跟踪系统参数
hat_z = [hat_C*hat_x];
ww = zeros(1,steps_test);
www = 0;

%% 系统模拟仿真参数处于LQR控制器
% 被控系统参数
x_LQR = [0 0 0]';
u_LQR = zeros(1);
z_LQR = [];
% 增广系统参数
tilde_x_LQR = [x_LQR' hat_x']';
tilde_z_LQR = [0];
% Monta Carlo 模拟变量
max_tilde_z_LQR = zeros(1,steps_test);
min_tilde_z_LQR = zeros(1,steps_test);
expect_tilde_z_LQR = zeros(1,steps_test);
% 噪声衰减
uu_LQR = zeros(1,steps_test);
uuu_LQR = 0;
tilde_zz_LQR = zeros(1,steps_test);
tilde_zzz_LQR = 0;
L2_gain_LQR = zeros(1,steps_test);

%% 系统模拟仿真参数处于H无穷控制器
% 被控系统参数
x_H = [0 0 0]';
u_H = zeros(1);
z_H = [];
% 增广系统参数
tilde_x_H = [x_H' hat_x']';
tilde_z_H = [0];
% Monta Carlo 模拟变量
max_tilde_z_H = zeros(1,steps_test);
min_tilde_z_H = zeros(1,steps_test);
expect_tilde_z_H = zeros(1,steps_test);
% 噪声衰减
uu_H = zeros(1,steps_test);
uuu_H = 0;
tilde_zz_H = zeros(1,steps_test);
tilde_zzz_H = 0;
L2_gain_H = zeros(1,steps_test);

%% 系统模拟仿真参数处于TD控制器
% 被控系统参数
x_TD = [0 0 0]';
u_TD = zeros(1);
z_TD = [];
% 增广系统参数
tilde_x_TD = [x_TD' hat_x']';
tilde_z_TD = [0];
% Monta Carlo 模拟变量
max_tilde_z_TD = zeros(1,steps_test);
min_tilde_z_TD = zeros(1,steps_test);
expect_tilde_z_TD = zeros(1,steps_test);
% 噪声衰减
uu_TD = zeros(1,steps_test);
uuu_TD = 0;
tilde_zz_TD = zeros(1,steps_test);
tilde_zzz_TD = 0;
L2_gain_TD = zeros(1,steps_test);

%% 仿真跟踪系统
for episode_test = 1:episodes_test
    mode_now = randsrc(1,1,[1 2 3;1/3 1/3 1/3]);
    for step_test = 1:steps_test
        tilde_A_now = tilde_A(:,:,mode_now); %由当前模态确定参数矩阵
        tilde_B_now = tilde_B(:,:,mode_now);
        tilde_F_now = tilde_F(:,:,mode_now);
        tilde_C_now = tilde_C(:,:,mode_now);
        tilde_H_now = tilde_H(:,:, mode_now);
        C_now = C(:,:, mode_now);
        D_now = D(mode_now);
        H_now = H(:,:, mode_now);
        
        Q_now = Q(:,:, mode_now); %由当前模态确定权重矩阵
        R_now = R(mode_now);
        
        K_u_LQR_now = K_u_LQR(:,:, mode_now); %由当前模态确定控制器矩阵
        K_u_H_now = K_u_H(:,:, mode_now);
        K_u_TD_now = K_u_TD(:,:, mode_now);
        K_w_now = K_w_H(:,:, mode_now);
        
        mode_old = mode_now; %储存当前模态
        mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %基于转移概率更新模态
        
        u_LQR(:,step_test) = K_u_LQR_now*tilde_x_LQR(:,step_test); %基于当前控制器给出反馈控制律
        u_H(:,step_test) = K_u_H_now*tilde_x_H(:,step_test);
        u_TD(:,step_test) = K_u_TD_now*tilde_x_TD(:,step_test);
        
%         w(:,step_test) = normrnd(0,0.05,[F_col+hat_F_col 1]); %给出随机高斯白噪声
        w(:,step_test) = K_w_now*tilde_x_H(:,step_test);
        %         w(:,step_test) = normrnd(0,0.1)*[sin(step_test)+cos(step_test)^(2);cos(step_test)+cos(step_test)^(2);sin(2*step_test)+cos(step_test)^(2);cos(step_test)+cos(step_test)^(3);sin(step_test)+cos(2*step_test)^(2);cos(step_test)+sin(step_test)^(2)];
        w_LQR(:,step_test) = w(:,step_test);
        w_H(:,step_test) = w(:,step_test);
        w_TD(:,step_test) = w(:,step_test);
        
        tilde_x_LQR(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_LQR_now)*tilde_x_LQR(:,step_test) + tilde_F_now*w_LQR(:,step_test); %基于当前模态以及反馈控制器更新状态
        tilde_x_H(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_H_now)*tilde_x_H(:,step_test) + tilde_F_now*w_H(:,step_test);
        tilde_x_TD(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_TD_now)*tilde_x_TD(:,step_test) + tilde_F_now*w_TD(:,step_test);
        
        tilde_z_LQR(:,step_test) = tilde_C_now*tilde_x_LQR(:,step_test) + D_now*u_LQR(:,step_test) + tilde_H_now*w_LQR(:,step_test); %基于当前模态以及反馈控制器、状态更新跟踪误差
        tilde_z_H(:,step_test) = tilde_C_now*tilde_x_H(:,step_test) + D_now*u_H(:,step_test) + tilde_H_now*w_H(:,step_test);
        tilde_z_TD(:,step_test) = tilde_C_now*tilde_x_TD(:,step_test) + D_now*u_TD(:,step_test) + tilde_H_now*w_TD(:,step_test);
        
        z_LQR(:,step_test) = C_now*tilde_x_LQR(1:A_row,step_test) + D_now*u_LQR(:,step_test) + H_now*w_LQR(1:F_col,step_test); %基于当前模态以及反馈控制器、状态更新输出
        z_H(:,step_test) = C_now*tilde_x_H(1:A_row,step_test) + D_now*u_H(:,step_test) + H_now*w_H(1:F_col,step_test);
        z_TD(:,step_test) = C_now*tilde_x_TD(1:A_row,step_test) + D_now*u_TD(:,step_test) + H_now*w_TD(1:F_col,step_test);
        hat_z(:,step_test) = hat_C*tilde_x_LQR(A_row+1:A_row+hat_A_row,step_test);
        
        ww(:,step_test) = ww(:,step_test) + gamma^(step_test)*w(:,step_test)'*w(:,step_test); %计算累计噪声
        
        uu_LQR(:,step_test) = uu_LQR(:,step_test) + gamma^(step_test)*u_LQR(:,step_test)'*R_now*u_LQR(:,step_test); %计算累计控制输入
        uu_H(:,step_test) = uu_H(:,step_test) + gamma^(step_test)*u_H(:,step_test)'*R_now*u_H(:,step_test);
        uu_TD(:,step_test) = uu_TD(:,step_test) + gamma^(step_test)*u_TD(:,step_test)'*R_now*u_TD(:,step_test);
        
        tilde_zz_LQR(:,step_test) = tilde_zz_LQR(:,step_test) + gamma^(step_test)*tilde_z_LQR(:,step_test)'*Q_now*tilde_z_LQR(:,step_test); %计算累计跟踪误差
        tilde_zz_H(:,step_test) = tilde_zz_H(:,step_test) + gamma^(step_test)*tilde_z_H(:,step_test)'*Q_now*tilde_z_H(:,step_test);
        tilde_zz_TD(:,step_test) = tilde_zz_TD(:,step_test) + gamma^(step_test)*tilde_z_TD(:,step_test)'*Q_now*tilde_z_TD(:,step_test);
        
        if episode_test == episodes_test
            www = www + ww(:,step_test); %计算累计噪声
            
            uuu_LQR = uuu_LQR+uu_LQR(:,step_test); %计算累计控制输入
            uuu_H = uuu_H+uu_H(:,step_test);
            uuu_TD = uuu_TD+uu_TD(:,step_test);
            
            tilde_zzz_LQR = tilde_zzz_LQR+tilde_zz_LQR(:,step_test); %计算累计跟踪误差
            tilde_zzz_H = tilde_zzz_H+tilde_zz_H(:,step_test);
            tilde_zzz_TD = tilde_zzz_TD+tilde_zz_TD(:,step_test);
            
            L2_gain_LQR(:,step_test) = (uuu_LQR+tilde_zzz_LQR)/www; %计算L2增益变化
            L2_gain_H(:,step_test) = (uuu_H+tilde_zzz_H)/www;
            L2_gain_TD(:,step_test) = (uuu_TD+tilde_zzz_TD)/www;
        end
        
        episode_test
        if episode_test == 1
            max_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
            min_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
            max_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
            min_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
            max_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
            min_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
        else
            if max_tilde_z_LQR(:,step_test) < tilde_z_LQR(:,step_test)
                max_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
            end
            if min_tilde_z_LQR(:,step_test) > tilde_z_LQR(:,step_test)
                min_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
            end
            if max_tilde_z_H(:,step_test) < tilde_z_H(:,step_test)
                max_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
            end
            if min_tilde_z_H(:,step_test) > tilde_z_H(:,step_test)
                min_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
            end
            if max_tilde_z_TD(:,step_test) < tilde_z_TD(:,step_test)
                max_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
            end
            if min_tilde_z_TD(:,step_test) > tilde_z_TD(:,step_test)
                min_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
            end
        end
        expect_tilde_z_LQR(:,step_test) = expect_tilde_z_LQR(:,step_test) + (tilde_z_LQR(:,step_test) - expect_tilde_z_LQR(:,step_test))/episode_test;
        expect_tilde_z_H(:,step_test) = expect_tilde_z_H(:,step_test) + (tilde_z_H(:,step_test) - expect_tilde_z_H(:,step_test))/episode_test;
        expect_tilde_z_TD(:,step_test) = expect_tilde_z_TD(:,step_test) + (tilde_z_TD(:,step_test) - expect_tilde_z_TD(:,step_test))/episode_test;
    end
end

% 直接迭代
figure(6)
fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_LQR(1,1:steps_test),fliplr(min_tilde_z_LQR(1,1:steps_test))],'b','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
hold on
fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_H(1,1:steps_test),fliplr(min_tilde_z_H(1,1:steps_test))],'r','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
hold on
plot(1:1:(steps_test),expect_tilde_z_LQR(1,1:steps_test),'Color','b',"LineWidth",0.5)
hold on
plot(1:1:(steps_test),expect_tilde_z_H(1,1:steps_test),'Color','r',"LineWidth",0.5)
legend('$z_{LQR}$','$z_H$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$x_{1}$ and $\check{x}_{1}$','interpreter','latex')
% axis([1 steps -2.5 2.5])
% xticks([0:steps/10:steps])
% yticks([-2.5:0.5:2.5])
% set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5);

figure(7)
fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_H(1,1:steps_test),fliplr(min_tilde_z_H(1,1:steps_test))],'r','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
hold on
fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_TD(1,1:steps_test),fliplr(min_tilde_z_TD(1,1:steps_test))],'g','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
hold on
plot(1:1:(steps_test),expect_tilde_z_H(1,1:steps_test),'Color','r',"LineWidth",0.5)
hold on
plot(1:1:(steps_test),expect_tilde_z_TD(1,1:steps_test),'Color','g',"LineWidth",0.5)
legend('$z_H$','$z_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$x_{1}$ and $\check{x}_{1}$','interpreter','latex')
% axis([1 steps -2.5 2.5])
% xticks([0:steps/10:steps])
% yticks([-2.5:0.5:2.5])
% set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5);

figure(8)
fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_LQR(1,1:steps_test),fliplr(min_tilde_z_LQR(1,1:steps_test))],'b','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
hold on
fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_H(1,1:steps_test),fliplr(min_tilde_z_H(1,1:steps_test))],'r','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
hold on
fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_TD(1,1:steps_test),fliplr(min_tilde_z_TD(1,1:steps_test))],'g','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
hold on
plot(1:1:(steps_test),expect_tilde_z_LQR(1,1:steps_test),'--','Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),expect_tilde_z_H(1,1:steps_test),'-.','Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),expect_tilde_z_TD(1,1:steps_test),'Color','g',"LineWidth",1)
legend('$\tilde{z}_{LQR}$','$\tilde{z}_H$','$\tilde{z}_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$\tilde{z}$ under differant controllers ','interpreter','latex')
axis([1 steps_test -1 1])
xticks([0:10:steps_test])
yticks([-1:0.5:1])
set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

figure(9)
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

figure(10)
plot(1:1:(steps_test),L2_gain_LQR(1,1:steps_test),'Color','b',"LineWidth",1)
hold on
plot(1:1:(steps_test),L2_gain_H(1,1:steps_test),'Color','r',"LineWidth",1)
hold on
plot(1:1:(steps_test),L2_gain_TD(1,1:steps_test),'Color','g',"LineWidth",1)
legend('$L_{2}-gain_{LQR}$','$L_{2}-gain_{H}$','$L_{2}-gain_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
xlabel('$Steps$','interpreter','latex')
ylabel('$L_{2}-gain$ under differant controller','interpreter','latex')
axis([1 steps_test 0 0.8])
xticks([0:steps_test/10:steps_test])
yticks([0:10:100])
set(gca,"FontName","Times New Roman","FontSize",32,"LineWidth",1);