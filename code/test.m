%% ���п���������
%% ϵͳģ��������
episodes_test = 1; %������������
steps_test = 150;
hat_x = [10 -10 0]'; % ������ϵͳ����
hat_z = [hat_C*hat_x];
ww = zeros(1,steps_test);
www = 0;

%% ϵͳģ������������LQR������
% ����ϵͳ����
x_LQR = [0 0 0]';
u_LQR = zeros(1);
z_LQR = [];
% ����ϵͳ����
tilde_x_LQR = [x_LQR' hat_x']';
tilde_z_LQR = [0];
% Monta Carlo ģ�����
max_tilde_z_LQR = zeros(1,steps_test);
min_tilde_z_LQR = zeros(1,steps_test);
expect_tilde_z_LQR = zeros(1,steps_test);
% ����˥��
uu_LQR = zeros(1,steps_test);
uuu_LQR = 0;
tilde_zz_LQR = zeros(1,steps_test);
tilde_zzz_LQR = 0;
L2_gain_LQR = zeros(1,steps_test);

%% ϵͳģ������������H���������
% ����ϵͳ����
x_H = [0 0 0]';
u_H = zeros(1);
z_H = [];
% ����ϵͳ����
tilde_x_H = [x_H' hat_x']';
tilde_z_H = [0];
% Monta Carlo ģ�����
max_tilde_z_H = zeros(1,steps_test);
min_tilde_z_H = zeros(1,steps_test);
expect_tilde_z_H = zeros(1,steps_test);
% ����˥��
uu_H = zeros(1,steps_test);
uuu_H = 0;
tilde_zz_H = zeros(1,steps_test);
tilde_zzz_H = 0;
L2_gain_H = zeros(1,steps_test);

%% ϵͳģ������������TD������
% ����ϵͳ����
x_TD = [0 0 0]';
u_TD = zeros(1);
z_TD = [];
% ����ϵͳ����
tilde_x_TD = [x_TD' hat_x']';
tilde_z_TD = [0];
% Monta Carlo ģ�����
max_tilde_z_TD = zeros(1,steps_test);
min_tilde_z_TD = zeros(1,steps_test);
expect_tilde_z_TD = zeros(1,steps_test);
% ����˥��
uu_TD = zeros(1,steps_test);
uuu_TD = 0;
tilde_zz_TD = zeros(1,steps_test);
tilde_zzz_TD = 0;
L2_gain_TD = zeros(1,steps_test);

%% �������ϵͳ
for episode_test = 1:episodes_test
    mode_now = randsrc(1,1,[1 2 3;1/3 1/3 1/3]);
    for step_test = 1:steps_test
        tilde_A_now = tilde_A(:,:,mode_now); %�ɵ�ǰģ̬ȷ����������
        tilde_B_now = tilde_B(:,:,mode_now);
        tilde_F_now = tilde_F(:,:,mode_now);
        tilde_C_now = tilde_C(:,:,mode_now);
        tilde_H_now = tilde_H(:,:, mode_now);
        C_now = C(:,:, mode_now);
        D_now = D(mode_now);
        H_now = H(:,:, mode_now);
        
        Q_now = Q(:,:, mode_now); %�ɵ�ǰģ̬ȷ��Ȩ�ؾ���
        R_now = R(mode_now);
        
        K_u_LQR_now = K_u_LQR(:,:, mode_now); %�ɵ�ǰģ̬ȷ������������
        K_u_H_now = K_u_H(:,:, mode_now);
        K_u_TD_now = K_u_TD(:,:, mode_now);
        K_w_now = K_w_H(:,:, mode_now);
        
        mode_old = mode_now; %���浱ǰģ̬
        mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %����ת�Ƹ��ʸ���ģ̬
        
        u_LQR(:,step_test) = K_u_LQR_now*tilde_x_LQR(:,step_test); %���ڵ�ǰ��������������������
        u_H(:,step_test) = K_u_H_now*tilde_x_H(:,step_test);
        u_TD(:,step_test) = K_u_TD_now*tilde_x_TD(:,step_test);
        
%         w(:,step_test) = normrnd(0,0.05,[F_col+hat_F_col 1]); %���������˹������
        w(:,step_test) = K_w_now*tilde_x_H(:,step_test);
        %         w(:,step_test) = normrnd(0,0.1)*[sin(step_test)+cos(step_test)^(2);cos(step_test)+cos(step_test)^(2);sin(2*step_test)+cos(step_test)^(2);cos(step_test)+cos(step_test)^(3);sin(step_test)+cos(2*step_test)^(2);cos(step_test)+sin(step_test)^(2)];
        w_LQR(:,step_test) = w(:,step_test);
        w_H(:,step_test) = w(:,step_test);
        w_TD(:,step_test) = w(:,step_test);
        
        tilde_x_LQR(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_LQR_now)*tilde_x_LQR(:,step_test) + tilde_F_now*w_LQR(:,step_test); %���ڵ�ǰģ̬�Լ���������������״̬
        tilde_x_H(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_H_now)*tilde_x_H(:,step_test) + tilde_F_now*w_H(:,step_test);
        tilde_x_TD(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_TD_now)*tilde_x_TD(:,step_test) + tilde_F_now*w_TD(:,step_test);
        
        tilde_z_LQR(:,step_test) = tilde_C_now*tilde_x_LQR(:,step_test) + D_now*u_LQR(:,step_test) + tilde_H_now*w_LQR(:,step_test); %���ڵ�ǰģ̬�Լ�������������״̬���¸������
        tilde_z_H(:,step_test) = tilde_C_now*tilde_x_H(:,step_test) + D_now*u_H(:,step_test) + tilde_H_now*w_H(:,step_test);
        tilde_z_TD(:,step_test) = tilde_C_now*tilde_x_TD(:,step_test) + D_now*u_TD(:,step_test) + tilde_H_now*w_TD(:,step_test);
        
        z_LQR(:,step_test) = C_now*tilde_x_LQR(1:A_row,step_test) + D_now*u_LQR(:,step_test) + H_now*w_LQR(1:F_col,step_test); %���ڵ�ǰģ̬�Լ�������������״̬�������
        z_H(:,step_test) = C_now*tilde_x_H(1:A_row,step_test) + D_now*u_H(:,step_test) + H_now*w_H(1:F_col,step_test);
        z_TD(:,step_test) = C_now*tilde_x_TD(1:A_row,step_test) + D_now*u_TD(:,step_test) + H_now*w_TD(1:F_col,step_test);
        hat_z(:,step_test) = hat_C*tilde_x_LQR(A_row+1:A_row+hat_A_row,step_test);
        
        ww(:,step_test) = ww(:,step_test) + gamma^(step_test)*w(:,step_test)'*w(:,step_test); %�����ۼ�����
        
        uu_LQR(:,step_test) = uu_LQR(:,step_test) + gamma^(step_test)*u_LQR(:,step_test)'*R_now*u_LQR(:,step_test); %�����ۼƿ�������
        uu_H(:,step_test) = uu_H(:,step_test) + gamma^(step_test)*u_H(:,step_test)'*R_now*u_H(:,step_test);
        uu_TD(:,step_test) = uu_TD(:,step_test) + gamma^(step_test)*u_TD(:,step_test)'*R_now*u_TD(:,step_test);
        
        tilde_zz_LQR(:,step_test) = tilde_zz_LQR(:,step_test) + gamma^(step_test)*tilde_z_LQR(:,step_test)'*Q_now*tilde_z_LQR(:,step_test); %�����ۼƸ������
        tilde_zz_H(:,step_test) = tilde_zz_H(:,step_test) + gamma^(step_test)*tilde_z_H(:,step_test)'*Q_now*tilde_z_H(:,step_test);
        tilde_zz_TD(:,step_test) = tilde_zz_TD(:,step_test) + gamma^(step_test)*tilde_z_TD(:,step_test)'*Q_now*tilde_z_TD(:,step_test);
        
        if episode_test == episodes_test
            www = www + ww(:,step_test); %�����ۼ�����
            
            uuu_LQR = uuu_LQR+uu_LQR(:,step_test); %�����ۼƿ�������
            uuu_H = uuu_H+uu_H(:,step_test);
            uuu_TD = uuu_TD+uu_TD(:,step_test);
            
            tilde_zzz_LQR = tilde_zzz_LQR+tilde_zz_LQR(:,step_test); %�����ۼƸ������
            tilde_zzz_H = tilde_zzz_H+tilde_zz_H(:,step_test);
            tilde_zzz_TD = tilde_zzz_TD+tilde_zz_TD(:,step_test);
            
            L2_gain_LQR(:,step_test) = (uuu_LQR+tilde_zzz_LQR)/www; %����L2����仯
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

% ֱ�ӵ���
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