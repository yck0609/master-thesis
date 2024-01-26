clc
clear all

%% ����ϵͳ����,�������� ��Mode-Independent H2-Control of a DC Motor Modeled as a Markov Jump Linear System��
%% ����ϵͳ����
modes = 3;
A(:,:,1) = [-0.4799 5.1546 0;-3.8162 14.4723 0;0.1399 0 -0.9925]; %����ϵͳ״̬����״̬ת�ƾ���
A(:,:,2) = [-1.6026 9.1632 0;-0.5918 3.0317 0;0.0740 0 -0.4338];
A(:,:,3) = [0.6436 0.9178 0;-0.5056 2.4811 0;0.3865 0 0.0982];

B(:,:,1) = [5.8705 15.5010 0]'; %����ϵͳ״̬���̿��������������
B(:,:,2) = [10.2851 2.2282 0]';
B(:,:,3) = [0.7874 1.5302 0]';

F(:,:,1) = [0.10 0 0;0.10 0 0;0.10 0 0]; %����ϵͳ״̬���������������
F(:,:,2) = [0.10 0 0;0.10 0 0;0.10 0 0];
F(:,:,3) = [0.10 0 0;0.10 0 0;0.10 0 0];

C(:,:,1) = [1.0230 2.1100 0.9500]; %����ϵͳ��������������
C(:,:,2) = [0.9800 2.0500 1.1000];
C(:,:,3) = [1.0000 2.0000 1.0500];

D = [1.000 0.5000 -0.5000]; %����ϵͳ������̿��������������

H(:,:,1) = [0 0.10 0.15]; %����ϵͳ������������������
H(:,:,2) = [0 0.10 0.15];
H(:,:,3) = [0 0.10 0.15];

Pr(:,:) = [0.95 0.05 0;0.36 0.6 0.04;0.1 0.1 0.8]; % ����ϵͳģ̬ת�Ƹ���

%% ����ϵͳ����
alpha = 0.15;
a = 0.4;
b = 0.5;
c = sqrt(1-a^2-b^2); %�����ϲ���������ά���Ҳ�
hat_A = [cos(alpha)+(1-cos(alpha))*a^2   (1-cos(alpha))*a*b-sin(alpha)*c  (1-cos(alpha))*a*c+sin(alpha)*b;
    (1-cos(alpha))*a*b+sin(alpha)*c  cos(alpha)+(1-cos(alpha))*b^2    (1-cos(alpha))*b*c-sin(alpha)*a;
    (1-cos(alpha))*a*c-sin(alpha)*b  (1-cos(alpha))*b*c+sin(alpha)*a  cos(alpha)+(1-cos(alpha))*c^2;];

hat_F = [0.2 0 0;0 0.2 0;0 0 0.2]; %������ϵͳ״̬���������������

hat_C = [1 1 1]; %������ϵͳ״̬���������������

%% ��������ָ�������������Q��R��˥������gamma
[A_row,~] = size(A(:,:,1));  %size����Ĭ�Ϸ��ؾ���ά��,��ȡ����ϵͳ״̬��ά��
[~,B_col] = size(B(:,:,1));  %��ȡ����ϵͳ���������ά��
[~,F_col] = size(F(:,:,1));  %��ȡ����ϵͳ������ά��
[C_row,~] = size(C(:,:,1));  %��ȡ����ϵͳ�����ά��
[hat_A_row,~] = size(hat_A); %��ȡ������ϵͳ״̬��ά��
[~,hat_F_col] = size(hat_F); %��ȡ������ϵͳ���������ά��

Q(:,:,1) = 10*eye(C_row); %ϵͳ�������Ȩ�ؾ���
Q(:,:,2) = 10*eye(C_row);
Q(:,:,3) = 10*eye(C_row);

R = 0.2*[1 1 1];  %ϵͳ��������Ȩ�ؾ���
% R = 0.01*[1 1 1];  %ϵͳ��������Ȩ�ؾ���

gamma = 0.99;  %˥������gamma

theta = 5.35; %L2����

%% ����[x r]����ϵͳ���������µ�ϵͳ����
for mode=1:modes
    tilde_A(:,:,mode) = blkdiag(A(:,:,mode),hat_A);  %blkdiag�����������о�������׼�ԽǾ���
    tilde_B(:,:,mode) = [B(:,:,mode)' zeros(B_col,hat_A_row)]';
    tilde_F(:,:,mode) = blkdiag(F(:,:,mode),hat_F);
    tilde_C(:,:,mode) = [C(:,:,mode) -hat_C];
    tilde_H(:,:,mode) = [H(:,:,mode) zeros(C_row,hat_F_col)];
end

%% LQR��������������
P_LQR = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % ������ĳ�ʼֵ
sigmma_P_LQR = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %���尴���ʼ�Ȩ��;���
K_u_LQR(:,:,:) = zeros(B_col,A_row+hat_A_row,modes);  % ������������K_u_LQR�ĳ�ʼֵ
norm_LQR = [];
episodes_LQR = 1000; %����LQR����������������

%% ���LQR������
for episode_LQR = 1:episodes_LQR
    for mode = 1:modes
        sigmma_P_LQR(:,:,mode) = Pr(mode,1)*P_LQR(:,:,1) + Pr(mode,2)*P_LQR(:,:,2) + Pr(mode,3)*P_LQR(:,:,3); %�ɵ�ǰP����sigmma_P
        P_LQR(:,:,mode) = tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_A(:,:,mode)-(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_B(:,:,mode))*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_A(:,:,mode));
        P_LQR(:,:,mode) = (P_LQR(:,:,mode)' + P_LQR(:,:,mode))/2;
    end
    norm(episode_LQR) = trace(P_LQR(:,:,1)'*P_LQR(:,:,1)+P_LQR(:,:,2)'*P_LQR(:,:,2)+P_LQR(:,:,3)'*P_LQR(:,:,3));
end

for mode = 1:modes
    sigmma_P_LQR(:,:,mode) = Pr(mode,1)*P_LQR(:,:,1)+Pr(mode,2)*P_LQR(:,:,2) + Pr(mode,3)*P_LQR(:,:,3);  %���㰴���ʼ�Ȩ��͵�P
    K_u_LQR(:,:,mode) = -inv(gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_B(:,:,mode) + R(mode))*(D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_LQR(:,:,mode)*tilde_A(:,:,mode));
end

figure(1)
plot(0:1:(episodes_LQR-1),norm)

%% H�����������������
P_H = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % ������ĳ�ʼֵ
sigmma_P_H = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %���尴���ʼ�Ȩ��;���
sigmma_V_H = zeros(A_row + hat_A_row,A_row + hat_A_row,modes);

% K_u_initial(:,:,1) = [0.2003  -0.9257   0.0299   0.0417   0.0298   0.0408];
% K_u_initial(:,:,2) = [0.1697  -1.0012   0.0303   0.0723   0.0588   0.0713];
% K_u_initial(:,:,3) = [0.0457  -1.5005   0.0150   0.2497   0.2302   0.2511];

K_u_initial(:,:,1) = [0.250  -0.850   0.050   0.050   0.050   0.050];
K_u_initial(:,:,2) = [0.201  -1.502   0.100   0.100   0.000   0.100];
K_u_initial(:,:,3) = [0.100  -0.998   0.100   0.100   0.100   0.100];

K_u_H = K_u_initial; % ��LQR��������Ϊ������ʼֵ
% K_u_H = K_u_LQR; % ��LQR��������Ϊ������ʼֵ
K_w_H(:,:,:) = zeros(F_col + hat_F_col,A_row+hat_A_row,modes);  % ������������K_w�ĳ�ʼֵ

K_H(:,:,1) = [eye(A_row+hat_A_row);K_u_H(:,:,1);K_w_H(:,:,1)];
K_H(:,:,2) = [eye(A_row+hat_A_row);K_u_H(:,:,2);K_w_H(:,:,2)];
K_H(:,:,3) = [eye(A_row+hat_A_row);K_u_H(:,:,3);K_w_H(:,:,3)];

P_H_episode = zeros(A_row+hat_A_row,A_row+hat_A_row,modes); % ��¼����������ÿһĻP��ֵ
K_H_episode = zeros(A_row+hat_A_row+B_col+F_col+hat_F_col,A_row+hat_A_row,modes);

norm_P_H_1_episode = [];  %��¼�����P��2����ÿһĻ�ı仯
norm_P_H_2_episode = [];
norm_P_H_3_episode = [];

norm_K_H_1_episode = [];
norm_K_H_2_episode = [];
norm_K_H_3_episode = [];

Gamma_H_episode = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
Upsilon_H_episode = zeros(A_row+hat_A_row+B_col+F_col+hat_F_col,A_row+hat_A_row+B_col+F_col+hat_F_col,modes);

episodes_H = 31; %������������

%% ���H���������
norm_P_H_1_episode(1) = trace(P_H(:,:,1)'*P_H(:,:,1));
norm_P_H_2_episode(1) = trace(P_H(:,:,2)'*P_H(:,:,2));
norm_P_H_3_episode(1) = trace(P_H(:,:,3)'*P_H(:,:,3));
norm_K_H_1_episode(1) = trace(K_H(:,:,1)'*K_H(:,:,1));
norm_K_H_2_episode(1) = trace(K_H(:,:,2)'*K_H(:,:,2));
norm_K_H_3_episode(1) = trace(K_H(:,:,3)'*K_H(:,:,3));

for episode_H = 1:episodes_H
    episode_H
    for mode = 1:modes %����K_u�Լ�K_w���P
        Gamma_H_episode(:,:,mode) = sqrt(gamma)*(tilde_A(:,:,mode) + tilde_B(:,:,mode)*K_u_H(:,:,mode) + tilde_F(:,:,mode)*K_w_H(:,:,mode));  %��¼��һĻ����ģ̬��A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
        Upsilon_H_episode(:,:,mode) = [tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode); ...
            D(mode)'*Q(:,:,mode)*tilde_C(:,:,mode) D(mode)'*Q(:,:,mode)*D(mode)+R(mode) D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode); ...
            tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode) tilde_H(:,:,mode)'*Q(:,:,mode)*D(mode) tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)-theta^(2)*eye(F_col + hat_F_col)];
        K_H_episode(:,:,mode) = [eye(A_row+hat_A_row);K_u_H(:,:,mode);K_w_H(:,:,mode)];
    end
    
    n = 1;
    V_H = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
    while n < 150  %�������Lyapunov���̵����
        for mode = 1:modes   %������ģ̬���е������
            sigmma_V_H(:,:,mode) = Pr(mode,1)*V_H(:,:,1)+Pr(mode,2)*V_H(:,:,2)+Pr(mode,3)*V_H(:,:,3); %���㰴���ʼ�Ȩ��͵�P
            V_H(:,:,mode) = Gamma_H_episode(:,:,mode)'*sigmma_V_H(:,:,mode)*Gamma_H_episode(:,:,mode) + K_H_episode(:,:,mode)'*Upsilon_H_episode(:,:,mode)*K_H_episode(:,:,mode);  % �������P
            V_H(:,:,mode) = (V_H(:,:,mode)' + V_H(:,:,mode))/2;
        end
        n = n + 1;
    end
    P_H = V_H;
    
    for mode = 1:modes %��������P����K_u�Լ�K_w
        sigmma_P_H(:,:,mode) = Pr(mode,1)*P_H(:,:,1)+Pr(mode,2)*P_H(:,:,2)+Pr(mode,3)*P_H(:,:,3);  %���㰴���ʼ�Ȩ��͵�P
        
        K_u_H(:,:,mode) = inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))') ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))');
        K_w_H(:,:,mode) = inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col+hat_F_col)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))) ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))');
        K_H(:,:,mode) = [eye(A_row+hat_A_row);K_u_H(:,:,mode);K_w_H(:,:,mode)];
    end
    P_H_episode(:,:,:,episode_H+1) = P_H;  %��������P
    K_H_episode(:,:,:,episode_H+1) = K_H;  %��������P
    
    norm_P_H_1_episode(episode_H+1) = log(trace(P_H(:,:,1)'*P_H(:,:,1)));
    norm_P_H_2_episode(episode_H+1) = log(trace(P_H(:,:,2)'*P_H(:,:,2)));
    norm_P_H_3_episode(episode_H+1) = log(trace(P_H(:,:,3)'*P_H(:,:,3)));
    norm_K_H_1_episode(episode_H+1) = log(trace(K_H(:,:,1)'*K_H(:,:,1)));
    norm_K_H_2_episode(episode_H+1) = log(trace(K_H(:,:,2)'*K_H(:,:,2)));
    norm_K_H_3_episode(episode_H+1) = log(trace(K_H(:,:,3)'*K_H(:,:,3)));
end

for mode = 1:modes %�õ����ս��
    sigmma_P_H(:,:,mode) = Pr(mode,1)*P_H(:,:,1)+Pr(mode,2)*P_H(:,:,2)+Pr(mode,3)*P_H(:,:,3);  %���㰴���ʼ�Ȩ��͵�P
    
    K_u_H(:,:,mode) = inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))') ...
        *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))');
    K_w_H(:,:,mode) = inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col+hat_F_col)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))) ...
        *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))');
    K_H(:,:,mode) = [eye(A_row+hat_A_row);K_u_H(:,:,mode);K_w_H(:,:,mode)];
    
    delta_H(:,:,mode) = P_H(:,:,mode)-tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_C(:,:,mode)-gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_A(:,:,mode)+[(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)) tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)]*inv([(gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)+R(mode)+D(mode)'*Q(:,:,mode)*D(mode)) (D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode));(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode))' (tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))])*[(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_B(:,:,mode)) tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_H(:,:,mode)*tilde_F(:,:,mode)]';
end

delta_H = trace(delta_H(:,:,1)'*delta_H(:,:,1)+delta_H(:,:,2)'*delta_H(:,:,2)+delta_H(:,:,3)'*delta_H(:,:,3))

sigmma_P_optimal = sigmma_P_H; %����sigmma_P������ֵ
K_u_optimal = K_u_H; %����������������ֵ
K_w_optimal = K_w_H; % ������������K_w������ֵ
K_optimal = K_H;

%% ��ͼ
figure(2)
plot(0:1:(episodes_H-1),norm_P_H_1_episode(1:episodes_H),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_P_H_2_episode(1:episodes_H),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_P_H_3_episode(1:episodes_H),'Color','g','LineWidth',1.5)
legend('$log(\left\|P_{1}\right\|_{2})$','$log(\left\|P_{2}\right\|_{2})$','$log(\left\|P_{3}\right\|_{2})$','Interpreter','latex'); %legend�������������ͼ��
axis([0 episodes_H-1 0 10.5]) %���������᷶Χ axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes_H-1)/10:episodes_H-1]) %���� x ��̶�ֵ
ylabel('$log(\left\|P_{i}\right\|_{2})$','interpreter','latex')
yticks([0:1:10])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %��������������ΪTimes New Roman����СΪ26���߿�0.5

figure(3)
plot(0:1:(episodes_H-1),norm_K_H_1_episode(1:episodes_H),'--','Color','b','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_K_H_2_episode(1:episodes_H),'-.','Color','r','LineWidth',1.5)
hold on
plot(0:1:(episodes_H-1),norm_K_H_3_episode(1:episodes_H),'Color','g','LineWidth',1.5)
legend('$log(\left\|K_{1}\right\|_{2})$','$log(\left\|K_{2}\right\|_{2})$','$log(\left\|K_{3}\right\|_{2})$','Interpreter','latex');
axis([0 episodes_H-1 0 9]) %���������᷶Χaxis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes_H-1)/10:(episodes_H-1)])
ylabel('$log(\left\|K_{i}\right\|_{2})$','interpreter','latex')
yticks([0:1:9])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %��������������ΪTimes New Roman����СΪ26���߿�0.5
annotation(figure(3),'ellipse',[0.7015625 0.311320754716981 0.0427083333333333 0.0451215932914054]); % ���� ellipse
annotation(figure(3),'arrow',[0.7 0.646354166666667],[0.362683438155136 0.419287211740042]); % ���� arrow

%% TPδ֪�����H�����������������
sigmma_P_TD(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
sigmma_V_TD(:,:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);

%������ʼ�򶨿�����
K_u_TD = K_u_initial;
K_w_TD(:,:,:) = zeros(F_col + hat_F_col,A_row+hat_A_row,modes);  % ������������K_w�ĳ�ʼֵ

%������ʼ����������
K_TD(:,:,1) = [eye(A_row+hat_A_row);K_u_TD(:,:,1);K_w_TD(:,:,1)];
K_TD(:,:,2) = [eye(A_row+hat_A_row);K_u_TD(:,:,2);K_w_TD(:,:,2)];
K_TD(:,:,3) = [eye(A_row+hat_A_row);K_u_TD(:,:,3);K_w_TD(:,:,3)];

Gamma_TD = zeros(A_row+hat_A_row,A_row+hat_A_row,modes);
Upsilon_TD = zeros(A_row+hat_A_row+B_col+F_col+hat_F_col,A_row+hat_A_row+B_col+F_col+hat_F_col,modes);

delta_sigmma_P_TD = [0 0 0]; %���ڴ���sigmma_P�������ŵĲ�ֵ
delta_K_TD = [0 0 0]; %���ڴ���K_u�������ŵĲ�ֵ

% episodes_TD = 16; %���������Ļ��
% steps_TD = 100; %����ÿһĻ�Ĳ���
% mu = 0.98; %����ر�Ȩ��

episodes_TD = 51; %���������Ļ��
steps_TD = 200; %����ÿһĻ�Ĳ���
mu = 0.10; %����ر�Ȩ��

%% TPδ֪�����H���������
for mode = 1:modes
    Gamma_TD(:,:,mode) = sqrt(gamma)*(tilde_A(:,:,mode)+tilde_B(:,:,mode)*K_u_TD(:,:,mode)+tilde_F(:,:,mode)*K_w_TD(:,:,mode));  %��¼��һĻ����ģ̬��A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
    delta_sigmma_P_TD(1,mode,1) = log(1+(trace((sigmma_P_optimal(:,:,mode)-sigmma_P_TD(:,:,mode))'*(sigmma_P_optimal(:,:,mode)-sigmma_P_TD(:,:,mode))))/(trace(sigmma_P_optimal(:,:,mode)'*sigmma_P_optimal(:,:,mode))));
    delta_K_TD(1,mode,1) = log(1+(trace((K_optimal(:,:,mode)-K_TD(:,:,mode))'*(K_optimal(:,:,mode)-K_TD(:,:,mode))))/(trace(K_optimal(:,:,mode)'*K_optimal(:,:,mode))));
end

for episode_TD = 1:episodes_TD
    episode_TD
    
    lambda = 1/(episode_TD); %���²���
    
    sigmma_V_TD = sigmma_P_TD;
    for mode = 1:modes
        mode_now = mode; %����ÿһĻ�ĵ�һ��ģ̬
        TD_sum(:,:) = zeros(A_row+hat_A_row,A_row+hat_A_row);
        for step_test=1:steps_TD
            mode_old = mode_now; %���浱ǰģ̬
            mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %����ת�Ƹ��ʸ���ģ̬
            
            Gamma_TD(:,:,mode_now) = sqrt(gamma)*(tilde_A(:,:,mode_now)+tilde_B(:,:,mode_now)*K_u_TD(:,:,mode_now)+tilde_F(:,:,mode_now)*K_w_TD(:,:,mode_now));  %��¼��һĻ����ģ̬��A_tilde(:,:,m)+B_tilde(:,:,m)*K(:,:,m);
            Upsilon_TD(:,:,mode_now) = [tilde_C(:,:,mode_now)'*Q(:,:,mode_now)*tilde_C(:,:,mode_now) tilde_C(:,:,mode_now)'*Q(:,:,mode_now)*D(mode_now) tilde_C(:,:,mode_now)'*Q(:,:,mode_now)*tilde_H(:,:,mode_now); ...
                D(mode_now)'*Q(:,:,mode_now)*tilde_C(:,:,mode_now) D(mode_now)'*Q(:,:,mode_now)*D(mode_now)+R(mode_now) D(mode_now)'*Q(:,:,mode_now)*tilde_H(:,:,mode_now); ...
                tilde_H(:,:,mode_now)'*Q(:,:,mode_now)*tilde_C(:,:,mode_now) tilde_H(:,:,mode_now)'*Q(:,:,mode_now)*D(mode_now) tilde_H(:,:,mode_now)'*Q(:,:,mode_now)*tilde_H(:,:,mode_now)-theta^(2)*eye(F_col + hat_F_col)];
            K_TD(:,:,mode_now) = [eye(A_row+hat_A_row);K_u_TD(:,:,mode_now);K_w_TD(:,:,mode_now)];
            
            if step_test == 1
                Lambda(:,:,mode) = eye(A_row+hat_A_row); %����Lambda����
            else
                Lambda(:,:,mode) = Gamma_TD(:,:,mode_old)*Lambda(:,:,mode); %����Lambda����
            end
            
            TD = Lambda(:,:,mode)'*(Gamma_TD(:,:,mode_now)'*sigmma_V_TD(:,:,mode_now)*Gamma_TD(:,:,mode_now)+K_TD(:,:,mode_now)'*Upsilon_TD(:,:,mode_now)*K_TD(:,:,mode_now)-sigmma_V_TD(:,:,mode_old))*Lambda(:,:,mode);
            
            TD_sum = TD_sum + mu^(step_test-1)*TD; %���е�step����sigmma_P����
        end
        sigmma_P_TD(:,:,mode) = sigmma_V_TD(:,:,mode) + lambda*TD_sum; %���е�step����sigmma_P����
        K_u_TD(:,:,mode) = inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode)+R(mode)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))') ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))*inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col + hat_F_col))*(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode))');
        K_w_TD(:,:,mode) = inv(tilde_H(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_F(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode)-theta^(2)*eye(F_col+hat_F_col)-(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode) + gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))) ...
            *((D(mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))'*inv(D(mode)'*Q(:,:,mode)*D(mode)+gamma*tilde_B(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode)+R(mode))*(tilde_C(:,:,mode)'*Q(:,:,mode)*D(mode) + gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_B(:,:,mode))'-(tilde_C(:,:,mode)'*Q(:,:,mode)*tilde_H(:,:,mode)+gamma*tilde_A(:,:,mode)'*sigmma_P_TD(:,:,mode)*tilde_F(:,:,mode))');
    end
    for mode = 1:modes
        delta_sigmma_P_TD(1,mode,episode_TD+1) = log(1+(trace((sigmma_P_optimal(:,:,mode)-sigmma_P_TD(:,:,mode))'*(sigmma_P_optimal(:,:,mode)-sigmma_P_TD(:,:,mode))))/(trace(sigmma_P_optimal(:,:,mode)'*sigmma_P_optimal(:,:,mode))));
        delta_K_TD(1,mode,episode_TD+1) = log(1+(trace((K_optimal(:,:,mode)-K_TD(:,:,mode))'*(K_optimal(:,:,mode)-K_TD(:,:,mode))))/(trace(K_optimal(:,:,mode)'*K_optimal(:,:,mode))));
    end
end

for i = 1:modes
    P_TD(:,:,i) = tilde_C(:,:,i)'*Q(:,:,i)*tilde_C(:,:,i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_A(:,:,i)-[(tilde_C(:,:,i)'*Q(:,:,i)*D(i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_B(:,:,i)) tilde_C(:,:,i)'*Q(:,:,i)*tilde_H(:,:,i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i)]*inv([(gamma*tilde_B(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_B(:,:,i)+R(i)+D(i)'*Q(:,:,i)*D(i)) (D(i)'*Q(:,:,i)*tilde_H(:,:,i) + gamma*tilde_B(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i));(D(i)'*Q(:,:,i)*tilde_H(:,:,i) + gamma*tilde_B(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i))' (tilde_H(:,:,i)'*Q(:,:,i)*tilde_H(:,:,i)+gamma*tilde_F(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i)-theta^(2)*eye(F_col + hat_F_col))])*[(tilde_C(:,:,i)'*Q(:,:,i)*D(i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_B(:,:,i)) tilde_C(:,:,i)'*Q(:,:,i)*tilde_H(:,:,i)+gamma*tilde_A(:,:,i)'*sigmma_P_TD(:,:,i)*tilde_F(:,:,i)]';
end
for i = 1:modes
    delta_TD(:,:,i) = sigmma_P_TD(:,:,i)-Pr(i,1)*P_TD(:,:,1)-Pr(i,2)*P_TD(:,:,2)-Pr(i,3)*P_TD(:,:,3);
end
delta_TD_norm = trace(delta_TD(:,:,1)'*delta_TD(:,:,1)+delta_TD(:,:,2)'*delta_TD(:,:,2)+delta_TD(:,:,3)'*delta_TD(:,:,3))

delta_sigmma_P1 = delta_sigmma_P_TD(1,1,:);
delta_sigmma_P1 = delta_sigmma_P1(:);  %ͨ������������ά�����Ϊһά����
delta_sigmma_P2 = delta_sigmma_P_TD(1,2,:);
delta_sigmma_P2 = delta_sigmma_P2(:);
delta_sigmma_P3 = delta_sigmma_P_TD(1,3,:);
delta_sigmma_P3 = delta_sigmma_P3(:);

delta_K_TD_1 = delta_K_TD(1,1,:);
delta_K_TD_1 = delta_K_TD_1(:);  %ͨ������������ά�����Ϊһά����
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
    [0.691145833333333 0.212788259958072 0.188194444444442 0.205170875242091]); %legend�������������ͼ��
xlabel('$Iteration$','interpreter','latex')
ylabel('error of value function','interpreter','latex')
axis([0 episodes_TD-1 -0.05 0.8])
xticks([0:1:episodes_TD-1])
% yticks([0:0.1:0.8])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %��������������ΪTimes New Roman����СΪ26���߿�0.5

figure(5)
plot(0:1:(episodes_TD-1),delta_K_TD_1(1:episodes_TD),'--','Color','b',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_K_TD_2(1:episodes_TD),'-.','Color','r',"LineWidth",1.5)
hold on
plot(0:1:(episodes_TD-1),delta_K_TD_3(1:episodes_TD),'Color','g',"LineWidth",1.5)
legend('$\delta^{(l)}_{1}$','$\delta^{(l)}_{2}$','$\delta^{(l)}_{3}$','Interpreter','latex','Position', ...
    [0.691145833333333 0.212788259958072 0.188194444444442 0.205170875242091]); %legend�������������ͼ��
xlabel('$Iteration$','interpreter','latex')
ylabel('error of controller','interpreter','latex')
axis([0 episodes_TD-1 -0.0005 0.0040])
xticks([0:1:episodes_TD-1])
% yticks([0:0.00005:0.0004])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %��������������ΪTimes New Roman����СΪ26���߿�0.5

% %% ���п���������
% %% ϵͳģ��������
% episodes_test = 1; %������������
% steps_test = 150;
% hat_x = [10 -10 0]'; % ������ϵͳ����
% hat_z = [hat_C*hat_x];
% ww = zeros(1,steps_test);
% www = 0;
% 
% %% ϵͳģ������������LQR������
% % ����ϵͳ����
% x_LQR = [0 0 0]';
% u_LQR = zeros(1);
% z_LQR = [];
% % ����ϵͳ����
% tilde_x_LQR = [x_LQR' hat_x']';
% tilde_z_LQR = [0];
% % Monta Carlo ģ�����
% max_tilde_z_LQR = zeros(1,steps_test);
% min_tilde_z_LQR = zeros(1,steps_test);
% expect_tilde_z_LQR = zeros(1,steps_test);
% % ����˥��
% uu_LQR = zeros(1,steps_test);
% uuu_LQR = 0;
% tilde_zz_LQR = zeros(1,steps_test);
% tilde_zzz_LQR = 0;
% L2_gain_LQR = zeros(1,steps_test);
% 
% %% ϵͳģ������������H���������
% % ����ϵͳ����
% x_H = [0 0 0]';
% u_H = zeros(1);
% z_H = [];
% % ����ϵͳ����
% tilde_x_H = [x_H' hat_x']';
% tilde_z_H = [0];
% % Monta Carlo ģ�����
% max_tilde_z_H = zeros(1,steps_test);
% min_tilde_z_H = zeros(1,steps_test);
% expect_tilde_z_H = zeros(1,steps_test);
% % ����˥��
% uu_H = zeros(1,steps_test);
% uuu_H = 0;
% tilde_zz_H = zeros(1,steps_test);
% tilde_zzz_H = 0;
% L2_gain_H = zeros(1,steps_test);
% 
% %% ϵͳģ������������TD������
% % ����ϵͳ����
% x_TD = [0 0 0]';
% u_TD = zeros(1);
% z_TD = [];
% % ����ϵͳ����
% tilde_x_TD = [x_TD' hat_x']';
% tilde_z_TD = [0];
% % Monta Carlo ģ�����
% max_tilde_z_TD = zeros(1,steps_test);
% min_tilde_z_TD = zeros(1,steps_test);
% expect_tilde_z_TD = zeros(1,steps_test);
% % ����˥��
% uu_TD = zeros(1,steps_test);
% uuu_TD = 0;
% tilde_zz_TD = zeros(1,steps_test);
% tilde_zzz_TD = 0;
% L2_gain_TD = zeros(1,steps_test);
% 
% %% �������ϵͳ
% for episode_test = 1:episodes_test
%     mode_now = randsrc(1,1,[1 2 3;1/3 1/3 1/3]);
%     for step_test = 1:steps_test
%         tilde_A_now = tilde_A(:,:,mode_now); %�ɵ�ǰģ̬ȷ����������
%         tilde_B_now = tilde_B(:,:,mode_now);
%         tilde_F_now = tilde_F(:,:,mode_now);
%         tilde_C_now = tilde_C(:,:,mode_now);
%         tilde_H_now = tilde_H(:,:, mode_now);
%         C_now = C(:,:, mode_now);
%         D_now = D(mode_now);
%         H_now = H(:,:, mode_now);
%         
%         Q_now = Q(:,:, mode_now); %�ɵ�ǰģ̬ȷ��Ȩ�ؾ���
%         R_now = R(mode_now);
%         
%         K_u_LQR_now = K_u_LQR(:,:, mode_now); %�ɵ�ǰģ̬ȷ������������
%         K_u_H_now = K_u_H(:,:, mode_now);
%         K_u_TD_now = K_u_TD(:,:, mode_now);
%         K_w_now = K_w_H(:,:, mode_now);
%         
%         mode_old = mode_now; %���浱ǰģ̬
%         mode_now = randsrc(1,1,[1 2 3;Pr(mode_old,:)]); %����ת�Ƹ��ʸ���ģ̬
%         
%         u_LQR(:,step_test) = K_u_LQR_now*tilde_x_LQR(:,step_test); %���ڵ�ǰ��������������������
%         u_H(:,step_test) = K_u_H_now*tilde_x_H(:,step_test);
%         u_TD(:,step_test) = K_u_TD_now*tilde_x_TD(:,step_test);
%         
% %         w(:,step_test) = normrnd(0,0.05,[F_col+hat_F_col 1]); %���������˹������
%         w(:,step_test) = K_w_now*tilde_x_H(:,step_test);
%         %         w(:,step_test) = normrnd(0,0.1)*[sin(step_test)+cos(step_test)^(2);cos(step_test)+cos(step_test)^(2);sin(2*step_test)+cos(step_test)^(2);cos(step_test)+cos(step_test)^(3);sin(step_test)+cos(2*step_test)^(2);cos(step_test)+sin(step_test)^(2)];
%         w_LQR(:,step_test) = w(:,step_test);
%         w_H(:,step_test) = w(:,step_test);
%         w_TD(:,step_test) = w(:,step_test);
%         
%         tilde_x_LQR(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_LQR_now)*tilde_x_LQR(:,step_test) + tilde_F_now*w_LQR(:,step_test); %���ڵ�ǰģ̬�Լ���������������״̬
%         tilde_x_H(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_H_now)*tilde_x_H(:,step_test) + tilde_F_now*w_H(:,step_test);
%         tilde_x_TD(:,step_test+1) = (tilde_A_now + tilde_B_now*K_u_TD_now)*tilde_x_TD(:,step_test) + tilde_F_now*w_TD(:,step_test);
%         
%         tilde_z_LQR(:,step_test) = tilde_C_now*tilde_x_LQR(:,step_test) + D_now*u_LQR(:,step_test) + tilde_H_now*w_LQR(:,step_test); %���ڵ�ǰģ̬�Լ�������������״̬���¸������
%         tilde_z_H(:,step_test) = tilde_C_now*tilde_x_H(:,step_test) + D_now*u_H(:,step_test) + tilde_H_now*w_H(:,step_test);
%         tilde_z_TD(:,step_test) = tilde_C_now*tilde_x_TD(:,step_test) + D_now*u_TD(:,step_test) + tilde_H_now*w_TD(:,step_test);
%         
%         z_LQR(:,step_test) = C_now*tilde_x_LQR(1:A_row,step_test) + D_now*u_LQR(:,step_test) + H_now*w_LQR(1:F_col,step_test); %���ڵ�ǰģ̬�Լ�������������״̬�������
%         z_H(:,step_test) = C_now*tilde_x_H(1:A_row,step_test) + D_now*u_H(:,step_test) + H_now*w_H(1:F_col,step_test);
%         z_TD(:,step_test) = C_now*tilde_x_TD(1:A_row,step_test) + D_now*u_TD(:,step_test) + H_now*w_TD(1:F_col,step_test);
%         hat_z(:,step_test) = hat_C*tilde_x_LQR(A_row+1:A_row+hat_A_row,step_test);
%         
%         ww(:,step_test) = ww(:,step_test) + gamma^(step_test)*w(:,step_test)'*w(:,step_test); %�����ۼ�����
%         
%         uu_LQR(:,step_test) = uu_LQR(:,step_test) + gamma^(step_test)*u_LQR(:,step_test)'*R_now*u_LQR(:,step_test); %�����ۼƿ�������
%         uu_H(:,step_test) = uu_H(:,step_test) + gamma^(step_test)*u_H(:,step_test)'*R_now*u_H(:,step_test);
%         uu_TD(:,step_test) = uu_TD(:,step_test) + gamma^(step_test)*u_TD(:,step_test)'*R_now*u_TD(:,step_test);
%         
%         tilde_zz_LQR(:,step_test) = tilde_zz_LQR(:,step_test) + gamma^(step_test)*tilde_z_LQR(:,step_test)'*Q_now*tilde_z_LQR(:,step_test); %�����ۼƸ������
%         tilde_zz_H(:,step_test) = tilde_zz_H(:,step_test) + gamma^(step_test)*tilde_z_H(:,step_test)'*Q_now*tilde_z_H(:,step_test);
%         tilde_zz_TD(:,step_test) = tilde_zz_TD(:,step_test) + gamma^(step_test)*tilde_z_TD(:,step_test)'*Q_now*tilde_z_TD(:,step_test);
%         
%         if episode_test == episodes_test
%             www = www + ww(:,step_test); %�����ۼ�����
%             
%             uuu_LQR = uuu_LQR+uu_LQR(:,step_test); %�����ۼƿ�������
%             uuu_H = uuu_H+uu_H(:,step_test);
%             uuu_TD = uuu_TD+uu_TD(:,step_test);
%             
%             tilde_zzz_LQR = tilde_zzz_LQR+tilde_zz_LQR(:,step_test); %�����ۼƸ������
%             tilde_zzz_H = tilde_zzz_H+tilde_zz_H(:,step_test);
%             tilde_zzz_TD = tilde_zzz_TD+tilde_zz_TD(:,step_test);
%             
%             L2_gain_LQR(:,step_test) = (uuu_LQR+tilde_zzz_LQR)/www; %����L2����仯
%             L2_gain_H(:,step_test) = (uuu_H+tilde_zzz_H)/www;
%             L2_gain_TD(:,step_test) = (uuu_TD+tilde_zzz_TD)/www;
%         end
%         
%         episode_test
%         if episode_test == 1
%             max_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
%             min_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
%             max_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
%             min_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
%             max_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
%             min_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
%         else
%             if max_tilde_z_LQR(:,step_test) < tilde_z_LQR(:,step_test)
%                 max_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
%             end
%             if min_tilde_z_LQR(:,step_test) > tilde_z_LQR(:,step_test)
%                 min_tilde_z_LQR(:,step_test) = tilde_z_LQR(:,step_test);
%             end
%             if max_tilde_z_H(:,step_test) < tilde_z_H(:,step_test)
%                 max_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
%             end
%             if min_tilde_z_H(:,step_test) > tilde_z_H(:,step_test)
%                 min_tilde_z_H(:,step_test) = tilde_z_H(:,step_test);
%             end
%             if max_tilde_z_TD(:,step_test) < tilde_z_TD(:,step_test)
%                 max_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
%             end
%             if min_tilde_z_TD(:,step_test) > tilde_z_TD(:,step_test)
%                 min_tilde_z_TD(:,step_test) = tilde_z_TD(:,step_test);
%             end
%         end
%         expect_tilde_z_LQR(:,step_test) = expect_tilde_z_LQR(:,step_test) + (tilde_z_LQR(:,step_test) - expect_tilde_z_LQR(:,step_test))/episode_test;
%         expect_tilde_z_H(:,step_test) = expect_tilde_z_H(:,step_test) + (tilde_z_H(:,step_test) - expect_tilde_z_H(:,step_test))/episode_test;
%         expect_tilde_z_TD(:,step_test) = expect_tilde_z_TD(:,step_test) + (tilde_z_TD(:,step_test) - expect_tilde_z_TD(:,step_test))/episode_test;
%     end
% end

%% ֱ�ӵ���
% figure(6)
% fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_LQR(1,1:steps_test),fliplr(min_tilde_z_LQR(1,1:steps_test))],'b','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
% hold on
% fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_H(1,1:steps_test),fliplr(min_tilde_z_H(1,1:steps_test))],'r','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
% hold on
% plot(1:1:(steps_test),expect_tilde_z_LQR(1,1:steps_test),'Color','b',"LineWidth",0.5)
% hold on
% plot(1:1:(steps_test),expect_tilde_z_H(1,1:steps_test),'Color','r',"LineWidth",0.5)
% legend('$z_{LQR}$','$z_H$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
% xlabel('$Steps$','interpreter','latex')
% ylabel('$x_{1}$ and $\check{x}_{1}$','interpreter','latex')
% % axis([1 steps -2.5 2.5])
% % xticks([0:steps/10:steps])
% % yticks([-2.5:0.5:2.5])
% % set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5);

% figure(7)
% fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_H(1,1:steps_test),fliplr(min_tilde_z_H(1,1:steps_test))],'r','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
% hold on
% fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_TD(1,1:steps_test),fliplr(min_tilde_z_TD(1,1:steps_test))],'g','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
% hold on
% plot(1:1:(steps_test),expect_tilde_z_H(1,1:steps_test),'Color','r',"LineWidth",0.5)
% hold on
% plot(1:1:(steps_test),expect_tilde_z_TD(1,1:steps_test),'Color','g',"LineWidth",0.5)
% legend('$z_H$','$z_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
% xlabel('$Steps$','interpreter','latex')
% ylabel('$x_{1}$ and $\check{x}_{1}$','interpreter','latex')
% % axis([1 steps -2.5 2.5])
% % xticks([0:steps/10:steps])
% % yticks([-2.5:0.5:2.5])
% % set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5);

% figure(8)
% fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_LQR(1,1:steps_test),fliplr(min_tilde_z_LQR(1,1:steps_test))],'b','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
% hold on
% fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_H(1,1:steps_test),fliplr(min_tilde_z_H(1,1:steps_test))],'r','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
% hold on
% fill([1:1:(steps_test),fliplr(1:1:(steps_test))],[max_tilde_z_TD(1,1:steps_test),fliplr(min_tilde_z_TD(1,1:steps_test))],'g','edgealpha',0,'facealpha',0.05,'HandleVisibility','off')
% hold on
% plot(1:1:(steps_test),expect_tilde_z_LQR(1,1:steps_test),'--','Color','b',"LineWidth",1)
% hold on
% plot(1:1:(steps_test),expect_tilde_z_H(1,1:steps_test),'-.','Color','r',"LineWidth",1)
% hold on
% plot(1:1:(steps_test),expect_tilde_z_TD(1,1:steps_test),'Color','g',"LineWidth",1)
% legend('$\tilde{z}_{LQR}$','$\tilde{z}_H$','$\tilde{z}_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
% xlabel('$Steps$','interpreter','latex')
% ylabel('$\tilde{z}$ under differant controllers ','interpreter','latex')
% axis([1 steps_test -1 1])
% xticks([0:10:steps_test])
% yticks([-1:0.5:1])
% set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

% figure(9)
% plot(1:1:(steps_test),z_H(1,1:steps_test),'--','Color','r',"LineWidth",1)
% hold on
% plot(1:1:(steps_test),z_TD(1,1:steps_test),'-.','Color','g',"LineWidth",1)
% hold on
% plot(1:1:(steps_test),hat_z(1,1:steps_test),'Color','k',"LineWidth",1)
% legend('$z_H$','$z_{TD}$','$\hat{z}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
% xlabel('$Steps$','interpreter','latex')
% ylabel('$z$ and $\hat{z}$','interpreter','latex')
% axis([1 steps_test -10 6])
% xticks([0:10:steps_test])
% yticks([-10:2:6])
% set(gca,"FontName","Times New Roman","FontSize",36,"LineWidth",0.5);

% figure(10)
% plot(1:1:(steps_test),L2_gain_LQR(1,1:steps_test),'Color','b',"LineWidth",1)
% hold on
% plot(1:1:(steps_test),L2_gain_H(1,1:steps_test),'Color','r',"LineWidth",1)
% hold on
% plot(1:1:(steps_test),L2_gain_TD(1,1:steps_test),'Color','g',"LineWidth",1)
% legend('$L_{2}-gain_{LQR}$','$L_{2}-gain_{H}$','$L_{2}-gain_{TD}$','Interpreter','latex','Position',[0.756235532407406 0.445595077077051 0.125014467592594 0.147696327535109]);
% xlabel('$Steps$','interpreter','latex')
% ylabel('$L_{2}-gain$ under differant controller','interpreter','latex')
% axis([1 steps_test 0 0.8])
% xticks([0:steps_test/10:steps_test])
% yticks([0:10:100])
% set(gca,"FontName","Times New Roman","FontSize",32,"LineWidth",1);