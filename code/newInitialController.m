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

episodes_H = 41; %������������

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
axis([0 episodes_H-1 0 200]) %���������᷶Χ axis([x_min x_max y_min y_max])
xlabel('$Iteration$','interpreter','latex')
xticks([0:(episodes_H-1)/10:episodes_H-1]) %���� x ��̶�ֵ
ylabel('$log(\left\|P_{i}\right\|_{2})$','interpreter','latex')
yticks([0:20:200])
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
yticks([0:10:100])
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

episodes_TD = 51; %���������Ļ��
steps_TD = 100; %����ÿһĻ�Ĳ���
mu = 0.09; %����ر�Ȩ��

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
xticks([0:10:episodes_TD-1])
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
axis([0 episodes_TD-1 -0.0005 0.040])
xticks([0:10:episodes_TD-1])
% yticks([0:0.00005:0.0004])
set(gca,"FontName","Times New Roman","FontSize",42,"LineWidth",0.5); %��������������ΪTimes New Roman����СΪ26���߿�0.5