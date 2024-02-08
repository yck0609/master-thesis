clc
clear all

%% ����ϵͳ����,�������� ��Mode-Independent H2-Control of a DC Motor Modeled as a Markov Jump Linear System��
%% ����ϵͳ����
modes = 3; %ϵͳģ̬����
A(:,:,1) = [-0.4799 5.1546 0;-3.8162 14.4723 0;0.1399 0 -0.9925]; %����ϵͳ״̬����״̬ת�ƾ���
A(:,:,2) = [-1.6026 9.1632 0;-0.5918 3.0317 0;0.0740 0 -0.4338];
A(:,:,3) = [0.6436 0.9178 0;-0.5056 2.4811 0;0.3865 0 0.0982];

B(:,:,1) = [5.8705 15.5010 0]'; %����ϵͳ״̬���̿��������������
B(:,:,2) = [10.2851 2.2282 0]';
B(:,:,3) = [0.7874 1.5302 0]';

F(:,:,1) = [0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0]; %����ϵͳ״̬���������������
F(:,:,2) = [0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0];
F(:,:,3) = [0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0;0.20 0 0 0 0 0 0 0 0];

C(:,:,1) = [1.0230 2.1100 0.9500]; %����ϵͳ��������������
C(:,:,2) = [0.9800 2.0500 1.1000];
C(:,:,3) = [1.0000 2.0000 1.0500];

D = [1.000 0.5000 -0.5000]; %����ϵͳ������̿��������������

H(:,:,1) = [0 0 0 0 0 0 0 0.10 0.15]; %����ϵͳ������������������
H(:,:,2) = [0 0 0 0 0 0 0 0.10 0.15]; 
H(:,:,3) = [0 0 0 0 0 0 0 0.10 0.15]; 

E(:,:,1) = [2 0 0;0 1.05 0]; %����ϵͳ���ⷽ���������
E(:,:,2) = [2 0 0;0 1.05 0];
E(:,:,3) = [2 0 0;0 1.05 0];

G(:,:,1) = [0 0 0 0.5 0 0 0 0 0;0 0 0 0 0.6 0 0 0 0]; %����ϵͳ���ⷽ�������������
G(:,:,2) = [0 0 0 0.5 0 0 0 0 0;0 0 0 0 0.6 0 0 0 0]; 
G(:,:,3) = [0 0 0 0.5 0 0 0 0 0;0 0 0 0 0.6 0 0 0 0]; 

Pr(:,:) = [0.95 0.05 0;0.36 0.6 0.04;0.1 0.1 0.8]; % ����ϵͳģ̬ת�Ƹ���

%% ����ϵͳ����
alpha = 0.15;
a = 0.4;
b = 0.5;
c = sqrt(1-a^2-b^2); %�����ϲ���������ά���Ҳ�
hat_A = [cos(alpha)+(1-cos(alpha))*a^2   (1-cos(alpha))*a*b-sin(alpha)*c  (1-cos(alpha))*a*c+sin(alpha)*b;
    (1-cos(alpha))*a*b+sin(alpha)*c  cos(alpha)+(1-cos(alpha))*b^2    (1-cos(alpha))*b*c-sin(alpha)*a;
    (1-cos(alpha))*a*c-sin(alpha)*b  (1-cos(alpha))*b*c+sin(alpha)*a  cos(alpha)+(1-cos(alpha))*c^2;]; %������ϵͳ״̬����״̬ת�ƾ���

hat_F = [0.2 0 0 0 0 0 0 0 0;0 0.2 0 0 0 0 0 0 0;0 0 0.2 0 0 0 0 0 0]; %������ϵͳ״̬���������������

hat_C = [1 1 1]; %������ϵͳ��������������

hat_E = [1.5 0 0;0 2 0]; %������ϵͳ���ⷽ���������

hat_G = [0 0 0 0.25 0 0 0 0 0;0 0 0 0 0.20 0 0 0 0];%������ϵͳ���ⷽ�������������

%% ��������ϵͳ����
[A_row,~] = size(A(:,:,1));  %size����Ĭ�Ϸ��ؾ���ά��,��ȡ����ϵͳ״̬��ά��
[~,B_col] = size(B(:,:,1));  %��ȡ����ϵͳ���������ά��
[~,F_col] = size(F(:,:,1));  %��ȡ����ϵͳ������ά��
[C_row,~] = size(C(:,:,1));  %��ȡ����ϵͳ�����ά��
[E_row,~] = size(E(:,:,1));  %��ȡ����ϵͳ�����ά��
[hat_A_row,~] = size(hat_A); %��ȡ������ϵͳ״̬��ά��
[~,hat_F_col] = size(hat_F); %��ȡ������ϵͳ��������ά��
[hat_E_row,~] = size(hat_E); %��ȡ������ϵͳ�����ά��

tilde_A = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %����ϵͳά�������������
tilde_B = zeros(A_row + hat_A_row,B_col,modes);
tilde_F = zeros(A_row + hat_A_row,F_col + hat_F_col,modes);
tilde_C = zeros(C_row,A_row + hat_A_row,modes);
tilde_H = zeros(C_row,F_col + hat_F_col,modes);
tilde_E = zeros(E_row + hat_E_row,A_row + hat_A_row,modes);
tilde_G = zeros(E_row + hat_E_row,F_col + hat_F_col,modes);

for mode=1:modes
    tilde_A(:,:,mode) = blkdiag(A(:,:,mode),hat_A);  %blkdiag�����������о�������׼�ԽǾ���
    tilde_B(:,:,mode) = [B(:,:,mode)' zeros(B_col,hat_A_row)]';
    tilde_F(:,:,mode) = blkdiag(F(:,:,mode),hat_F);
    tilde_C(:,:,mode) = [C(:,:,mode) -hat_C];
    tilde_H(:,:,mode) = [H(:,:,mode) zeros(C_row,hat_F_col)];
    tilde_E(:,:,mode) = blkdiag(E(:,:,mode),hat_E);
    tilde_G(:,:,mode) = blkdiag(G(:,:,mode),hat_G);
end

%% ��������ָ�������������Q��R��˥������gamma
Q(:,:,1) = 10*eye(C_row); %ϵͳ�������Ȩ�ؾ���
Q(:,:,2) = 10*eye(C_row);
Q(:,:,3) = 10*eye(C_row);

R = 0.5*[1 1 1];  %ϵͳ��������Ȩ�ؾ���

gamma = 0.99;  %˥������gamma
theta = 22.35; %H�������L2����
theta_filtering = 3.60; %H�����˲�L2����

%% ���H�����˲����˲���
P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); % ������ĳ�ʼֵ
sigmma_P = zeros(A_row + hat_A_row,A_row + hat_A_row,modes); %���尴���ʼ�Ȩ��;���
L_y(:,:,:) = zeros(A_row + hat_A_row,E_row + hat_E_row,modes);  % ������������K_u_LQR�ĳ�ʼֵ
norm_P = [];
episodes = 21; %����LQR����������������

for episode = 1:episodes
    for mode = 1:modes
        sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3); %�ɵ�ǰP����sigmma_P
        P(:,:,mode) = tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_A(:,:,mode)' + tilde_F(:,:,mode)*tilde_F(:,:,mode)' - (tilde_F(:,:,mode)*tilde_G(:,:,mode)' ...
            + tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)')*inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)') ...
            *(tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)')';
        P(:,:,mode) = (P(:,:,mode) + P(:,:,mode)')/2;
    end
    norm_P(episode) = trace(P(:,:,1)'*P(:,:,1) + P(:,:,2)'*P(:,:,2) + P(:,:,3)'*P(:,:,3));
end

for mode = 1:modes
    sigmma_P(:,:,mode) = Pr(mode,1)*P(:,:,1) + Pr(mode,2)*P(:,:,2) + Pr(mode,3)*P(:,:,3); %���㰴���ʼ�Ȩ��͵�P
    L_y(:,:,mode) = -(tilde_F(:,:,mode)*tilde_G(:,:,mode)' + tilde_A(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)' - tilde_A(:,:,mode)*sigmma_P(:,:,mode)*inv(sigmma_P(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P(:,:,mode)*tilde_E(:,:,mode)') ...
        *inv(tilde_G(:,:,mode)*tilde_G(:,:,mode)' + tilde_E(:,:,mode)*sigmma_P(:,:,mode)*tilde_E(:,:,mode)' - tilde_E(:,:,mode)*sigmma_P(:,:,mode)*inv(sigmma_P(:,:,mode) - theta_filtering^(2)*eye(A_row + hat_A_row))*sigmma_P(:,:,mode)*tilde_E(:,:,mode)');
end

figure(6)
plot(0:1:(episodes-1),norm_P)
