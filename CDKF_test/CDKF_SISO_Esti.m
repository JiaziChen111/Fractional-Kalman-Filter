%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�һ�׷��������Ĳ�ֿ������˲����㷨����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й��ƣ���ֵϵͳ������ֵ����
%
%   ��ע�����������Ĳ�ֿ������˲������㷨����
%         ԭϵͳ����:  D^{0.7} x_k = 3*sin(2*x_{k-1}) + w_k
%                             y_k = x_k + v_k
%         ��ԭϵͳ������ֵ���й��Ʋ�������ԭ����ԭϵͳ��������ϵͳ״̬��ʱ���һ
%         ֱ���󣩡�
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 1000;
LineWidth = 2;
h_con = sqrt(3);

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
X_esti = zeros(1,N);         %״̬���Ź���ֵ
P_esti = zeros(1,N);         %����������
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ
q_con = 1;                   %ϵͳ������ֵ
r_con = 0.5;                   %����������ֵ
Q_con = 1;                 %ϵͳ�����������
R_con = 0.49;                 %���������������

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %���ɵ�λ��
%I(3,3) = 0;

%����
q = q_con;                   %ϵͳ������ֵ
r = r_con;                   %����������ֵ
Q = Q_con;                   %ϵͳ�����������
R = R_con;                   %���������������
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

%Q_test = 1;

%��ʼֵ���ã���ʼ������Ϊ�㣩
q_esti = zeros(1,N);         %ϵͳ������ֵ����ֵ
r_esti = zeros(1,N);         %����������ֵ����ֵ
Q_esti = zeros(1,N);         %ϵͳ��������������ֵ
R_esti = zeros(1,N);         %������������������ֵ

R_esti(1,1) = 1;

P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_esti(1,1) = P_pred_0;      %��ʼ���Ʒ�����
x_0  = 0;                    %��ʼ״̬     
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
X_esti(1,1) = 0;             %״̬���Ƴ�ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

%GL�����¶̼���ԭ��ĳ���
L = 101;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = 0;

for k=2:1:N
    %����ʵ��״̬
    X_real(1,k) = 3*sin(2*X_real(1,k-1)) -X_real(1,k-1) + W_noise(1,k-1);

    %ʵ�ʹ۲�ֵ
    Z_meas(1,k) = X_real(1,k) + V_noise(1,k);
    %�������˲�
        %״̬Ԥ��:X_pre
        X_pre = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1) + q_esti(k-1);    %һ��״̬Ԥ��

        %Ԥ�����Э�������:P_pred
            %�����������cholsky�ֽ�
            S_chol = chol(P_esti(1,k-1))';

        %��ʱ���� temp_fun : ������ֵ,����Ϊ������
        temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol) ) - (X_esti(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti(1,k-1)-h_con*S_chol) ) + (X_esti(1,k-1)-h_con*S_chol);
        P_pred = 1/(4*h_con^2) * temp_fun^2 + Q;

        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r;

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = P_pred + R;

        %���㿨��������:Kk(2*1)
        Kk = P_pred/P_zpred;

        %r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;

        %״̬����
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %���Ʒ���������
        P_esti(1,k) = P_pred - Kk*P_zpred*Kk';

        q_esti(k) = ( (k-1)*q_esti(k-1) + X_esti(1,k) - X_pre  + q_esti(k-1) )/k;

        %R_esti(k) = ( (k-1)*R_esti(k-1) + (Z_meas(1,k) - Z_esti)^2 - P_pred)/k; %

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;
figure;
plot(k,Z_meas,':r',k,X_real(1,:),'b','linewidth',LineWidth);
legend('�������','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(1,:),'b',k,X_esti(1,:),':r','linewidth',LineWidth);
legend('ʵ��״̬1','����״̬1','Location','best');

%���ͼ
%������
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti;
mse = err_state.^2;
plot(k,mse(1,:),'b','linewidth',LineWidth);
legend('�������1','Location','best');
% plot(k,err_state(1,:),'b','linewidth',LineWidth);
% legend('�������1','Location','best');

figure;
plot(k,q_esti(1,:),'b','linewidth',LineWidth);
legend('ϵͳ������������������','Location','best');
line([0,N],[R_con,R_con],'lineStyle','--')




