%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�һ�׷��������Ĳ�ֿ������˲����㷨����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й��ƣ���ֵϵͳ������ֵ����
%
%   ��ע����������չ�������˲������㷨����
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 100;
LineWidth = 1.5;

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
X_esti = zeros(1,N);         %״̬���Ź���ֵ
P_esti = zeros(1,N);         %����������
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %���ɵ�λ��
%I(3,3) = 0;

%����
q = 10;                     %ϵͳ������ֵ
r = 20;                     %����������ֵ
Q = 0.3;                    %ϵͳ�����������
R = 0.3;                    %���������������
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

%��ʼֵ����
P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_esti(1,1) = P_pred_0;      %��ʼ���Ʒ�����
x_0  = 0;                    %��ʼ״̬
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
X_esti(1,1) = 0;             %״̬���Ƴ�ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

% %����GL΢�ֶ�����ϵ������
% gamma = cell(1,N);
% temp_matrx = zeros(1,1);
% for i = 1:1:N
%     temp_matrx(1,1) = bino_fir(1,i);
%     temp_matrx(2,2) = bino_sec(1,i);
%     temp_matrx(3,3) = bino_thi(1,i);
%     gamma{1,i} = temp_matrx;
% end

% %ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ-1
% U_input = ones(1,N);
% for i = N/2+1:1:N
%     U_input(:,i) = -1;
% end

%GL�����¶̼���ԭ��ĳ���
L = 101;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = 0;

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = 3*sin(2*X_real(1,k-1)) -X_real(1,k-1) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(1,k) = X_real(1,k) + V_noise(1,k);
  %�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1);
            %��������
            rema = 0;
            if k>L
                for i = 2:1:L+1
                   rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred
            %��������
            rema_P = 0;
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            end
            F = 6*cos(2*X_esti(1,k-1)) - 1;

        P_pred = (F-bino_fir(1,2))*P_esti(1,k-1)*(F-bino_fir(1,2))'+ Q + rema_P;

        %���㿨��������:Kk(2*1)
        H = 1;
        Kk = P_pred*H'/(H*P_pred*H' + R);

        %״̬����
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - X_pre -r);

        %���Ʒ���������
        P_esti(1,k) = (I-Kk*H)*P_pred;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;
figure;
plot(k,Z_meas,':r',k,X_real(1,:),'b','linewidth',LineWidth);
legend('�������','ʵ��״̬','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(1,:),'b',k,X_esti(1,:),':r','linewidth',LineWidth);
legend('ʵ��״̬','����״̬','Location','best');

%���ͼ
%������
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti;
mse = err_state.^2;
plot(k,mse(1,:),'b','linewidth',LineWidth);
legend('�������','Location','best');






