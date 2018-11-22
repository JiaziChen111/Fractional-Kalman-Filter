%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   ��ע������������ϵͳ�������˲����������� ͼ��
%         Ӧ���ǳɹ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 100;
LineWidth = 1.5;

%״̬������ʼ��
X_real = zeros(2,N);            %��ʵ״̬
X_esti = zeros(2,N);            %״̬���Ź���ֵ
P_esti = cell(1,N);             %����������
Z_meas = zeros(1,N);            %ʵ�ʹ۲�ֵ

%ϵͳ��������
A = [0,1; -0.1,-0.2];           %ϵͳ����
B = [0; 1];                     %
C = [0.1,0.3];                  %
I = eye(2,2);                   %���ɵ�λ��

%����
Q = [0.3,0; 0,0.3];             %ϵͳ�����������
R = 0.3;                        %���������������
W_noise = Q*randn(2,N);         %ϵͳ����
V_noise = R*randn(1,N);         %��������

%��ʼֵ����
P_pred_0 = [100,0; 0,100];      %��ʼԤ�ⷽ����
P_esti{:,1} = P_pred_0;         %��ʼ���Ʒ�����
x_0  = [0; 0];                  %��ʼ״̬
X_real(:,1) = x_0;              %��ʵ״̬��ʼֵ
X_esti(:,1) = X_real(:,1);      %״̬���Ƴ�ʼֵ
Z_meas(:,1) = C*X_real(:,1);    %�������ݳ�ʼֵ

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);          %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
alpha = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);  
end

%����GL΢�ֶ�����ϵ������
gamma = cell(1,N);
temp_matrx = zeros(2,2);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    gamma{1,i} = temp_matrx;
end

%ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ-1
U_input = ones(1,N);            
for i = N/2+1:1:N
   U_input(:,i) = -1;
end

%GL�����¶̼���ԭ��ĳ���
L = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = [0;0];

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = A*X_real(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);

    rema = [0;0];
    for i = 2:1:k
        rema = rema + gamma{1,i}*X_real(:,k+1-i);
    end
    X_real(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(:,k) = C*X_real(:,k) + V_noise(1,k);

    %�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = A*X_esti(:,k-1) + B*U_input(k-1);
            %��������
            rema = [0;0];
            if k>L
                for i = 2:1:L+1
                   rema = rema + gamma{1,i}*X_esti(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + gamma{1,i}*X_esti(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred
            %��������
            rema_P = [0,0;0,0];
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_esti{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + gamma{1,i}*P_esti{1,k+1-i}*gamma{1,i}';
                end
            end
        P_pred = (A-gamma{1,2})*P_esti{1,k-1}*(A-gamma{1,2})'+ Q + rema_P;
        %���㿨��������:Kk(2*1)
        Kk = P_pred*C'/(C*P_pred*C' + R);
        %״̬����
        X_esti(:,k) = X_pre + Kk*(Z_meas(k)-C*X_pre);
        %���Ʒ���������
        P_esti{1,k} = (I-Kk*C)*P_pred;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;
figure;
plot(k,U_input,'b',k,Z_meas,':r','linewidth',LineWidth);
legend('����','�������','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(1,:),'b',k,X_esti(1,:),':r','linewidth',LineWidth);
legend('ʵ��״̬1','����״̬1','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(2,:),'b',k,X_esti(2,:),':r','linewidth',LineWidth);
legend('ʵ��״̬2','����״̬2','Location','best');

%���ͼ
%������
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti;
plot(k,err_state(1,:),'b',k,err_state(2,:),':r','linewidth',LineWidth);
legend('�������1','�������2','Location','best');






