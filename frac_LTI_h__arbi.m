%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   ��ע������������ϵͳ�������˲���
%         ���h��Ϊ1��ʱ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 200;
%GL�����¶̼���ԭ��ĳ���
L = 6;

LineWidth = 1.2;

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
Q = [0,0; 0,0];                %ϵͳ�����������
R = 0;                         %���������������
W_noise = Q*randn(2,N);        %ϵͳ�������˴���Ϊ�˺����ķ���һ�£�ʵ��Ӧ��sqrt
V_noise = R*randn(1,N);        %�����������˴���Ϊ�˺����ķ���һ�£�ʵ��Ӧ��sqrt

%ϵͳ�״�
alpha_1 = 0.7;
alpha_2 = 1.2;

%���h��Ϊ1
h = 1;
H_alpha = [1/(h^(alpha_1)),0; 0,1/(h^(alpha_2))];

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);         %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha_1 = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha_1+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);         %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
alpha_2 = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha_2+1)/(i-1))*bino_sec(1,i-1);  
end

%����GL΢�ֶ�����ϵ������
gamma = cell(1,N);
temp_matrx = zeros(2,2);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    gamma{1,i} = temp_matrx;
end

%��ʼֵ����
P_pred_0 = [100,0; 0,100];         %��ʼԤ�ⷽ����
P_esti{:,1} = P_pred_0;            %��ʼ���Ʒ�����
x_0  = [0; 0];                     %��ʼ״̬
X_real(:,1) = x_0;                 %��ʵ״̬��ʼֵ
X_esti(:,1) = X_real(:,1);         %״̬���Ƴ�ʼֵ
Z_meas(:,1) = C*X_real(:,1);       %�������ݳ�ʼֵ

%ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ-1
U_input = ones(1,N);            
for i = 51:1:100
   U_input(:,i) = -1;
end
for i = 101:1:200
   U_input(:,i) = U_input(:,i-100);
end

diff_X_real = [0;0];
new_gamma = cell(1,N);
for i = 1:1:N
    new_gamma{1,i} = H_alpha*gamma{1,i};
end

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = A*X_real(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);

    rema = [0;0];
    for i = 2:1:k
        rema = rema + new_gamma{1,i}*X_real(:,k+1-i);
    end
    X_real(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(:,k) = C*X_real(:,k) + V_noise(1,k);
end

% for k=2:1:N
%     %�������˲�
%         %״̬Ԥ��:X_pre
%         diff_X_esti = A*X_esti(:,k-1) + B*U_input(k-1);
%             %��������
%             rema = [0;0];
%             if k>L
%                 for i = 2:1:L+1
%                     rema = rema + H_alpha*gamma{1,i}*X_esti(:,k+1-i);
%                 end
%             else
%                 for i = 2:1:k
%                     rema = rema + H_alpha*gamma{1,i}*X_esti(:,k+1-i);
%                 end
%             end
%         X_pre = diff_X_esti - rema;     %һ��״̬Ԥ��
%         %Ԥ�����Э�������:P_pred
%             %��������
%             rema_P = [0,0;0,0];
%             if k>L+1
%                 for i = 3:1:L+2
%                     rema_P = rema_P + H_alpha*gamma{1,i}*P_esti{1,k+1-i}*gamma{1,i}'*H_alpha';
%                 end
%             else
%                 for i = 3:1:k
%                     rema_P = rema_P + H_alpha*gamma{1,i}*P_esti{1,k+1-i}*gamma{1,i}'*H_alpha';
%                 end
%             end
%         P_pred = (A-H_alpha*gamma{1,2})*P_esti{1,k-1}*(A-gamma{1,2}*H_alpha)'+ Q + rema_P;
%         %���㿨��������:Kk(2*1)
%         Kk = P_pred*C'/(C*P_pred*C' + R);
%         %״̬����
%         X_esti(:,k) = X_pre + Kk*(Z_meas(k)-C*X_pre);
%         %���Ʒ���������
%         P_esti{1,k} = (I-Kk*C)*P_pred;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;
figure;
plot(k,U_input,'b',k,Z_meas,':r','linewidth',LineWidth);
legend('����','�������','Location','best');

% %״̬����ͼ
% figure;
% plot(k,X_real(1,:),'-r',k,X_esti(1,:),':b');
% legend('ʵ��״̬1','����״̬1','Location','best');
% 
% %״̬����ͼ
% figure;
% plot(k,X_real(2,:),'-r',k,X_esti(2,:),':b','linewidth',LineWidth);
% legend('ʵ��״̬2','����״̬2','Location','best');





