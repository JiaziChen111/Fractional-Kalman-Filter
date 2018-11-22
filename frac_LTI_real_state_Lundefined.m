%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   ��ע������������ϵͳ�������˲�����������
%  
%   ʵ��״̬����
%   ���ڶ̼���ԭ��L = 3, 6, 50, 200 
%           
%   ��ע����Ϊ����������ϵͳ�ĸ��֣�����������µ�״̬��
%         ͼ���ĸ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 200;

LineWidth = 1.5;

%״̬������ʼ��
X_real_3 = zeros(2,N);            %��ʵ״̬
X_esti = zeros(2,N);            %״̬���Ź���ֵ
P_esti = cell(1,N);             %����������
Z_meas_3 = zeros(1,N);            %ʵ�ʹ۲�ֵ

%ϵͳ��������
A = [0,1; -0.1,-0.2];           %ϵͳ����
B = [0; 1];                     %
C = [0.1,0.3];                  %
I = eye(2,2);                   %���ɵ�λ��

%����
Q = [0,0; 0,0];             %ϵͳ�����������
R = 0;                        %���������������
W_noise = sqrt(Q)*randn(2,N);   %ϵͳ����
V_noise = sqrt(R)*randn(1,N);   %��������

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N+1);          %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N+1
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N+1);          %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
alpha = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:N+1
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);
end
%����GL΢�ֶ�����ϵ������
gamma = cell(1,N+1);
temp_matrx = zeros(2,2);
for i = 1:1:N+1
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    gamma{1,i} = temp_matrx;
end

%ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ-1
U_input = ones(1,N);            
for i = 51:1:100
   U_input(:,i) = -1;
end
for i = 101:1:200
   U_input(:,i) = U_input(:,i-100);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��ⲻͬ�̼���ԭ���µ�״̬��L = 3
L = 3;

%��ʼֵ����
x_0  = [0; 0];                    %��ʼ״̬
X_real_3(:,1) = x_0;              %��ʵ״̬��ʼֵ
Z_meas_3(:,1) = V_noise(1,1);     %�������ݳ�ʼֵ

diff_X_real = [0;0];
for k=2:1:N
    %����ʵ��״̬
    diff_X_real = A*X_real_3(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);

    rema = [0;0];
    if k>L
        for i = 2:1:L+1
            rema = rema + gamma{1,i}*X_real_3(:,k+1-i);
        end
    else
        for i = 2:1:k
            rema = rema + gamma{1,i}*X_real_3(:,k+1-i);
        end
    end
    X_real_3(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas_3(:,k) = C*X_real_3(:,k) + V_noise(1,k);
end



%��ⲻͬ�̼���ԭ���µ�״̬��L = 6
L = 6;

%��ʼֵ����
x_0  = [0; 0];                    %��ʼ״̬
X_real_6(:,1) = x_0;              %��ʵ״̬��ʼֵ
Z_meas_6(:,1) = V_noise(1,1);     %�������ݳ�ʼֵ

diff_X_real = [0;0];
for k=2:1:N
    %����ʵ��״̬
    diff_X_real = A*X_real_6(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);

    rema = [0;0];
    if k>L
        for i = 2:1:L+1
            rema = rema + gamma{1,i}*X_real_6(:,k+1-i);
        end
    else
        for i = 2:1:k
            rema = rema + gamma{1,i}*X_real_6(:,k+1-i);
        end
    end
    X_real_6(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas_6(:,k) = C*X_real_6(:,k) + V_noise(1,k);
end

%��ⲻͬ�̼���ԭ���µ�״̬��L = 50
L = 50;

%��ʼֵ����
x_0  = [0; 0];                    %��ʼ״̬
X_real_50(:,1) = x_0;              %��ʵ״̬��ʼֵ
Z_meas_50(:,1) = V_noise(1,1);     %�������ݳ�ʼֵ

diff_X_real = [0;0];
for k=2:1:N
    %����ʵ��״̬
    diff_X_real = A*X_real_50(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);

    rema = [0;0];
    if k>L
        for i = 2:1:L+1
            rema = rema + gamma{1,i}*X_real_50(:,k+1-i);
        end
    else
        for i = 2:1:k
            rema = rema + gamma{1,i}*X_real_50(:,k+1-i);
        end
    end
    X_real_50(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas_50(:,k) = C*X_real_50(:,k) + V_noise(1,k);
end

%��ⲻͬ�̼���ԭ���µ�״̬��L = 10
L = 200;

%��ʼֵ����
x_0  = [0; 0];                    %��ʼ״̬
X_real_200(:,1) = x_0;              %��ʵ״̬��ʼֵ
Z_meas_200(:,1) = V_noise(1,1);     %�������ݳ�ʼֵ

diff_X_real = [0;0];
for k=2:1:N
    %����ʵ��״̬
    diff_X_real = A*X_real_200(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);

    rema = [0;0];
    if k>L
        for i = 2:1:L+1
            rema = rema + gamma{1,i}*X_real_200(:,k+1-i);
        end
    else
        for i = 2:1:k
            rema = rema + gamma{1,i}*X_real_200(:,k+1-i);
        end
    end
    X_real_200(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas_200(:,k) = C*X_real_200(:,k) + V_noise(1,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;
figure;
plot(k,Z_meas_3,':',k,Z_meas_6,'-.',k,Z_meas_50,'--',k,Z_meas_200,'-');
legend('L=3','L=6','L=50','L=200');

% %״̬����ͼ
% figure;
% k = 1:1:N;
% plot(k,X_real_3(1,:),'b',k,X_real_3(2,:),':r','linewidth',LineWidth);
% legend('ʵ��״̬1','ʵ��״̬2','Location','best');







