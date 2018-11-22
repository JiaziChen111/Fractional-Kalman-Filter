%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   备注：分数阶线性系统卡尔曼滤波器算例复现
%  
%   实际状态测试
%   基于短记忆原理：L = 3, 6, 50, 200 
%           
%   备注：此为论文中线性系统的复现（无噪声情况下的状态）
%         图三的复现
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N = 200;

LineWidth = 1.5;

%状态测量初始化
X_real_3 = zeros(2,N);            %真实状态
X_esti = zeros(2,N);            %状态最优估计值
P_esti = cell(1,N);             %估计误差方差阵
Z_meas_3 = zeros(1,N);            %实际观测值

%系统矩阵设置
A = [0,1; -0.1,-0.2];           %系统矩阵
B = [0; 1];                     %
C = [0.1,0.3];                  %
I = eye(2,2);                   %生成单位阵

%噪声
Q = [0,0; 0,0];             %系统噪声方差矩阵
R = 0;                        %测量噪声方差矩阵
W_noise = sqrt(Q)*randn(2,N);   %系统噪声
V_noise = sqrt(R)*randn(1,N);   %测量噪声

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N+1);          %微分阶次为0.7时GL定义下的系数
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N+1
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N+1);          %微分阶次为1.2时GL定义下的系数
alpha = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:N+1
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);
end
%计算GL微分定义下系数矩阵
gamma = cell(1,N+1);
temp_matrx = zeros(2,2);
for i = 1:1:N+1
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    gamma{1,i} = temp_matrx;
end

%系统输入设置:前一半输入为1，后一半输入为-1
U_input = ones(1,N);            
for i = 51:1:100
   U_input(:,i) = -1;
end
for i = 101:1:200
   U_input(:,i) = U_input(:,i-100);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    表示k时刻状态的微分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%求解不同短记忆原理下的状态：L = 3
L = 3;

%初始值设置
x_0  = [0; 0];                    %初始状态
X_real_3(:,1) = x_0;              %真实状态初始值
Z_meas_3(:,1) = V_noise(1,1);     %测量数据初始值

diff_X_real = [0;0];
for k=2:1:N
    %计算实际状态
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
    %实际观测值
    Z_meas_3(:,k) = C*X_real_3(:,k) + V_noise(1,k);
end



%求解不同短记忆原理下的状态：L = 6
L = 6;

%初始值设置
x_0  = [0; 0];                    %初始状态
X_real_6(:,1) = x_0;              %真实状态初始值
Z_meas_6(:,1) = V_noise(1,1);     %测量数据初始值

diff_X_real = [0;0];
for k=2:1:N
    %计算实际状态
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
    %实际观测值
    Z_meas_6(:,k) = C*X_real_6(:,k) + V_noise(1,k);
end

%求解不同短记忆原理下的状态：L = 50
L = 50;

%初始值设置
x_0  = [0; 0];                    %初始状态
X_real_50(:,1) = x_0;              %真实状态初始值
Z_meas_50(:,1) = V_noise(1,1);     %测量数据初始值

diff_X_real = [0;0];
for k=2:1:N
    %计算实际状态
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
    %实际观测值
    Z_meas_50(:,k) = C*X_real_50(:,k) + V_noise(1,k);
end

%求解不同短记忆原理下的状态：L = 10
L = 200;

%初始值设置
x_0  = [0; 0];                    %初始状态
X_real_200(:,1) = x_0;              %真实状态初始值
Z_meas_200(:,1) = V_noise(1,1);     %测量数据初始值

diff_X_real = [0;0];
for k=2:1:N
    %计算实际状态
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
    %实际观测值
    Z_meas_200(:,k) = C*X_real_200(:,k) + V_noise(1,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入与测量输出图
k = 1:1:N;
figure;
plot(k,Z_meas_3,':',k,Z_meas_6,'-.',k,Z_meas_50,'--',k,Z_meas_200,'-');
legend('L=3','L=6','L=50','L=200');

% %状态估计图
% figure;
% k = 1:1:N;
% plot(k,X_real_3(1,:),'b',k,X_real_3(2,:),':r','linewidth',LineWidth);
% legend('实际状态1','实际状态2','Location','best');







