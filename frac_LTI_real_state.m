%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   备注：分数阶线性系统卡尔曼滤波器算例复现
%
%   实际状态测试
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N = 100;
%GL定义下短记忆原理的长度
L = 50;
LineWidth = 1.5;

%状态测量初始化
X_real = zeros(2,N);            %真实状态
X_esti = zeros(2,N);            %状态最优估计值
P_esti = cell(1,N);             %估计误差方差阵
Z_meas = zeros(1,N);            %实际观测值

%系统矩阵设置
A = [0,1; -0.1,-0.2];           %系统矩阵
B = [0; 1];                     %
C = [0.1,0.3];                  %
I = eye(2,2);                   %生成单位阵

%噪声
Q = [0.3,0; 0,0.3];             %系统噪声方差矩阵
R = 0.3;                        %测量噪声方差矩阵
W_noise = Q*randn(2,N);   %系统噪声
V_noise = R*randn(1,N);   %测量噪声

%初始值设置
x_0  = [0; 0];                  %初始状态
X_real(:,1) = x_0;              %真实状态初始值
Z_meas(:,1) = C*X_real(:,1);   %测量数据初始值
%Z_meas(:,1) = V_noise(1,1);    %测量数据初始值

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
for i = N/2+1:1:N
   U_input(:,i) = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    表示k时刻状态的微分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = [0;0];

for k=2:1:N
    %计算实际状态
    diff_X_real = A*X_real(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);
    %
    rema = [0;0];
    if k>L
        for i = 2:1:L+1
            rema = rema + gamma{1,i}*X_real(:,k+1-i);
        end
    else
        for i = 2:1:k
            rema = rema + gamma{1,i}*X_real(:,k+1-i);
        end
    end
    X_real(:,k) = diff_X_real - rema;
    %实际观测值
    Z_meas(:,k) = C*X_real(:,k) + V_noise(1,k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入与测量输出图
k = 1:1:N;
figure;
plot(k,U_input,'b',k,Z_meas,':r','linewidth',LineWidth);
legend('输入','测量输出','Location','best');

%状态估计图
figure;
k = 1:1:N;
plot(k,X_real(1,:),'b',k,X_real(2,:),':r','linewidth',LineWidth);
legend('实际状态1','实际状态2','Location','best');







