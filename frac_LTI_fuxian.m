%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   备注：分数阶线性系统卡尔曼滤波器算例复现 图二
%         应该是成功了
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N = 100;
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
W_noise = Q*randn(2,N);         %系统噪声
V_noise = R*randn(1,N);         %测量噪声

%初始值设置
P_pred_0 = [100,0; 0,100];      %初始预测方差阵
P_esti{:,1} = P_pred_0;         %初始估计方差阵
x_0  = [0; 0];                  %初始状态
X_real(:,1) = x_0;              %真实状态初始值
X_esti(:,1) = X_real(:,1);      %状态估计初始值
Z_meas(:,1) = C*X_real(:,1);    %测量数据初始值

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);          %微分阶次为0.7时GL定义下的系数
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %微分阶次为1.2时GL定义下的系数
alpha = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);  
end

%计算GL微分定义下系数矩阵
gamma = cell(1,N);
temp_matrx = zeros(2,2);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    gamma{1,i} = temp_matrx;
end

%系统输入设置:前一半输入为1，后一半输入为-1
U_input = ones(1,N);            
for i = N/2+1:1:N
   U_input(:,i) = -1;
end

%GL定义下短记忆原理的长度
L = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    表示k时刻状态的微分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = [0;0];

for k=2:1:N
    %计算实际状态
    diff_X_real = A*X_real(:,k-1) + B*U_input(k-1) + W_noise(:,k-1);

    rema = [0;0];
    for i = 2:1:k
        rema = rema + gamma{1,i}*X_real(:,k+1-i);
    end
    X_real(:,k) = diff_X_real - rema;
    %实际观测值
    Z_meas(:,k) = C*X_real(:,k) + V_noise(1,k);

    %卡尔曼滤波
        %状态预测:X_pre
        diff_X_esti = A*X_esti(:,k-1) + B*U_input(k-1);
            %计算余项
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
        X_pre = diff_X_esti - rema;     %一步状态预测
        %预测误差协方差矩阵:P_pred
            %计算余项
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
        %计算卡尔曼增益:Kk(2*1)
        Kk = P_pred*C'/(C*P_pred*C' + R);
        %状态更新
        X_esti(:,k) = X_pre + Kk*(Z_meas(k)-C*X_pre);
        %估计方差矩阵更新
        P_esti{1,k} = (I-Kk*C)*P_pred;
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
plot(k,X_real(1,:),'b',k,X_esti(1,:),':r','linewidth',LineWidth);
legend('实际状态1','估计状态1','Location','best');

%状态估计图
figure;
plot(k,X_real(2,:),'b',k,X_esti(2,:),':r','linewidth',LineWidth);
legend('实际状态2','估计状态2','Location','best');

%误差图
%误差计算
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti;
plot(k,err_state(1,:),'b',k,err_state(2,:),':r','linewidth',LineWidth);
legend('估计误差1','估计误差2','Location','best');






