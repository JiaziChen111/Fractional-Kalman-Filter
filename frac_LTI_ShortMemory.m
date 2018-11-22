%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   备注：分数阶线性系统卡尔曼滤波器短记忆原理验证
%         图四的复现
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N = 200;

LineWidth = 1.2;

%状态测量初始化
X_real = zeros(2,N);              %真实状态
X_esti_6 = zeros(2,N);            %状态最优估计值
P_esti_6 = cell(1,N);             %估计误差方差阵
X_esti_50 = zeros(2,N);           %状态最优估计值
P_esti_50 = cell(1,N);            %估计误差方差阵
Z_meas = zeros(1,N);              %实际观测值

%系统矩阵设置
A = [0,1; -0.1,-0.2];           %系统矩阵
B = [0; 1];                     %
C = [0.1,0.3];                  %
I = eye(2,2);                   %生成单位阵

%噪声
Q = [0,0; 0,0];             %系统噪声方差矩阵
R = 0;                        %测量噪声方差矩阵
W_noise = Q*randn(2,N);         %系统噪声，此处是为了和论文仿真一致，实际应加sqrt
V_noise = R*randn(1,N);         %测量噪声，此处是为了和论文仿真一致，实际应加sqrt

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

%初始值设置
P_pred_0 = [100,0; 0,100];         %初始预测方差阵
P_esti_6{:,1} = P_pred_0;          %初始估计方差阵
P_esti_50{:,1} = P_pred_0;         %初始估计方差阵
x_0  = [0; 0];                     %初始状态
X_real(:,1) = x_0;                 %真实状态初始值
X_esti_6(:,1) = X_real(:,1);       %状态估计初始值
X_esti_50(:,1) = X_real(:,1);      %状态估计初始值
Z_meas(:,1) = C*X_real(:,1);       %测量数据初始值

%系统输入设置:前一半输入为1，后一半输入为-1
U_input = ones(1,N);            
for i = 51:1:100
   U_input(:,i) = -1;
end
for i = 101:1:200
   U_input(:,i) = U_input(:,i-100);
end

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
end

% 短记忆原理 L = 6
%GL定义下短记忆原理的长度
L = 6;
for k=2:1:N
    %卡尔曼滤波
        %状态预测:X_pre
        diff_X_esti = A*X_esti_6(:,k-1) + B*U_input(k-1);
            %计算余项
            rema = [0;0];
            if k>L
                for i = 2:1:L+1
                   rema = rema + gamma{1,i}*X_esti_6(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + gamma{1,i}*X_esti_6(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema;     %一步状态预测
        %预测误差协方差矩阵:P_pred
            %计算余项
            rema_P = [0,0;0,0];
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_esti_6{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + gamma{1,i}*P_esti_6{1,k+1-i}*gamma{1,i}';
                end
            end
        P_pred = (A-gamma{1,2})*P_esti_6{1,k-1}*(A-gamma{1,2})'+ Q + rema_P;
        %计算卡尔曼增益:Kk(2*1)
        Kk = P_pred*C'/(C*P_pred*C' + R);
        %状态更新
        X_esti_6(:,k) = X_pre + Kk*(Z_meas(k)-C*X_pre);
        %估计方差矩阵更新
        P_esti_6{1,k} = (I-Kk*C)*P_pred;
end

% 短记忆原理 L = 50
L = 50;
for k=2:1:N
    %卡尔曼滤波
        %状态预测:X_pre
        diff_X_esti = A*X_esti_50(:,k-1) + B*U_input(k-1);
            %计算余项
            rema = [0;0];
            if k>L
                for i = 2:1:L+1
                   rema = rema + gamma{1,i}*X_esti_50(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + gamma{1,i}*X_esti_50(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema;     %一步状态预测
        %预测误差协方差矩阵:P_pred
            %计算余项
            rema_P = [0,0;0,0];
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_esti_50{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + gamma{1,i}*P_esti_50{1,k+1-i}*gamma{1,i}';
                end
            end
        P_pred = (A-gamma{1,2})*P_esti_50{1,k-1}*(A-gamma{1,2})'+ Q + rema_P;
        %计算卡尔曼增益:Kk(2*1)
        Kk = P_pred*C'/(C*P_pred*C' + R);
        %状态更新
        X_esti_50(:,k) = X_pre + Kk*(Z_meas(k)-C*X_pre);
        %估计方差矩阵更新
        P_esti_50{1,k} = (I-Kk*C)*P_pred;
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
plot(k,X_real(1,:),'-r',k,X_esti_6(1,:),':b',k,X_esti_50(1,:),'-.');
legend('实际状态1','L6状态1','L50状态1','Location','best');

%状态估计图
figure;
plot(k,X_real(2,:),'-r',k,X_esti_6(2,:),':b',k,X_esti_50(2,:),'-.','linewidth',LineWidth);
legend('实际状态2','L6状态2','L50状态2','Location','best');





