%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：     fractional order CDKF
%   目的：一阶分数阶中心差分卡尔曼滤波器算法测试
%         对系统噪声均值进行估计
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   结果：较好的对状态进行估计，常值系统噪声均值收敛
%
%   备注：分数阶中心差分卡尔曼滤波器的算法测试
%         原系统函数:  D^{0.7} x_k = 3*sin(2*x_{k-1}) + w_k
%                             y_k = x_k + v_k
%         对原系统噪声均值进行估计不收敛的原因是原系统不收敛（系统状态随时间而一
%         直增大）。
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N = 1000;
LineWidth = 2;
h_con = sqrt(3);

%状态测量初始化
X_real = zeros(1,N);         %真实状态
X_esti = zeros(1,N);         %状态最优估计值
P_esti = zeros(1,N);         %估计误差方差阵
Z_meas = zeros(1,N);         %实际观测值
q_con = 1;                   %系统噪声均值
r_con = 0.5;                   %测量噪声均值
Q_con = 1;                 %系统噪声方差矩阵
R_con = 0.49;                 %测量噪声方差矩阵

%系统矩阵设置
% A = [0,1; -0.1,-0.2];      %系统矩阵
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %生成单位阵
%I(3,3) = 0;

%噪声
q = q_con;                   %系统噪声均值
r = r_con;                   %测量噪声均值
Q = Q_con;                   %系统噪声方差矩阵
R = R_con;                   %测量噪声方差矩阵
W_noise = sqrt(Q)*randn(1,N) + q;  %系统噪声
V_noise = sqrt(R)*randn(1,N) + r;  %测量噪声

%Q_test = 1;

%初始值设置（初始矩阵不能为零）
q_esti = zeros(1,N);         %系统噪声均值估计值
r_esti = zeros(1,N);         %测量噪声均值估计值
Q_esti = zeros(1,N);         %系统噪声方差矩阵估计值
R_esti = zeros(1,N);         %测量噪声方差矩阵估计值

R_esti(1,1) = 1;

P_pred_0 = 100;              %初始预测方差阵
P_esti(1,1) = P_pred_0;      %初始估计方差阵
x_0  = 0;                    %初始状态     
X_real(1,1) = x_0;           %真实状态初始值
X_esti(1,1) = 0;             %状态估计初始值
Z_meas(1,1) = V_noise(1,1);  %测量数据初始值

%GL定义下短记忆原理的长度
L = 101;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    表示k时刻状态的微分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = 0;

for k=2:1:N
    %计算实际状态
    X_real(1,k) = 3*sin(2*X_real(1,k-1)) -X_real(1,k-1) + W_noise(1,k-1);

    %实际观测值
    Z_meas(1,k) = X_real(1,k) + V_noise(1,k);
    %卡尔曼滤波
        %状态预测:X_pre
        X_pre = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1) + q_esti(k-1);    %一步状态预测

        %预测误差协方差矩阵:P_pred
            %对误差矩阵进行cholsky分解
            S_chol = chol(P_esti(1,k-1))';

        %临时变量 temp_fun : 函数差值,函数为单变量
        temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol) ) - (X_esti(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti(1,k-1)-h_con*S_chol) ) + (X_esti(1,k-1)-h_con*S_chol);
        P_pred = 1/(4*h_con^2) * temp_fun^2 + Q;

        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r;

        %测量预测误差协方差矩阵:P_zpred ---- P_z_k|k-1
        P_zpred = P_pred + R;

        %计算卡尔曼增益:Kk(2*1)
        Kk = P_pred/P_zpred;

        %r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;

        %状态更新
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %估计方差矩阵更新
        P_esti(1,k) = P_pred - Kk*P_zpred*Kk';

        q_esti(k) = ( (k-1)*q_esti(k-1) + X_esti(1,k) - X_pre  + q_esti(k-1) )/k;

        %R_esti(k) = ( (k-1)*R_esti(k-1) + (Z_meas(1,k) - Z_esti)^2 - P_pred)/k; %

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入与测量输出图
k = 1:1:N;
figure;
plot(k,Z_meas,':r',k,X_real(1,:),'b','linewidth',LineWidth);
legend('测量输出','Location','best');

%状态估计图
figure;
plot(k,X_real(1,:),'b',k,X_esti(1,:),':r','linewidth',LineWidth);
legend('实际状态1','估计状态1','Location','best');

%误差图
%误差计算
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti;
mse = err_state.^2;
plot(k,mse(1,:),'b','linewidth',LineWidth);
legend('均方误差1','Location','best');
% plot(k,err_state(1,:),'b','linewidth',LineWidth);
% legend('估计误差1','Location','best');

figure;
plot(k,q_esti(1,:),'b','linewidth',LineWidth);
legend('系统噪声方差收敛趋势线','Location','best');
line([0,N],[R_con,R_con],'lineStyle','--')




