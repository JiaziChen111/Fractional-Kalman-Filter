%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%
%   论文：fractional order CDKF
%
%   目的：一阶分数阶中心差分卡尔曼滤波器算法测试-
%
%   函数实验: MIMO system （credit） 
%
%   结果：较好的对状态进行估计，常值系统噪声均值收敛
%
%   备注：分数阶中心差分卡尔曼滤波器的算法测试
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N = 50;
LineWidth = 2;
h_con = sqrt(3);

%状态测量初始化
X_real = zeros(2,N);         %真实状态
X_esti = zeros(2,N);         %状态最优估计值
%P_esti = zeros(1,N);        %估计误差方差阵
P_esti = cell(1,N);          %估计误差方差阵
Z_meas = zeros(1,N);         %实际观测值

%系统矩阵设置
% A = [0,1; -0.1,-0.2];      %系统矩阵
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(2,2);                %生成单位阵
%I(3,3) = 0;
 
%噪声
q = [1,1]';                 %系统噪声均值
r = 1;                    %测量噪声均值
Q =  0.09*eye(2,2);         %系统噪声方差矩阵
R =  0.09;                  %测量噪声方差矩阵
W_noise = sqrt(Q)*randn(2,N) + q;   %系统噪声
V_noise = sqrt(R)*randn(1,N) + r;   %测量噪声

%Q_test = 100;

%初始值设置（初始矩阵不能为零）
bino = cell(1,N);            %二项式系数矩阵
r_esti = zeros(2,N);         %测量噪声均值估计值
P_pred_0 = 100*eye(2,2);     %初始预测方差阵
P_esti{1,1} = P_pred_0;      %初始估计方差阵
x_0  = [0,0]';               %初始状态     
X_real(:,1) = x_0;           %真实状态初始值
X_esti(:,1) = [0,0]';        %状态估计初始值
Z_meas(1,1) = V_noise(1,1);  %测量数据初始值

%计算alpha阶次对应的GL定义系数 binomial coefficient
bino_fir = zeros(1,N);       %微分阶次为0.03时GL定义下的系数
alpha = 0.5;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);
end

bino_sec = zeros(1,N);       %微分阶次为1.2时GL定义下的系数
alpha = 0.7;
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);
end

temp_matrx = zeros(2,2);
for k = 1:1:N
    temp_matrx(1,1) = bino_fir(1,k);
    temp_matrx(2,2) = bino_sec(1,k);
    bino{1,k} = temp_matrx;
end
% %计算GL微分定义下系数矩阵
% gamma = cell(1,N);
% temp_matrx = zeros(1,1);
% for i = 1:1:N 
%     temp_matrx(1,1) = bino_fir(1,i);
%     temp_matrx(2,2) = bino_sec(1,i);
%     temp_matrx(3,3) = bino_thi(1,i);
%     gamma{1,i} = temp_matrx;
% end

%系统输入设置:前一半输入为0.04，后一半输入为0.04          
% for i = N/2+1:1:N
%     U_input(:,i) = -1;
% end

%GL定义下短记忆原理的长度
L = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% diff_X_real   表示k时刻状态的微分 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real= [0,0]';

r_esti_A = zeros(1,N);            %测量噪声均值估计值
r_esti_A(1,1) = 0.5; 

for k=2:1:N
    %计算实际状态
    diff_X_real(1,1) = 3*cos(X_real(2,k-1)) - X_real(1,k-1) + W_noise(1,k-1);
    diff_X_real(2,1) = 2*cos(X_real(1,k-1)) - X_real(2,k-1) + W_noise(2,k-1);
    
    %diff_X_real = 3*sin(2*X_real(1,k-1)) -X_real(1,k-1) + W_noise(1,k-1);
    %rema_fir = 0; rema_sec = 0;
    rema = [0,0]';
    for i = 2:1:k
        %rema_fir = rema_fir + bino_fir(1,i)*X_real(1,k+1-i);
        %rema_sec = rema_sec + bino_sec(1,i)*X_real(1,k+1-i);
        rema = rema + bino{1,i}*X_real(:,k+1-i);
    end
    X_real(:,k) = diff_X_real - rema;

    %实际观测值
    Z_meas(1,k) = 3*X_real(1,k)*X_real(2,k)+ V_noise(1,k);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%卡尔曼滤波
        %状态预测:X_pre
        %diff_X_esti = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1);
        diff_X_esti = [3*cos(X_esti(2,k-1)) - X_esti(1,k-1);
                       2*cos(X_esti(1,k-1)) - X_esti(2,k-1)];

            %计算余项
            rema = [0;0];
            if k>L
                for i = 2:1:L+1
                   rema = rema + bino{1,i}*X_esti(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + bino{1,i}*X_esti(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q;     %一步状态预测
        
%             %计算余项
%             rema_fir = 0; rema_sec = 0;
%             
%             if k>L
%                 for i = 2:1:L+1
%                    %rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
%                    rema_fir = rema_fir + bino_fir(1,i)*X_esti(1,k+1-i);
%                    rema_sec = rema_sec + bino_sec(1,i)*X_esti(1,k+1-i);
%                 end
%             else
%                 for i = 2:1:k
%                    %rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
%                    rema_fir = rema_fir + bino_fir(1,i)*X_esti(1,k+1-i);
%                    rema_sec = rema_sec + bino_sec(1,i)*X_esti(1,k+1-i);
%                 end
%             end
%         %一步状态预测
%         %X_pre = diff_X_esti - rema + q;
%         X_pre(1,1) = diff_X_esti(1,1) - rema_fir + q(1,1);
%         X_pre(2,1) = diff_X_esti(2,1) - rema_sec + q(2,1);
        
        %预测误差协方差矩阵:P_pred

            %对误差矩阵进行cholsky分解
            S_chol = chol(P_esti{1,k-1})';

            %计算余项
            rema_P = zeros(2,2);
            if k>L+1
                for i = 2:1:L+2
                    %rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                    rema_P = rema_P + bino{1,i}*P_esti{1,k+1-i}*bino{1,i}';
                end
            else
                for i = 2:1:k
                    %rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                    rema_P = rema_P + bino{1,i}*P_esti{1,k+1-i}*bino{1,i}';
                end
            end

%         %临时变量 temp_fun : 函数差值,函数为单变量
%         temp_fun = zeros(2,2);
%         %X_esti_new_plus = X_esti(1,k-1) + h_con*S_chol(:,1);
%         for i=1:1:2
%             temp_fun(1,i) = 0.3*cos( X_esti(2,k-1) + h_con*S_chol(2,i) ) - ...
%                             ( X_esti(1,k-1) + h_con*S_chol(1,i) ) - ...
%                             0.3*cos( X_esti(2,k-1) - h_con*S_chol(2,i) ) ...
%                             + ( X_esti(1,k-1) - h_con*S_chol(1,i) );
%             temp_fun(2,i) = 0.2*cos( X_esti(1,k-1) + h_con*S_chol(1,i) ) -...
%                             0.2*cos( X_esti(1,k-1) - h_con*S_chol(1,i) );
%         end
%         temp_matrx = temp_fun(:,1)*temp_fun(:,1)' + temp_fun(:,2)*temp_fun(:,2)';

        temp_f_plus_hs1 = [3*cos( X_esti(2,k-1) + h_con*S_chol(2,1) ) - ...
                          ( X_esti(1,k-1) + h_con*S_chol(1,1) ); ...
                          2*cos( X_esti(1,k-1) + h_con*S_chol(1,1) ); ...
                          ];
        temp_f_subtraction_hs1 = [3*cos( X_esti(2,k-1) - h_con*S_chol(2,1) ) - ...
                          ( X_esti(1,k-1) - h_con*S_chol(1,1) ); ...
                          2*cos( X_esti(1,k-1) - h_con*S_chol(1,1) ); ...
                          ];

        temp_f_plus_hs2 = [3*cos( X_esti(2,k-1) + h_con*S_chol(2,2) ) - ...
                          ( X_esti(1,k-1) + h_con*S_chol(1,2) ); ...
                          2*cos( X_esti(1,k-1) + h_con*S_chol(1,2) ); ...
                          ];
        temp_f_subtraction_hs2 = [3*cos( X_esti(2,k-1) - h_con*S_chol(2,2) ) - ...
                          ( X_esti(1,k-1) - h_con*S_chol(1,2) ); ...
                          2*cos( X_esti(1,k-1) - h_con*S_chol(1,2) ); ...
                          ];

        %论文算法1中，状态预测误差协方差矩阵中第二项和第三项中的G^{k-1}_f为temp_fun的一半
        temp_initial =  ( 0.5*(temp_f_plus_hs1-temp_f_subtraction_hs1)*S_chol(:,1)'+ ...
                         0.5*(temp_f_plus_hs2-temp_f_subtraction_hs2)*S_chol(:,2)' ...
                         )*bino{1,2}'; 

        %P_pred = 1/(4*h_con^2) * temp_matrx + rema_P - temp_initial + Q; %

        P_pred = 1/(4*h_con^2)*( ...
                  (temp_f_plus_hs1 - temp_f_subtraction_hs1)*(temp_f_plus_hs1 - temp_f_subtraction_hs1)'  ...
                 +(temp_f_plus_hs2 - temp_f_subtraction_hs2)*(temp_f_plus_hs2 - temp_f_subtraction_hs2)'  ...
                ) + ...
                - 1/h_con * temp_initial - 1/h_con * temp_initial' + ...
                + Q + rema_P;
        
        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = 3*X_pre(1,1)*X_pre(2,1) + r_esti_A;

        %对状态预测误差协方差矩阵进行cholsky分解
        S_chol_pred = chol(P_pred)';

        %测量预测误差协方差矩阵:P_zpred ---- P_z_k|k-1
        P_zpred = 1/(4*h_con^2) * ( ... 
                  ( 3* (X_pre(1,1)+h_con*S_chol_pred(1,1)) * (X_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                   - ...
                    3* (X_pre(1,1)-h_con*S_chol_pred(1,1)) * (X_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                   )^2 + ...
                  ...
                  ( 3* (X_pre(1,1)+h_con*S_chol_pred(1,2)) * (X_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                  - ...
                    3* (X_pre(1,1)-h_con*S_chol_pred(1,2)) * (X_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                   )^2 ) ...
                  + R;

        %预测误差互协方差矩阵:P_xzpred
        %X_pre =  X_pre; % inv(S_chol_pred) * 
        P_xzpred = 1/(2*h_con) * ( ... 
                  ( 3*(X_pre(1,1)+h_con*S_chol_pred(1,1)) * (X_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                  - ...
                  ( 3* (X_pre(1,1)-h_con*S_chol_pred(1,1)) * (X_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                  ) ) * S_chol_pred(:,1)  ...
                  + ...
                  (3*(X_pre(1,1)+h_con*S_chol_pred(1,2)) * (X_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                   - ...
                  3* (X_pre(1,1)-h_con*S_chol_pred(1,2)) * (X_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                  ) * S_chol_pred(:,2) ...
                  );

        %计算卡尔曼增益:Kk(2*1)
        Kk = P_xzpred/P_zpred;

        r_esti_A(k) = ( (k-1)*r_esti_A(k-1) + Z_meas(1,k) - (3*X_pre(1,1)*X_pre(2,1)) )/k;

        %状态更新
        X_esti(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );
        
        %估计方差矩阵更新
        P_esti{1,k} = P_pred - Kk*P_zpred*Kk';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入与测量输出图
 k = 1:1:N;

%  figure;
% plot(k,X_real(1,k),'b',k,X_real(2,k),'r',k,X_real(3,k),'g','linewidth',LineWidth);
% set(gcf,'Position',[200 200 400 300]); 
% %axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
% axis normal
% ylabel('state estimate','FontSize',8)
% xlabel('iteration times','FontSize',8)
% %设置坐标轴刻度字体名称，大小
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
% %print 2.eps -depsc2 -r600
 
 
%状态估计图
figure;
plot(k,X_real(1,k),'b',k,X_esti(1,k),'--r','linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
axis normal
ylabel('state estimate','FontSize',8)
xlabel('iteration times','FontSize',8)
%设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
%print 2.eps -depsc2 -r600

figure;
plot(k,X_real(2,k),'b',k,X_esti(2,k),'--r','linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
axis normal
ylabel('state estimate','FontSize',8)
xlabel('iteration times','FontSize',8)
%设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
%print 2.eps -depsc2 -r600

%  k = 1:1:N;
% %误差图
% %误差计算
% % err_state = zeros(2,N);
% figure;
% err_state = X_real - X_esti;
% mse = err_state.^2;
% plot(k,mse(1,:),'b','linewidth',LineWidth);
% set(gcf,'Position',[200 200 400 300]); 
% %axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
% axis normal
% ylabel('square error','FontSize',8)
% xlabel('time(sec)','FontSize',8)
% %设置坐标轴刻度字体名称，大小
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('square error','Location','best');




