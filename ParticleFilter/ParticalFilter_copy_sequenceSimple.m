clc;
clear all;

LineWidth = 2;
h_con = sqrt(3);
tf = 100; % 仿真时长
N = 100;  % 粒子个数

%系统矩阵设置
% A = [0,1; -0.1,-0.2];      %系统矩阵
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %生成单位阵
%I(3,3) = 0;

R = 1; % 测量协方差矩阵
Q = 1; % 系统协方差矩阵
x = zeros(1,tf); % 系统状态真实值 初始值
y = zeros(1,tf); % 系统状态真实值 初始值
y(1,1) = x(1,1)^2 / 20 + sqrt(R) * randn;

P = zeros(1,tf); % 采样方差
P(1,1) = 2;      % 初始采样分布的方差
xhatPart = zeros(1,tf);%状态估计值

for i = 1 : N    
    xpart(i) = x(1,1) + sqrt(P(1,1)) * randn;%初始状态服从x=0均值，方差为sqrt(P)的高斯分布
end
% xArr = [x];
% yArr = [];
% xhatArr = [x];
% PArr = [P];
%xhatPartArr = [xhatPart]; %

for k = 2 : tf    

x(1,k) = 0.5 * x(1,k-1) + 25 * x(1,k-1) / (1 + x(1,k-1)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    %k时刻真实值
    y(1,k) = x(1,k)^2 / 20 + sqrt(R) * randn;  %k时刻观测值
    
 %% 采样N个粒子   
 for i = 1 : N 
     xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) ...
         + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;    %采样获得N个粒子
     ypart = xpartminus(i)^2 / 20;%每个粒子对应观测值
     vhat = y(1,k) - ypart;%与真实观测之间的似然
     q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R); 
     %每个粒子的似然即相似度
 end
 
 %%
%权值归一化
qsum = sum(q);
for i = 1 : N
    q(i) = q(i) / qsum; %归一化后的权值 q
end

%%
 %根据权值重新采样
  for i = 1 : N 
      u = rand;
      qtempsum = 0; 
      for j = 1 : N
          qtempsum = qtempsum + q(j); 
          if qtempsum >= u 
              xpart(i) = xpartminus(j);
              break;
          end 
      end
  end
xhatPart(1,k) = mean(xpart);

%%
%最后的状态估计值即为N个粒子的平均值，这里经过重新采样后各个粒子的权值相同
% xArr = [xArr x];   
% yArr = [yArr y];  
% % xhatArr = [xhatArr xhat]; 
% PArr = [PArr P]; 
% xhatPartArr = [xhatPartArr xhatPart];
end

t = 1 : tf;
figure;
plot(t, x, 'b-.', t, xhatPart, 'k-');
legend('Real Value','Estimated Value');
set(gca,'FontSize',10); 
xlabel('time step'); 
ylabel('state');
title('Particle filter')
%xhatRMS = sqrt((norm(x - xhat))^2 / tf);
%xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
figure;
plot(t,abs(x-xhatPart),'b');
title('The error of PF')

% t = 0 : tf;
% figure;
% plot(t, xArr, 'b-.', t, xhatPartArr, 'k-');
% legend('Real Value','Estimated Value');
% set(gca,'FontSize',10); 
% xlabel('time step'); 
% ylabel('state');
% title('Particle filter')
% xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
% xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
% figure;
% plot(t,abs(xArr-xhatPartArr),'b');
% title('The error of PF')

