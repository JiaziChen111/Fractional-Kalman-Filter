clc;
clear all;

LineWidth = 2;
h_con = sqrt(3);
tf = 100; % ����ʱ��
N = 100;  % ���Ӹ���

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %���ɵ�λ��
%I(3,3) = 0;

R = 1; % ����Э�������
Q = 1; % ϵͳЭ�������
x = zeros(1,tf); % ϵͳ״̬��ʵֵ ��ʼֵ
y = zeros(1,tf); % ϵͳ״̬��ʵֵ ��ʼֵ
y(1,1) = x(1,1)^2 / 20 + sqrt(R) * randn;

P = zeros(1,tf); % ��������
P(1,1) = 2;      % ��ʼ�����ֲ��ķ���
xhatPart = zeros(1,tf);%״̬����ֵ

for i = 1 : N    
    xpart(i) = x(1,1) + sqrt(P(1,1)) * randn;%��ʼ״̬����x=0��ֵ������Ϊsqrt(P)�ĸ�˹�ֲ�
end
% xArr = [x];
% yArr = [];
% xhatArr = [x];
% PArr = [P];
%xhatPartArr = [xhatPart]; %

for k = 2 : tf    

x(1,k) = 0.5 * x(1,k-1) + 25 * x(1,k-1) / (1 + x(1,k-1)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    %kʱ����ʵֵ
    y(1,k) = x(1,k)^2 / 20 + sqrt(R) * randn;  %kʱ�̹۲�ֵ
    
 %% ����N������   
 for i = 1 : N 
     xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) ...
         + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;    %�������N������
     ypart = xpartminus(i)^2 / 20;%ÿ�����Ӷ�Ӧ�۲�ֵ
     vhat = y(1,k) - ypart;%����ʵ�۲�֮�����Ȼ
     q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R); 
     %ÿ�����ӵ���Ȼ�����ƶ�
 end
 
 %%
%Ȩֵ��һ��
qsum = sum(q);
for i = 1 : N
    q(i) = q(i) / qsum; %��һ�����Ȩֵ q
end

%%
 %����Ȩֵ���²���
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
%����״̬����ֵ��ΪN�����ӵ�ƽ��ֵ�����ﾭ�����²�����������ӵ�Ȩֵ��ͬ
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

