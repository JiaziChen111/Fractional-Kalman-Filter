clc;
clear;
%close all;

x = 0; % ϵͳ״̬��ʵֵ ��ʼֵ
R = 1; % ����Э�������
Q = 1; % ϵͳЭ�������
tf = 100; % ����ʱ��
N = 100;  % ���Ӹ���
P = 2;    % ��ʼ�����ֲ��ķ���
xhatPart = x;

% ��ʼ������������ʲ�����ʼ����Ⱥ
for i = 1 : N    
    xpart(i) = x + sqrt(P) * randn;%��ʼ״̬����x=0��ֵ������Ϊsqrt(P)�ĸ�˹�ֲ�
end

xArr = [x];
yArr = [x^2 / 20 + sqrt(R) * randn];
xhatArr = [x];
PArr = [P];
xhatPartArr = [xhatPart];

for k = 1 : tf    

 %% ϵͳʵ��״̬���۲�ֵ
    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    %kʱ����ʵֵ
    y = x^2 / 20 + sqrt(R) * randn;  %kʱ�̹۲�ֵ
    
 %% ����N������   
 for i = 1 : N 
     xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) ...
         + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;    %�������N������
     ypart = xpartminus(i)^2 / 20;%ÿ�����Ӷ�Ӧ�۲�ֵ
     
     vhat = y - ypart;%����ʵ�۲�֮�����Ȼ
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

%%

xhatPart = mean(xpart);
P = 0;
for i = 1 : N
    P = P + u*(xpart(i) - xhatPart)*(xpart(i) - xhatPart)';
end 
  
%%
%����״̬����ֵ��ΪN�����ӵ�ƽ��ֵ�����ﾭ�����²�����������ӵ�Ȩֵ��ͬ

xArr = [xArr x];       %ʵ��״ֵ̬
yArr = [yArr y];       %ʵ�ʲ���ֵ
% xhatArr = [xhatArr xhat]; 
PArr = [PArr P];       %
xhatPartArr = [xhatPartArr xhatPart];

end


t = 0 : tf;
figure;
plot(t, xArr, 'b-.', t, xhatPartArr, 'k-');
legend('Real Value','Estimated Value');
set(gca,'FontSize',10); 
xlabel('time step'); 
ylabel('state');
title('Particle filter')
xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
figure;
plot(t,abs(xArr-xhatPartArr),'b');
title('The error of PF')

