%*************************************************************************%
%                      �޼��������˲������渴��                           %
%_________________________________________________________________________%
%   ���� : ������ϵͳ�˲����� P108 5.4.2 ǿ������ϵͳģ�ͷ���
%   Ŀ�� : �޼��������˲������渴��
%   ����ʵ�� ::
%                   | 3*sin(2*x_2)                 |   | 1 |
%             x_k = | x_1 + exp(-0.05x_3) + 10     | + | 1 |w_k
%                   | = x_1(x_2 + x_3)/5 + |x_1|/2 |   | 1 |
%
%   ��� : �˲�Ч������
%
%   ��ע : 
%
%*************************************************************************%

clc
clear

%���沽��
N = 50;

q = 0.3;                %ϵͳ������ֵ
r = 0.5;                %����������ֵ
Q = 0.7;                %ϵͳ�����������
R = 1.0;                %���������������

% %GL�����¶̼���ԭ��ĳ���
% L = N+1;

% %����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
% bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
% alpha = 0.7;
% bino_fir(1,1) = 1;
% for i = 2:1:N
%     bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
% end

I = eye(3,3);                %���ɵ�λ��

%err_state_FEKF  = zeros(kk_N,N);

X_state_real = zeros(3,N);         %��ʵ״̬
Z_state_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%����
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

x_0  = [-0.7; 1; 1];               %��ʼ״̬     
X_state_real(:,1) = x_0;           %��ʵ״̬��ʼֵ
Z_state_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

f=@(x)[3*sin(2*x(2)); ...
       x(1) + exp(-0.05*x(3)) + 10; ...
       0.2 * x(1) * (x(2) + x(3)) + 0.5*abs(x(1)) ];

h=@(x)x(1) + x(2) * x(3);
   
for k=2:1:N
    %����ʵ��״̬
    X_state_real(:,k) = f(X_state_real(:,k-1)) + [1; 1; 1] * W_noise(1,k-1);
        
    %ʵ�ʹ۲�ֵ
    Z_state_meas(1,k) = h(X_state_real(:,k)) + V_noise(1,k); 
end

%*************************************************************************%
%-------------------------�޼��������˲������ܲ���------------------------%
%*************************************************************************%

X_state_esti = zeros(3,N);      %״̬���Ź���ֵ
P_xesti      = cell(1,N);      %����������

%��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0     = eye(3,3);        %��ʼԤ�ⷽ����
P_xesti{1,1} = P_pred_0;        %��ʼ���Ʒ�����

state_dim = 3;
L_sample  = 2 * state_dim +1;

SigmaPoints = zeros(3, L_sample);
SigmaWeight = zeros(1, L_sample);
GammaPoints = zeros(3, L_sample);
ChiPoints   = zeros(1, L_sample);

for k=2:1:N
    %�ԳƲ���
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_state_esti(:,k-1), P_xesti{1,k-1}); 
    
    for i = 1 : 1 : L_sample
        GammaPoints(:,i) = f(SigmaPoints(:,i)) + [1; 1; 1]*q;
    end
    % Predicted state
    X_pre = [0; 0; 0];
    for i = 1 : 1 : L_sample
        X_pre = X_pre +  SigmaWeight(:,i)*GammaPoints(:,i);
    end
    
    % Predicted state error covariance 
    P_xpre = Q*I;
    for i = 1 : 1 : L_sample 
        P_xpre = P_xpre + SigmaWeight(:,i)*(GammaPoints(:,i)-X_pre)* ...
                 (GammaPoints(:,i)-X_pre)';
    end
    
    % measurement update
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_pre, P_xpre); 

    for i = 1 : 1 : L_sample
        ChiPoints(:,i) = h(SigmaPoints(:,i)) + r;
    end

    % Predicted measurement
    Z_pre = 0;
    for i = 1 : 1 : L_sample
        Z_pre = Z_pre +  SigmaWeight(:,i)*ChiPoints(:,i);
    end

    % Predicted measurement error covariance 
    P_zpre = R;
    for i = 1 : 1 : L_sample 
        P_zpre = P_zpre + SigmaWeight(:,i)*(ChiPoints(:,i)-Z_pre)* ...
                 (ChiPoints(:,i)-Z_pre)';
    end

    % cross-variance 
    P_xzpre = [0; 0; 0];
    for i = 1 : 1 : L_sample 
        P_xzpre = P_xzpre + SigmaWeight(:,i)*(SigmaPoints(:,i)-X_pre)* ...
                 (ChiPoints(:,i)-Z_pre)';
    end
    
    % Kalman gain
    Kk = P_xzpre/P_zpre;
    
    % estimated state
    X_state_esti(:,k) = X_pre + Kk*( Z_state_meas(1,k) - Z_pre );

    % estimation error covariance
    P_xesti{1,k} = P_xpre - Kk*P_zpre*Kk';
    
end


%*******************************%
%   ��ͼ��� ��ֵ�������ɢ��ͼ
%*******************************%
%������������ͼ
k = 1:1:N;

LineWidth = 2;

%square error
figure;
plot(k,X_state_real(3,:),'r',k,X_state_esti(3,:),'b--','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
% axis normal
% axis([ -10 N 0 6 ])
ylabel('$x_3$','FontSize',8)
xlabel('time(sec)','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
legend('Real state 3','Estimation state 3','Location','best');

figure;
plot(k,X_state_real(2,:),'r',k,X_state_esti(2,:),'b--','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
% axis normal
% axis([ -10 N 0 6 ])
ylabel('$x_3$','FontSize',8)
xlabel('time(sec)','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
legend('Real state 2','Estimation state 2','Location','best');


figure;
plot(k,X_state_real(1,:),'r',k,X_state_esti(1,:),'b--','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
% axis normal
% axis([ -10 N 0 6 ])
ylabel('$x_3$','FontSize',8)
xlabel('time(sec)','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
legend('Real state 1','Estimation state 1','Location','best');






