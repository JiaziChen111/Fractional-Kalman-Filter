%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�FCDKF��FEKF�����ܱȽ�
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   �����
%
%   ��ע��
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 50;

h_con = sqrt(4);

q = 0.5;                   %ϵͳ������ֵ
r = 0.5;                   %����������ֵ
Q = 0.81;                %ϵͳ�����������
R = 0.81;                %���������������

%GL�����¶̼���ԭ��ĳ���
L = N+1;

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%ϵͳ��������
I = eye(1,1);                %���ɵ�λ��

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%����
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

x_0  = 0;                    %��ʼ״̬     
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

% ϵͳ�������������
f=@(x)3*sin(2*x)-x;
h=@(x)x;

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = f(X_real(1,k-1)) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(1,k) = h(X_real(1,k)) + V_noise(1,k); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------�����׿������˲������ܲ���------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti = zeros(1,N);         %״̬���Ź���ֵ
P_xesti = zeros(1,N);        %����������

P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_xesti(1,1) = P_pred_0;     %��ʼ���Ʒ�����

for k=2:1:N
      %�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = f(X_esti(1,k-1));
            %��������
            rema = 0;
            if k>L
                for i = 2:1:L+1
                   rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred

            %�����������cholsky�ֽ�
            S_chol = chol(P_xesti(1,k-1))';

            %��������
            rema_P = 0;
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            end

        %��ʱ���� temp_fun : ������ֵ,����Ϊ������
        temp_fun = f(X_esti(1,k-1)+h_con*S_chol) - f(X_esti(1,k-1)-h_con*S_chol);
        temp = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                  1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                  1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
        P_xpred = temp + Q;
        
        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = h(X_pre) + r;

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = P_xpred + R;

        %���㿨��������:Kk(2*1)
        Kk = P_xpred/P_zpred;

        %״̬����
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %���Ʒ���������
        P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';
end

%������������ͼ
k = 1:1:N;

LineWidth = 1.5;

figure;
plot(k,X_real(1,:),'r',k,X_esti(1,:),'b--','linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 N -6 6 ])
ylabel('x','FontSize',8)
xlabel('time(sec)','FontSize',8)
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');







