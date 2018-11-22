%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%
%   ���ģ�fractional order CDKF
%
%   Ŀ�ģ�һ�׷��������Ĳ�ֿ������˲����㷨����-
%
%   ����ʵ��: MIMO system ��credit�� 
%
%   ������ϺõĶ�״̬���й��ƣ���ֵϵͳ������ֵ����
%
%   ��ע�����������Ĳ�ֿ������˲������㷨����
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 200;
LineWidth = 2;
h_con = sqrt(1.5);

%״̬������ʼ��
X_real = zeros(2,N);         %��ʵ״̬
X_esti = zeros(2,N);         %״̬���Ź���ֵ
%P_esti = zeros(1,N);        %����������
P_esti = cell(1,N);          %����������
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(2,2);                %���ɵ�λ��
%I(3,3) = 0;
 
%����
q = [0.2 0.2]';                 %ϵͳ������ֵ
r = 0.2;                      %����������ֵ
Q =  0.09*eye(2,2);         %ϵͳ�����������
R =  0.09;                  %���������������
W_noise = sqrt(Q)*randn(2,N) + q;   %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;   %��������

%Q_test = 100;

%��ʼֵ���ã���ʼ������Ϊ�㣩
bino = cell(1,N);            %����ʽϵ������
r_esti = zeros(2,N);         %����������ֵ����ֵ
P_pred_0 = 100*eye(2,2);      %��ʼԤ�ⷽ����
P_esti{1,1} = P_pred_0;      %��ʼ���Ʒ�����
x_0  = [0 0]';               %��ʼ״̬     
X_real(:,1) = x_0;           %��ʵ״̬��ʼֵ
X_esti(:,1) = [0 0]';        %״̬���Ƴ�ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.03ʱGL�����µ�ϵ��
alpha = 0.5;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);
end

bino_sec = zeros(1,N);       %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
alpha = 0.3;
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
% %����GL΢�ֶ�����ϵ������
% gamma = cell(1,N);
% temp_matrx = zeros(1,1);
% for i = 1:1:N 
%     temp_matrx(1,1) = bino_fir(1,i);
%     temp_matrx(2,2) = bino_sec(1,i);
%     temp_matrx(3,3) = bino_thi(1,i);
%     gamma{1,i} = temp_matrx;
% end

%ϵͳ��������:ǰһ������Ϊ0.04����һ������Ϊ0.04          
% for i = N/2+1:1:N
%     U_input(:,i) = -1;
% end

%GL�����¶̼���ԭ��ĳ���
L = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% diff_X_real   ��ʾkʱ��״̬��΢�� %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real= [0 0]';

for k=2:1:N
    %����ʵ��״̬
    diff_X_real(1,1) = 0.3*cos(X_real(2,k-1)) + W_noise(1,k-1);
    diff_X_real(2,1) = 0.2*cos(X_real(1,k-1)) + W_noise(2,k-1);
    
    %diff_X_real = 3*sin(2*X_real(1,k-1)) -X_real(1,k-1) + W_noise(1,k-1);
    rema_fir = 0; rema_sec = 0;
    for i = 2:1:k
        rema_fir = rema_fir + bino_fir(1,i)*X_real(1,k+1-i);
        rema_sec = rema_sec + bino_sec(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real(1,1) - rema_fir;
    X_real(2,k) = diff_X_real(2,1) - rema_sec;

    %ʵ�ʹ۲�ֵ
    Z_meas(1,k) = 0.3*X_real(1,k)*X_real(2,k)+ V_noise(1,k);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%�������˲�
        %״̬Ԥ��:X_pre
        %diff_X_esti = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1);
        diff_X_esti(1,1) = 0.3*cos(X_esti(2,k-1));
        diff_X_esti(2,1) = 0.2*cos(X_esti(1,k-1));

            %��������
            rema_fir = 0; rema_sec = 0;
            
            if k>L
                for i = 2:1:L+1
                   %rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                   rema_fir = rema_fir + bino_fir(1,i)*X_esti(1,k+1-i);
                   rema_sec = rema_sec + bino_sec(1,i)*X_esti(1,k+1-i);
                end
            else
                for i = 2:1:k
                   %rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                   rema_fir = rema_fir + bino_fir(1,i)*X_esti(1,k+1-i);
                   rema_sec = rema_sec + bino_sec(1,i)*X_esti(1,k+1-i);
                end
            end
        %һ��״̬Ԥ��
        %X_pre = diff_X_esti - rema + q;
        X_pre(1,1) = diff_X_esti(1,1) - rema_fir + q(1,1);
        X_pre(2,1) = diff_X_esti(2,1) - rema_sec + q(2,1);
        X_pre
        
        %Ԥ�����Э�������:P_pred

            %�����������cholsky�ֽ�
            S_chol = chol(P_esti{1,k-1})';

            %��������
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

        %��ʱ���� temp_fun : ������ֵ,����Ϊ������
        temp_fun = zeros(2,2);
        %X_esti_new_plus = X_esti(1,k-1) + h_con*S_chol(:,1);
        for i=1:1:2
            temp_fun(1,i) = 0.3*cos( X_esti(2,k-1) + h_con*S_chol(2,i) ) - ...
                            0.3*cos( X_esti(2,k-1) - h_con*S_chol(2,i) ) ;
            temp_fun(2,i) = 0.2*cos( X_esti(1,k-1) + h_con*S_chol(1,i) ) -...
                            0.2*cos( X_esti(1,k-1) - h_con*S_chol(1,i) );
        end
        temp_matrx = temp_fun(:,1)*temp_fun(:,1)' + temp_fun(:,2)*temp_fun(:,2)';

        %�����㷨1�У�״̬Ԥ�����Э��������еڶ���͵������е�G^{k-1}_fΪtemp_fun��һ��
        temp_initial = + 1/h_con*0.5*temp_fun*S_chol'*(-bino{1,2})' + ...
                         1/h_con*(-bino{1,2})*S_chol*0.5*temp_fun';

        P_pred = 1/(4*h_con^2) * temp_matrx + rema_P + temp_initial + Q; %

        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = 0.3*X_pre(1,1)*X_pre(2,1) + r;

        %��״̬Ԥ�����Э����������cholsky�ֽ�
        S_chol_pred = chol(P_pred)';

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = 1/(4*h_con^2) * ( ... 
                  ( 0.3* (X_pre(1,1)+h_con*S_chol_pred(1,1)) * (X_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                   - ...
                  0.3* (X_pre(1,1)-h_con*S_chol_pred(1,1)) * (X_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                   )^2 + ...
                  ...
                  ( 0.3* (X_pre(1,1)+h_con*S_chol_pred(1,2)) * (X_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                  - ...
                  0.3* (X_pre(1,1)-h_con*S_chol_pred(1,2)) * (X_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                   )^2 ) ...
                  + R;

        %Ԥ����Э�������:P_xzpred
        Y_pre =  X_pre; % inv(S_chol_pred) * 
        P_xzpred = 1/(2*h_con) * ( ... 
                  S_chol_pred(:,1)* ( 0.3*(Y_pre(1,1)+h_con*S_chol_pred(1,1)) * (Y_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                  - ...
                  ( 0.3* (Y_pre(1,1)-h_con*S_chol_pred(1,1)) * (Y_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                  ) ) ...
                  + ...
                  S_chol_pred(:,2)* (0.3*(Y_pre(1,1)+h_con*S_chol_pred(1,2)) * (Y_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                   - ...
                  0.3* (Y_pre(1,1)-h_con*S_chol_pred(1,2)) * (Y_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                  ) ...
                  );

        %���㿨��������:Kk(2*1)
        Kk = P_xzpred/P_zpred

        %r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;

        %״̬����
        X_esti(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );
         Z_meas(1,k) - Z_esti 
        
        %���Ʒ���������
        P_esti{1,k} = P_pred - Kk*P_zpred*Kk';
        k
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
 k = 1:1:N;

%  figure;
% plot(k,X_real(1,k),'b',k,X_real(2,k),'r',k,X_real(3,k),'g','linewidth',LineWidth);
% set(gcf,'Position',[200 200 400 300]); 
% %axis([xmin xmax ymin ymax])������������ָ��������
% axis normal
% ylabel('state estimate','FontSize',8)
% xlabel('iteration times','FontSize',8)
% %����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
% %print 2.eps -depsc2 -r600
 
 
%״̬����ͼ
figure;
plot(k,X_real(1,k),'b',k,X_esti(1,k),'--r','linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('state estimate','FontSize',8)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
%print 2.eps -depsc2 -r600

figure;
plot(k,X_real(2,k),'b',k,X_esti(2,k),'--r','linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('state estimate','FontSize',8)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
%print 2.eps -depsc2 -r600

%  k = 1:1:N;
% %���ͼ
% %������
% % err_state = zeros(2,N);
% figure;
% err_state = X_real - X_esti;
% mse = err_state.^2;
% plot(k,mse(1,:),'b','linewidth',LineWidth);
% set(gcf,'Position',[200 200 400 300]); 
% %axis([xmin xmax ymin ymax])������������ָ��������
% axis normal
% ylabel('square error','FontSize',8)
% xlabel('time(sec)','FontSize',8)
% %����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('square error','Location','best');




