%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   ��ע�����������Ĳ�ֿ������˲�����������
%         ָ������ʵ��: D^{0.7} x_k = 3*sin(2*x_{k-1}) + w_k
%                              y_k = x_k^3 + v_k
%   
%   ��ע�����������Ĳ�ֿ������˲���
%         ���׽���Ч������
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 100;
LineWidth = 1.5;
h_con = sqrt(3);

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
X_esti = zeros(1,N);         %״̬���Ź���ֵ
P_esti = zeros(1,N);         %����������
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %���ɵ�λ��
%I(3,3) = 0;
 
%����  
Q = 0.01;                       %ϵͳ�����������
R = 0.01;                       %���������������
W_noise = Q*randn(1,N);      %ϵͳ����
V_noise = R*randn(1,N);      %��������

%��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_esti(1,1) = P_pred_0;      %��ʼ���Ʒ�����
x_0  = 0;                    %��ʼ״̬
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
X_esti(1,1) = 0.5;           %״̬���Ƴ�ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
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

% %ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ-1
% U_input = ones(1,N);            
% for i = N/2+1:1:N
%     U_input(:,i) = -1;
% end

%GL�����¶̼���ԭ��ĳ���
L = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = 0;

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = 3*sin(2*X_real(1,k-1)) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(1,k) = X_real(1,k)^3 + V_noise(1,k);

    %�������˲�
        %�Թ�������������cholsky�ֽ�
        S_chol_esti = chol(P_esti(1,k-1))';

        %*************Step 1��״̬Ԥ��:X_pre*******************************%
        diff_X_esti = 3*sin(2*X_esti(1,k-1));
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

        add_temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol_esti) ) + ...
                       3*sin( 2*(X_esti(1,k-1)-h_con*S_chol_esti) ) - ...
                       2*3*sin( 2*X_esti(1,k-1) ) ;
        add_sec = 1/(2*h_con^2) * add_temp_fun;

        X_pre(k) = diff_X_esti + add_sec - rema;     %һ��״̬Ԥ��

        %*************Step 2��Ԥ�����Э�������:P_pred*********************%

            %��������
            rema_P = 0;
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            end
            %��ʱ���� temp_fun : ������ֵ,����Ϊ������
        temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol_esti) ) - ...
                   3*sin( 2*(X_esti(1,k-1)-h_con*S_chol_esti) );

        %Ԥ�����Э�������:P_pred
        P_pred = 1/(4*h_con^2) * temp_fun^2 + Q + rema_P + ...
                 1/h_con*0.5*temp_fun*S_chol_esti'*(-bino_fir(1,2))' + ...
                 1/h_con*(-bino_fir(1,2))*S_chol_esti*0.5*temp_fun'  + ...
                 (h_con^2-1)/(4*h_con^2) * add_temp_fun^2 ;

        %�����������cholsky�ֽ�
        S_chol_pre = chol(P_pred)';

        %*************Step 3�������������ֵZk��Z_esti**********************%

        temp_fun_esti = ( X_pre(k)+h_con*S_chol_pre )^3 + ...
                        ( X_pre(k)-h_con*S_chol_pre )^3 - 2*X_pre(k)^3;
        Z_esti = X_pre(k)^3 + 1/(2*h_con^2)*temp_fun_esti;

        %�������Ԥ�����Э����: P_z_pre
        P_z_pre = 1/(4*h_con^2)*( (X_pre(k)+h_con*S_chol_pre)^3 - ...
                  (X_pre(k)-h_con*S_chol_pre)^3 )^2 + R + ...
                  (4*h_con^2 -1)/(4*h_con^2)* temp_fun_esti^2;

        %�������Ԥ����� ��Э����: P_xz_pre
        P_xz_pre = 1/(2*h_con)*S_chol_pre*( (X_pre(k)+h_con*S_chol_pre)^3 - ...
                  (X_pre(k)-h_con*S_chol_pre)^3 );

        %���㿨��������:Kk(2*1)
        Kk = P_xz_pre/P_z_pre;

        %״̬����
        X_esti(1,k) = X_pre(k) + Kk*( Z_meas(1,k) - X_pre(k)^3 );

        %���Ʒ���������
        P_esti(1,k) = P_pred - Kk*P_z_pre*Kk';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;
figure;
plot(k,Z_meas,':r',k,X_real(1,:),'b','linewidth',LineWidth);
legend('�������','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(1,:),'b',k,X_esti(1,:),':r','linewidth',LineWidth);
legend('ʵ��״̬1','����״̬1','Location','best');

%���ͼ
%������
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti;
mse = err_state.^2;
plot(k,mse(1,:),'b','linewidth',LineWidth);
legend('���ƾ������','Location','best');






