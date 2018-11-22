%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�Fractional Kalman filter algorithm for the states, parameters and
%        order of fractional system estimation
%   ��ע����������չ�������˲�����������
%  
%   ʵ��״̬����
%           
%   ��ע����Ϊ�����з����׿������˲����ĸ���
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 100;
LineWidth = 1.5;
h_con = sqrt(3);

%״̬������ʼ��
X_real = zeros(3,N);            %��ʵ״̬
X_esti = zeros(3,N);            %״̬���Ź���ֵ
P_esti = cell(1,N);             %����������
Z_meas = zeros(1,N);            %ʵ�ʹ۲�ֵ

%ϵͳ��������
% A = [0,1; -0.1,-0.2];           %ϵͳ����
% B = [0; 1];                     %
% C = [0.1,0.3];                  %
I = eye(3,3);                     %���ɵ�λ��
%I(3,3) = 0;

%����
q = [0 0 0]';                      %ϵͳ������ֵ
r = 1;                             %����������ֵ
Q = [0.3,0,0; 0,0.3,0; 0,0,0.001];%ϵͳ�����������
R = sqrt(0.5);                           %���������������
W_noise = Q*randn(3,N) + q;        %ϵͳ����
V_noise = R*randn(1,N) + r;        %��������

%��ʼֵ����
P_pred_0 = [10,0,0; 0,10,0;0,0,10];%��ʼԤ�ⷽ����
P_esti{1,1} = P_pred_0;               %��ʼ���Ʒ�����
x_0  = [0; 0; 0.2];                   %��ʼ״̬
X_real(:,1) = x_0;                    %��ʵ״̬��ʼֵ
X_esti(:,1) = [0;0;0];                %״̬���Ƴ�ʼֵ
Z_meas(:,1) = V_noise(:,1);           %�������ݳ�ʼֵ

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);          %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
alpha = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,N);          %΢�ֽ״�Ϊ1ʱGL�����µ�ϵ��
alpha = 1;
bino_thi(1,1) = 1;
for i = 2:1:N
    bino_thi(1,i) = (1-(alpha+1)/(i-1))*bino_thi(1,i-1);  
end

%����GL΢�ֶ�����ϵ������
gamma = cell(1,N);
temp_matrx = zeros(3,3);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

%ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ
U_input = 10*ones(1,N);            
for i = N/2+1:1:N
    U_input(:,i) = -10;
end

%GL�����¶̼���ԭ��ĳ���
L = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = [0;0;0];

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = [X_real(2,k-1); ...
                  -0.1*X_real(1,k-1)-X_real(2,k-1)*X_real(3,k-1)+ ...
                  U_input(1,k-1); 0] + W_noise(:,k-1);
    rema = [0;0;0];
    for i = 2:1:k
        rema = rema + gamma{1,i}*X_real(:,k+1-i);
    end
    X_real(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(:,k) = 0.1*X_real(1,k) + 0.3*X_real(2,k) + V_noise(1,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------����Ӧ�����׿������˲������ܲ���-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��ʼֵ���ã���ʼ������Ϊ�㣩
r_esti_A = zeros(1,N);      %����������ֵ����ֵ
r_esti_A(1,1) = 1.5;  

q_esti_A = zeros(3,N);      %����������ֵ����ֵ
q_esti_A(:,1) = [-1,-5,1]';  

% R_esti_A(1,1) = 0.3;        
% R_esti_A = zeros(1,N);      %������������������ֵ

for k=2:1:N
%�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = [X_esti(2,k-1); ...
                       -0.1*X_esti(1,k-1)-X_esti(2,k-1)*X_esti(3,k-1)+ ...
                       U_input(1,k-1); 0];
            %��������
            rema = [0;0;0];
            if k>L
                for i = 2:1:L+1
                   rema = rema + gamma{1,i}*X_esti(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + gamma{1,i}*X_esti(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q_esti_A(:,k-1);     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred

            %�����������cholsky�ֽ�
            S_chol = chol(P_esti{1,k-1})';

            %��������
            rema_P = [0,0,0;0,0,0;0,0,0];
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_esti{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + gamma{1,i}*P_esti{1,k+1-i}*gamma{1,i}';
                end
            end

        temp_f_plus_hs1 = [X_esti(2,k-1)+ h_con*S_chol(2,1); ...
                           -0.1*(X_esti(1,k-1)+h_con*S_chol(1,1))- ...
                           (X_esti(2,k-1)+h_con*S_chol(2,1))*( X_esti(3,k-1) ...
                           +h_con*S_chol(3,1))+U_input(1,k-1); 0];
        temp_f_subtraction_hs1 = [X_esti(2,k-1)-h_con*S_chol(2,1); ...
                         -0.1*(X_esti(1,k-1)-h_con*S_chol(1,1))- ...
                         (X_esti(2,k-1)-h_con*S_chol(2,1))* ...
                         (X_esti(3,k-1)-h_con*S_chol(3,1))+U_input(1,k-1); 0];   

        temp_f_plus_hs2 = [X_esti(2,k-1)+ h_con*S_chol(2,2); ...
                         -0.1*(X_esti(1,k-1)+h_con*S_chol(1,2))- ...
                         (X_esti(2,k-1)+h_con*S_chol(2,2))* ...
                         (X_esti(3,k-1)+h_con*S_chol(3,2))+U_input(1,k-1); 0];
        temp_f_subtraction_hs2 = [X_esti(2,k-1)-h_con*S_chol(2,2); ...
                         -0.1*(X_esti(1,k-1)-h_con*S_chol(1,2))- ...
                         (X_esti(2,k-1)-h_con*S_chol(2,2)) * ...
                         (X_esti(3,k-1)-h_con*S_chol(3,2))+U_input(1,k-1); 0]; 

        temp_f_plus_hs3 = [X_esti(2,k-1)+h_con*S_chol(2,3); ...
                         -0.1*(X_esti(1,k-1)+h_con*S_chol(1,3))- ...
                         (X_esti(2,k-1)+h_con*S_chol(2,3))* ...
                         (X_esti(3,k-1)+h_con*S_chol(3,3))+U_input(1,k-1); 0];
        temp_f_subtraction_hs3 = [X_esti(2,k-1)-h_con*S_chol(2,3); ...
                         -0.1*(X_esti(1,k-1)-h_con*S_chol(1,3))- ...
                         (X_esti(2,k-1)-h_con*S_chol(2,3)) * ...
                         (X_esti(3,k-1)-h_con*S_chol(3,3))+U_input(1,k-1); 0]; 

        temp_initial = ( 0.5*(temp_f_plus_hs1-temp_f_subtraction_hs1)*S_chol(:,1)'+ ...
                         0.5*(temp_f_plus_hs2-temp_f_subtraction_hs2)*S_chol(:,2)'+ ...
                         0.5*(temp_f_plus_hs3-temp_f_subtraction_hs3)*S_chol(:,3)' ...
                         )*gamma{1,2}'; 

        %F = [0,1,0; -0.1,-X_esti(3,k-1),-X_esti(2,k-1); 0,0,0];
        %P_pred = (F-gamma{1,2})*P_esti{1,k-1}*(F-gamma{1,2})'+ Q + rema_P;
        P_pred = 1/(4*h_con^2)*( ...
                  (temp_f_plus_hs1 - temp_f_subtraction_hs1)*(temp_f_plus_hs1 - temp_f_subtraction_hs1)'  ...
                 +(temp_f_plus_hs2 - temp_f_subtraction_hs2)*(temp_f_plus_hs2 - temp_f_subtraction_hs2)'  ...
                 +(temp_f_plus_hs3 - temp_f_subtraction_hs3)*(temp_f_plus_hs3 - temp_f_subtraction_hs3)'  ...
                ) + ...
                - 1/h_con * temp_initial - 1/h_con * temp_initial' + ...
                + Q + rema_P;

        %��״̬Ԥ�����Э����������cholsky�ֽ�
        S_chol_pred = chol(P_pred)';

        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = 0.1*X_pre(1,1) + 0.3*X_pre(2,1) + r; %_esti_A(1,k-1)

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = 1/(4*h_con^2) * ( ... 
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,1)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,1)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                   )^2 + ...
                  ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,2)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,2)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                   )^2 + ...
                   ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,3)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,3)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,3)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,3)) ...
                   )^2 ...
                  ) ...
                  + R;%_esti_A(1,k-1)

        %Ԥ����Э�������:P_xzpred
        Y_pre =  X_pre; % inv(S_chol_pred) * 
        P_xzpred = 1/(2*h_con) * ( ... 
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,1)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,1)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                   ) * S_chol_pred(:,1) + ...
                  ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,2)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,2)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                   ) * S_chol_pred(:,2) + ...
                   ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,3)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,3)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,3)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,3)) ...
                   ) * S_chol_pred(:,3)  ...
                  );

        %���㿨��������:Kk(2*1)
        Kk = P_xzpred/P_zpred;

        %״̬����
        X_esti(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %����������ֵ����
        %r_esti_A(k) = ( (k-1)*r_esti_A(k-1) + Z_meas(1,k) - (0.1*X_pre(1,1) + 0.3*X_pre(2,1)) )/k;
        %q_esti_A(:,k) = ( (k-1)*q_esti_A(:,k-1) + X_esti(:,k) - (diff_X_esti - rema) )/k; %(diff_X_esti - rema)
        %q_esti_A(:,k) = q_esti_A(:,k-1) - ( X_esti(:,k) - X_pre )/k; %(diff_X_esti - rema)
        q_esti_A(:,k) = ( (k-1)*q_esti_A(:,k-1) + X_esti(:,k) - (X_pre - q_esti_A(:,k-1)) )/k; %(diff_X_esti - rema)
        
        %���Ʒ���������
        P_esti{1,k} = P_pred - Kk*P_zpred*Kk';

        %���������������
        %R_esti_A(1,k) = ( (k-1)*R_esti_A(1,k-1) - P_zpred + R_esti_A(1,k-1) + ...
        %             ( Z_meas(:,k) - Z_esti )^2 )/k;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;
figure;
plot(k,U_input,'b',k,Z_meas,':r','linewidth',LineWidth);
legend('����','�������','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(1,:),'b',k,X_esti(1,:),':r','linewidth',LineWidth);
legend('ʵ��״̬1','����״̬1','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(2,:),'b',k,X_esti(2,:),':r','linewidth',LineWidth);
legend('ʵ��״̬2','����״̬2','Location','best');

%״̬����ͼ
figure;
plot(k,X_real(3,:),'b',k,X_esti(3,:),':r','linewidth',LineWidth);
legend('ʵ��״̬3','����״̬3','Location','best');

%���ͼ
%������
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti;
plot(k,r_esti_A(1,k),':r','linewidth',LineWidth);
legend('���Ʋ���������ֵ','Location','best');
line([0,N],[r,r]);






