clc
clear
N = 100;
r = 1;
LineWidth = 1.5;

load FEKF_AFCDKF_MIMO_COMPARISON_DATA

%״̬����ͼ
k = 1:1:N;
figure;
plot(k,X_real(3,:),'b',k,X_esti_A(3,:),'--r',k,X_esti_B(3,:),':g','linewidth',LineWidth);
legend('ʵ��״̬3','FCDKF����״̬3','FEKF����״̬3','Location','best');

figure;
plot(k,X_real(2,:),'b',k,X_esti_A(2,:),'--r',k,X_esti_B(2,:),':g','linewidth',LineWidth);
legend('ʵ��״̬2','FCDKF����״̬2','FEKF����״̬2','Location','best');

hold on
plot(k,X_real(1,:),'b',k,X_esti_A(1,:),'--r',k,X_esti_B(1,:),':g','linewidth',LineWidth);
legend('ʵ��״̬1','FCDKF����״̬1','FEKF����״̬1','Location','best');

figure;
err_state = X_real - X_esti_A;
plot(k,r_esti_A(1,k),':r','linewidth',LineWidth);
legend('���Ʋ���������ֵ','Location','best');
line([0,N],[r,r]);





