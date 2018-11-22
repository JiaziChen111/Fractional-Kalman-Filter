clc
clear
N = 100;
r = 1;
LineWidth = 1.5;

load FEKF_AFCDKF_MIMO_COMPARISON_DATA

%×´Ì¬¹À¼ÆÍ¼
k = 1:1:N;
figure;
plot(k,X_real(3,:),'b',k,X_esti_A(3,:),'--r',k,X_esti_B(3,:),':g','linewidth',LineWidth);
legend('Êµ¼Ê×´Ì¬3','FCDKF¹À¼Æ×´Ì¬3','FEKF¹À¼Æ×´Ì¬3','Location','best');

figure;
plot(k,X_real(2,:),'b',k,X_esti_A(2,:),'--r',k,X_esti_B(2,:),':g','linewidth',LineWidth);
legend('Êµ¼Ê×´Ì¬2','FCDKF¹À¼Æ×´Ì¬2','FEKF¹À¼Æ×´Ì¬2','Location','best');

hold on
plot(k,X_real(1,:),'b',k,X_esti_A(1,:),'--r',k,X_esti_B(1,:),':g','linewidth',LineWidth);
legend('Êµ¼Ê×´Ì¬1','FCDKF¹À¼Æ×´Ì¬1','FEKF¹À¼Æ×´Ì¬1','Location','best');

figure;
err_state = X_real - X_esti_A;
plot(k,r_esti_A(1,k),':r','linewidth',LineWidth);
legend('¹À¼Æ²âÁ¿ÔëÉù¾ùÖµ','Location','best');
line([0,N],[r,r]);





