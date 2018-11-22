%˵��������PF�㷨������˸��������µ��״�Ŀ��������⡣
%     ��Ŀ��������ֱ���˶����״�λ��(x0,y0)��
%     ״̬����Ϊ:X(k)=PHI*X(k-1)+G*w(k-1);
%     �״�۲ⷽ��Ϊ:Z(k)=h(X(k))+v(k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main
 
M=100; %��������
T=1;   %�������
N=100;  %������
number=10; %���ؿ���������
x0=50000;y0=50000;vx=300;vy=-100;  %Ŀ���˶���ʼ״̬
delta_w=0.1;      %����������׼��
delta_r=50;       %��˸�����¹۲�����׼��
delta_theta1=1*pi/180;   %��������Ӧ��λ�Ǳ�׼��
delta_theta2=5*pi/180;   %��˸ЧӦ��Ӧ��λ�Ǳ�׼��
eta=0.3;                %�˲�������������ʽ��=0Ϊ��˹����������Ϊ��˸����
Q=delta_w^2*eye(2);     %��������������
R1=diag([delta_r^2,delta_theta1^2]);
R2=diag([delta_r^2,delta_theta2^2]);
R=(1-eta)*R1+eta*R2;    %��������������
G=[T^2/2,0;T,0;0,T^2/2;0,T];
%������ʵ����&����
X=zeros(4,M);
Z=zeros(2,M);
Xn=zeros(2,M);
w=sqrtm(Q)*randn(2,M);
v=sqrtm(R)*randn(2,M);
X(:,1)=[x0,vx,y0,vy]'; %��ʼ״̬
Z(:,1)=feval('hfun',X(:,1),x0,y0)+v(:,1);
Xn(:,1)=ffun(Z(:,1),x0,y0);
for t=2:M
    X(:,t)=feval('sfun',X(:,t-1),T)+G*w(:,t); %��ʵ����
    Z(:,t)=feval('hfun',X(:,t),x0,y0)+v(:,t);
    Xn(:,t)=ffun(Z(:,t),x0,y0); %����
end
%�����˲����Ƴ�ʼ��
Xmean_pf=zeros(number,4,M);
for i=1:number
    Xmean_pf(i,:,1)=X(:,1)+randn(4,1);
end
%��ʼ���棨number�Σ� 
for j=1:number
    %���Ӽ���ʼ��
    Xparticle_pf=zeros(4,M,N);
    XparticlePred_pf=zeros(4,M,N);
    zPred_pf=zeros(2,M,N);
    weight=zeros(M,N); %����Ȩֵ
    %��ʼ��
    for i=1:N
        Xparticle_pf(:,1,i)=[x0,vx,y0,vy]'+20*randn(4,1);
    end
    ww=randn(2,M);
    for t=2:M
        %����
        for i=1:N
            XparticlePred_pf(:,t,i)=feval('sfun',Xparticle_pf(:,t-1,i),T)...
                +G*sqrtm(Q)*ww(:,t-1);
        end
        %��Ҫ��Ȩֵ����
        for i=1:N
            zPred_pf(:,t,i)=feval('hfun',XparticlePred_pf(:,t,i),x0,y0);
            weight(t,i)=(1-eta)*inv(sqrt(2*pi*det(R1)))*exp(-.5*(Z(:,t)...
                -zPred_pf(:,t,i))'*inv(R1)*(Z(:,t)-zPred_pf(:,t,i)))...
                +eta*inv(sqrt(2*pi*det(R2)))*exp(-.5*(Z(:,t)-...
                zPred_pf(:,t,i))'*inv(R2)*(Z(:,t)-zPred_pf(:,t,i)))...
                + 1e-99; %Ȩֵ���㣬Ϊ����ȨֵΪ0����1e-99
        end
        weight(t,:)=weight(t,:)./sum(weight(t,:));%��һ��Ȩֵ
        outIndex = randomR(1:N,weight(t,:)');     %�������
        Xparticle_pf(:,t,:) = XparticlePred_pf(:,t,outIndex);%��ȡ�²���ֵ
        %״̬����
        mx=mean(Xparticle_pf(1,t,:));
        my=mean(Xparticle_pf(3,t,:));
        mvx=mean(Xparticle_pf(2,t,:));
        mvy=mean(Xparticle_pf(4,t,:));
        Xmean_pf(j,:,t)=[mx,mvx,my,mvy]';
    end
end
%��numberci���ؿ�����������վ�ֵ
Xpf=zeros(4,M);
for k=1:M
    Xpf(:,k)=[mean(Xmean_pf(:,1,k)),mean(Xmean_pf(:,2,k)),...
        mean(Xmean_pf(:,3,k)),mean(Xmean_pf(:,4,k))]';
end
%�������˲�����״̬����ʵ״̬֮���ƫ��
Div_Of_Xpf_X=Xpf-X;
%���������׼���RMSE 
for k=1:M
    sumX=zeros(4,1);
    for j=1:number
        sumX=sumX+(Xmean_pf(j,:,k)'-X(:,k)).^2;
    end
    RMSE(:,k)=sumX/number;
    Div_Std_Xpf(:,k)=sqrt(RMSE(:,k)-Div_Of_Xpf_X(:,k).^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); %���ٹ켣ͼ
plot(X(1,:),X(3,:),'b',Xn(1,:),Xn(2,:),'g',Xpf(1,:),Xpf(3,:),'r');
legend('��ʵ�켣','�۲�켣','���ƹ켣');
xlabel('X/m');ylabel('X/m');
figure(2);
subplot(2,2,1);plot(Div_Of_Xpf_X(1,:),'b');
ylabel('value/m');xlabel('(a) x����λ�ù�������ֵ����');
subplot(2,2,2);plot(Div_Of_Xpf_X(2,:),'b');
ylabel('value');xlabel('(b) x�����ٶȹ�������ֵ����');
subplot(2,2,3);plot(Div_Of_Xpf_X(3,:),'b');
ylabel('value/m');xlabel('(c) y����λ�ù�������ֵ����');
subplot(2,2,4);plot(Div_Of_Xpf_X(4,:),'b');
ylabel('value');xlabel('(d) y�����ٶȹ�������ֵ����');
figure(3);
subplot(2,2,1);plot(Div_Std_Xpf(1,:),'b');
ylabel('value');xlabel('(a) x����λ�ù�������׼������');
subplot(2,2,2);plot(Div_Std_Xpf(2,:),'b');
ylabel('value');xlabel('(b) x�����ٶȹ�������׼������');
subplot(2,2,3);plot(Div_Std_Xpf(3,:),'b');
ylabel('value');xlabel('(c) y����λ�ù�������׼������');
subplot(2,2,4);plot(Div_Std_Xpf(4,:),'b');
ylabel('value');xlabel('(d) y�����ٶȹ�������׼������');
figure(4);
subplot(2,2,1);plot(RMSE(1,:),'b');
ylabel('value');xlabel('(a) x����λ�ù���������������');
subplot(2,2,2);plot(RMSE(2,:),'b');
ylabel('value');xlabel('(b) x�����ٶȹ���������������');
subplot(2,2,3);plot(RMSE(3,:),'b');
ylabel('value');xlabel('(c) y����λ�ù���������������');
subplot(2,2,4);plot(RMSE(4,:),'b');
ylabel('value');xlabel('(d) y�����ٶȹ���������������');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%