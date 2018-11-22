clear;
echo off; 
%-----------------Parameter--------------------------
xMEAN=4.5;xSIGMA=0.2;              %��ʼ���������ӵ�λ�ú��ٶ�                
yMEAN=4;ySIGMA=0.3;
VxMEAN=0.2;VxSIGMA=0.01;             
VyMEAN=0.15;VySIGMA=0.02;

UxMEAN=0;UxSIGMA=0.015625;          % ��ʼ��x�����y������������Ϊ���ֵ��˹������
UyMEAN=0;UySIGMA=0.015625;          % ��׼�� Ϊ 0.015625 
UzMEAN=0;UzSIGMA=0.001;             

rSTEP=30;                         % ����30��ʱ��
wk=-1/(2*UzSIGMA);                % ����Ȩ�س���
wj=1/sqrt(2*pi*UzSIGMA);          % ����Ȩ�س���1/sqr(2*pi)
rN=8000;                          % ������Ϊ1024

w_buffer=zeros(rN,1);             % �洢Ȩֵ�Ŀռ�
r_buffer=zeros(rN,1);             % �洢���ƴ����Ŀռ�
i_buffer=zeros(rN,1);             % �洢����ָ��Ŀռ�

tX=zeros(rSTEP,1);tY=zeros(rSTEP,1);tVx=zeros(rSTEP,1);tVy=zeros(rSTEP,1);tZ=zeros(rSTEP,1);  %ÿһʱ��״̬�Ĺ���ֵ
rX=zeros(rSTEP,1);rY=zeros(rSTEP,1);rVx=zeros(rSTEP,1);rVy=zeros(rSTEP,1);                    %ÿһʱ��״̬����ʵֵ
x0=4.5;y0=4.5;Vx0=0.2;Vy0=0.2;           % ��ʼ����ʵ״̬�ĳ�ֵ


%-----------------��ʵ״ֵ̬�͹۲�Ƕ�ֵ-----------------------------
tX(1,1)=x0;tY(1,1)=y0;tVx(1,1)=Vx0;tVy(1,1)=Vy0;
w0=normrnd(UxMEAN,UxSIGMA,rSTEP);                      %����
w1=normrnd(UyMEAN,UySIGMA,rSTEP);
w2=normrnd(UzMEAN,UzSIGMA,rSTEP);
tZ(1,1)=atan(tY(1,1)./tX(1,1))+w2(1,1);            %��ʼ������          
for t=2:rSTEP
    tX(t,1)=tX(t-1,1)+tVx(t-1,1)+w0(t,1);          % tX(t,1)�����洢Ŀ����x�����ϵ���ʵ����ֵ
    tY(t,1)=tY(t-1,1)+tVy(t-1,1)+w1(t,1);          % tY(t,1)�����洢Ŀ����y�����ϵ���ʵ����ֵ
    tVx(t,1)=tVx(t-1,1)+0.5.*w0(t,1);              % tVx(t,1)�����洢Ŀ����x�����ϵ���ʵ�ٶ�ֵ
    tVy(t,1)=tVy(t-1,1)+0.5.*w1(t,1);              % tVy(t,1)�����洢Ŀ����y�����ϵ���ʵ�ٶ�ֵ
    tZ(t,1)=atan(tY(t,1)./tX(t,1))+w2(t,1);        % tZ(t,1)�����洢ÿһʱ�̵Ĳ�����
end;

%---------------��ʼ��----------------------------
x_buffer=normrnd(xMEAN,xSIGMA,rN);       
Vx_buffer=normrnd(VxMEAN,VxSIGMA,rN);    
y_buffer=normrnd(yMEAN,ySIGMA,rN);     
Vy_buffer=normrnd(VyMEAN,VySIGMA,rN);
for i=1:rN
    w_buffer(i,1)=1/rN;                % ��ʼȨ��ֵ
    r_buffer(i,1)=1;                   % ��ʼ���ƴ���
    i_buffer(i,1)=i;                   % ��ʼ����ָ��
end;
iR=rN;                                 % ���ز�����������ʼΪrN
    
%------------------�����˲�-------------------------------
for t=1:rSTEP  
   %----------------����----------------------------
    indr=1;indd=rN;
    w0=normrnd(UxMEAN,UxSIGMA,rN);
    w1=normrnd(UyMEAN,UySIGMA,rN);
    reg=zeros(1,4);
    for indr=1:iR
       x=i_buffer(indr,1);           % ����ָ��
       reg(1,1)=x_buffer(x,1);       % reg�����洢ÿһʱ��Ŀ��x,y�����ϵĹ�������ֵ��x,y�����ϵĵĹ����ٶ�ֵ
       reg(1,2)=Vx_buffer(x,1);      % 
       reg(1,3)=y_buffer(x,1);
       reg(1,4)=Vy_buffer(x,1);
     
       x_buffer(x,1)=reg(1,1)+reg(1,2)+w0(indr,1); %
       Vx_buffer(x,1)=reg(1,2)+0.5.*w0(indr,1); 
       y_buffer(x,1)=reg(1,3)+reg(1,4)+w1(indr,1); 
       Vy_buffer(x,1)=reg(1,4)+0.5.*w1(indr,1); 
     
       for k=r_buffer(indr,1)-1:-1:1
           x_buffer(i_buffer(indd,1),1)=reg(1,1)+reg(1,2)+w0(indd,1);   %��������
           Vx_buffer(i_buffer(indd,1),1)=reg(1,2)+0.5.*w0(indd,1);  
           y_buffer(i_buffer(indd,1),1)=reg(1,3)+reg(1,4)+w1(indd,1); 
           Vy_buffer(i_buffer(indd,1),1)=reg(1,4)+0.5.*w1(indd,1); 
           indd=indd-1;
      end;
   end;    

  %----------------Ȩֵ����---------------------

   W=0;
   for i=1:rN;
      x=tZ(t,1)-atan(y_buffer(i,1)./x_buffer(i,1));     %���㹫ʽ��������ʽ1.5
      w_buffer(i,1)=exp(wk*x*x);                       %��i�����ӵ�Ȩֵ
      W=W+w_buffer(i,1);                                %Ȩֵ���ۼӺ�
  end;
%----------------״̬���---------------------------
   X=0;Vx=0;Y=0;Vy=0;
   for i=1:rN;
      X=X+w_buffer(i,1).*x_buffer(i,1);            % ����x��λֵ
      Vx=Vx+w_buffer(i,1).*Vx_buffer(i,1);         % ����x�����ٶ�ֵ
      Y=Y+w_buffer(i,1).*y_buffer(i,1);            % ����y��λֵ
      Vy=Vy+w_buffer(i,1).*Vy_buffer(i,1);         % ����y�����ٶ�ֵ
  end; 
   rX(t,1)=X/W;                                    % ���ֵ��һ��
   rVx(t,1)=Vx/W;
   rY(t,1)=Y/W;
   rVy(t,1)=Vy/W;
 
  %----------------�ز���-------------------------
  indr=1;indd=rN;                       %����ָ��ĳ�ʼֵ
  K=rN/W;                               %�����м����K
  U=rand(1,1);                          %����һ�������ֵ
  for i=1:rN;                           %��ѭ��
      temp=K.*w_buffer(i,1)-U;          %���һ���м����temp
      r_buffer(indr,1)=quzheng(temp);   %�洢���ƴ���
      U=r_buffer(indr,1)-temp;          %������ֵ
      if r_buffer(indr,1)>0             %
         i_buffer(indr,1)=i;            %�洢���������ӵĵ�ַ
         indr=indr+1;
      else
         i_buffer(indd,1)=i;            % �洢���������ӵĵ�ַ
         indd=indd-1;                   
     end;
  end;
  iR=indr-1;
 %--------------------------------------------- 
end;

%------------------��ͼ--------------------------
figure(1);
plot(tX,tY,'-.g*',rX,rY,'-ro');
legend('true','estimate');
