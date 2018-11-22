clear;
echo off; 
%-----------------Parameter--------------------------
xMEAN=4.5;xSIGMA=0.2;              %初始化采样粒子的位置和速度                
yMEAN=4;ySIGMA=0.3;
VxMEAN=0.2;VxSIGMA=0.01;             
VyMEAN=0.15;VySIGMA=0.02;

UxMEAN=0;UxSIGMA=0.015625;          % 初始化x方向和y方向噪声，均为零均值高斯白噪声
UyMEAN=0;UySIGMA=0.015625;          % 标准差 为 0.015625 
UzMEAN=0;UzSIGMA=0.001;             

rSTEP=30;                         % 采样30个时刻
wk=-1/(2*UzSIGMA);                % 粒子权重常数
wj=1/sqrt(2*pi*UzSIGMA);          % 粒子权重常数1/sqr(2*pi)
rN=8000;                          % 粒子数为1024

w_buffer=zeros(rN,1);             % 存储权值的空间
r_buffer=zeros(rN,1);             % 存储复制次数的空间
i_buffer=zeros(rN,1);             % 存储粒子指针的空间

tX=zeros(rSTEP,1);tY=zeros(rSTEP,1);tVx=zeros(rSTEP,1);tVy=zeros(rSTEP,1);tZ=zeros(rSTEP,1);  %每一时刻状态的估计值
rX=zeros(rSTEP,1);rY=zeros(rSTEP,1);rVx=zeros(rSTEP,1);rVy=zeros(rSTEP,1);                    %每一时刻状态的真实值
x0=4.5;y0=4.5;Vx0=0.2;Vy0=0.2;           % 初始化真实状态的初值


%-----------------真实状态值和观测角度值-----------------------------
tX(1,1)=x0;tY(1,1)=y0;tVx(1,1)=Vx0;tVy(1,1)=Vy0;
w0=normrnd(UxMEAN,UxSIGMA,rSTEP);                      %噪声
w1=normrnd(UyMEAN,UySIGMA,rSTEP);
w2=normrnd(UzMEAN,UzSIGMA,rSTEP);
tZ(1,1)=atan(tY(1,1)./tX(1,1))+w2(1,1);            %初始测量角          
for t=2:rSTEP
    tX(t,1)=tX(t-1,1)+tVx(t-1,1)+w0(t,1);          % tX(t,1)用来存储目标在x方向上的真实坐标值
    tY(t,1)=tY(t-1,1)+tVy(t-1,1)+w1(t,1);          % tY(t,1)用来存储目标在y方向上的真实坐标值
    tVx(t,1)=tVx(t-1,1)+0.5.*w0(t,1);              % tVx(t,1)用来存储目标在x方向上的真实速度值
    tVy(t,1)=tVy(t-1,1)+0.5.*w1(t,1);              % tVy(t,1)用来存储目标在y方向上的真实速度值
    tZ(t,1)=atan(tY(t,1)./tX(t,1))+w2(t,1);        % tZ(t,1)用来存储每一时刻的测量角
end;

%---------------初始化----------------------------
x_buffer=normrnd(xMEAN,xSIGMA,rN);       
Vx_buffer=normrnd(VxMEAN,VxSIGMA,rN);    
y_buffer=normrnd(yMEAN,ySIGMA,rN);     
Vy_buffer=normrnd(VyMEAN,VySIGMA,rN);
for i=1:rN
    w_buffer(i,1)=1/rN;                % 初始权重值
    r_buffer(i,1)=1;                   % 初始复制次数
    i_buffer(i,1)=i;                   % 初始粒子指针
end;
iR=rN;                                 % 需重采样粒子数初始为rN
    
%------------------粒子滤波-------------------------------
for t=1:rSTEP  
   %----------------采样----------------------------
    indr=1;indd=rN;
    w0=normrnd(UxMEAN,UxSIGMA,rN);
    w1=normrnd(UyMEAN,UySIGMA,rN);
    reg=zeros(1,4);
    for indr=1:iR
       x=i_buffer(indr,1);           % 粒子指针
       reg(1,1)=x_buffer(x,1);       % reg用来存储每一时刻目标x,y方向上的估计坐标值和x,y方向上的的估计速度值
       reg(1,2)=Vx_buffer(x,1);      % 
       reg(1,3)=y_buffer(x,1);
       reg(1,4)=Vy_buffer(x,1);
     
       x_buffer(x,1)=reg(1,1)+reg(1,2)+w0(indr,1); %
       Vx_buffer(x,1)=reg(1,2)+0.5.*w0(indr,1); 
       y_buffer(x,1)=reg(1,3)+reg(1,4)+w1(indr,1); 
       Vy_buffer(x,1)=reg(1,4)+0.5.*w1(indr,1); 
     
       for k=r_buffer(indr,1)-1:-1:1
           x_buffer(i_buffer(indd,1),1)=reg(1,1)+reg(1,2)+w0(indd,1);   %复制粒子
           Vx_buffer(i_buffer(indd,1),1)=reg(1,2)+0.5.*w0(indd,1);  
           y_buffer(i_buffer(indd,1),1)=reg(1,3)+reg(1,4)+w1(indd,1); 
           Vy_buffer(i_buffer(indd,1),1)=reg(1,4)+0.5.*w1(indd,1); 
           indd=indd-1;
      end;
   end;    

  %----------------权值计算---------------------

   W=0;
   for i=1:rN;
      x=tZ(t,1)-atan(y_buffer(i,1)./x_buffer(i,1));     %计算公式见报告中式1.5
      w_buffer(i,1)=exp(wk*x*x);                       %第i个粒子的权值
      W=W+w_buffer(i,1);                                %权值的累加和
  end;
%----------------状态输出---------------------------
   X=0;Vx=0;Y=0;Vy=0;
   for i=1:rN;
      X=X+w_buffer(i,1).*x_buffer(i,1);            % 计算x方位值
      Vx=Vx+w_buffer(i,1).*Vx_buffer(i,1);         % 计算x方向速度值
      Y=Y+w_buffer(i,1).*y_buffer(i,1);            % 计算y方位值
      Vy=Vy+w_buffer(i,1).*Vy_buffer(i,1);         % 计算y方向速度值
  end; 
   rX(t,1)=X/W;                                    % 输出值归一化
   rVx(t,1)=Vx/W;
   rY(t,1)=Y/W;
   rVy(t,1)=Vy/W;
 
  %----------------重采样-------------------------
  indr=1;indd=rN;                       %设置指针的初始值
  K=rN/W;                               %计算中间变量K
  U=rand(1,1);                          %生成一个随机阈值
  for i=1:rN;                           %主循环
      temp=K.*w_buffer(i,1)-U;          %添加一个中间变量temp
      r_buffer(indr,1)=quzheng(temp);   %存储复制次数
      U=r_buffer(indr,1)-temp;          %更新阈值
      if r_buffer(indr,1)>0             %
         i_buffer(indr,1)=i;            %存储被复制粒子的地址
         indr=indr+1;
      else
         i_buffer(indd,1)=i;            % 存储被抛弃粒子的地址
         indd=indd-1;                   
     end;
  end;
  iR=indr-1;
 %--------------------------------------------- 
end;

%------------------作图--------------------------
figure(1);
plot(tX,tY,'-.g*',rX,rY,'-ro');
legend('true','estimate');
