% a test for the SOR solver
clear all;close all;
clc;

tic; % 开始计时
dx=1;  % 网格间距 
N=64;  % 方形网格，每维181个点
ix=2:N-1;  % 格点
jy=2:N-1;

r_d2x=1/(2*dx);
r_dx2=1/(dx^2);
t=10000;  % 时间步
Ah=1.e-5;  % 粘性系数
Amp=1;

coeffi_Jacobian= -1/(12*dx^2);  % Arakawa 1966 Eq 45

% 涡度场，3个时间层
zeta0= zeros(N,N); 
zeta1= zeros(N,N); 
zeta2= zeros(N,N); 

% 初始化原始相对涡度场
for i=1:N
    for j=1:N
        r2=(j/N-0.35)^2*60+(i/N-0.25)^2*100;
        r3=(j/N-0.65)^2*30+(i/N-0.75)^2*100;
        zeta0(i,j)=exp(-1*r2)-exp(-1*r3);
    end
end
zeta=Amp*zeta0 ;
[xx,yy]=meshgrid(1:N,1:N);

% figure;
% imagesc(zeta0);
% title('initial vorticity'); 
% colorbar

% SOR求解初始场的psi
SOR2d;

% psi=SOR2d_check(zeta0,psi_guess );
% figure
% imagesc(psi);
% title('SOR solution of \psi for the initial \zeta field')
% colorbar

u=max(max(abs(diff(psi)/dx)));
dt=(dx/u)*0.2;

jac
zeta=zeta0+J0*dt;  % Euler 时间向前预测下一步
SOR2d;  
jac;
J1=J0;
zeta1=zeta;  

u=max(max(abs(diff(psi)/dx)));
dt=min(dx/u)*0.1;

%% 绘图
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = '2_2_涡旋.gif'; % 定义将要创建的动画文件的名称

count=0;T = 0;plt = 200;
while count<t
    count=count+1;
    T=T+dt;

    zeta2=zeta1+dt*(1.5*J1-0.5*J0) ;  % 非线性项的Adam-Bath二阶时间格式
    zeta2(jy,ix)=zeta2(jy,ix)+(zeta2(jy+1,ix)+zeta2(jy-1,ix)+zeta2(jy,ix+1)+zeta2(jy,ix-1))*Ah*r_dx2*dt; 
    diffusion_periobc;

    zeta2=zeta2/(1+4*Ah*dt*r_dx2);  
    
    zeta1=zeta2;
    zeta=zeta2;

    SOR2d;  
    jac;        
    J2=J0;
    
    J0=J1;
    J1=J2;
    
    diag_ts(count,1)=sum(sum(zeta2.^2));
    
    if mod(count,plt)==0%画图的部分
        subplot(1,2,1);
        contourf(zeta2,25,'LineStyle','none');
        colorbar;
        caxis([-1 1]*Amp*0.8);
        set(gca,'ydir','reverse')
   
        subplot(1,2,2)
        plot(diag_ts,'b*-')
        xlabel('time step')
        ylabel('enstrophy')
        pause(0.1);
        u=max(max(abs(diff(psi)/dx)));
        dt=min(dx/u)*0.1;
        
        sgtitle(['第',num2str(count),'步'],'Fontsize',15)
        drawnow 
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if count == plt
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
        end
    end
end

time1 = toc; % 停止计时并获取时间
disp(['执行时间: ', num2str(time1), ' 秒']);
