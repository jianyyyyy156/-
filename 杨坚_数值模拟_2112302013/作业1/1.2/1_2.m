%% 1.2
% FTCS
clear all;close all;clc
L=1;%棒子的长度
Nx=180;%水平方向节点数
x=linspace(0,L,Nx);%水平网格点所在位置
dx=L/(Nx-1);%网格的大小
K=1.e-6;%热传导系数
A=5;%初始温度分布的振幅
dt=10;%时间步长，即每一步10s
mu=(K*dt)/(dx^2);%u

%-------------给定初值t=0时刻温度分布------------
T0=exp(-(x-0.5*L).^2/(0.1*L).^2);
% area(x,T0);
T1=zeros(size(T0));%用于存储t=1时刻温度计算值

%%CTS:central time central space:
xi=2:Nx-1;%计算从第二个网格到倒数第二个网格，第一个网格和最后一个网格的数值由边界条件给定

area(x,T0);
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = 'FTCS1_2.gif'; % 定义将要创建的动画文件的名称

count=0;plt = 50;
while count<5000
    count=count+1;
    % FTCS
    %T1(xi)=(1-2*mu)*T0(xi)+mu*(T0(xi+1)+T0(xi-1));
    %T1(xi)=mu.*(T0(xi+1)-2*T0(xi)+T0(xi-1)) + T0(xi);
	T1(xi)=mu.*(T0(xi+1)+T0(xi-1))+(1-2*mu).*T0(xi);
    %将计算好的T值赋值到未计算部分
    T1(1)=T1(2);
    T1(end)=T1(end-1);
    
    T0=T1;%用于迭代。t=1时刻的温度成了下一次计算需要的t=0时刻的温度
    % T1=T2;%最新算的t=2时刻的温度成了下次计算需要的t=1时刻的温度。。。如此往复
    
    if mod(count,plt)==0%画图的部分
        area(x,T0)
        title(['第',num2str(count),'步']);
        axis([0 1 0 1]);
        drawnow % 强制渲染画布
        % 强制MATLAB立即渲染当前的画布，以便可以捕获当前状态。
        % 保存每一帧为 gif 图像
        frame = getframe(fig);
        % 捕获当前图形窗口的一帧。
        im = frame2im(frame);
        % 将捕获的帧转换为图像对象。
        [imind,cm] = rgb2ind(im,256);
        % 将RGB图像转换为索引图像，并生成一个颜色映射表。
        if count == plt
            % 如果是第一帧
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'Loopcount',inf);
            % 写入第一帧到GIF文件，设置延迟时间为0.1秒，动画无限循环。
        else
            % 如果不是第一帧
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
            % 将后续帧追加到GIF文件中，同样设置延迟时间和写入模式为追加。
        end
    end
end

% 1.2.2解析解
clear;clc
Nx=180;%水平方向节点数
K=1.e-6;%热传导系数
A=5;%初始温度分布的振幅
L=1;%棒子的长度 
r=pi/L;
dx=L/(Nx-1);%网格的大小
dt=10;  % 时间间隔
% dt=0.5*dx^2/K*20; % time step required by the CFD stability condition:
mu=(K*dt)/(dx^2);%u
x=linspace(0,L,Nx);% 水平网格点所在位置

T0=exp(-(x-0.5*L).^2/(0.1*L).^2);
an=1/Nx*fft(T0)';% 可以调整振幅
kn=[0:Nx/2  , -Nx/2+1:-1]'.*2*pi/L; %L=Nx*dx
expikx= exp(1j.*(kn*x)); % 1j返回虚数单位

fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = '解析解1_2.gif'; % 定义将要创建的动画文件的名称
plt = 50;
for t = 1:5000  % 时间
%     t = n;
    % Tx = an*exp(-K*Kn.^2*t)*exp(i*Kn*x);
    T_exact = (an.*exp(-K*kn.^2*t*dt))'*expikx ;
    if mod(t,plt)==0  % 画图的部分
        area(x,real(T_exact))
        title(['第',num2str(t),'步']);
        axis([0 1 0 1]);
        drawnow % 强制渲染画布
        % 强制MATLAB立即渲染当前的画布，以便可以捕获当前状态。
        % 保存每一帧为 gif 图像
        frame = getframe(fig);
        % 捕获当前图形窗口的一帧。
        im = frame2im(frame);
        % 将捕获的帧转换为图像对象。
        [imind,cm] = rgb2ind(im,256);
        % 将RGB图像转换为索引图像，并生成一个颜色映射表。
        if t == plt
            % 如果是第一帧
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'Loopcount',inf);
            % 写入第一帧到GIF文件，设置延迟时间为0.1秒，动画无限循环。
        else
            % 如果不是第一帧
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
            % 将后续帧追加到GIF文件中，同样设置延迟时间和写入模式为追加。
        end
    end
end

