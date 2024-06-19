%% 1.1
% 解析解
clear;clc
Nx=180;%水平方向节点数
K=1.e-6;%热传导系数
A=5;%初始温度分布的振幅
L=1;%棒子的长度 
r=pi/L;
dx=L/(Nx-1);%网格的大小
x=linspace(0,L,Nx);% 水平网格点所在位置
T0=A*exp(-1*K*r^2*1)*cos(r*x);
t = linspace(0,1,100000); % 时间

fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = '解析解1_1.gif'; % 定义将要创建的动画文件的名称
for n = 1:numel(t)
    t = n;
    Txt=A*exp(-1*K*r^2*t)*cos(r*x);
    if mod(t,5000)==0  % 画图的部分
        area(x,Txt)
        title(['第',num2str(t),'步']);
        axis([0 L -A A]);
        drawnow % 强制渲染画布
        % 强制MATLAB立即渲染当前的画布，以便可以捕获当前状态。
        % 保存每一帧为 gif 图像
        frame = getframe(fig);
        % 捕获当前图形窗口的一帧。
        im = frame2im(frame);
        % 将捕获的帧转换为图像对象。
        [imind,cm] = rgb2ind(im,256);
        % 将RGB图像转换为索引图像，并生成一个颜色映射表。
        if n == 5000
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

% FTCS
clear all;close all;clc
L=1;%棒子的长度
Nx=180;%水平方向节点数
x=linspace(0,L,Nx);%水平网格点所在位置
dx=L/(Nx-1);%网格的大小
K=1.e-7;%热传导系数
A=5;%初始温度分布的振幅
dt=10;%时间步长，即每一步10s
mu=(K*dt)/(dx^2);%u

%-------------给定初值t=0时刻温度分布------------
% T0=0.5*A*(1+cos(1*x/L*pi));
% T0=A*exp(-(x-0.5*L).^4*1.e5);
T0=A*cos(pi/L*x);
% area(x,T0);
T1=zeros(size(T0));%用于存储t=1时刻温度计算值

%%CTS:central time central space:
xi=2:Nx-1;%计算从第二个网格到倒数第二个网格，第一个网格和最后一个网格的数值由边界条件给定

area(x,T0);
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = 'FTCS1_1.gif'; % 定义将要创建的动画文件的名称

count=0;
while count<100000
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
    
    if mod(count,5000)==0%画图的部分
        area(x,T0)
        title(['第',num2str(count),'步']);
        axis([0 L -A A]);
        drawnow % 强制渲染画布
        % 强制MATLAB立即渲染当前的画布，以便可以捕获当前状态。
        % 保存每一帧为 gif 图像
        frame = getframe(fig);
        % 捕获当前图形窗口的一帧。
        im = frame2im(frame);
        % 将捕获的帧转换为图像对象。
        [imind,cm] = rgb2ind(im,256);
        % 将RGB图像转换为索引图像，并生成一个颜色映射表。
        if count == 5000
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
