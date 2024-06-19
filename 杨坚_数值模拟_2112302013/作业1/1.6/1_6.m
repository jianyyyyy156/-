%% 1.6
%% a FTCS
clear all;close all;clc

tic; % 开始计时
L=1;% 棒子的长度
Nx=180;%水平方向节点数
x=linspace(0,L,Nx);%水平网格点所在位置
dx=L/(Nx-1);%网格的大小
K=1.e-2;%热传导系数
u=1;%u=1cm/s
A=5;%初始温度分布的振幅
dt=0.2*dx/u; %时间步长
% dt = 0.5*dt;
mu=(K*dt)/(dx^2);%u
lambda=u*dt/dx; % 0.2<1

T0=A*exp(-(x-0.5*L).^2/(0.1*L).^2);
% area(x,T0);
T1=zeros(size(T0));%用于存储t=1时刻温度计算值

%%CTS:central time central space:
xi=2:Nx-1;%计算从第二个网格到倒数第二个网格，第一个网格和最后一个网格的数值由边界条件给定

TO(1) = 0;
TO(end) = 0;

% area(x,T0);
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = 'FTCS_1_6_a.gif'; % 定义将要创建的动画文件的名称

count=0;plt = 50;
while count<5000
    count=count+1;
    % FTCS
	T1(xi)=(1-2*mu)*T0(xi)+(mu-lambda/2)*T0(xi+1)+(mu+lambda/2)*T0(xi-1);
    %将计算好的T值赋值到未计算部分
%     T1(1)=TO(1);
%     T1(end)=TO(end);
    %将计算好的T值赋值到未计算部分
    T1(1)=T1(2);
    T1(end)=T1(end-1);
    
    T0=T1;%用于迭代。t=1时刻的温度成了下一次计算需要的t=0时刻的温度
    
    if mod(count,plt)==0%画图的部分
        area(x,T0,-A)
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
time1 = toc; % 停止计时并获取时间
disp(['FTCS执行时间: ', num2str(time1), ' 秒']);

%% b BTCS
clear all;close all;clc

tic; % 开始计时
L=1;%棒子的长度
Nx=180;%水平方向节点数
x=linspace(0,L,Nx);%水平网格点所在位置
dx=L/(Nx-1);%网格的大小
u=1;%u=1cm/s
K=1.e-2;%热传导系数
A=5;%初始温度分布的振幅
dt=0.2*dx/u; %时间步长
dt = 10*dt;%放大10倍
mu=(K*dt)/(dx^2);%u
lambda=u*dt/dx; % 2>1

%-------------给定初值t=0时刻温度分布------------
T0=A*exp(-(x-0.5*L).^2/(0.1*L).^2);
% area(x,T0);
T1=zeros(size(T0));%用于存储t=1时刻温度计算值

%%CTS:central time central space:
xi=2:Nx-1;%计算从第二个网格到倒数第二个网格，第一个网格和最后一个网格的数值由边界条件给定

D=zeros(Nx-2,Nx-2);
%差分矩阵
for i=1:Nx-2
    D(i,i)= - (1+2*mu);
    if i<Nx-2
        D(i,i+1)=mu-lambda/2;
    end
    if i>1
        D(i,i-1)=mu+lambda/2;
    end
end
% D(1,1)=-(1+mu);
% D(end,end)=-(1+mu);
D_inv= D^-1;

% area(x,T0);
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = 'BTCS_1_6_b.gif'; % 定义将要创建的动画文件的名称

count=0;plt = 10;
while count<500
    count=count+1;
    % BTCS
    %  再迭代，用T0计算 T1 (显式)
    T1(xi)= - D_inv*T0(xi)';
    % 赋值边界条件: dT/dx(x=0,x=L)=0，棒子的两端始终没有温度梯度
    T1(1)=T0(1)+dt/dx^2*K*T1(2) ;
    
    T1(1)=T1(2);
    T1(end)=T1(end-1);

    T0=T1; %最新计算的t=1时刻温度，用于下次计算，作为t=0时刻的温度
    
    if mod(count,plt)==0%画图的部分
        area(x,T0,-A)
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
time1 = toc; % 停止计时并获取时间
disp(['BTCS执行时间: ', num2str(time1), ' 秒']);
