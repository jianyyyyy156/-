%% 2.1
%% 蓝色格式
clear all;close all;clc
tic; % 开始计时
L=1;% 长度
Nx=180;%水平方向节点数
x=linspace(0,L,Nx);%水平网格点所在位置
dx=L/(Nx-1);%网格的大小
u=1;%u=1cm/s
A=5;%初始温度分布的振幅
dt=0.2*dx/u; %时间步长
dt = 0.02*dt;
r=u*dt/dx; 

%初始条件
T0=sin(2*pi*x);%2.5.12a
% T0=1.5+sin(2*pi*x);%2.5.12b

% area(x,T0);
T1=zeros(size(T0));%用于存储t=1时刻计算值

%%CTS:central time central space:
xi=2:Nx-1;%计算从第二个网格到倒数第二个网格，第一个网格和最后一个网格的数值由边界条件给定

% leap-frog第一步启动时候，需要从0->1
T1(xi)= (T1(xi-1)+T1(xi)).^2/4*r-(T1(xi+1)+T1(xi)).^2/4*r+T0(xi-1);
% % 边界条件
T1(1)= T0(end);
T1(end)= T0(1);
T2=T1;

% area(x,T0);
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = '2_1_蓝色格式_a初始条件.gif'; % 定义将要创建的动画文件的名称
um = zeros(500,1);
count=0;plt = 5;
while count<500
    count=count+1;
    % CTCS
    T2(xi)=(T1(xi-1)+T1(xi)).^2/4*r-(T1(xi+1)+T1(xi)).^2/4*r+T0(xi-1);
    %boundary condition
    T2(1)= T1(end);
    T2(end)= T2(1);
    
    T0=T1;
    T1=T2;
    um(count) = 1/2*sum(T2)^2;
    
    if mod(count,plt)==0%画图的部分
        area(x,T2,-A)
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
figure
plot(um)
string = {'a初始条件动能变化曲线'};
title(string)
disp(['执行时间: ', num2str(time1), ' 秒']);

%% 绿色格式
clear all;close all;clc
tic; % 开始计时
L=1;% 长度
Nx=180;%水平方向节点数
x=linspace(0,L,Nx);%水平网格点所在位置
dx=L/(Nx-1);%网格的大小
u=1;%u=1cm/s
A=5;%初始温度分布的振幅
dt=0.2*dx/u; %时间步长
dt = 0.02*dt;
r=u*dt/dx; 

%初始条件
T0=sin(2*pi*x);%2.5.12a
% T0=1.5+sin(2*pi*x);%2.5.12b

% area(x,T0);
T1=zeros(size(T0));%用于存储t=1时刻计算值
T2=T1;
%%CTS:central time central space:
xi=2:Nx-1;%计算从第二个网格到倒数第二个网格，第一个网格和最后一个网格的数值由边界条件给定

% leap-frog第一步启动时候，需要从0->1
T2(xi)= T0(xi)-1/12.*r.*((T1(xi)+T2(xi)).*((T1(xi+1)+T2(xi+1))-(T1(xi-1)+T2(xi-1)))+(T1(xi+1)+T2(xi+1)).^2-(T1(xi-1)+T2(xi-1)).^2);
% % 边界条件
T1(1)= T0(end);
T1(end)= T0(1);


% area(x,T0);
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = '2_1_绿色格式_a初始条件.gif'; % 定义将要创建的动画文件的名称

count=0;plt = 5;
while count<500
    count=count+1;
    % CTCS
    T2(xi)=T0(xi)-1/12.*r.*((T1(xi)+T2(xi)).*((T1(xi+1)+T2(xi+1))-(T1(xi-1)+T2(xi-1)))+(T1(xi+1)+T2(xi+1)).^2-(T1(xi-1)+T2(xi-1)).^2);
    %boundary condition
    T2(1)= T1(end);
    T2(end)= T2(1);
    
    T0=T1;
    T1=T2;
    um(count) = 1/2*sum(T2)^2;
    
    if mod(count,plt)==0%画图的部分
        area(x,T2,-A)
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
figure
plot(um)
string = {'b初始条件动能变化曲线'};
title(string)
disp(['执行时间: ', num2str(time1), ' 秒']);
