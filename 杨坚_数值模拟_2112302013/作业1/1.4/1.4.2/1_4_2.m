%% 1.4.2
% CTCS
clear all;close all;clc

tic; % 开始计时
L=1; %长度
Nx=196; %水平方向节点数
x=linspace(0,L,Nx);
dx=L/(Nx-1);
u=1;%U=1m/s
A=5;
dt=0.2*dx/u; 
lambda=u*dt/dx;

%-------------给定初值t=0时刻分布------------
T0=A*exp(-(x-0.5*L).^2/(0.1*L).^2);
% area(x,T0);

T1=zeros(size(T0));
T2=T1;
%%CTS:central time central space:
xi=2:Nx-1;%计算从第二个网格到倒数第二个网格，第一个网格和最后一个网格的数值由边界条件给定

% leap-frog第一步启动时候，需要从0->1
T1(xi)= T0(xi)- 0.5*lambda*( T0(xi+1)-T0(xi-1)  );
% % 边界条件
T1(1)= T0(end);
T1(end)= T0(1);

% area(x,T0);
fig = figure; % 创建一个新的 figure
% 创建一个新的图形窗口。
axis tight manual % 设置坐标轴
% 设置坐标轴属性，使得坐标轴的边界紧贴图形，并且手动设置，不自动调整。
filename = 'CTCS_1_4_2.gif'; % 定义将要创建的动画文件的名称

count=0;plt = 50;
while count<5000
    count=count+1;
    % CTCS
    T2(xi)= T0(xi)- lambda*( T1(xi+1)-T1(xi-1)  );
    % %boundary condition
    T2(1)= T0(1)- lambda*( T1(2)-T1(end)  );
    T2(end)= T0(end)- lambda*( T1(1)-T1(end-1)  );

    T0=T1;
    T1=T2;
    
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
disp(['CTCS执行时间: ', num2str(time1), ' 秒']);

