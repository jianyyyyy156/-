clc;clear;close all
% 参数设置
f=5.e-3;     % Coriolis parameter s^-1
bta=5.e-3;  % Ekman friction, s^-1
g=9.8;

L=1.e6;% dimension of the domain in meter
Nx=100;
Ny=Nx;
dx=L/Nx;
xc=0.5*Nx;

% initialize the matrix
u=zeros(Nx,Ny);
v=zeros(Nx,Ny);
u1=u;
v1=v;

h=zeros(Nx,Ny);
[X,Y]=meshgrid(1:Nx,1:Ny);
H0=600; %mean water depth,
h=H0-500*exp( (-(X-xc*1.6).^2-(Y-xc*1.6).^2)/15^2); %add bathmetry
% figure;imagesc(h);colorbar;title('bathmetry, m')


% dt=min([bta/f^2,dx/sqrt(g*max(h(:)))])*0.5; %time step, s
dt=bta/f^2*0.5; %time step, s

Coe_uv=1./(1+bta*dt); %coefficient before the u v equations

eta=zeros(Nx,Ny);
eta1=eta;
% eta=0.2* sin(X/Nx*2*pi*3+Y/Nx*2*pi*2).*exp( (-(X-xc).^2-(Y-xc).^2)/12^2); %initial surf ssh
eta=0.2*exp( (-(X-xc).^2-(Y-xc).^2)/4^2); %initial surf ssh
% eta(:,end)=0;%continent
% eta(1,:)=0;%continent
% figure;mesh(eta);colorbar;title('initial vortex ssh, m')

ts=sum(eta(:)) ;%used for check mass conservation;
disp(['initial mass is ',num2str(ts),'!']);

fig = figure; % 创建一个新的 figure
axis tight manual 
filename = '线性H-N.gif'; 

count=0;plt = 5;
while count<400
    count=count+1;
    j=2:Nx-1; %eastward
    k=2:Ny; %southward
    u1(k,j)=u(k,j)+f*0.25*dt*(v(k-1,j)+v(k-1,j-1)+v(k,j)+v(k,j-1))-g*dt/dx*(eta(k,j)-eta(k,j-1));
    u1(k,j)=u1(k,j)*Coe_uv;
    % boundary conditions for u is implied in the above j k indices so no need
    % to explicitly define the boundary u1 

    j=1:Nx-1; %eastward
    k=2:Ny-1; %southward
    v1(k,j)=v(k,j)-f*0.25*dt*(u(k,j)+u(k,j+1)+u(k+1,j)+u(k+1,j+1))-g*dt/dx*(eta(k,j)-eta(k+1,j));
    v1(k,j)=v1(k,j)*Coe_uv;


    j=1:Nx-1; %eastward
    k=2:Ny; %southward

    % consider a varying topology
    eta1(k,j)=eta(k,j)-dt/dx*( ......
                                        0.5*(h(k-1,j+1)+h(k,j+1)).*u1(k,j+1)- 0.5*(h(k-1,j)+h(k,j)).*u1(k,j) +.....
                                        0.5*(h(k-1,j)+h(k-1,j+1)).*v1(k-1,j)-  0.5*(h(k,j)+h(k,j+1)).*v1(k,j)   .....
         );

    u=u1;
    v=v1;
    % Now ,u1 v1 could be used to restore other vars, 

    % lets include nonlinearity in the mass conservation eq.
    %  use u1 v1 to restore eta*u1 eta*v1
    j=2:Nx-1; %eastward
    k=2:Ny; %southward
    u1(k,j)=u1(k,j).*(eta(k,j)+eta(k,j-1))*0.5; % u1=u1*eta

    j=1:Nx-1; %eastward
    k=2:Ny-1; %southward
    v1(k,j)=v1(k,j).*(eta(k,j)+eta(k+1,j))*0.5; % v1=v1*eta

    j=1:Nx-1; %eastward
    k=2:Ny;
    eta1(k,j)=eta1(k,j)-dt/dx*( ......
        u1(k,j+1)- u1(k,j) +.....
        v1(k-1,j)- v1(k,j)   .....
        );
    ts=[ts sum(eta1(:))];
    eta=eta1;
    if mod(count,plt)==0%画图的部分
        imagesc(eta);
        title(['第',num2str(count),'步'])
        colorbar;
        
        drawnow % 强制渲染画布
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


figure;
plot(0:count,ts,'bo-')
title('time series of mass')
xlabel('steps')

figure('position',[10,10,800,400])
subplot(1,2,1)
imagesc(u)
title('u')
subplot(1,2,2)
imagesc(v)
title('v')
