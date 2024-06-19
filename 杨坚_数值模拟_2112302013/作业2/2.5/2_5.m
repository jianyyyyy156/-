clc;clear;close all
% 参数设置
f=1.e-4;   % Coriolis parameter s^-1
bta=1.e-7; % Ekman friction, s^-1
g=9.8;     % gravitational acceleration
L=1.e6;     % dimension of the domain in meter
Nx=200;    % grid resolution
Ny=Nx;
dx=L/Nx;   % grid interval
xc=0.5*Nx; % location of initial perturbation

% initialize the matrix
u=zeros(Nx,Ny); % u(n)
v=zeros(Nx,Ny); % v(n)
u1=u;           % u(n+1)
v1=v;           % v(n+1)
eta=zeros(Nx,Ny); % eta(n)
eta1=eta;         %eta(n+1)
h=zeros(Nx,Ny);   %bathmetry

[X,Y]=meshgrid(1:Nx,1:Ny);
H0=2000; %        mean water depth,
% h=H0*(1-0.8*exp( -((X-xc*1.3).^2+(Y-xc*1.3).^2 )/200 )  ); % bathmetry
h=ones(Nx,Ny)*H0;
% figure;mesh(-h);colorbar;title('bathmetry, m')

% eta=0.2* sin(X/Nx*2*pi*3+Y/Nx*2*pi*2).*exp( (-(X-xc).^2-(Y-xc).^2)/12^2); %initial surf ssh
eta=0.05*exp( -(X-xc).^2/4^2-(Y-xc).^2/8^2);
% eta(:,end)=0;%continent
% eta(1,:)=0;%continent

dt=min([dx/sqrt(g*max(h(:))),bta/f^2]); %time step, s
% dt=10;%seconds
ns=max([10/dt, 10]); %output frequency
Coe_uv=1./(1+bta*dt); %coefficient before the u v equations

mu=sqrt(9.8*H0)*dt/dx;

ts=sum(eta(:)) ;%used for check mass conservation;
disp(['initial mass is ',num2str(ts),'!']);

fig = figure; % 创建一个新的 figure
axis tight manual 
filename = '线性H-N开边界.gif'; 

count=0;T=0;
while count*dt<=3*3600
    count=count+1;
    T=T+dt;
        
    j=2:Nx-1; %eastward
    k=1:Ny-1; %northward
    u1(k,j)=u(k,j)+f*0.25*dt*(v(k,j)+v(k+1,j)+v(k,j-1)+v(k+1,j-1))-g*dt/dx*(eta(k,j)-eta(k,j-1));
    u1(k,j)=u1(k,j)*Coe_uv;
    
    j=Nx;    
    u1(k,j)=(1-mu).*u(k,j)+mu.*u(k,j-1);

    j=1:Nx-1;
    k=2:Ny-1;
    v1(k,j)=v(k,j)-f*0.25*dt*(u(k,j)+u(k,j+1)+u(k-1,j)+u(k-1,j+1))-g*dt/dx*(eta(k,j)-eta(k-1,j));
    v1(k,j)=v1(k,j)*Coe_uv;
    
    j=Nx;
    v1(k,j)=(1-mu).*v(k,j)+mu.*v(k,j-1);
    
    j=1:Nx-1;
    k=1:Ny-1;
    eta1(k,j)=eta(k,j)-dt/dx*( ......
        0.5*(h(k+1,j+1)+h(k,j+1)).*u1(k,j+1)-0.5*(h(k+1,j)+h(k,j)).*u1(k,j) +.....
        0.5*(h(k+1,j)+h(k+1,j+1)).*v1(k+1,j)-0.5*(h(k,j)+h(k,j+1)).*v1(k,j));
    
    j=Nx;
    eta1(k,j)=(1-mu).*eta(k,j)+mu.*eta(k,j-1);
    
    u=u1;
    v=v1;
    eta=eta1;
   
    ts=[ts sum(eta1(:))];
    
    if mod(count,ns)==0%画图的部分
        mesh((1:Nx)*dx*1.e-3,(1:Nx)*dx*1.e-3,eta);
        set(gca,'zlim',[-0.02  0.05]);colorbar
        zlabel('\eta, m')
        title([num2str(round(count*dt)),' s'])
        xlabel('x, km')
        ylabel('y, km')
        zlabel('\eta, m')
        
        drawnow % 强制渲染画布
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if count == ns
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
        end
    end
end

figure;
plot( ts,'b-')
title('time series of mass')
xlabel('steps')

figure('position',[10,10,800,400])
subplot(121)
imagesc(u)
title('u')
subplot(122)
imagesc(v)
title('v')
