%% 1.3
% clear;clc
% 用函数
tic; % 开始计时
% 直接数值积分
n = 8;
% 定义积分函数
f = @(x) (x.^n)./(x+5); % 定义匿名函数

% 定义积分区间
a = 0; % 积分下限
b = 1; % 积分上限

% 使用MATLAB内置的数值积分函数进行积分
integral_result = integral(f, a, b);
time1 = toc; % 停止计时并获取时间
% 打印第一段代码的执行时间
disp(['函数的执行时间: ', num2str(time1), ' 秒']);
integral_result

% 计算梯形法则的数值积分
tic;
dt = 0.0001;
x=(0:dt:1);n=8;
% sn1 = zeros(1,101);
sn1 = 0;
for i = 1:length(x)
    sn1 = x(i)^n/(x(i)+5)*dt + sn1;
end
time1 = toc; % 停止计时并获取时间
% 打印第一段代码的执行时间
disp(['计算梯形法则的数值积分执行时间: ', num2str(time1), ' 秒']);
sn1

%% 稳定数值算法
% clear;clc
tic; % 开始计时
n = 9;

sn_ul = 1/(5*(n+1));   % 上限
sn_ll = 1/(6*(n+1));   % 下限

% 对SN上下限平均，估计一个 SN 值
sn = 1/2*(sn_ul+sn_ll);
sn_1 = 1/(5*n)-1/5*sn;
time2 = toc; % 停止计时并获取时间
% 打印第二段代码的执行时间
disp(['稳定数值算法的执行时间: ', num2str(time2), ' 秒']);
sn_1
