close all
clear all
clc


a = [2,2^0.5,2^0.25];
u = [0.25,0.5,1,1.5,1.75];
for a_i = 1:length(a)
    aa = a(a_i);
    disp('--------------------------------------------------')
    for u_i = 1:length(u)
        uu = u(u_i);
        [A0,B0] = get_A0_B0(aa,uu);
        fprintf('a=%.2f, u0=%.2f, A0=%.3f, B0=%.3f, B0/A0=%.3f \n',aa,uu,A0,B0,B0/A0);
    end
end
        



% 墨西哥草帽小波函数
function result = Psi(t) 
result = 2/sqrt(3)*power(pi,-1/4)*(t.^2-1).*exp(-t.^2/2);
end

% 墨西哥草帽小波傅里叶变换
function result = Psi_hat(w)
result = -sqrt(8)*power(pi,1/4)*w.^2/sqrt(3).*exp(-w.^2/2);
end

% 定义theta()函数，对应书中的（5.78）式
function result = theta(a,ksi)
SUM = zeros(1,numel(1:0.001:a));
w_i = 0;
for w=1:0.001:a
    w_i = w_i + 1;
    sum = 0;
    for j=-50:50
        sum = sum + abs(Psi_hat(a^j*w))*abs(Psi_hat(a^j*w+ksi));
    end
    SUM(w_i) = sum;
end
result = max(SUM);
end

% 定义delta，对应书中的（5.78）式
function result = delta(a,u0)
sum = 0;
for k =-50:50
    if k==0
        continue;
    end
    sum = sum + power(theta(a,2*pi*k/u0)*theta(a,-2*k*pi/u0),1/2);
end
result = sum;
end

% 计算A0和B0，对应书中的（5.79）和（5.80）式
function [A0,B0] = get_A0_B0(a,u0)
SUM = zeros(1,numel(1:0.001:a));
w_i = 0;
for w=1:0.001:a
    w_i = w_i + 1;
    sum = 0;
    for j=-10:10
        sum = sum + abs(Psi_hat(a^j*w))^2;
    end
    SUM(w_i) = sum;
end
A = min(SUM);
B = max(SUM);
deltaa = delta(a,u0);
A0 = (A-deltaa)/u0;
B0 = (B+deltaa)/u0;
end

