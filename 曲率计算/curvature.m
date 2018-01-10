function [ c,start,final ] = curvature( Lattitude,Longitude )
%% 变量声明
%curvature 计算离散点 U弦长曲率（见文献计算方法）
%参数说明：
% Lattitude：纬度数组
% Longitude：经度数组
% Heading：朝向数组
% x:以起始点为原点的坐标系下的x坐标数组
% y:以起始点为原点的坐标系下的y坐标数组
% n:经纬度数组纬度
% U:支持领域欧氏距离，输入参数，根据实际情况调节
% start:可计算曲率的起始点
% fianl:可计算曲率的终点
% c:计算曲率数组
%% 计算参数设定
c=[];
U=0.4;  %根据实际情况调节
n=size(Longitude,1);
%% 以起始点为原点计算所有点的 x , y 坐标
x = -(Longitude- Longitude(1)) * (111413*cos(Lattitude(1)) - 94*cos(3 * Lattitude(1)));
y = (Lattitude - Lattitude(1)) * 3600 * 30.864;
%% 计算离散点 U弦长曲率
for start=1:n    %找到可计算曲率的起始点
    dist=sqrt((x(start)-x(1))^2+(y(start)-y(1))^2);
    if(dist>U)
        break;
    end
end

for final=n:-1:1 %找到可计算曲率的终点
    dist=sqrt((x(final)-x(n))^2+(y(final)-y(n))^2);
    if(dist>U)
        break;
    end
end
for i=start:final
    for b=i:-1:1 %向前遍历，找到距离大于U的第一个点
        dist=sqrt((x(i)-x(b))^2+(y(i)-y(b))^2);
        if(dist>U)
            break;
        end
    end
    for f=i:n %向后遍历，找到距离大于U的第一个点
       dist=sqrt((x(i)-x(f))^2+(y(i)-y(f))^2);
       if(dist>U)
        break;
       end
    end
    b=b+1;
    f=f-1;
    % 由于采样点密度较大，这里不进行隐式精化数字曲线策略优化修正，直接计算曲率
    ci=sign((x(i) - x(b))*(y(f)- y(b)) - (x(f)- x(b))*(y(i) - y(b)))*sqrt(1-(((sqrt((x(b)-x(f))^2+(y(b)-y(f))^2))/(2*U))^2));
    c=[c ci];
end
c=[zeros(1,start-1) c zeros(1,n-final)];
c=c';
%% 处理异常
for i=start:final
    if(c(i)==0)
        c(i)=c(i-1);
    end
end
end