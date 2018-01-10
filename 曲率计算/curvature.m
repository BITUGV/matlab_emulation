function [ c,start,final ] = curvature( Lattitude,Longitude )
%% ��������
%curvature ������ɢ�� U�ҳ����ʣ������׼��㷽����
%����˵����
% Lattitude��γ������
% Longitude����������
% Heading����������
% x:����ʼ��Ϊԭ�������ϵ�µ�x��������
% y:����ʼ��Ϊԭ�������ϵ�µ�y��������
% n:��γ������γ��
% U:֧������ŷ�Ͼ��룬�������������ʵ���������
% start:�ɼ������ʵ���ʼ��
% fianl:�ɼ������ʵ��յ�
% c:������������
%% ��������趨
c=[];
U=0.4;  %����ʵ���������
n=size(Longitude,1);
%% ����ʼ��Ϊԭ��������е�� x , y ����
x = -(Longitude- Longitude(1)) * (111413*cos(Lattitude(1)) - 94*cos(3 * Lattitude(1)));
y = (Lattitude - Lattitude(1)) * 3600 * 30.864;
%% ������ɢ�� U�ҳ�����
for start=1:n    %�ҵ��ɼ������ʵ���ʼ��
    dist=sqrt((x(start)-x(1))^2+(y(start)-y(1))^2);
    if(dist>U)
        break;
    end
end

for final=n:-1:1 %�ҵ��ɼ������ʵ��յ�
    dist=sqrt((x(final)-x(n))^2+(y(final)-y(n))^2);
    if(dist>U)
        break;
    end
end
for i=start:final
    for b=i:-1:1 %��ǰ�������ҵ��������U�ĵ�һ����
        dist=sqrt((x(i)-x(b))^2+(y(i)-y(b))^2);
        if(dist>U)
            break;
        end
    end
    for f=i:n %���������ҵ��������U�ĵ�һ����
       dist=sqrt((x(i)-x(f))^2+(y(i)-y(f))^2);
       if(dist>U)
        break;
       end
    end
    b=b+1;
    f=f-1;
    % ���ڲ������ܶȽϴ����ﲻ������ʽ�����������߲����Ż�������ֱ�Ӽ�������
    ci=sign((x(i) - x(b))*(y(f)- y(b)) - (x(f)- x(b))*(y(i) - y(b)))*sqrt(1-(((sqrt((x(b)-x(f))^2+(y(b)-y(f))^2))/(2*U))^2));
    c=[c ci];
end
c=[zeros(1,start-1) c zeros(1,n-final)];
c=c';
%% �����쳣
for i=start:final
    if(c(i)==0)
        c(i)=c(i-1);
    end
end
end