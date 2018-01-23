function [y_k,y] = kalman_filter(Lattitude,Longitude,Heading)
%kalman_filter ��GPS·�����п������˲�
%   Lattitude:ά��
%   Longitude:����
%   Heading:�����
%   x,y:����任�������,����ֵ
%   y_k:k ʱ�̵����Ž�
%   Kg:����������
%����任
x = (Lattitude - Lattitude(1)) * 3600 * 30.864;
y = -(Longitude- Longitude(1)) * (111413*cos(Lattitude(1)) - 94*cos(3 * Lattitude(1)));
%�����ʼ��
w=randn(1,length(Heading));   
v=randn(1,length(Heading));
R=(std(v)).^2; 
Q=(std(w)).^2; 
P_kk(1)=0;
P_k(1)=0;
%��ʼ��״̬��
y_k=y(1);
H=1;
Kg(1)=P_kk(1)*H'/(H*P_kk(1)*H'+R);
P_k(1)=Q;
A=1;
%����
for i=2:length(Heading)
    g=y_k(i-1)+tan(Heading(i-1)/180*pi)*(x(i)-x(i-1));  %ע�����복��ת���������
    A=1;
    y_k(i)=g+Kg(i-1)*(y(i-1)-H*A*y_k(i-1));
    Kg(i)=P_k(i-1)*H'\((H*P_k(i-1)*H'+R));
    P_k(i)=A*P_kk(i-1)*A'+Q;
    P_kk(i)=P_k(i-1)-Kg(i)*H*P_kk(i-1);
end

%�쳣����
temp=y_k'-y;
for i=1:length(temp)
    if abs(temp(i))>0.3
        y_k(i)=y(i);
    end
end
%��ͼ��
plot(x,y_k,'k');
hold on
plot(x,y,'r')
end

