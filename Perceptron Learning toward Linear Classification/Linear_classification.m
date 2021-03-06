close all;
% 数据加载
x1 = load('F:\学习文件\航线 大三下\模式识别与机器学习\15031145崔雪建Experiment 1\x1.txt');
x2 = load('F:\学习文件\航线 大三下\模式识别与机器学习\15031145崔雪建Experiment 1\x2.txt');
for i = 1:100
    r1(i) = x1(i,1);
    r2(i) = x1(i,2);
    r3(i) = x2(i,1);
    r4(i) = x2(i,2);
end
figure(1);
plot(r1,r2,'+',r3,r4,'*'); %原始数据点绘图
hold on;

x1(:,3) = 1; %数据扩维
x2(:,3) = 1;
% Initialize
w1 = rand(1,3);
w = [1,1,1];
x1 = x1.';
x2 = x2.';
m1 = 1;
m2 = -1;
p = 0.1; %学习率
n = 0; %迭代次数
s = 1; %循环标志符
x1 = x1*m1;
x2 = x2*m2;
J = 0; %代价函数

while s
    w1 = w;
    j = [0;0;0];
    J = 0;
    for i = 1:100
        y1 = w1*x1(:,i);
        if y1>0
            w = w1;
        else
            j = j+x1(:,i);
            J = J-y1;
        end
    end
    for i = 1:100
        y2 = w1*x2(:,i);
        if y2>0
            w = w1;
        else
            j = j+x2(:,i);
            J = J-y2;
        end
    end
    if J==0
        s = 0; %代价函数为0时终止迭代
    else
        j = j.';
        w1 = w+p*j;
        p = p+1;
        n = n+1;
    end
    w = w1; %更新权向量
end

x = linspace(0,10,5000);
y = (-w(1)/w(2)*x-w(3)/w(2));
plot(x,y,'r');
disp(n); %显示迭代次数
axis([0,8,0,5]);
