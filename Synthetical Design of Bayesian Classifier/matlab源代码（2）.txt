N = 80;
mua = [0,2];
mub = [0,0];
muc = [3,0];
sigmaa = [1,0;0,1];
sigmab = [2,0;0,2];
sigmac = [1.5,0;0,1.5];
a = mvnrnd(mua,sigmaa,N);
b = mvnrnd(mub,sigmab,N);
c = mvnrnd(muc,sigmac,N);
figure(1)
scatter(a(:,1),a(:,2),5,'r'); %散点图
hold on
scatter(b(:,1),b(:,2),5,'y');
hold on
scatter(c(:,1),c(:,2),5,'b');

mu1 = mean(a); %a的均值
sigma1 = cov(a(:,1),a(:,2)); %a的方差
P1 = 0.2; %a的先验概率
mu2 = mean(b);
sigma2 = cov(b(:,1),b(:,2));
P2 = 0.3;
mu3 = mean(c);
sigma3 = cov(c(:,1),c(:,2));
P3 = 0.5;

%x = unifrnd(-10,10,1000,2);
x1 = -5:0.1:10;
x2 = -5:0.1:10;
d = size(x1);
N = d(2);
d = d(1)*2;

p1 = zeros(N,N);
for i = 1:N
    for j = 1:N
        p1(i,j) = (1/((2*pi)^(d/2)*det(sigma1)^0.5))*exp(-0.5*([x1(i),x2(j)]-mu1)*inv(sigma1)*([x1(i),x2(j)]-mu1)');
    end
end

p2 = zeros(N,N);
for i = 1:N
    for j = 1:N
        p2(i,j) = (1/((2*pi)^(d/2)*det(sigma2)^0.5))*exp(-0.5*([x1(i),x2(j)]-mu2)*inv(sigma2)*([x1(i),x2(j)]-mu2)');
    end
end

p3 = zeros(N,N);
for i = 1:N
    for j = 1:N
        p3(i,j) = (1/((2*pi)^(d/2)*det(sigma3)^0.5))*exp(-0.5*([x1(i),x2(j)]-mu3)*inv(sigma3)*([x1(i),x2(j)]-mu3)');
    end
end

figure(2)
subplot(311),p = meshz(x1,x2,p1);
set(p,'FaceColor','white','EdgeColor','red');
hold on
subplot(311),p = meshz(x1,x2,p2);
set(p,'FaceColor','white','EdgeColor','yellow');
hold on
subplot(311),p = meshz(x1,x2,p3);
set(p,'FaceColor','white','EdgeColor','blue');
legend('P(X|ω1)','P(X|ω2)','P(X|ω3)'),xlabel('x1'),ylabel('x2'),zlabel('P')
%----------------类条件概率----------------

p_1 = p1*P1./(p1*P1+p2*P2+p3*P3); 
p_2 = p2*P2./(p1*P1+p2*P2+p3*P3);
p_3 = p3*P3./(p1*P1+p2*P2+p3*P3);

subplot(312),p = meshz(x1,x2,p_1);
set(p,'FaceColor','white','EdgeColor','red');
hold on
subplot(312),p = meshz(x1,x2,p_2);
set(p,'FaceColor','white','EdgeColor','yellow');
hold on
subplot(312),p = meshz(x1,x2,p_3);
set(p,'FaceColor','white','EdgeColor','blue');
legend('P(ω1|X)','P(ω2|X)','P(ω3|X)'),xlabel('x1'),ylabel('x2'),zlabel('P');
%----------------后验概率----------------

cost = [0,10,4;10,0,20;8,6,0]; %风险设置
cost1 = cost(1,2)*p_2+cost(1,3)*p_3;
cost2 = cost(2,1)*p_1+cost(2,3)*p_3;
cost3 = cost(3,1)*p_1+cost(3,2)*p_2;

subplot(313),p = meshz(x1,x2,1./(cost1+1));
set(p,'FaceColor','white','EdgeColor','red');
hold on
subplot(313),p = meshz(x1,x2,1./(cost2+1));
set(p,'FaceColor','white','EdgeColor','yellow');
hold on
subplot(313),p = meshz(x1,x2,1./(cost3+1));
set(p,'FaceColor','white','EdgeColor','blue');
legend('1/(cost1+1)','1/(cost2+1)','1/(cost3+1)'),xlabel('x1'),ylabel('x2');

x = unifrnd(-5,10,500,2);
y1 = zeros(500,1);
y2 = zeros(500,1);
v = [1,2,3];
for i = 1:500
    px1 = (1/((2*pi)^(d/2)*det(sigma1)^0.5))*exp(-0.5*(x(i,:)-mu1)*inv(sigma1)*(x(i,:)-mu1)');
    px2 = (1/((2*pi)^(d/2)*det(sigma2)^0.5))*exp(-0.5*(x(i,:)-mu2)*inv(sigma2)*(x(i,:)-mu2)');
    px3 = (1/((2*pi)^(d/2)*det(sigma3)^0.5))*exp(-0.5*(x(i,:)-mu3)*inv(sigma3)*(x(i,:)-mu3)');
    px = [px1,px2,px3];
    costx1 = cost(1,2)*px2+cost(1,3)*px3;
    costx2 = cost(2,1)*px1+cost(2,3)*px3;
    costx3 = cost(3,1)*px1+cost(3,2)*px2;
    costx = [costx1,costx2,costx3];
    for j = 1:3
        if px(j) == max(px)
            y1(i) = v(j);
        end
        if costx(j) == min(costx)
            y2(i) = v(j);
        end
    end
end
    figure(3)
    for i = 1:500 
        switch y1(i)
            case 1
                subplot(121),scatter(x(i,1),x(i,2),5,'r');
            case 2
                subplot(121),scatter(x(i,1),x(i,2),5,'y');
            case 3
                subplot(121),scatter(x(i,1),x(i,2),5,'b');
            otherwise
        end
        hold on 
    end
    for i = 1:500 
        switch y2(i)
            case 1
                subplot(122),scatter(x(i,1),x(i,2),5,'r');
            case 2
                subplot(122),scatter(x(i,1),x(i,2),5,'y');
            case 3
                subplot(122),scatter(x(i,1),x(i,2),5,'b');
            otherwise
        end
        hold on 
    end
    




