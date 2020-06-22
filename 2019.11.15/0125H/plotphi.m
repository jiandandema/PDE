clear;
H=1;
L=1;
W=0.1;
x = 0:0.01:1;
y=x;
phi = zeros(101,101);
for i = 1:101
    for j = 1:101
        bool = sqrt((x(i) - 1/2)^2+(y(j) - 1/2)^2);
        phi(i,j) = -tanh(2.4*(bool - 1/8)/W) + 1;
    end
end
phi = phi./2;
[XX,YY]=meshgrid(x,y);
C = XX.*YY;
figure(1)
surf(XX,YY,phi,C)
figure(2)
plot(x,phi(:,50))