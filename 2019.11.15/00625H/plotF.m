clear;
H=1;
L=1;
W=0.1;
m=0.5;
x = 0:1/199:1;
y = x;
[XX,YY] = meshgrid(x,y);
for i = 1:200
    for j = 1:200
        bool = sqrt((x(i)-1/2)^2+(y(j)-1/2)^2);
        F2_1= 58.8*(y(j)-1/2);
        F2_2 = (1-(tanh(24*(bool-1/16)))^2)/bool;
        F(i,j) = F2_1*F2_2;
    end
end
figure(1)
surf(XX,YY,F)
figure(2)
plot(x,F(100,:))
hold on
plot(x,F(:,100))
legend('x=0.5','y=0.5')