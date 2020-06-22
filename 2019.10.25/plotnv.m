clear;
H=1;
L=1;
rate = 999;
W=0.01;
m = 0.5;
x = 0:1/199:1;
y=x;
phi = zeros(101,101);
for i = 1:99
    k = tanh((x(i) - 1/3)/W);
    y(i) = (rate-1)/2*(k + 1) + 1;
end
for i = 100:200
    k = tanh((x(i) - 2/3)/W);
    y(i) = (rate-1)/2*(-k + 1) + 1;
end
plot(x,y)