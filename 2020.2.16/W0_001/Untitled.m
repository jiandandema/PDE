clear;
syms x C1 C2
W = 0.01;
rate = 50;
phi = 0.5*(tanh((x - 1/2)/W) + 1);
eta = 1 + rate * phi - phi;
u = int(C1/eta) + C2
x = 0:0.01:1;
u1 = 49*log(tanh(100*x - 50) + 1);u1 = u1 / 10000;
u2 = 49*log(49*tanh(100*x - 50) + 51); u2 = u2 / 10000;
u3 = x/50;
u = u3 + u1 - u2 + 0.1;
plot(x,u)