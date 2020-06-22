clear;
H=1;
step = 1/199;
rate = 999;
lidvelocity = 0.1;
Re = 10000;
L=1;
W=0.01;
m=0.5;
x = 0:1/199:1;
y = x;
[XX,YY] = meshgrid(x,y);
for i = 1:200
    for j = 1:200
     if (j < 100)
      n = tanh((j * step - 1.0 / 3.0) / W);
      n = 1 - n*n;
      grad_nu =  rate/(2.0 * W) * n;
      F(j) = -grad_nu * lidvelocity;
      F(j) = F(j) / Re;
     else
      n = tanh((j * step - 2.0 / 3.0) / W);
      n = 1 - n*n;
      grad_nu =  -rate/(2.0 * W) * n;
      F(j) = -grad_nu * lidvelocity;
      F(j) = F(j) / Re;
      end
    end
end
plot(y,F)