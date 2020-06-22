clear;
H=1;
rate = 50;
lidvelocity = 0.1;
Re = 500;
L=1;
W=0.1;
m=0.5;
x = 0:1/199:1;
for i = 1:99
    k = tanh((x(i) - H/3)/W);
    y(i) = (rate-1)/2*(k + 1) + 1;
end
for i = 100:200
    k = tanh((x(i) - 2*H/3)/W);
    y(i) = (rate-1)/2*(-k + 1) + 1;
end
y = 8.*lidvelocity.*y;
for i = 1:100
    n = tanh((x(i) - H/3)/W);
    n = 1 - n*n;
    grad_nu = (rate-1)/2.0 / W * n;
    grad_u = 4 * lidvelocity * (1 - 2.0 * x(i));
    y(i) = y(i) - grad_nu*grad_u;
end
for i = 101:200
    n = tanh((x(i) - 2*H/3) / W);
    n = 1 - n*n;
    grad_nu = -(rate-1)/2.0 / W * n;
    grad_u = 4 * lidvelocity * (1 - 2.0 * x(i));
    y(i) = y(i) - grad_nu*grad_u;
end
y = y/Re;
plot(x,y)


fid1 = fopen('./exF.dat', 'r');
for j = 1 : 200
    for i = 1 : 200
        for k = 1 : 9
            u(j, i, k) = fscanf(fid1, '%g', 1);
        end
        u1(j, i) = u(j, i, 1);
    end
end
fclose(fid1);
hold on
plot(x,u1(:,1)');
ee = u1(:,1)' - y;
ee = ee';
e = norm(ee);