clear;
H=1;
L=1;
rate = 999;
W=0.1;
m = 0.5;
x = 0:1/199:1;
y=x;
for i = 1:100
    k = tanh((x(i) - 1/3)/W);
    y(i) = (rate-1)/2*(k + 1) + 1;
end
for i = 101:200
    k = tanh((x(i) - 2/3)/W);
    y(i) = (rate-1)/2*(-k + 1) + 1;
end
y = 3*y/10000;
plot(x,y)

fid1 = fopen('./extao.dat', 'r');
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