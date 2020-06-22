W = 0.01;
rate = 50;
eta_l = 0.002;
l = 0:1/199:1;

for j = 1 : 200
      if (j <= 100)
        n = tanh((l(j)  - 1.0 / 3.0) / W);
        phi = 0.5 * (-n + 1);
        nu = ((rate - 1) * phi + 1) * eta_l;
        x(j) = nu;
      else
        n = tanh((l(j) - 2.0 / 3.0) / W);
        phi = 0.5 * (n + 1);
        nu = ((rate - 1) * phi + 1) * eta_l;
        x(j) = nu;
      end
end
fid1 = fopen('./exnu.dat', 'r');
for j = 1 : 200
    for i = 1 : 200
        for k = 1 : 9
            u(j, i, k) = fscanf(fid1, '%g', 1);
        end
        u1(j, i) = u(j, i, 1);
    end
end
fclose(fid1);
figure(1)
plot(l,x)
hold on
plot(l', u1(:,1));
legend('精确','数值')
title('tanh图像');


for j = 1 : 200
      if (j <= 100)
          n = tanh((l(j) - 1.0 / 3.0) / W);
          y(j) = -0.5 * (rate - 1) * eta_l * (1 - n * n) / W;
      else
          n = tanh((l(j) - 2.0 / 3.0) / W);
          y(j) = 0.5 * (rate - 1) * eta_l * (1 - n * n) / W;
      end
end
fid1 = fopen('./exgrad_nu.dat', 'r');
for j = 1 : 200
    for i = 1 : 200
        for k = 1 : 9
            u(j, i, k) = fscanf(fid1, '%g', 1);
        end
        u1(j, i) = u(j, i, 2);
    end
end
fclose(fid1);
figure(2)
plot(l,y)
hold on
plot(l', u1(:,1));
legend('精确','数值')
title('tanh导数图像');