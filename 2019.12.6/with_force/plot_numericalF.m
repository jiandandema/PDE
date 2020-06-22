clc, clear, close all;
Nx = 200;
Ny = 200;
AQ = 9;
BQ = 6;
rate = 50;
Re=500;
W=0.1;
lidvelocity = 0.1;

xl = 0 : 1 / (Nx - 1) : 1;
yl = 0 : 1 / (Ny - 1) : 1;

fid1 = fopen('./exF.dat', 'r');
%fid1 = fopen('./Vorticity2001out.dat', 'r');

for j = 1 : Ny
    for i = 1 : Nx
        for k = 1 : AQ
            u(j, i, k) = fscanf(fid1, '%g', 1);
        end
        u1(j, i) = u(j, i, 1);
    end
end
figure(1)
surf(xl, yl, u1');
for i = 1:200
    for j = 1:200
        if(j < 100)
            n = tanh((j * 1/199 - 1.0 / 3.0) / W);
            n = n*n - 1;
            grad_nu =  (rate-1)/(2.0 * W) * n;
            F(j,i) = -grad_nu * lidvelocity;
            F(j,i) = F(j,i) / Re;
        else
            n = tanh((j * 1/199 - 2.0 / 3.0) / W);
            n = 1 - n*n;
            grad_nu = rate/(2.0 * W) * n;
            F(j,i) = -grad_nu * lidvelocity;
            F(j,i) = F(j,i) / Re;
        end
    end
end
fclose(fid1);
figure(2)
surf(xl,yl,F);