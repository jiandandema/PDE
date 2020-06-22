clc, clear, close all;
Nx = 200;
Ny = 200;
AQ = 9;
BQ = 6;
rate = 50;

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
max = max(max(u1));
min = min(min(u1));