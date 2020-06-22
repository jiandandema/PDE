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
for i = 1:200
    for j = 1:200
        bool = sqrt((xl(i)-1/2)^2+(yl(j)-1/2)^2);
        F2_1= 58.8*(yl(j)-1/2);
        F2_2 = (1-(tanh(24*(bool-1/4)))^2)/bool;
        F(i,j) = F2_1*F2_2;
    end
end
fclose(fid1); 
figure(2)
surf(xl, yl, F/500);
e = norm(F/500 - u1');