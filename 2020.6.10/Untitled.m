clc, clear, close all;                
Nx = 200;
Ny = 200;
AQ = 9;
BQ = 6;

xl = 0 : 1 / (Nx - 1) : 1;
yl = 0 : 1 / (Ny - 1) : 1;

fid1 = fopen('./ex741F.dat', 'r'); 
%fid1 = fopen('./Vorticity2001out.dat', 'r');

for j = 1 : Ny
  for i = 1 : Nx
    for k = 1 : AQ

      u(j, i, k) = fscanf(fid1, '%g', 1);
    end
         u1(j, i) = u(j, i, 1);
         u2(j, i) = u(j, i, 2);
  end
end
fclose(fid1);
surf(xl, yl, u2); 
%三维曲面�?%contourf(xl, yl, u1, 120); %等高线图
%quiver(xl, yl, f1, f2); %离散箭头矢量�?%streamline(xl, yl, f1, f2, 0.0, 0.52);