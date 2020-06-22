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
surf(xl, yl, u1);
% for j = 1 : Ny
%     for i = 1 : Nx
%         r = sqrt((xl(i) - 1/2)^2 + (yl(j) - 1/2)^2);
%         phi(j,i) = 0.5*(1 - tanh(24*(r - 1/4)));
%         nu(j,i) = (rate - 1)*phi(j,i) + 1;
%         nu(j,i) = nu(j,i);
%     end
% end
% 
% tao = 3 * nu / 500;
% figure(2)
% surf(xl, yl, tao);
for j = 1 : Ny
    for i = 1 : Nx
        r = sqrt((xl(i) - 1/2)^2 + (yl(j) - 1/2)^2);
        phi(j,i) = 0.5*(1 - tanh(24*(r - 1/4)));
        nu(j,i) = (rate - 1)*phi(j,i) + 1;
        F(j,i) = 0.8*nu(j,i);
        grad_r = (yl(j) - 1/2)/r;
        n = 1 - (tanh(24*(r - 1/4)))^2;
        grad_nu = -(rate - 1)/2 * n * grad_r * 24;
        grad_u = 0.4*(1 - 2*yl(j));
        F(j,i) = F(j,i) - grad_nu*grad_u;
        F(j,i) = F(j,i) / 500;
    end
end
fclose(fid1); 

figure(3)
surf(xl, yl, F);
norm(F - u1)