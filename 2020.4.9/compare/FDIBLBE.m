clc, clear, close all;                
Nx = 200;
Ny = 200;
AQ = 9;
BQ = 6;

xl = 0 : 1 / (Nx - 1) : 1;
yl = 0 : 1 / (Ny - 1) : 1;

fid1 = fopen('./ex0_001out.dat', 'r'); 
%fid1 = fopen('./Vorticity2001out.dat', 'r');

for j = 1 : Ny
  for i = 1 : Nx
    for k = 1 : AQ

      u(j, i, k) = fscanf(fid1, '%g', 1);
    end
         u1(j, i) = u(j, i, 1);
         u2(j, i) = u(j, i, 2);
         u3(j, i) = u(j, i, 3);
         u4(j, i) = u(j, i, 4);
         u5(j, i) = u(j, i, 5);
         u6(j, i) = u(j, i, 6);
         u7(j, i) = u(j, i, 7);
         u8(j, i) = u(j, i, 8);
         u9(j, i) = u(j, i, 9); 
         
         p(j, i) = u1(j, i) + u2(j, i) + u3(j, i) + u4(j, i) + u5(j, i) + u6(j, i) + u7(j, i) + u8(j, i) + u9(j, i);
         fp1(j, i) = u2(j, i) - u4(j, i) + u6(j, i) - u7(j, i) - u8(j, i) + u9(j, i);
         fp2(j, i) = u3(j, i) - u5(j, i) + u6(j, i) + u7(j, i) - u8(j, i) - u9(j, i);
         f1(j, i) = fp1(j, i) / p(j, i);
         f2(j, i) = fp2(j, i) / p(j, i);
         f(j, i) = f1(j, i)^2 + f2(j, i)^2;
  end
end
fclose(fid1);
plot(xl,f1(:,200))
hold on
fid1 = fopen('./ex0_01out.dat', 'r'); 
%fid1 = fopen('./Vorticity2001out.dat', 'r');

for j = 1 : Ny
  for i = 1 : Nx
    for k = 1 : AQ

      u(j, i, k) = fscanf(fid1, '%g', 1);
    end
         u1(j, i) = u(j, i, 1);
         u2(j, i) = u(j, i, 2);
         u3(j, i) = u(j, i, 3);
         u4(j, i) = u(j, i, 4);
         u5(j, i) = u(j, i, 5);
         u6(j, i) = u(j, i, 6);
         u7(j, i) = u(j, i, 7);
         u8(j, i) = u(j, i, 8);
         u9(j, i) = u(j, i, 9); 
         
         p(j, i) = u1(j, i) + u2(j, i) + u3(j, i) + u4(j, i) + u5(j, i) + u6(j, i) + u7(j, i) + u8(j, i) + u9(j, i);
         fp1(j, i) = u2(j, i) - u4(j, i) + u6(j, i) - u7(j, i) - u8(j, i) + u9(j, i);
         fp2(j, i) = u3(j, i) - u5(j, i) + u6(j, i) + u7(j, i) - u8(j, i) - u9(j, i);
         f1(j, i) = fp1(j, i) / p(j, i);
         f2(j, i) = fp2(j, i) / p(j, i);
         f(j, i) = f1(j, i)^2 + f2(j, i)^2;
  end
end
fclose(fid1);
plot(xl,f1(:,200))
hold on

fid1 = fopen('./ex0_1out.dat', 'r'); 
%fid1 = fopen('./Vorticity2001out.dat', 'r');

for j = 1 : Ny
  for i = 1 : Nx
    for k = 1 : AQ

      u(j, i, k) = fscanf(fid1, '%g', 1);
    end
         u1(j, i) = u(j, i, 1);
         u2(j, i) = u(j, i, 2);
         u3(j, i) = u(j, i, 3);
         u4(j, i) = u(j, i, 4);
         u5(j, i) = u(j, i, 5);
         u6(j, i) = u(j, i, 6);
         u7(j, i) = u(j, i, 7);
         u8(j, i) = u(j, i, 8);
         u9(j, i) = u(j, i, 9); 
         
         p(j, i) = u1(j, i) + u2(j, i) + u3(j, i) + u4(j, i) + u5(j, i) + u6(j, i) + u7(j, i) + u8(j, i) + u9(j, i);
         fp1(j, i) = u2(j, i) - u4(j, i) + u6(j, i) - u7(j, i) - u8(j, i) + u9(j, i);
         fp2(j, i) = u3(j, i) - u5(j, i) + u6(j, i) + u7(j, i) - u8(j, i) - u9(j, i);
         f1(j, i) = fp1(j, i) / p(j, i);
         f2(j, i) = fp2(j, i) / p(j, i);
         f(j, i) = f1(j, i)^2 + f2(j, i)^2;
  end
end
fclose(fid1);
plot(xl,f1(:,200))
hold on
fid1 = fopen('./ex0_5out.dat', 'r'); 
%fid1 = fopen('./Vorticity2001out.dat', 'r');

for j = 1 : Ny
  for i = 1 : Nx
    for k = 1 : AQ

      u(j, i, k) = fscanf(fid1, '%g', 1);
    end
         u1(j, i) = u(j, i, 1);
         u2(j, i) = u(j, i, 2);
         u3(j, i) = u(j, i, 3);
         u4(j, i) = u(j, i, 4);
         u5(j, i) = u(j, i, 5);
         u6(j, i) = u(j, i, 6);
         u7(j, i) = u(j, i, 7);
         u8(j, i) = u(j, i, 8);
         u9(j, i) = u(j, i, 9); 
         
         p(j, i) = u1(j, i) + u2(j, i) + u3(j, i) + u4(j, i) + u5(j, i) + u6(j, i) + u7(j, i) + u8(j, i) + u9(j, i);
         fp1(j, i) = u2(j, i) - u4(j, i) + u6(j, i) - u7(j, i) - u8(j, i) + u9(j, i);
         fp2(j, i) = u3(j, i) - u5(j, i) + u6(j, i) + u7(j, i) - u8(j, i) - u9(j, i);
         f1(j, i) = fp1(j, i) / p(j, i);
         f2(j, i) = fp2(j, i) / p(j, i);
         f(j, i) = f1(j, i)^2 + f2(j, i)^2;
  end
end
fclose(fid1);
plot(xl,f1(:,200))
hold on
plot(xl,0.01*xl)
legend('W=0.001','W=0.01','W=0.1','W=0.5')