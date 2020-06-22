clear;
H=1;
L=1;
W=0.1;
rate = 50;
m = 0.5;
Re = 500;
x = 0:1/199:1;
y=x;
phi = zeros(100,100);
for i = 1:200
    for j = 1:200
        bool = sqrt((x(i)-L/2)^2+(y(j)-H/2)^2);
        phi = 1 - tanh(2.4*(bool - 1/4)/W);
        phi = phi / 2;
        F1(i,j) = 1+(rate-1)*phi;
    end
end
F1 = 3*F1/Re;
[XX,YY]=meshgrid(x,y);
figure(1);
surf(XX,YY,F1)

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
figure(2)
surf(x, y, u1');

e = norm(u1' - F1);