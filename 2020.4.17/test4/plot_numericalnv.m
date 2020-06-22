clear;
H=1;
L=1;
W=0.5;
rate = 50;
Re = 500;
m = 0.5;
x = 0:1/199:1;
y=x;
phi = zeros(101,101);
% for i = 1:200
%     for j = 1:200
%         if (j < 100)
%             n = tanh((j * 1/199 - 1.0 / 3.0) / W);
%             phi(i,j) = 0.5 * (-n + 1);
%             nu =  (rate-1) * phi(i,j) + 1;
%             tao(i,j) = nu * 3.0 / Re;
%         else
%             n = tanh((j * 1/199 - 2.0 / 3.0) / W);
%             phi(i,j) = 0.5 * (n + 1);
%             nu = (rate-1) * phi(i,j) + 1;
%             tao(i,j) = nu * 3.0 / Re;
%         end
%     end
% end
% [XX,YY]=meshgrid(x,y);
% figure(1);
% surf(XX,YY,tao')

fid1 = fopen('./ex5001out.dat', 'r');
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
surf(x, y, u1);
min_=min(min(u1));
max_ = max(max(u1));
% ee = u1-tao';
% e = norm(u1 - tao')
% figure(3)
% surf(x, y, phi)