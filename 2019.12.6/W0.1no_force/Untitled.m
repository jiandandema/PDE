clear
W = 0.1;
Re = 500;
rate = 50;
x = 0:1/199:1;
y = x;
for i = 1:200
    for j = 1:200
        if (j < 100)
            n = tanh((j * 1/199 - 1.0 / 3.0) / W);
            phi(j,i) = 0.5 * (-n + 1);
        else
            n = tanh((j * 1/199 - 2.0 / 3.0) / W);
            phi(j,i) = 0.5 * (n + 1);
        end
    end
end
[XX,YY]=meshgrid(x,y);
surf(XX,YY,phi)
a = min(min(phi))