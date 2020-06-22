function dydt = twoode(t,y)
W = 0.001;rate = 99;Re = 500;lidvelocity = 0.1;


n_1 = tanh((t - 1.0 / 3.0) / W);
phi_1 = 0.5 * (-n_1 + 1);
nu_1 =  rate * phi_1 + 1;
nu_1 = nu_1/Re;
grad_nu_1 = -rate/(2*W)*(1-n_1^2)/Re;


n_2 = tanh((t - 2.0 / 3.0) / W);
phi_2 = 0.5 * (n_2 + 1);
nu_2 =  rate * phi_2 + 1;
nu_2 = nu_2/Re;
grad_nu_2 = rate/(2*W)*(1-n_2^2)/Re;

F = -8 * lidvelocity * 50.0/Re;

dydt = [y(2); (F - grad_nu_1*y(2))/nu_1*(t<0.5)+(F - grad_nu_2*y(2))/nu_2*(t>=0.5)];