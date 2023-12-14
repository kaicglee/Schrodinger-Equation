% initialize
idtype = 1;
vtype = 1;
idpar = [0.8, 0.8, 0.05, 0.05, -5, -5];
vpar = [0.1, 0.3, 0.1, 0.3, 100000];
tmax = 0.05;
lambda = 0.05;
level = 6;
name = "Barrier";

[x, y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

figure(1);
hold on;
animation_func(psire, x, y, level, lambda, tmax, name);

% initialize
idtype = 1;
vtype = 1;
idpar = [0.5, 0.5, 0.05, 0.05, -4, -4];
vpar = [0.1, 0.3, 0.1, 0.3, -100000];
tmax = 0.05;
lambda = 0.05;
level = 6;
name = "Well";

[x, y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

figure(2);
hold on;
animation_func(psire, x, y, level, lambda, tmax, name);

% initialize
idtype = 1;
vtype = 2;
idpar = [0.5, 0.8, 0.15, 0.06, 0, -40];
vpar = [0.25, 0.3, 0.65, 0.7, 100000];
tmax = 0.05;
lambda = 0.05;
level = 6;
name = "Double Slit";

[x, y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

figure(3);
hold on;
animation_func(psire, x, y, level, lambda, tmax, name);