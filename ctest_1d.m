% initialization
idpar = 3;
tmax = 0.25;
lambda = 0.1;
lmin = 6;
lmax = 9;
idtype = 0;
vtype = 0;
vpar = 0;

% call function
[x6, t6, psi6, psire6, psiim6, psimod6, prob6, v6] = sch_1d_cn(tmax, lmin, lambda, idtype, idpar, vtype, vpar);
[x7, t7, psi7, psire7, psiim7, psimod7, prob7, v7] = sch_1d_cn(tmax, lmin+1, lambda, idtype, idpar, vtype, vpar);
[x8, t8, psi8, psire8, psiim8, psimod8, prob8, v8] = sch_1d_cn(tmax, lmin+2, lambda, idtype, idpar, vtype, vpar);
[x9, t9, psi9, psire9, psiim9, psimod9, prob9, v9] = sch_1d_cn(tmax, lmax, lambda, idtype, idpar, vtype, vpar);

% reshape
rpsi7 = psire7(1:2:end, 1:2:end);
rpsi8 = psire8(1:4:end, 1:4:end);
rpsi9 = psire9(1:8:end, 1:8:end);

% find dpsi
dpsi6 = rpsi7 - psire6;
dpsi7 = rpsi8 - rpsi7;
dpsi8 = rpsi9 - rpsi8;

% find mod(dpsi)
dpsi6 = rms(dpsi6, 2);
dpsi7 = rms(dpsi7, 2);
dpsi8 = rms(dpsi8, 2);

% scale dpsi
dpsi7 = 4*dpsi7;
dpsi8 = 16*dpsi8;

% plot scaled dpsi convergence
figure(1);
hold on; 
plot(t6, dpsi6, 'r-.o'); 
plot(t6, dpsi7, 'g-.+');
plot(t6, dpsi8, 'b-.*'); 
xlabel("Time");
ylabel("Scaled dpsi");
legend('dpsi6', '4 * dpsi7', '16 * dpsi8');
title("Scaled differences between dpsi for Exact Family");

% calculate error for level 6
[nt6, nx6] = size(psire6);
psixct6 = zeros(nt6, nx6);
for n = 1:nt6
    psixct6(n, :) = exp(-1i * idpar^2 * pi^2 * t6(n)) * sin(idpar * pi * x6);
end
psierr6 = psire6 - psixct6;
psierr6 = real(psierr6);

% calculate error for level 7
[nt7, nx7] = size(psire7);
psixct7 = zeros(nt7, nx7);
for n = 1:nt7
    psixct7(n, :) = exp(-1i * idpar^2 * pi^2 * t7(n)) * sin(idpar * pi * x7);
end
psierr7 = psire7 - psixct7;
psierr7 = real(psierr7);

% calculate error for level 8
[nt8, nx8] = size(psire8);
psixct8 = zeros(nt8, nx8);
for n = 1:nt8
    psixct8(n, :) = exp(-1i * idpar^2 * pi^2 * t8(n)) * sin(idpar * pi * x8);
end
psierr8 = psire8 - psixct8;
psierr8 = real(psierr8);

% calculate error for level 9
[nt9, nx9] = size(psire9);
psixct9 = zeros(nt9, nx9);
for n = 1:nt9
    psixct9(n, :) = exp(-1i * idpar^2 * pi^2 * t9(n)) * sin(idpar * pi * x9);
end
psierr9 = psire9 - psixct9;
psierr9 = real(psierr9);

% reshape
psierr7 = psierr7(1:2:end, 1:2:end);
psierr8 = psierr8(1:4:end, 1:4:end);
psierr9 = psierr9(1:8:end, 1:8:end);

% find mod(psierr)
psierr6 = rms(psierr6, 2);
psierr7 = rms(psierr7, 2);
psierr8 = rms(psierr8, 2);
psierr9 = rms(psierr9, 2);

% scale dpsi
psierr7 = 4 * psierr7;
psierr8 = 16 * psierr8;
psierr9 = 64 * psierr9;

% plot errors against exact solution
figure(2);
hold on; 
plot(t6, psierr6, 'r-.o'); 
plot(t6, psierr7, 'g-.+');
plot(t6, psierr8, 'b-.*'); 
plot(t6, psierr9, 'k-.x');
xlabel("Time");
ylabel("Psierr");
legend('psierr6', '4*psierr7', '16*psierr8', '64*psierr9');
title("Scaled Differences between Exact Solution and Numerical Solution");

% initialization
idpar = [0.5, 0.075, 0.0];
tmax = 0.01;
lambda = 0.01;
idtype = 1;
lmin = 6;
lmax = 9;
vtype = 0;
vpar = 0;

% call function
[x6, t6, psi6, psire6, psiim6, psimod6, prob6, v6] = sch_1d_cn(tmax, lmin, lambda, idtype, idpar, vtype, vpar);
[x7, t7, psi7, psire7, psiim7, psimod7, prob7, v7] = sch_1d_cn(tmax, lmin+1, lambda, idtype, idpar, vtype, vpar);
[x8, t8, psi8, psire8, psiim8, psimod8, prob8, v8] = sch_1d_cn(tmax, lmin+2, lambda, idtype, idpar, vtype, vpar);
[x9, t9, psi9, psire9, psiim9, psimod9, prob9, v9] = sch_1d_cn(tmax, lmax, lambda, idtype, idpar, vtype, vpar);

% reshape
psire7 = psire7(1:2:end, 1:2:end);
psire8 = psire8(1:4:end, 1:4:end);
psire9 = psire9(1:8:end, 1:8:end);

% find dpsi
dpsi6 = psire7 - psire6;
dpsi7 = psire8 - psire7;
dpsi8 = psire9 - psire8;

% find mod(dpsi)
dpsi6 = rms(dpsi6, 2);
dpsi7 = rms(dpsi7, 2);
dpsi8 = rms(dpsi8, 2);

% scale dpsi
dpsi7 = 4*dpsi7;
dpsi8 = 16*dpsi8;

% plot scaled dpsi convergence
figure(3);
hold on; 
plot(t6, dpsi6, 'r-.o'); 
plot(t6, dpsi7, 'g-.+');
plot(t6, dpsi8, 'b-.*'); 
xlabel("Time");
ylabel("Scaled dpsi");
legend('dpsi6', '4*dpsi7', '16*dpsi8');
title("Scaled differences between dpsi for Boosted Gaussian");