% initialization
idpar = [2, 3];
tmax = 0.05;
lambda = 0.05;
lmin = 6;
lmax = 9;
idtype = 0;
vtype = 0;
vpar = 0;

% call function
[x6, y6, t6, psi6, psire6, psiim6, psimod6, v6] = sch_2d_adi(tmax, lmin, lambda, idtype, idpar, vtype, vpar);
[x7, y7, t7, psi7, psire7, psiim7, psimod7, v7] = sch_2d_adi(tmax, lmin+1, lambda, idtype, idpar, vtype, vpar);
[x8, y8, t8, psi8, psire8, psiim8, psimod8, v8] = sch_2d_adi(tmax, lmin+2, lambda, idtype, idpar, vtype, vpar);
[x9, y9, t9, psi9, psire9, psiim9, psimod9, v9] = sch_2d_adi(tmax, lmax, lambda, idtype, idpar, vtype, vpar);

% reshape
rpsi7 = psire7(1:2:end, 1:2:end, 1:2:end);
rpsi8 = psire8(1:4:end, 1:4:end, 1:4:end);
rpsi9 = psire9(1:8:end, 1:8:end, 1:8:end);

% find dpsi
dpsi6 = rpsi7 - psire6;
dpsi7 = rpsi8 - rpsi7;
dpsi8 = rpsi9 - rpsi8;

% find spatial norms of dpsi
dpsi6 = norm2d(dpsi6, lmin, lambda, tmax);
dpsi7 = norm2d(dpsi7, lmin, lambda, tmax);
dpsi8 = norm2d(dpsi8, lmin, lambda, tmax);

% scale dpsi
dpsi7 = 4 * dpsi7;
dpsi8 = 16 * dpsi8;

% plot scaled dpsi convergence
figure(1);
hold on; 
plot(t6, dpsi6, 'r-.o'); 
plot(t6, dpsi7, 'g-.+');
plot(t6, dpsi8, 'b-.*'); 
xlabel("Time");
ylabel("Scaled dpsi");
legend('dpsi6', '4*dpsi7', '16*dpsi8');
title("Scaled differences between dpsi for Exact Family");

% calculate error for level 6
[nt6, nx6, ny6] = size(psire6);
psixct6 = zeros(nt6, nx6, ny6);
for n = 1:nt6
    for i = 1:nx6
        for j = 1:ny6
            psixct6(n, i, j) = exp(-1i * (idpar(1)^2 + idpar(2)^2) * pi^2 * t6(n)) .* sin(idpar(1) * pi * x6(i)) .* sin(idpar(2) * pi * y6(j)).'; 
        end
    end
end
psierr6 = psire6 - real(psixct6);

% calculate error for level 7
[nt7, nx7, ny7] = size(psire7);
psixct7 = zeros(nt7, nx7, ny7);
for n = 1:nt7
    for i = 1:nx7
        for j = 1:ny7
            psixct7(n, i, j) = exp(-1i * (idpar(1)^2 + idpar(2)^2) * pi^2 * t7(n)) .* sin(idpar(1) * pi * x7(i)) .* sin(idpar(2) * pi * y7(j)).'; 
        end
    end
end
psierr7 = psire7 - real(psixct7);

% calculate error for level 8
[nt8, nx8, ny8] = size(psire8);
psixct8 = zeros(nt8, nx8, ny8);
for n = 1:nt8
    for i = 1:nx8
        for j = 1:ny8
            psixct8(n, i, j) = exp(-1i * (idpar(1)^2 + idpar(2)^2) * pi^2 * t8(n)) .* sin(idpar(1) * pi * x8(i)) .* sin(idpar(2) * pi * y8(j)).'; 
        end
    end
end
psierr8 = psire8 - real(psixct8);


%calculate error for level 9
[nt9, nx9, ny9] = size(psire9);
psixct9 = zeros(nt9, nx9, ny9);
for n = 1:nt9
    for i = 1:nx9
        for j = 1:ny9
            psixct9(n, i, j) = exp(-1i * (idpar(1)^2 + idpar(2)^2) * pi^2 * t9(n)) .* sin(idpar(1) * pi * x9(i)) .* sin(idpar(2) * pi * y9(j)).'; 
        end
    end
end
psierr9 = psire9 - real(psixct9);


% reshape
psierr7 = psierr7(1:2:end, 1:2:end, 1:2:end);
psierr8 = psierr8(1:4:end, 1:4:end, 1:4:end);
psierr9 = psierr9(1:8:end, 1:8:end, 1:8:end);

% find spatial norms of psierr
psierr6 = norm2d(psierr6, lmin, lambda, tmax);
psierr7 = norm2d(psierr7, lmin, lambda, tmax);
psierr8 = norm2d(psierr8, lmin, lambda, tmax);
psierr9 = norm2d(psierr9, lmin, lambda, tmax);

% scale psierr
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

% function for spatial norm
function norm = norm2d(psi, level, lambda, tmax)
    nx = 2^level + 1;
    ny = nx;
    dx = 2^(-level);
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;
    
    norm = zeros(1, nt);
    
    for n = 1:nt
        psi_t = reshape(psi(n, :, :), [nx, ny]);
        norm(n) = sqrt(sum(sum(abs(psi_t)^2)) ./ (nx * ny));
    end
end