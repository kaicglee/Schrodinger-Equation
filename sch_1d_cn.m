function [x, t, psi, psire, psiim, psimod, prob, v] = ...
   sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Inputs
%
%     tmax:    (real scalar) Maximum integration time
%     level:   (integer scalar) Discretization level
%     lambda:  (real scalar) dt/dx
%     idtype:  (0 or 1) Selects initial data type
%     idpar:   (real vector (1 or 3)) Initial data parameters
%     vtype:   (0 or 1) Selects potential type
%     vpar:    (real vector (1 or 3)) Potential parameters
%
% Outputs
%
%     x:       (nx vector) x coordinates
%     t:       (nt vector) t coordinates
%     psi:     (nt x nx array) Computed psi values
%     psire:   (nt x nx array) Computed psi_re values
%     psiim:   (nt x nx array) Computed psi_im values
%     psimod:  (nt x nx array) Computed sqrt(psi psi*) values
%     prob:    (nt x nx array) Computed running integral values
%     v:       (nx vector) Potential values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % increments in space and time
    dx = 2^(-level);
    dt = lambda * dx;

    % number of increments
    nx = 2^level + 1;
    nt = round(tmax / dt) + 1;
    
    % set up grids
    x = linspace(0.0, 1.0, nx);
    t = [0 : nt-1] * dt;

    % set up solutions
    psi = zeros(nt, nx);
    psire = zeros(nt, nx);
    psiim = zeros(nt, nx);
    psimod = zeros(nt, nx);
    prob = zeros(nt, nx);
    
    % initial conditions
    if idtype == 0
        m = idpar;
        psi(1, :) = sin(m * pi * x);
    elseif idtype == 1
        x0 = idpar(1);
        sigma = idpar(2);
        p = idpar(3);
        psi(1, :) = exp(1i * p * x) .* exp(-((x - x0) / sigma).^2);
    else
      fprintf('sch_1d_cn: Invalid idtype=%d\n', idtype);
      return;
    end
    
    % initialize potential
    v = zeros(1, nx);
    if vtype == 0
        v = zeros(1, nx);
    elseif vtype == 1
        xmin = vpar(1);
        xmax = vpar(2);
        vc = vpar(3);
        a = find(x >= xmin & x <= xmax);
        v(a) = vc;
    else
      fprintf('sch_1d_cn: Invalid vtype=%d\n', vtype);
      return
    end
    
    % initialize storage for sparse matrix and RHS ...
    dl = zeros(nx, 1);
    d  = zeros(nx, 1);
    du = zeros(nx, 1);
    f  = zeros(nx, 1);
    
    % set up tridiagonal system
    dl = 0.5 / dx^2 * ones(nx, 1);
    d = 1i / dt - 1 / dx^2 - v' / 2;
    du = dl;
    
    % boundary conditions
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx - 1) = 0.0;
    d(nx) = 1.0;
    
    % define sparse matrix
    A = spdiags([dl d du], -1:1, nx, nx);
    
    % compute solution using CN scheme
    for n = 1:nt-1

       % RHS
       f(2:nx-1) = -0.5 * psi(n, 1:nx-2) / dx^2 + (1i / dt + 1 / dx^2 + ...
           0.5 * v(2:nx-1)) .* psi(n, 2:nx-1) - 0.5 * psi(n, 3:nx) / dx^2;
       f(1) = 0.0;
       f(nx) = 0.0;

       % next time step
       psi(n + 1, :) = A \ f;
    end

    % calculate real, imaginary, mod, and probability density for psi
    psire = real(psi);
    psiim = imag(psi);
    psimod = abs(psi);
    rho = psimod.^2;
    for i = 1:nx-1
        next_prob = 0.5 * ((rho(:, i) + rho(:, i + 1)) .* (x(i + 1) - x(i)));
        prob(:, i + 1) = prob(:, i) + next_prob;
    end
end
