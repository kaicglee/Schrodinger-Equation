function [x, y, t, psi, psire, psiim, psimod, v] = ...
         sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Inputs
%
%     tmax:    (real scalar) Maximum integration time
%     level:   (integer scalar) Discretization level
%     lambda:  (real scalar) dt/dx
%     idtype:  (0 or 1) Selects initial data type
%     idpar:   (real vector (2 or 6)) Vector of initial data parameters
%     vtype:   (0, 1, or 2) Selects potential type
%     vpar:    (real vector (1 or 5)) Vector of potential parameters
%
% Outputs
%
%     x:       (nx vector) x coordinates
%     y:       (ny vector) y coordinates
%     t:       (nt vector) t coordinates
%     psi:     (nt x nx x ny array) computed psi values
%     psire:   (nt x nx x ny array) computed psi_re values
%     psiim:   (nt x nx x ny array) computed psi_im values
%     psimod:  (nt x nx x ny array) computed sqrt(psi psi*) values
%     v:       (nx x ny array) Potential values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % increments in space and time
    dx = 2^(-level);
    dy = dx;
    dt = lambda * dx;

    % number of increments
    nx = 2^level + 1;
    ny = nx;
    nt = round(tmax / dt) + 1;
    
    % set up grids
    x = linspace(0.0, 1.0, nx);
    y = linspace(0.0, 1.0, ny);
    t = [0 : nt-1] * dt;

    % set up solutions
    psi = zeros(nt, nx, ny);
    psire = zeros(nt, nx, ny);
    psiim = zeros(nt, nx, ny);
    psimod = zeros(nt, nx, ny);

    % initial conditions
    if idtype == 0
        mx = idpar(1);
        my = idpar(2);
        psi_init = zeros(nx, ny);
        for n = 1:nx
            for m = 1:ny
                psi_init(n, m) = sin(mx * pi * x(n)) .* sin(my * pi * y(m)).';
            end
        end
        psi(1, :, :) = psi_init;
    elseif idtype == 1
        x0 = idpar(1);
        y0 = idpar(2);
        sigmax = idpar(3);
        sigmay = idpar(4);
        px = idpar(5);
        py = idpar(6);
        psi_init = zeros(nx, ny);
        for n = 1:nx
            for m = 1:ny
                psi_init(n, m) = exp(1i * px * x(n)) .* exp(1i * py * y(m)) ...
                    .* exp(-((x(n) - x0) / sigmax).^2 ...
                    + ((y(m) - y0) / sigmay).^2);
            end
        end
        psi(1, :, :) = psi_init;
    else
      fprintf('sch_2d_adi: Invalid idtype=%d\n', idtype);
      return;
    end

    % initialize potential
    v = zeros(nx, ny);
    if vtype == 0
        v = zeros(nx, ny);
    elseif vtype == 1
        xmin = vpar(1);
        xmax = vpar(2);
        ymin = vpar(3);
        ymax = vpar(4);
        vc = vpar(5);
        a = find(x >= xmin & x <= xmax);
        b = find(y >= ymin & y <= ymax);
        v(a, b) = vc;
    elseif vtype == 2
        x1 = vpar(1);
        x2 = vpar(2);
        x3 = vpar(3);
        x4 = vpar(4);
        vc = vpar(5);
        jprime = (ny - 1) / 4 + 1;
        c = find((x < x1) | (x > x2 & x < x3) | x > x4);
        v(c, jprime) = vc;
        v(c, jprime + 1) = vc;
    else
      fprintf('sch_2d_adi: Invalid vtype=%d\n', vtype);
      return
    end

    % set up base differential tridiagonal system
    dl = (1 / dx^2) * (1i * dt / 2) * ones(nx, 1);
    d = (-2 / dx^2) * (1i * dt / 2) * ones(nx, 1);
    du = (1 / dx^2) * (1i * dt / 2) * ones(nx, 1);

    % boundary conditions
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx - 1) = 0.0;
    d(nx) = 1.0;

    % define base differential matrix
    D_base = spdiags([dl d du], -1:1, nx, nx);
    xd = ones(nx, 1);
    xd(1) = 0.0;
    xd(nx) = 0.0;
    xd_matrix = spdiags(xd, 0, nx, nx);

    % build matrices for RHS and LHS
    D_right = D_base + xd_matrix;
    D_left = xd_matrix - D_base;

    % create V matrix with i * dt / 2 factor
    V = v * 1i * dt / 2;

    % Compute solution using ADI
    for n = 1:nt-1

        % boundary conditions
        psi(n, 1, :) = 0.0;
        psi(n, end ,:) = 0.0;
        psi(n,:,1) = 0.0;
        psi(n,:,end) = 0.0;
        
        % set up intermediate psi value
        psi_int = zeros(nx,ny);

        % solve for psi_int
        for i = 2:nx-1 
            psi_i = reshape(psi(n, i, :), [ny, 1]);
            yd = ones(ny, 1) - V(i, :).';
            yd(1) = 0.0;
            yd(ny) = 0.0;
            yd_matrix = spdiags(yd, 0, ny, ny);
            D_right = D_base + yd_matrix;
            psi_int(i, :)= D_right * psi_i;
        end
    
        % set up psi_half
        psi_half = zeros(nx,ny);

        % solve for psi(n+1/2)
        for j = 2:ny-1 
            psi_j = psi_int(:, j);
            f = D_right * psi_j;
            psi_half(:, j) = D_left \ f;
        end

        % solve for psi(n+1)
        for i = 2:nx-1 
            yd = ones(ny, 1) +  V(i, :).';
            yd(1) = 0.0;
            yd(ny) = 0.0;
            yd_matrix = spdiags(yd, 0, ny, ny);
            D_left =  yd_matrix - D_base;
            psi(n + 1, i, :) = D_left \ (psi_half(i, :).');
        end 
    end

    % calculate real, imaginary, and mod for psi
    psire = real(psi);
    psiim = imag(psi);
    psimod = abs(psi);
end