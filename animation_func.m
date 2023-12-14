function animation_func(psire, x, y, level, lambda, tmax, name)

    nx = 2^level + 1;
    ny = nx;
    dx = 2^(-level);
    dy = dx;
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;

    set(gca,"NextPlot", "replacechildren")
    video = VideoWriter("Animation: " + name, 'MPEG-4');
    video.Quality = 95;
    open(video);

    for n = 1:nt
        psi_t = reshape(psire(n, :, :), [nx, ny]);
        [X, Y] = meshgrid(x, y);
        contourf(X, Y, abs(psi_t));
        xlabel("x");
        ylabel("y");
        title("2D Schrodinger Equation: " + name);
        frame = getframe(gcf);
        writeVideo(video, frame);

    end

    close(video);
end