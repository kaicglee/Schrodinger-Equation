% initialization
tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.40, 0.075, 20.0];
vtype = 1;
x1 = 0.8;
x2 = 1.0;
n = 251;

% set up n uniformly spaced values for vc
vc = linspace(exp(-2), exp(5), n);

% set up fe
fe = zeros(1, n);

% calculate fe for each value in vc
for k = 1:n
    vpar = [x1, x2, vc(k)];
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

    % get temporal average and normalize
    tempavg = mean(prob, 1);
    tempavg = tempavg / prob(end);

    % find ln(fe)
    [val1, idx1] = min(abs(x - x1));
    [val2, idx2] = min(abs(x - x2));
    num = tempavg(idx2) - tempavg(idx1);
    denom = x2 - x1;
    fe(k) = num / denom;
end

clf;
hold on; 
plot(log(vc), log(fe)); 
xlabel("ln(V0)");
ylabel("ln(f(x1, x2))");
title("Barrier Survey plot");