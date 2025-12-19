function [delta, V, N, cpu_time] = fdwlsse( ...
    nfrom, nto, r, xline, bsh, Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)

% FDWLSSE  Weighted least–squares state estimator (fast–decoupled style)
%   [delta, V, N, cpu_time] = fdwlsse(nfrom, nto, r, xline, bsh, ...
%       Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)

%% Build Ybus
nb = max(max(nfrom), max(nto));
nl = length(nfrom);

Ybus = zeros(nb, nb);
for ell = 1:nl
    i = nfrom(ell);
    j = nto(ell);
    z = r(ell) + 1i * xline(ell);
    y = 1 / z;
    ysh = 1i * bsh(ell) / 2;
    Ybus(i,i) = Ybus(i,i) + y + ysh;
    Ybus(j,j) = Ybus(j,j) + y + ysh;
    Ybus(i,j) = Ybus(i,j) - y;
    Ybus(j,i) = Ybus(j,i) - y;
end
G = real(Ybus);
B = imag(Ybus);

%% Measurement bookkeeping
npinj = size(Pinj,1);
nqinj = size(Qinj,1);
npf   = size(Pflow,1);
nqf   = size(Qflow,1);
nv    = size(Vnode,1);

z = [Pinj(:,2);
     Qinj(:,2);
     Pflow(:,3);
     Qflow(:,3);
     Vnode(:,2)];

sigma = [Pinj(:,3);
         Qinj(:,3);
         Pflow(:,4);
         Qflow(:,4);
         Vnode(:,3)];

W = diag(1./(sigma.^2));          % weight matrix
nz = length(z);

%% State vector: x = [delta2..delta_nb, V1..Vnb]'
nstate = (nb-1) + nb;
x = zeros(nstate,1);
x(nb:end) = 1.0;   % flat V=1 p.u.

%% Helper: h(x)  (predicted measurements)
    function h = calc_h(x)

        % unpack state
        delta_loc = zeros(nb,1);
        V_loc     = ones(nb,1);
        delta_loc(2:nb) = x(1:nb-1);
        V_loc(:)       = x(nb:end);
        Vc = V_loc .* exp(1i*delta_loc);    % complex voltages

        % active injections
        hP = zeros(npinj,1);
        for k = 1:npinj
            i = Pinj(k,1);
            P = 0;
            for j = 1:nb
                P = P + V_loc(i)*V_loc(j) * ...
                    ( G(i,j)*cos(delta_loc(i)-delta_loc(j)) + ...
                      B(i,j)*sin(delta_loc(i)-delta_loc(j)) );
            end
            hP(k) = P;
        end

        % reactive injections
        hQ = zeros(nqinj,1);
        for k = 1:nqinj
            i = Qinj(k,1);
            Q = 0;
            for j = 1:nb
                Q = Q + V_loc(i)*V_loc(j) * ...
                    ( G(i,j)*sin(delta_loc(i)-delta_loc(j)) - ...
                      B(i,j)*cos(delta_loc(i)-delta_loc(j)) );
            end
            hQ(k) = Q;
        end

        % active line flows
        hPf = zeros(npf,1);
        for k = 1:npf
            m = Pflow(k,1);
            n = Pflow(k,2);
            ell = find(nfrom == m & nto == n, 1);
            zline = r(ell) + 1i*xline(ell);
            yline = 1/zline;
            ysh   = 1i*bsh(ell)/2;
            Vm = Vc(m); Vn = Vc(n);
            Imn = (Vm - Vn)*yline + Vm*ysh;
            Smn = Vm * conj(Imn);
            hPf(k) = real(Smn);
        end

        % reactive line flows
        hQf = zeros(nqf,1);
        for k = 1:nqf
            m = Qflow(k,1);
            n = Qflow(k,2);
            ell = find(nfrom == m & nto == n, 1);
            zline = r(ell) + 1i*xline(ell);
            yline = 1/zline;
            ysh   = 1i*bsh(ell)/2;
            Vm = Vc(m); Vn = Vc(n);
            Imn = (Vm - Vn)*yline + Vm*ysh;
            Smn = Vm * conj(Imn);
            hQf(k) = imag(Smn);
        end

        % voltage magnitudes
        hV = zeros(nv,1);
        for k = 1:nv
            i = Vnode(k,1);
            hV(k) = V_loc(i);
        end

        h = [hP; hQ; hPf; hQf; hV];
    end

%% Helper: numerical Jacobian H(x)
    function H = calc_H(x)
        epsFD = 1e-5;
        h0 = calc_h(x);
        H  = zeros(nz, nstate);
        for p = 1:nstate
            xp = x;
            xp(p) = xp(p) + epsFD;
            hp = calc_h(xp);
            H(:,p) = (hp - h0) / epsFD;
        end
    end

%% WLS iteration
tic;
for k = 1:maxiter
    h  = calc_h(x);
    rvec = z - h;
    H  = calc_H(x);
    Gmat = H' * W * H;
    rhs  = H' * W * rvec;

    dx = Gmat \ rhs;
    x  = x + dx;

    if norm(dx, inf) < toler
        break;
    end
end
cpu_time = toc;
N = k;

%% unpack final state
delta = zeros(nb,1);
V     = zeros(nb,1);
delta(1)    = 0.0;              % reference
delta(2:nb) = x(1:nb-1);
V(:)        = x(nb:end);

end
