function out = Tabc_k_fast_2D(Ham, kpoint, opts)
% BandQuantitiesSolver_point_2D_final
% 2D, single-k, multi-band, optimized final version
%
% Internal storage:
%   V(n,m,a)
%   DD(n,m,a,b)
%   Delta(n,m,a)
%   dvv_off_diag(n,m,a,b,c)
%   K(n,m,a,b,c)
%
% Returned storage:
%   out.vel(a,n,m)
%   out.ddH_mel(a,b,n,m)
%   out.v_diag_diff(a,n,m)
%   out.dvv_off_diag(a,b,c,n,m)
%   out.Berry_curvature(a,b,n)
%   out.BCP(a,b,n)
%   out.T_tensor(a,b,c,n)
%   out.dOmega(a,b,c,n)
%
% Notes:
%   - 2D only: a,b,c in {1,2}
%   - single k point only
%   - uses G_nm = Re[1/(E_n-E_m+i*eta)] to match your Python convention
%   - avoids squeeze in the computation stage

arguments
    Ham
    kpoint double
    opts.eta double = 1e-5
    opts.hermitize logical = false
end

dimk   = 2;
Nbands = Ham.Basis_num;

%% ==================== diagonalization ====================
[WAV_ki, EIG_all, dH_dk_xyz, dH_dk_dk_xyz] = Ham.fft_2(kpoint);

eigvals = reshape(EIG_all, [], 1);
eigvecs = WAV_ki;

%% ==================== energy denominators ====================
dEnm = eigvals - eigvals.';                  % E_n - E_m
G    = real(1 ./ (dEnm + 1i*opts.eta));      % Python convention
G(1:Nbands+1:end) = 0;                       % remove n = m

G2 = G.^2;
G3 = G.^3;
G5 = G.^5;

%% ==================== velocity matrices V(n,m,a) ====================
V  = complex(zeros(Nbands, Nbands, dimk));
Vt = complex(zeros(Nbands, Nbands, dimk));

for a = 1:dimk
    Va = WAV_ki' * dH_dk_xyz(:,:,a) * WAV_ki;
    if opts.hermitize
        Va = 0.5 * (Va + Va');
    end
    V(:,:,a)  = Va;
    Vt(:,:,a) = Va.';
end

%% ==================== second derivatives DD(n,m,a,b) ====================
DD = complex(zeros(Nbands, Nbands, dimk, dimk));

for a = 1:dimk
    for b = 1:dimk
        M = WAV_ki' * dH_dk_dk_xyz(:,:,a,b) * WAV_ki;
        if opts.hermitize
            M = 0.5 * (M + M');
        end
        DD(:,:,a,b) = M;
    end
end

%% ==================== diagonal velocities and Delta(n,m,a) ====================
diagV = zeros(Nbands, dimk);
Delta = zeros(Nbands, Nbands, dimk);

for a = 1:dimk
    dv = real(diag(V(:,:,a)));
    diagV(:,a)   = dv;
    Delta(:,:,a) = dv - dv.';
end

%% ==================== dvel_diag(a,b,n) ====================
% ∂_a v_b^(nn) = <n|∂_a∂_b H|n> + 2 Re Σ_m [ v_a^(nm) v_b^(mn) G_nm ]
dvel_diag = zeros(dimk, dimk, Nbands);

for a = 1:dimk
    Va = V(:,:,a);
    for b = 1:dimk
        Vb   = V(:,:,b);
        DDab = DD(:,:,a,b);
        tmp  = 2 * real(sum(Va .* Vb.' .* G, 2)) + real(diag(DDab));
        dvel_diag(a,b,:) = tmp;
    end
end

%% ==================== Berry curvature(a,b,n) ====================
Berry_curvature = zeros(dimk, dimk, Nbands);

for a = 1:dimk
    Va = V(:,:,a);
    for b = 1:dimk
        Vb = V(:,:,b);
        Berry_curvature(a,b,:) = -2 * imag(sum(Va .* Vb.' .* G2, 2));
    end
end

%% ==================== BCP(a,b,n) ====================
BCP = zeros(dimk, dimk, Nbands);

for a = 1:dimk
    Va = V(:,:,a);
    for b = 1:dimk
        Vb = V(:,:,b);
        BCP(a,b,:) = 2 * real(sum(Va .* Vb.' .* G3, 2));
    end
end

%% ==================== dvv_off_diag(n,m,a,b,c) ====================
% Complete multi-band formula:
% ∂_c(v_a^{nm} v_b^{mn})
dvv_off_diag = complex(zeros(Nbands, Nbands, dimk, dimk, dimk));

for a = 1:dimk
    Va = V(:,:,a);

    for b = 1:dimk
        Vb  = V(:,:,b);
        Vbt = Vt(:,:,b);

        for c = 1:dimk
            Vc = V(:,:,c);

            % C1(n,l) = v_c(n,l) G(n,l)
            C1 = Vc .* G;

            % C2(l,m) = v_c(l,m) G(m,l)
            C2 = Vc .* G.';

            % regrouped exact multi-band sums
            % Ablk = C1 * Va + Va * C2;
            % Bblk = Vb * C1.' + C2 * Vb;

            Ablk = C1 * Va + Va * C2;
            Bblk = C1 * Vb + Vb * C2;


            Tblk  = Ablk .* Vbt + Va .* Bblk.';
            Tdd   = DD(:,:,a,c) .* Vbt + Va .* DD(:,:,b,c).';
            M     = Tblk + Tdd;
            M(1:Nbands+1:end) = 0;

            dvv_off_diag(:,:,a,b,c) = M;
        end
    end
end

%% ==================== covD_A(n,m,a,b) ====================
% A_{nm;a}^b in the notation used previously
%
% part 1:
%   i [ v_b^{nm} Δ_a^{nm} + v_a^{nm} Δ_b^{nm} ] / ε_nm^2
%
% part 2A:
%   i / ε_nm * sum_{l!=n,m} v_a^{nl} v_b^{lm} / ε_lm
%
% part 2B:
%  -i / ε_nm * sum_{l!=n,m} v_b^{nl} v_a^{lm} / ε_nl
%
% part 3:
%  -i <n|∂_a∂_b H|m> / ε_nm
covD_A = complex(zeros(Nbands, Nbands, dimk, dimk));

oneRow = ones(1, Nbands);
oneCol = ones(Nbands, 1);

for a = 1:dimk
    Va    = V(:,:,a);
    diagA = diagV(:,a);              % N x 1

    for b = 1:dimk
        Vb    = V(:,:,b);
        diagB = diagV(:,b);          % N x 1

        % part 1
        part1 = 1i * (Vb .* Delta(:,:,a) + Va .* Delta(:,:,b)) .* G2;

        % part 2A:
        % full matmul includes l=m (killed by G_lm=0 already) and l=n extra term
        % remove l=n term explicitly
        M2A_full  = Va * (Vb .* G);
        M2A_sub_n = (diagA * oneRow) .* Vb .* G;
        part2A    = 1i * G .* (M2A_full - M2A_sub_n);

        % part 2B:
        % full matmul includes l=n (killed by G_nl with l=n) and l=m extra term
        % remove l=m term explicitly
        M2B_full  = (Vb .* G) * Va;
        M2B_sub_m = (Vb .* G) .* (oneCol * diagA.');
        part2B    = -1i * G .* (M2B_full - M2B_sub_m);

        % part 3
        part3 = -1i * DD(:,:,a,b) .* G;

        Aab = part1 + part2A + part2B + part3;
        Aab(1:Nbands+1:end) = 0;

        covD_A(:,:,a,b) = Aab;
    end
end

%% ==================== K(n,m,a,b,c) ====================
% K_{nm}^{abc} = -(v_a^{nm}/ε_nm^3) A_{mn;b}^c + i v_a^{nm} Δ_b^{nm} v_c^{mn}/ε_nm^5
K = complex(zeros(Nbands, Nbands, dimk, dimk, dimk));

for a = 1:dimk
    Va = V(:,:,a);

    for b = 1:dimk
        Deltab = Delta(:,:,b);

        for c = 1:dimk
            Abc_mn = covD_A(:,:,b,c).';   % A_{mn;b}^c as matrix in (n,m)
            Vc_t   = V(:,:,c).';

            Kmat = -(Va .* Abc_mn) .* G3 ...
                 + 1i * (Va .* Deltab .* Vc_t) .* G5;

            Kmat(1:Nbands+1:end) = 0;
            K(:,:,a,b,c) = Kmat;
        end
    end
end

%% ==================== Rsum(n,a,b,c) ====================
% R_{nlm}^{abc} = -i * v_a^{nl} v_b^{lm} v_c^{mn} / (ε_ln^2 ε_lm ε_mn^2)
Rsum = complex(zeros(Nbands, dimk, dimk, dimk));

for n = 1:Nbands
    for l = 1:Nbands
        if l == n
            continue;
        end
        Gln2 = G(l,n)^2;

        for m = 1:Nbands
            if m == n || m == l
                continue;
            end

            coeff = -1i * Gln2 * G(l,m) * G(m,n)^2;

            for a = 1:dimk
                val_a = V(n,l,a);
                for b = 1:dimk
                    val_ab = val_a * V(l,m,b);
                    for c = 1:dimk
                        Rsum(n,a,b,c) = Rsum(n,a,b,c) + coeff * val_ab * V(m,n,c);
                    end
                end
            end
        end
    end
end

%% ==================== T_tensor(a,b,c,n) ====================
% Use the same 3-term combination as in your Tabc-style implementation:
%   K^{abc} + K^{bac} + K^{cba}
%   R^{abc} + R^{bac} + R^{cba}
T_tensor = zeros(dimk, dimk, dimk, Nbands);

for a = 1:dimk
    for b = 1:dimk
        for c = 1:dimk
            K1 = reshape(K(:,:,a,b,c), [Nbands, Nbands]);
            K2 = reshape(K(:,:,b,a,c), [Nbands, Nbands]);
            K3 = reshape(K(:,:,c,b,a), [Nbands, Nbands]);

            rowSumK = sum(K1, 2) + sum(K2, 2) + sum(K3, 2);
            rowSumR = Rsum(:,a,b,c) + Rsum(:,b,a,c) + Rsum(:,c,b,a);

            T_tensor(a,b,c,:) = -real(rowSumK + rowSumR);
        end
    end
end

%% ==================== dOmega(a,b,c,n) ====================
% ∂_c Ω_ab^n = -2 Im Σ_m [ d_c(v_a^{nm}v_b^{mn})/ε_nm^2 - 2 v_a^{nm}v_b^{mn}Δ_c^{nm}/ε_nm^3 ]
dOmega = zeros(dimk, dimk, dimk, Nbands);

for a = 1:dimk
    Va = V(:,:,a);
    for b = 1:dimk
        Vb_t = V(:,:,b).';
        for c = 1:dimk
            M1 = reshape(dvv_off_diag(:,:,a,b,c), [Nbands, Nbands]) .* G2;
            M2 = -2 * (Va .* Vb_t .* Delta(:,:,c)) .* G3;
            dOmega(a,b,c,:) = -2 * imag(sum(M1 + M2, 2));
        end
    end
end

%% ==================== pack output ====================
out = struct();

out.kpoint = kpoint(:).';

out.eigvals = eigvals;
out.eigvecs = eigvecs;

out.e_diff_inv = G;

% return in the same external order you used before
out.vel           = permute(V,  [3 1 2]);      % (a,n,m)
out.ddH_mel       = permute(DD, [3 4 1 2]);    % (a,b,n,m)
out.v_diag_diff   = permute(Delta, [3 1 2]);   % (a,n,m)

out.dvel_diag       = dvel_diag;               % (a,b,n)
out.Berry_curvature = Berry_curvature;         % (a,b,n)
out.BCP             = BCP;                     % (a,b,n)
out.dvv_off_diag    = permute(dvv_off_diag, [3 4 5 1 2]); % (a,b,c,n,m)
out.T_tensor        = T_tensor;                % (a,b,c,n)
out.dOmega          = dOmega;                  % (a,b,c,n)

end