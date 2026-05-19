function [Aab,VEC_ki,dEnm,Delta_abc] = shift_vector(Ham, kpoint, eps, checkSym, hermitize)
arguments
    Ham TBkit
    kpoint double
    eps = 1e-6;
    checkSym logical = true;   % 是否检查对称性
    hermitize logical = false; % 是否强制厄米化
end
%-------------------------------------------------------------
% \mathcal{A}_{n m ; a}^b
% = (i/ε_nm) [ (v^b_nm Δ^a_nm + v^a_nm Δ^b_nm)/ε_nm
%   + Σ_{p≠n,m}( v^b_np v^a_pm/ε_pm - v^a_np v^b_pm/ε_np ) ]
%   + (1/(i ε_nm)) <u_n| ∂_a ∂_b H_k |u_m>
%-------------------------------------------------------------

Nbands = Ham.Basis_num;

%---------------------------------------------
% 1. FFT 得到本征态、本征值及导数
%---------------------------------------------
[WAV_ki,EIG_ki,dH_dk_xyz,dH_dk_dk_xyz] = Ham.fft_2(kpoint);

% 能量差
dEnm = EIG_ki - EIG_ki';   % (n,m)

%---------------------------------------------
% 2. 速度矩阵 v^a = <u| dH/dk_a |u>
%---------------------------------------------
VEC_ki = zeros(Nbands,Nbands,3);
for a = 1:3
    VEC_ki(:,:,a) = WAV_ki' * dH_dk_xyz(:,:,a) * WAV_ki;
end

%---------------------------------------------
% 3. Δ_nm^a = v_nn^a - v_mm^a
%---------------------------------------------
Delta_abc = zeros(Nbands,Nbands,3);
for a = 1:3
    v_diag = diag(VEC_ki(:,:,a));
    Delta_abc(:,:,a) = v_diag.' - v_diag;  % Δ_nm^a
end

%---------------------------------------------
% 4. 计算 A_{nm;a}^b
%---------------------------------------------
Aab = zeros(Nbands,Nbands,3,3);

for a = 1:3
    for b = 1:3
        Va = VEC_ki(:,:,a);
        Vb = VEC_ki(:,:,b);
        Mab = WAV_ki' * dH_dk_dk_xyz(:,:,a,b) * WAV_ki;  % <u| ∂_a ∂_b H |u>

        for n = 1:Nbands
            for m = 1:Nbands
                if n == m || abs(dEnm(n,m)) < eps
                    continue
                end

                Enm = dEnm(n,m);

                % Δ_nm^a, Δ_nm^b
                Delta_a = Delta_abc(n,m,a);
                Delta_b = Delta_abc(n,m,b);

                % 第一部分: (v^b_nm Δ^a_nm + v^a_nm Δ^b_nm)/Enm
                term1 = ( Vb(n,m)*Delta_a + Va(n,m)*Delta_b ) / Enm;

                % 第二部分: sum_{p≠n,m}
                term2 = 0;
                for p = 1:Nbands
                    if p ~= n && p ~= m
                        if abs(dEnm(p,m)) > eps
                            term2 = term2 + (Vb(n,p)*Va(p,m)) / dEnm(p,m);
                        end
                        if abs(dEnm(n,p)) > eps
                            term2 = term2 - (Va(n,p)*Vb(p,m)) / dEnm(n,p);
                        end
                    end
                end

                % 第三部分: <u_n| d^2H/dk_a dk_b |u_m>
                term3 = Mab(n,m);

                % 组合公式
                Aab(n,m,a,b) = (1i / Enm) * ( (term1 + term2)/Enm ) ...
                             + (1/(1i*Enm)) * term3;
            end
        end
    end
end

%---------------------------------------------
% 5. 可选：检查对称性
%---------------------------------------------
if checkSym
    tol = 1e-8;
    for a = 1:3
        for b = 1:3
            diffMat = Aab(:,:,a,b) - conj(Aab(:,:,b,a));
            maxErr = max(abs(diffMat(:)));
            if maxErr > tol
                warning('Symmetry check failed: (a=%d,b=%d), maxErr=%.3e', a,b,maxErr);
            end
        end
    end
end

%---------------------------------------------
% 6. 可选：强制厄米化
%---------------------------------------------
if hermitize
    for a = 1:3
        for b = 1:3
            Aab(:,:,a,b) = 0.5*( Aab(:,:,a,b) + conj(Aab(:,:,b,a)) );
        end
    end
end

end
