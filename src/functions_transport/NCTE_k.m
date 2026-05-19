function alpha_mu = NCTE_k(Ham, tensor_index, kpoint, mu_list, T, eps, Coeffs,eta)
% --------------------------------------------------------
% NCTE_k
%   Compute nonlinear chiral thermoelectric kernel α_{abc}(k, μ, T)
%   Single-index mode:   tensor_index = [a,b,c] → alpha_mu (1 × Nmu)
%   All-27 mode:         tensor_index = []      → alpha_mu (27 × Nmu)
% --------------------------------------------------------
arguments
    Ham TBkit
    tensor_index double = []
    kpoint double = []
    mu_list double = []
    T double = 50
    eps double = 1e-4
    Coeffs double = [1 1 -1]
    eta = 0.001;
end

Nbands = Ham.Basis_num;

% default coeffs
if isempty(Coeffs)
    c1 = 1; c2 = 1; c3 = -1;
else
    c1 = Coeffs(1); c2 = Coeffs(2); c3 = Coeffs(3);
end

% ---------------- Eigen + dH/dk ----------------
[WAV,EIG,dH] = Ham.fft(kpoint);
EIG = EIG(:);

dEnm = EIG - EIG.';
inv_dEnm_sq = zeros(Nbands);
valid = abs(dEnm)+eta > eps;
inv_dEnm_sq(valid) = 1 ./ (dEnm(valid)+1i*eta).^2;

% ---------------- velocity matrices ----------------
V   = zeros(Nbands, Nbands, 3);
Vnn = zeros(Nbands, 3);
for i = 1:3
    tmp      = WAV' * dH(:,:,i) * WAV;
    V(:,:,i) = tmp;
    Vnn(:,i) = real(diag(tmp));
end

% ---------------- Energy factors ----------------
mu_list    = mu_list(:).';
E_minus_mu = EIG - mu_list;    % Nbands × Nmu

% --------------------------------------------------------
%  Case A: single index (a,b,c)
% --------------------------------------------------------
if ~isempty(tensor_index)
    if numel(tensor_index) ~= 3
        error("NCTE_k:InvalidIndex", "tensor_index must be empty or a 1x3 vector.");
    end
    a = tensor_index(1);
    b = tensor_index(2);
    c = tensor_index(3);
    
    alpha_n = zeros(Nbands,1);

    if b ~= c
        ImV   = imag( V(:,:,b) .* conj(V(:,:,c)) );
        Omega = -2 * sum(ImV .* inv_dEnm_sq, 2);
        alpha_n = alpha_n + c1 * Vnn(:,a) .* Omega;
    end

    if c ~= a
        ImV   = imag( V(:,:,c) .* conj(V(:,:,a)) );
        Omega = -2 * sum(ImV .* inv_dEnm_sq, 2);
        alpha_n = alpha_n + c2 * Vnn(:,b) .* Omega;
    end

    if a ~= b
        ImV   = imag( V(:,:,a) .* conj(V(:,:,b)) );
        Omega = -2 * sum(ImV .* inv_dEnm_sq, 2);
        alpha_n = alpha_n + c3 * Vnn(:,c) .* Omega;
    end

    % --------- temperature loop (vectorized in μ) ----------
    if isscalar(T)
        f1 = Fermi_1(E_minus_mu, T);       % Nbands × Nmu
        W  = f1 .* E_minus_mu;            % Nbands × Nmu
        alpha_mu = alpha_n.' * W;         % 1 × Nmu   (只转置小向量 alpha_n)
    else
        NT = numel(T);
        alpha_mu = cell(1,NT);
        for it = 1:NT
            f1 = Fermi_1(E_minus_mu, T(it));
            W  = f1 .* E_minus_mu;
            alpha_mu{it} = alpha_n.' * W; % 1 × Nmu
        end
    end
    return;
end

% --------------------------------------------------------
%  Case B: all 27 tensors
% --------------------------------------------------------
alpha_n_all = zeros(Nbands,27);
idx = 0;

for a = 1:3
    for b = 1:3
        for c = 1:3
            idx = idx + 1;
            col = zeros(Nbands,1);

            if b ~= c
                ImV   = imag( V(:,:,b) .* conj(V(:,:,c)) );
                Omega = -2 * sum(ImV .* inv_dEnm_sq, 2);
                col   = col + c1 * Vnn(:,a) .* Omega;
            end
            if c ~= a
                ImV   = imag( V(:,:,c) .* conj(V(:,:,a)) );
                Omega = -2 * sum(ImV .* inv_dEnm_sq, 2);
                col   = col + c2 * Vnn(:,b) .* Omega;
            end
            if a ~= b
                ImV   = imag( V(:,:,a) .* conj(V(:,:,b)) );
                Omega = -2 * sum(ImV .* inv_dEnm_sq, 2);
                col   = col + c3 * Vnn(:,c) .* Omega;
            end

            alpha_n_all(:,idx) = col;
        end
    end
end

% ---------- final μ,T combination ----------
if isscalar(T)
    f1 = Fermi_1(E_minus_mu, T);      % Nbands × Nmu
    W  = f1 .* E_minus_mu;           % Nbands × Nmu
    alpha_mu = alpha_n_all.' * W;    % 27 × Nmu   (只转置 Nbands×27，小矩阵)
else
    NT = numel(T);
    alpha_mu = cell(1,NT);
    for it = 1:NT
        f1 = Fermi_1(E_minus_mu, T(it));
        W  = f1 .* E_minus_mu;
        alpha_mu{it} = alpha_n_all.' * W;   % 27 × Nmu
    end
end

end
