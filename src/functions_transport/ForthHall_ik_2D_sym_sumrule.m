function result = ForthHall_ik_2D_sym_sumrule(Ham, kpoint, mu_list, T, tau, selectbands, varargin)
% ForthHall_ik_2D_sym_sumrule
%
% Single-k 2D local fourth-Hall response using the sum-rule /
% integration-by-parts form, with optional chemical-potential correction:
%
%   delta chi_{abcde}^{tau1} =
%       - tau * Q2_{de}  * sum_n d_c Omega_{ab,n} f0'
%       + tau * Q3_{cde} * sum_n v_{a,n} v_{b,n} f0''.
%
% The correction tensors Q2 and Q3 must be computed from the full k mesh
% first, then supplied to every k point through the name-value option
%       'mu_correction', muCorr
% where muCorr.Q2 and muCorr.Q3 have either scalar-T shapes
%       Q2: [Nmu,2,2],       Q3: [Nmu,2,2,2]
% or vector-T shapes
%       Q2: [NT,Nmu,2,2],    Q3: [NT,Nmu,2,2,2].
%
% Other varargin entries are passed directly to Tabc_k_fast_2D.
%
% Scalar T output, backward-compatible:
%   result.vec_tau1  : [Nmu,10]
%   result.vec_tau3  : [Nmu,10] legacy sampled tensor components
%   result.vec_tau3_ab_cde : [Nmu,4] independent Hall-sector components
%                            [xyxxx, xyxxy, xyxyy, xyyyy]
%   result.chi_tau1  : [Nmu,6]
%   result.sym_tau1  : [Nmu,2,2,2,2,2]
%
% Vector T output:
%   result.vec_tau1  : [NT,Nmu,10]
%   result.vec_tau3  : [NT,Nmu,10] legacy sampled tensor components
%   result.vec_tau3_ab_cde : [NT,Nmu,4] independent Hall-sector components
%   result.chi_tau1  : [NT,Nmu,6]
%   result.sym_tau1  : [NT,Nmu,2,2,2,2,2]

if nargin < 6 || isempty(selectbands)
    selectbands = 1:Ham.Basis_num;
end

% Strip local-only options before forwarding the remaining varargin to
% Tabc_k_fast_2D. This keeps Tabc_k_fast_2D's original interface intact.
[local_opts, tabc_args] = parse_local_options(varargin{:});

mu_col  = mu_list(:);       % [Nmu,1]
T_vec   = T(:);             % [NT,1]
Nmu     = numel(mu_col);
NT      = numel(T_vec);
scalarT = (NT == 1);

if any(~isfinite(T_vec)) || any(T_vec <= 0)
    error('This sum-rule implementation assumes finite temperature T > 0.');
end

% Optional chemical-potential correction tensors.  Internally everything is
% stored as [NT,Nmu,...] even for scalar T.
use_mu_correction = local_opts.include_mu_correction && ~isempty(local_opts.mu_correction);
if use_mu_correction
    [Q2_mu, Q3_mu] = normalize_mu_correction_local(local_opts.mu_correction, NT, Nmu, scalarT);
else
    Q2_mu = zeros(NT, Nmu, 2, 2);
    Q3_mu = zeros(NT, Nmu, 2, 2, 2);
end

% ------------------------------------------------------------
% Evaluate single-point band quantities from Tabc_k_fast_2D
% ------------------------------------------------------------
out = Tabc_k_fast_2D(Ham, kpoint, tabc_args{:});

sel  = selectbands(:).';
nsel = numel(sel);

% eigvals: force row vector for stable contractions
E = out.eigvals(sel);
E = real(E(:).');                         % [1,nsel]

% diagonal velocities v_a^{nn}
Vdiag = zeros(2, nsel);
for a = 1:2
    Va_full = reshape(out.vel(a,:,:), size(out.vel,2), size(out.vel,3));
    dv_full = diag(Va_full);
    Vdiag(a,:) = real(dv_full(sel)).';
end

% dvel_diag(a,b,n)
dvel = real(out.dvel_diag(:,:,sel));      % [2,2,nsel]

% T_tensor(c,d,e,n)
Tten = real(out.T_tensor(:,:,:,sel));     % [2,2,2,nsel]

% BCP / G_{de}(n)
Gten = real(out.BCP(:,:,sel));            % [2,2,nsel]

% dOmega(a,b,c,n) = partial_c Omega_ab
% The uploaded/original code uses dOmega(a,b,c,n); keep the same convention.
dOm  = real(out.dOmega(:,:,:,sel));       % [2,2,2,nsel]

% ------------------------------------------------------------
% Raw local integrands, vectorized over (T,mu)
% ------------------------------------------------------------
raw_tau1      = zeros(NT, Nmu, 2,2,2,2,2);
raw_tau3      = zeros(NT, Nmu, 2,2,2,2,2);
raw_tau1_base = zeros(NT, Nmu, 2,2,2,2,2);
raw_tau1_mu   = zeros(NT, Nmu, 2,2,2,2,2);

for iT = 1:NT
    Tnow = T_vec(iT);

    % Emu: [Nmu,nsel], Fermi derivatives: [Nmu,nsel]
    Emu  = E - mu_col;
    Fp   = Fermi_1_local(Emu, Tnow);
    Fpp  = Fermi_2_local(Emu, Tnow);

    for a = 1:2
        va = Vdiag(a,:);            % [1,nsel]
        for b = 1:2
            vb = Vdiag(b,:);
            vavb_sum = Fpp * (va .* vb).';    % [Nmu,1]

            for c = 1:2
                dOmega_abc = reshape(dOm(a,b,c,:), 1, nsel);
                dOmega_sum = Fp * dOmega_abc.';      % [Nmu,1]

                for d = 1:2
                    vd = Vdiag(d,:);
                    for e = 1:2
                        ve = Vdiag(e,:);

                        % ---- pieces used repeatedly, all [1,nsel] ----
                        T_cde = reshape(Tten(c,d,e,:), 1, nsel);
                        T_bcd = reshape(Tten(b,c,d,:), 1, nsel);
                        T_acd = reshape(Tten(a,c,d,:), 1, nsel);

                        dvel_ba = reshape(dvel(b,a,:), 1, nsel);
                        dvel_ae = reshape(dvel(a,e,:), 1, nsel);
                        dvel_be = reshape(dvel(b,e,:), 1, nsel);
                        dvel_ed = reshape(dvel(e,d,:), 1, nsel);

                        G_de = reshape(Gten(d,e,:), 1, nsel);

                        % ==================================================
                        % tau1 base part. Each contraction returns [Nmu,1].
                        % ==================================================
                        core_tau1_base = ...
                              (-2/3) * (Fpp * (va .* vb .* T_cde).') ...
                            + (-4/3) * (Fp  * (dvel_ba .* T_cde).') ...
                            + (   1) * (Fp  * (dvel_ae .* T_bcd).') ...
                            - (   1) * (Fp  * (dvel_be .* T_acd).') ...
                            + (   1) * (Fpp * (va .* ve .* T_bcd).') ...
                            - (   1) * (Fpp * (vb .* ve .* T_acd).') ...
                            - ( 1/2) * (Fp  * (dOmega_abc .* G_de).');

                        % ==================================================
                        % Chemical-potential correction to tau1:
                        %   - Q2_de  * sum_n d_c Omega_ab f0'
                        %   + Q3_cde * sum_n v_a v_b f0''
                        % Then the whole correction is multiplied by tau.
                        % ==================================================
                        if use_mu_correction
                            Q2_de  = reshape(Q2_mu(iT,:,d,e),  Nmu, 1);
                            Q3_cde = reshape(Q3_mu(iT,:,c,d,e), Nmu, 1);
                            core_tau1_mu = - Q2_de .* dOmega_sum + Q3_cde .* vavb_sum;
                        else
                            core_tau1_mu = zeros(Nmu,1);
                        end

                        raw_tau1_base(iT,:,a,b,c,d,e) = reshape(tau * core_tau1_base, 1, Nmu);
                        raw_tau1_mu(iT,:,a,b,c,d,e)   = reshape(tau * core_tau1_mu,   1, Nmu);
                        raw_tau1(iT,:,a,b,c,d,e)      = reshape(tau * (core_tau1_base + core_tau1_mu), 1, Nmu);

                        % ==================================================
                        % tau3 part. Each contraction returns [Nmu,1].
                        % ==================================================
                        core_tau3 = ...
                              (Fp  * (dOmega_abc .* dvel_ed).') ...
                            + (Fpp * (dOmega_abc .* vd .* ve).');

                        raw_tau3(iT,:,a,b,c,d,e) = reshape(tau^3 * core_tau3, 1, Nmu);
                    end
                end
            end
        end
    end
end

% ------------------------------------------------------------
% Tensor projection / physical index symmetries
%
% tau1: chi^{tau1}_{a b c d e} is symmetrized over the four
%       electric-field indices (b,c,d,e).
%
% tau3: chi^{tau3}_{a b c d e} is first antisymmetrized over the Hall
%       index pair (a,b), then symmetrized over the remaining three
%       electric-field indices (c,d,e).  Do NOT symmetrize tau3 over
%       (b,c,d,e); that would erase the Hall antisymmetry structure.
% ------------------------------------------------------------
ordered_bcdes = [
    1 1 1 1
    1 1 1 2
    1 1 2 2
    1 2 2 2
    2 2 2 2
];

ordered_cdes = [
    1 1 1
    1 1 2
    1 2 2
    2 2 2
];

sym_tau1      = zeros(size(raw_tau1));
sym_tau1_base = zeros(size(raw_tau1_base));
sym_tau1_mu   = zeros(size(raw_tau1_mu));

% ---- tau1: full symmetrization over (b,c,d,e) ----
for a = 1:2
    for iord = 1:size(ordered_bcdes,1)
        bcde = ordered_bcdes(iord,:);
        plist = unique(perms(bcde), 'rows');
        nperm = size(plist,1);

        acc1      = zeros(NT, Nmu);
        acc1_base = zeros(NT, Nmu);
        acc1_mu   = zeros(NT, Nmu);

        for ip = 1:nperm
            b = plist(ip,1);
            c = plist(ip,2);
            d = plist(ip,3);
            e = plist(ip,4);

            acc1      = acc1      + raw_tau1(:,:,a,b,c,d,e);
            acc1_base = acc1_base + raw_tau1_base(:,:,a,b,c,d,e);
            acc1_mu   = acc1_mu   + raw_tau1_mu(:,:,a,b,c,d,e);
        end

        acc1      = acc1      / nperm;
        acc1_base = acc1_base / nperm;
        acc1_mu   = acc1_mu   / nperm;

        for ip = 1:nperm
            b = plist(ip,1);
            c = plist(ip,2);
            d = plist(ip,3);
            e = plist(ip,4);

            sym_tau1(:,:,a,b,c,d,e)      = acc1;
            sym_tau1_base(:,:,a,b,c,d,e) = acc1_base;
            sym_tau1_mu(:,:,a,b,c,d,e)   = acc1_mu;
        end
    end
end

% ---- tau3: antisymmetrization over (a,b), then symmetrization over (c,d,e) ----
anti_tau3 = zeros(size(raw_tau3));
for a = 1:2
    for b = 1:2
        anti_tau3(:,:,a,b,:,:,:) = 0.5 .* ( ...
            raw_tau3(:,:,a,b,:,:,:) - raw_tau3(:,:,b,a,:,:,:));
    end
end

sym_tau3 = zeros(size(raw_tau3));
for a = 1:2
    for b = 1:2
        for iord = 1:size(ordered_cdes,1)
            cde = ordered_cdes(iord,:);
            plist = unique(perms(cde), 'rows');
            nperm = size(plist,1);

            acc3 = zeros(NT, Nmu);
            for ip = 1:nperm
                c = plist(ip,1);
                d = plist(ip,2);
                e = plist(ip,3);
                acc3 = acc3 + anti_tau3(:,:,a,b,c,d,e);
            end
            acc3 = acc3 / nperm;

            for ip = 1:nperm
                c = plist(ip,1);
                d = plist(ip,2);
                e = plist(ip,3);
                sym_tau3(:,:,a,b,c,d,e) = acc3;
            end
        end
    end
end

sym_total = sym_tau1 + sym_tau3;

% ------------------------------------------------------------
% Flatten 10 independent components and chi6 post-processing
% ------------------------------------------------------------
vec_tau1       = zeros(NT, Nmu, 10);
vec_tau3       = zeros(NT, Nmu, 10);
vec_tau3_ab_cde = zeros(NT, Nmu, 4);
vec_total      = zeros(NT, Nmu, 10);
vec_tau1_base  = zeros(NT, Nmu, 10);
vec_tau1_mu    = zeros(NT, Nmu, 10);

chi_tau1       = zeros(NT, Nmu, 6);
chi_tau3       = zeros(NT, Nmu, 6);
chi_total      = zeros(NT, Nmu, 6);
chi_tau1_base  = zeros(NT, Nmu, 6);
chi_tau1_mu    = zeros(NT, Nmu, 6);

for iT = 1:NT
    for imu = 1:Nmu
        sym1      = squeeze(sym_tau1(iT,imu,:,:,:,:,:));
        sym3      = squeeze(sym_tau3(iT,imu,:,:,:,:,:));
        symt      = squeeze(sym_total(iT,imu,:,:,:,:,:));
        sym1_base = squeeze(sym_tau1_base(iT,imu,:,:,:,:,:));
        sym1_mu   = squeeze(sym_tau1_mu(iT,imu,:,:,:,:,:));

        vec1      = flatten_10_local(sym1);
        vec3      = flatten_10_local(sym3);
        vec3_ab_cde = flatten_tau3_ab_cde_4_local(sym3);
        vect      = flatten_10_local(symt);
        vec1_base = flatten_10_local(sym1_base);
        vec1_mu   = flatten_10_local(sym1_mu);

        vec_tau1(iT,imu,:)      = vec1;
        vec_tau3(iT,imu,:)      = vec3;
        vec_tau3_ab_cde(iT,imu,:) = vec3_ab_cde;
        vec_total(iT,imu,:)     = vect;
        vec_tau1_base(iT,imu,:) = vec1_base;
        vec_tau1_mu(iT,imu,:)   = vec1_mu;

        chi_tau1(iT,imu,:)      = vec10_to_chi6_local(vec1);
        chi_tau3(iT,imu,:)      = vec10_to_chi6_local(vec3);
        chi_total(iT,imu,:)     = vec10_to_chi6_local(vect);
        chi_tau1_base(iT,imu,:) = vec10_to_chi6_local(vec1_base);
        chi_tau1_mu(iT,imu,:)   = vec10_to_chi6_local(vec1_mu);
    end
end

% ------------------------------------------------------------
% Backward-compatible scalar-T output shape
% ------------------------------------------------------------
if scalarT
    vec_tau1_out      = reshape(vec_tau1(1,:,:),      Nmu, 10);
    vec_tau3_out      = reshape(vec_tau3(1,:,:),      Nmu, 10);
    vec_tau3_ab_cde_out = reshape(vec_tau3_ab_cde(1,:,:), Nmu, 4);
    vec_total_out     = reshape(vec_total(1,:,:),     Nmu, 10);
    vec_tau1_base_out = reshape(vec_tau1_base(1,:,:), Nmu, 10);
    vec_tau1_mu_out   = reshape(vec_tau1_mu(1,:,:),   Nmu, 10);

    chi_tau1_out      = reshape(chi_tau1(1,:,:),      Nmu, 6);
    chi_tau3_out      = reshape(chi_tau3(1,:,:),      Nmu, 6);
    chi_total_out     = reshape(chi_total(1,:,:),     Nmu, 6);
    chi_tau1_base_out = reshape(chi_tau1_base(1,:,:), Nmu, 6);
    chi_tau1_mu_out   = reshape(chi_tau1_mu(1,:,:),   Nmu, 6);

    sym_tau1_out      = reshape(sym_tau1(1,:,:,:,:,:,:),      Nmu,2,2,2,2,2);
    sym_tau3_out      = reshape(sym_tau3(1,:,:,:,:,:,:),      Nmu,2,2,2,2,2);
    sym_total_out     = reshape(sym_total(1,:,:,:,:,:,:),     Nmu,2,2,2,2,2);
    sym_tau1_base_out = reshape(sym_tau1_base(1,:,:,:,:,:,:), Nmu,2,2,2,2,2);
    sym_tau1_mu_out   = reshape(sym_tau1_mu(1,:,:,:,:,:,:),   Nmu,2,2,2,2,2);
else
    vec_tau1_out      = vec_tau1;
    vec_tau3_out      = vec_tau3;
    vec_tau3_ab_cde_out = vec_tau3_ab_cde;
    vec_total_out     = vec_total;
    vec_tau1_base_out = vec_tau1_base;
    vec_tau1_mu_out   = vec_tau1_mu;

    chi_tau1_out      = chi_tau1;
    chi_tau3_out      = chi_tau3;
    chi_total_out     = chi_total;
    chi_tau1_base_out = chi_tau1_base;
    chi_tau1_mu_out   = chi_tau1_mu;

    sym_tau1_out      = sym_tau1;
    sym_tau3_out      = sym_tau3;
    sym_total_out     = sym_total;
    sym_tau1_base_out = sym_tau1_base;
    sym_tau1_mu_out   = sym_tau1_mu;
end

% ------------------------------------------------------------
% Pack output
% ------------------------------------------------------------
result = struct();
result.kpoint      = kpoint(:).';
result.mu_list     = mu_col;
result.temperature = T_vec;
result.tau         = tau;
result.selectbands = sel;
result.include_mu_correction = use_mu_correction;

result.chi_tau1  = chi_tau1_out;
result.chi_tau3  = chi_tau3_out;
result.chi_total = chi_total_out;

result.vec_tau1  = vec_tau1_out;
result.vec_tau3  = vec_tau3_out;
result.vec_tau3_ab_cde = vec_tau3_ab_cde_out;
result.vec_total = vec_total_out;

result.sym_tau1  = sym_tau1_out;
result.sym_tau3  = sym_tau3_out;
result.sym_total = sym_total_out;

% Separate bookkeeping for debugging: tau1 = tau1_base + tau1_mu.
result.chi_tau1_base = chi_tau1_base_out;
result.chi_tau1_mu   = chi_tau1_mu_out;
result.vec_tau1_base = vec_tau1_base_out;
result.vec_tau1_mu   = vec_tau1_mu_out;
result.sym_tau1_base = sym_tau1_base_out;
result.sym_tau1_mu   = sym_tau1_mu_out;

end

% ============================================================
% local helpers
% ============================================================

function [local_opts, tabc_args] = parse_local_options(varargin)
local_opts = struct();
local_opts.mu_correction = [];
local_opts.include_mu_correction = true;

keep = true(1, numel(varargin));
i = 1;
while i <= numel(varargin)
    if ischar(varargin{i}) || isstring(varargin{i})
        key = lower(char(varargin{i}));
        if any(strcmp(key, {'mu_correction','mucorrection','qcorr','mu_corr'}))
            if i == numel(varargin)
                error('Missing value after %s.', char(varargin{i}));
            end
            local_opts.mu_correction = varargin{i+1};
            keep(i:i+1) = false;
            i = i + 2;
            continue;
        elseif any(strcmp(key, {'include_mu_correction','use_mu_correction'}))
            if i == numel(varargin)
                error('Missing value after %s.', char(varargin{i}));
            end
            local_opts.include_mu_correction = logical(varargin{i+1});
            keep(i:i+1) = false;
            i = i + 2;
            continue;
        end
    end
    i = i + 1;
end

tabc_args = varargin(keep);
end

function [Q2, Q3] = normalize_mu_correction_local(muCorr, NT, Nmu, scalarT)
if ~isstruct(muCorr) || ~isfield(muCorr, 'Q2') || ~isfield(muCorr, 'Q3')
    error('mu_correction must be a struct with fields Q2 and Q3.');
end

Q2in = muCorr.Q2;
Q3in = muCorr.Q3;

% Accept scalar-T squeezed shapes [Nmu,2,2] and [Nmu,2,2,2].
if ndims(Q2in) == 3
    if size(Q2in,1) ~= Nmu || size(Q2in,2) ~= 2 || size(Q2in,3) ~= 2
        error('mu_correction.Q2 scalar-T shape must be [Nmu,2,2].');
    end
    Q2 = reshape(Q2in, 1, Nmu, 2, 2);
elseif ndims(Q2in) == 4
    if size(Q2in,1) ~= NT || size(Q2in,2) ~= Nmu || size(Q2in,3) ~= 2 || size(Q2in,4) ~= 2
        error('mu_correction.Q2 vector-T shape must be [NT,Nmu,2,2].');
    end
    Q2 = Q2in;
else
    error('mu_correction.Q2 must have shape [Nmu,2,2] or [NT,Nmu,2,2].');
end

if ndims(Q3in) == 4
    if size(Q3in,1) ~= Nmu || size(Q3in,2) ~= 2 || size(Q3in,3) ~= 2 || size(Q3in,4) ~= 2
        error('mu_correction.Q3 scalar-T shape must be [Nmu,2,2,2].');
    end
    Q3 = reshape(Q3in, 1, Nmu, 2, 2, 2);
elseif ndims(Q3in) == 5
    if size(Q3in,1) ~= NT || size(Q3in,2) ~= Nmu || size(Q3in,3) ~= 2 || size(Q3in,4) ~= 2 || size(Q3in,5) ~= 2
        error('mu_correction.Q3 vector-T shape must be [NT,Nmu,2,2,2].');
    end
    Q3 = Q3in;
else
    error('mu_correction.Q3 must have shape [Nmu,2,2,2] or [NT,Nmu,2,2,2].');
end

% If scalar T is requested but user supplied [1,Nmu,...], this is already OK.
% If vector T is requested, scalar-T Q is not enough because Q depends on T.
if ~scalarT && size(Q2,1) ~= NT
    error('For vector T, mu_correction.Q2/Q3 must include the NT dimension.');
end
end

function f1 = Fermi_1_local(x, T)
% f'(x) for f = 1/(exp(x/T)+1)
u  = x ./ (2*T);
cu = cosh(u);
f1 = -(1 ./ (4*T)) ./ (cu.^2);
end

function f2 = Fermi_2_local(x, T)
% f''(x) for f = 1/(exp(x/T)+1)
u  = x ./ (2*T);
cu = cosh(u);
tu = tanh(u);
f2 = (1 ./ (4*T^2)) .* (1 ./ (cu.^2)) .* tu;
end

function vec = flatten_10_local(sym_raw)
% same ordering as the old MATLAB code:
% a = 1,2 ; b<=c<=d<=e
vec = zeros(1,10);
idx = 1;
for a = 1:2
    for b = 1:2
        for c = b:2
            for d = c:2
                for e = d:2
                    vec(idx) = sym_raw(a,b,c,d,e);
                    idx = idx + 1;
                end
            end
        end
    end
end
end


function vec = flatten_tau3_ab_cde_4_local(sym_raw)
% Independent tau3 Hall-sector components after antisym(a,b) and sym(c,d,e):
%   [xyxxx, xyxxy, xyxyy, xyyyy]
% The yx components are the negative partners by construction.
vec = zeros(1,4);
vec(1) = sym_raw(1,2,1,1,1);
vec(2) = sym_raw(1,2,1,1,2);
vec(3) = sym_raw(1,2,1,2,2);
vec(4) = sym_raw(1,2,2,2,2);
end

function chi6 = vec10_to_chi6_local(vec)
% Same legacy post-processing as the old MATLAB code.
chi1 = 4*vec(8)  - vec(1);
chi2 = 6*vec(9)  - 4*vec(2);
chi3 = 4*vec(10) - 6*vec(3);
chi4 = vec(6)    - 4*vec(4);
chi5 = -vec(5);
chi6_val = -vec(7);
chi6 = [chi1, chi2, chi3, chi4, chi5, chi6_val];
end
