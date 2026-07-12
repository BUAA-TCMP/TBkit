function result = ForthHall_sumrule(Ham, klist, mu_list, optsParallel, opts, varargin)
% ForthHall_sumrule
%   Corrected all-k main loop for ForthHall_ik_2D_sym_sumrule.
%
%   Correct single-k call:
%       r = ForthHall_ik_2D_sym_sumrule( ...
%               Ham, kpoint, mu_list, T, tau, selectbands, varargin{:});
%
% Inputs
%   Ham
%   klist                 [Nk,dim]
%   mu_list               chemical-potential list
%   optsParallel.ncore    default 4
%   opts.T                k_B*T in the same energy unit as eigvals; scalar or vector
%   opts.tau              relaxation time, scalar, default 1
%   opts.selectbands      selected band indices, default [] means all bands
%   opts.kweight          [Nk,1] optional k weights, default ones(Nk,1)
%   opts.Nk_effect        normalization count, default sum(kweight)
%   opts.include_mu_correction  run the two-stage chemical-potential correction,
%                         default true
%   opts.const_factor     if given, use this directly. If empty, use the old convention:
%                         e^3 * cell_volume/Nk_effect / hbar_eV_s^2
%                         Set opts.const_factor = 1 if you only want the raw BZ sum.
%   varargin              passed directly to Tabc_k_fast_2D through ForthHall_ik_2D_sym_sumrule
%
% Outputs
%   result.vec_tau1, result.vec_tau3, result.vec_tau3_ab_cde, result.vec_total
%   result.chi_tau1, result.chi_tau3, result.chi_total
%   result.sym_tau1, result.sym_tau3, result.sym_total

if nargin < 4 || isempty(optsParallel)
    optsParallel = struct();
end
if nargin < 5 || isempty(opts)
    opts = struct();
end

nk         = size(klist, 1);
mu_col     = mu_list(:);
Nmu        = numel(mu_col);
ncore      = get_opt_local(optsParallel, 'ncore', 4);
T_vec      = get_opt_local(opts, 'T', 0.025851991);
T_vec      = T_vec(:);
NT         = numel(T_vec);
scalarT    = (NT == 1);
tau        = get_opt_local(opts, 'tau', 1);
selectbands = get_opt_local(opts, 'selectbands', []);
kweight    = get_opt_local(opts, 'kweight', ones(nk,1));
Nk_effect  = get_opt_local(opts, 'Nk_effect', []);
const_factor = get_opt_local(opts, 'const_factor', []);
include_mu_correction = get_opt_local(opts, 'include_mu_correction', true);

tabc_args = varargin;

% Stage 1: Q2_de and Q3_cde are global quantities and must be accumulated
% over the complete k mesh before evaluating any corrected tau1 integrand.
if include_mu_correction
    mu_correction = MuCorrection_2D_sumrule( ...
        Ham, klist, mu_col, optsParallel, opts, tabc_args{:});
else
    mu_correction = [];
end

if numel(kweight) ~= nk
    error('opts.kweight must have length size(klist,1).');
end
kweight = kweight(:);

if isempty(Nk_effect)
    Nk_effect = sum(kweight);
end
if Nk_effect <= 0
    error('Nk_effect must be positive.');
end

% Keep the old prefactor convention unless opts.const_factor is provided.
% If your klist/units use a different BZ normalization, set opts.const_factor explicitly.
if isempty(const_factor)
    cell_volume = dot(cross(Ham.Rm(1,:), Ham.Rm(2,:)), Ham.Rm(3,:));
    const_factor = (constants.charge_C^3) ...
        * (cell_volume / Nk_effect) ...
        / (constants.hbar_eV_s^2);
end

if scalarT
    vec_tau1_sum  = zeros(Nmu,10);
    vec_tau3_sum  = zeros(Nmu,10);
    vec_tau3_ab_cde_sum = zeros(Nmu,4);
    vec_total_sum = zeros(Nmu,10);

    chi_tau1_sum  = zeros(Nmu,6);
    chi_tau3_sum  = zeros(Nmu,6);
    chi_total_sum = zeros(Nmu,6);

    sym_tau1_sum  = zeros(Nmu,2,2,2,2,2);
    sym_tau3_sum  = zeros(Nmu,2,2,2,2,2);
    sym_total_sum = zeros(Nmu,2,2,2,2,2);
    vec_tau1_base_sum = zeros(Nmu,10);
    vec_tau1_mu_sum   = zeros(Nmu,10);
    chi_tau1_base_sum = zeros(Nmu,6);
    chi_tau1_mu_sum   = zeros(Nmu,6);
    sym_tau1_base_sum = zeros(Nmu,2,2,2,2,2);
    sym_tau1_mu_sum   = zeros(Nmu,2,2,2,2,2);
else
    vec_tau1_sum  = zeros(NT,Nmu,10);
    vec_tau3_sum  = zeros(NT,Nmu,10);
    vec_tau3_ab_cde_sum = zeros(NT,Nmu,4);
    vec_total_sum = zeros(NT,Nmu,10);

    chi_tau1_sum  = zeros(NT,Nmu,6);
    chi_tau3_sum  = zeros(NT,Nmu,6);
    chi_total_sum = zeros(NT,Nmu,6);

    sym_tau1_sum  = zeros(NT,Nmu,2,2,2,2,2);
    sym_tau3_sum  = zeros(NT,Nmu,2,2,2,2,2);
    sym_total_sum = zeros(NT,Nmu,2,2,2,2,2);
    vec_tau1_base_sum = zeros(NT,Nmu,10);
    vec_tau1_mu_sum   = zeros(NT,Nmu,10);
    chi_tau1_base_sum = zeros(NT,Nmu,6);
    chi_tau1_mu_sum   = zeros(NT,Nmu,6);
    sym_tau1_base_sum = zeros(NT,Nmu,2,2,2,2,2);
    sym_tau1_mu_sum   = zeros(NT,Nmu,2,2,2,2,2);
end

use_parallel = (ncore > 1);
if use_parallel
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= ncore
        if ~isempty(pool); delete(pool); end
        parpool(ncore);
    end
end

if use_parallel
    parfor ki = 1:nk
        r = ForthHall_ik_2D_sym_sumrule( ...
            Ham, klist(ki,:), mu_col, T_vec, tau, selectbands, tabc_args{:}, ...
            'mu_correction', mu_correction, ...
            'include_mu_correction', include_mu_correction);
        wk = kweight(ki);

        vec_tau1_sum  = vec_tau1_sum  + wk .* r.vec_tau1;
        vec_tau3_sum  = vec_tau3_sum  + wk .* r.vec_tau3;
        vec_tau3_ab_cde_sum = vec_tau3_ab_cde_sum + wk .* r.vec_tau3_ab_cde;
        vec_total_sum = vec_total_sum + wk .* r.vec_total;

        chi_tau1_sum  = chi_tau1_sum  + wk .* r.chi_tau1;
        chi_tau3_sum  = chi_tau3_sum  + wk .* r.chi_tau3;
        chi_total_sum = chi_total_sum + wk .* r.chi_total;

        sym_tau1_sum  = sym_tau1_sum  + wk .* r.sym_tau1;
        sym_tau3_sum  = sym_tau3_sum  + wk .* r.sym_tau3;
        sym_total_sum = sym_total_sum + wk .* r.sym_total;
        vec_tau1_base_sum = vec_tau1_base_sum + wk .* r.vec_tau1_base;
        vec_tau1_mu_sum   = vec_tau1_mu_sum   + wk .* r.vec_tau1_mu;
        chi_tau1_base_sum = chi_tau1_base_sum + wk .* r.chi_tau1_base;
        chi_tau1_mu_sum   = chi_tau1_mu_sum   + wk .* r.chi_tau1_mu;
        sym_tau1_base_sum = sym_tau1_base_sum + wk .* r.sym_tau1_base;
        sym_tau1_mu_sum   = sym_tau1_mu_sum   + wk .* r.sym_tau1_mu;
    end
else
    for ki = 1:nk
        r = ForthHall_ik_2D_sym_sumrule( ...
            Ham, klist(ki,:), mu_col, T_vec, tau, selectbands, tabc_args{:}, ...
            'mu_correction', mu_correction, ...
            'include_mu_correction', include_mu_correction);
        wk = kweight(ki);

        vec_tau1_sum  = vec_tau1_sum  + wk .* r.vec_tau1;
        vec_tau3_sum  = vec_tau3_sum  + wk .* r.vec_tau3;
        vec_tau3_ab_cde_sum = vec_tau3_ab_cde_sum + wk .* r.vec_tau3_ab_cde;
        vec_total_sum = vec_total_sum + wk .* r.vec_total;

        chi_tau1_sum  = chi_tau1_sum  + wk .* r.chi_tau1;
        chi_tau3_sum  = chi_tau3_sum  + wk .* r.chi_tau3;
        chi_total_sum = chi_total_sum + wk .* r.chi_total;

        sym_tau1_sum  = sym_tau1_sum  + wk .* r.sym_tau1;
        sym_tau3_sum  = sym_tau3_sum  + wk .* r.sym_tau3;
        sym_total_sum = sym_total_sum + wk .* r.sym_total;
        vec_tau1_base_sum = vec_tau1_base_sum + wk .* r.vec_tau1_base;
        vec_tau1_mu_sum   = vec_tau1_mu_sum   + wk .* r.vec_tau1_mu;
        chi_tau1_base_sum = chi_tau1_base_sum + wk .* r.chi_tau1_base;
        chi_tau1_mu_sum   = chi_tau1_mu_sum   + wk .* r.chi_tau1_mu;
        sym_tau1_base_sum = sym_tau1_base_sum + wk .* r.sym_tau1_base;
        sym_tau1_mu_sum   = sym_tau1_mu_sum   + wk .* r.sym_tau1_mu;
    end
end

result = struct();
result.mu_list      = mu_col;
result.temperature  = T_vec;
result.tau          = tau;
result.selectbands  = selectbands;
result.kweight_sum  = sum(kweight);
result.Nk_effect    = Nk_effect;
result.const_factor = const_factor;
result.include_mu_correction = include_mu_correction;
result.mu_correction = mu_correction;

result.vec_tau1  = const_factor .* vec_tau1_sum;
result.vec_tau3  = const_factor .* vec_tau3_sum;
result.vec_tau3_ab_cde = const_factor .* vec_tau3_ab_cde_sum;
result.vec_total = const_factor .* vec_total_sum;

result.chi_tau1  = const_factor .* chi_tau1_sum;
result.chi_tau3  = const_factor .* chi_tau3_sum;
result.chi_total = const_factor .* chi_total_sum;

result.sym_tau1  = const_factor .* sym_tau1_sum;
result.sym_tau3  = const_factor .* sym_tau3_sum;
result.sym_total = const_factor .* sym_total_sum;
result.vec_tau1_base = const_factor .* vec_tau1_base_sum;
result.vec_tau1_mu   = const_factor .* vec_tau1_mu_sum;
result.chi_tau1_base = const_factor .* chi_tau1_base_sum;
result.chi_tau1_mu   = const_factor .* chi_tau1_mu_sum;
result.sym_tau1_base = const_factor .* sym_tau1_base_sum;
result.sym_tau1_mu   = const_factor .* sym_tau1_mu_sum;

% Raw sums before prefactor are kept for debugging and unit checks.
result.raw.vec_tau1  = vec_tau1_sum;
result.raw.vec_tau3  = vec_tau3_sum;
result.raw.vec_tau3_ab_cde = vec_tau3_ab_cde_sum;
result.raw.vec_total = vec_total_sum;
result.raw.chi_tau1  = chi_tau1_sum;
result.raw.chi_tau3  = chi_tau3_sum;
result.raw.chi_total = chi_total_sum;
result.raw.sym_tau1  = sym_tau1_sum;
result.raw.sym_tau3  = sym_tau3_sum;
result.raw.sym_total = sym_total_sum;
result.raw.vec_tau1_base = vec_tau1_base_sum;
result.raw.vec_tau1_mu   = vec_tau1_mu_sum;
result.raw.chi_tau1_base = chi_tau1_base_sum;
result.raw.chi_tau1_mu   = chi_tau1_mu_sum;
result.raw.sym_tau1_base = sym_tau1_base_sum;
result.raw.sym_tau1_mu   = sym_tau1_mu_sum;

end

function val = get_opt_local(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end
