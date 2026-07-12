function result = MuCorrection_2D_sumrule(Ham, klist, mu_list, optsParallel, opts, varargin)
% MuCorrection_2D_sumrule
%   Chemical-potential-shift correction tensors after summing all k points:
%
%       Q2_ab(mu)  = -(1/2) * int[dk] sum_n G_ab,n(k) f0'(eps_n-mu)
%                            / int[dk] sum_n f0'(eps_n-mu)
%
%       Q3_abc(mu) = -(2/3) * int[dk] sum_n T_abc,n(k) f0'(eps_n-mu)
%                            / int[dk] sum_n f0'(eps_n-mu)
%
%   Single-k term is evaluated by
%       out = Tabc_k_fast_2D(Ham, kpoint, varargin{:});
%
% Inputs
%   Ham
%   klist                 [Nk,dim]
%   mu_list               [Nmu,1] or [1,Nmu], in the same energy unit as eigvals
%   optsParallel.ncore    number of parallel workers, default 4
%   opts.T                temperature broadening k_B*T, same energy unit as eigvals, default 0.025851991
%                         If your input is Kelvin, convert before calling this function.
%   opts.selectbands      selected band indices, default [] means all bands
%   opts.kweight          [Nk,1] optional k weights, default ones(Nk,1)
%   opts.denom_tol        small denominator threshold, default 1e-14
%   varargin              passed directly to Tabc_k_fast_2D
%
% Outputs
%   result.Q2             scalar T: [Nmu,2,2]; vector T: [NT,Nmu,2,2]
%   result.Q3             scalar T: [Nmu,2,2,2]; vector T: [NT,Nmu,2,2,2]
%   result.denom          summed denominator int f0'
%   result.numG           summed numerator int G_ab f0'
%   result.numT           summed numerator int T_abc f0'
%
% Vectorized flattened outputs:
%   result.Q2_vec: [xx, xy, yx, yy]
%   result.Q3_vec: [xxx, xxy, xyx, xyy, yxx, yxy, yyx, yyy]

if nargin < 4 || isempty(optsParallel)
    optsParallel = struct();
end
if nargin < 5 || isempty(opts)
    opts = struct();
end

nk     = size(klist, 1);
mu_col = mu_list(:);
Nmu    = numel(mu_col);

ncore       = get_opt_local(optsParallel, 'ncore', 4);
T_vec       = get_opt_local(opts, 'T', 0.025851991);
T_vec       = T_vec(:);
NT          = numel(T_vec);
scalarT     = (NT == 1);
selectbands = get_opt_local(opts, 'selectbands', []);
kweight     = get_opt_local(opts, 'kweight', ones(nk,1));
denom_tol   = get_opt_local(opts, 'denom_tol', 1e-14);

tabc_args = varargin;

if numel(kweight) ~= nk
    error('opts.kweight must have length size(klist,1).');
end
kweight = kweight(:);

if any(~isfinite(T_vec)) || any(T_vec <= 0)
    error('opts.T must be finite and positive. Use k_B*T in the same energy unit as eigvals.');
end

% Sums before the final ratio. The common BZ weight cancels in Q, but
% kweight is kept for non-uniform/adaptive meshes.
denom_sum = zeros(NT, Nmu);
numG_sum  = zeros(NT, Nmu, 2, 2);
numT_sum  = zeros(NT, Nmu, 2, 2, 2);

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
        [den_loc, numG_loc, numT_loc] = MuCorrection_ik_2D_local( ...
            Ham, klist(ki,:), mu_col, T_vec, selectbands, tabc_args{:});
        wk = kweight(ki);
        denom_sum = denom_sum + wk .* den_loc;
        numG_sum  = numG_sum  + wk .* numG_loc;
        numT_sum  = numT_sum  + wk .* numT_loc;
    end
else
    for ki = 1:nk
        [den_loc, numG_loc, numT_loc] = MuCorrection_ik_2D_local( ...
            Ham, klist(ki,:), mu_col, T_vec, selectbands, tabc_args{:});
        wk = kweight(ki);
        denom_sum = denom_sum + wk .* den_loc;
        numG_sum  = numG_sum  + wk .* numG_loc;
        numT_sum  = numT_sum  + wk .* numT_loc;
    end
end

Q2 = nan(NT, Nmu, 2, 2);
Q3 = nan(NT, Nmu, 2, 2, 2);

for iT = 1:NT
    for imu = 1:Nmu
        den = denom_sum(iT, imu);
        if isfinite(den) && abs(den) > denom_tol
            Q2(iT, imu, :, :)    = -0.5  .* numG_sum(iT, imu, :, :)    ./ den;
            Q3(iT, imu, :, :, :) = -(2/3) .* numT_sum(iT, imu, :, :, :) ./ den;
        end
    end
end

% Q2 and Q3 are chemical-potential shift tensors with electric-field
% indices only, so enforce the physical electric-field symmetries before
% they are passed into the tau1 correction:
%   Q2_{ab}     = Q2_{ba}
%   Q3_{abc}    = Q3_{(abc)}
Q2_raw = Q2;
Q3_raw = Q3;
Q2 = symmetrize_Q2_local(Q2);
Q3 = symmetrize_Q3_local(Q3);

Q2_vec = flatten_Q2_local(Q2);
Q3_vec = flatten_Q3_local(Q3);

result = struct();
result.mu_list      = mu_col;
result.temperature  = T_vec;
result.selectbands  = selectbands;
result.kweight_sum  = sum(kweight);
result.denom        = denom_sum;
result.numG         = numG_sum;
result.numT         = numT_sum;
result.Q2           = Q2;
result.Q3           = Q3;
result.Q2_raw       = Q2_raw;
result.Q3_raw       = Q3_raw;
result.Q2_vec       = Q2_vec;
result.Q3_vec       = Q3_vec;

% Backward-compatible scalar-T shape.
if scalarT
    result.denom  = reshape(result.denom(1,:), Nmu, 1);        % [Nmu,1]
    result.numG   = reshape(result.numG(1,:,:,:), Nmu, 2, 2);  % [Nmu,2,2]
    result.numT   = reshape(result.numT(1,:,:,:,:), Nmu, 2, 2, 2); % [Nmu,2,2,2]
    result.Q2     = reshape(result.Q2(1,:,:,:), Nmu, 2, 2);    % [Nmu,2,2]
    result.Q3     = reshape(result.Q3(1,:,:,:,:), Nmu, 2, 2, 2); % [Nmu,2,2,2]
    result.Q2_raw = reshape(result.Q2_raw(1,:,:,:), Nmu, 2, 2); % [Nmu,2,2]
    result.Q3_raw = reshape(result.Q3_raw(1,:,:,:,:), Nmu, 2, 2, 2); % [Nmu,2,2,2]
    result.Q2_vec = reshape(result.Q2_vec(1,:,:), Nmu, 4);     % [Nmu,4]
    result.Q3_vec = reshape(result.Q3_vec(1,:,:), Nmu, 8);     % [Nmu,8]
end

end

% ============================================================
% local helpers
% ============================================================

function [den_loc, numG_loc, numT_loc] = MuCorrection_ik_2D_local( ...
    Ham, kpoint, mu_col, T_vec, selectbands, varargin)

out = Tabc_k_fast_2D(Ham, kpoint, varargin{:});

if isempty(selectbands)
    sel = 1:numel(out.eigvals);
else
    sel = selectbands(:).';
end

E = out.eigvals(sel);
E = real(E(:).');                     % [1,nsel]
nsel = numel(sel);

Gten = real(out.BCP(:,:,sel));         % [2,2,nsel]
Tten = real(out.T_tensor(:,:,:,sel));  % [2,2,2,nsel]

Nmu = numel(mu_col);
NT  = numel(T_vec);

den_loc  = zeros(NT, Nmu);
numG_loc = zeros(NT, Nmu, 2, 2);
numT_loc = zeros(NT, Nmu, 2, 2, 2);

for iT = 1:NT
    Tnow = T_vec(iT);
    Emu  = E - mu_col;                 % [Nmu,nsel]
    Fp   = Fermi_1_local(Emu, Tnow);   % [Nmu,nsel], negative

    den_loc(iT,:) = reshape(sum(Fp, 2), 1, Nmu);

    for a = 1:2
        for b = 1:2
            Gab = reshape(Gten(a,b,:), nsel, 1);
            numG_loc(iT,:,a,b) = reshape(Fp * Gab, 1, Nmu);

            for c = 1:2
                Tabc = reshape(Tten(a,b,c,:), nsel, 1);
                numT_loc(iT,:,a,b,c) = reshape(Fp * Tabc, 1, Nmu);
            end
        end
    end
end

end

function f1 = Fermi_1_local(x, T)
% f'(x), f = 1/(exp(x/T)+1), T must be k_B*T in energy units.
u  = x ./ (2*T);
cu = cosh(u);
f1 = -(1 ./ (4*T)) ./ (cu.^2);
end

function val = get_opt_local(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end


function Q2s = symmetrize_Q2_local(Q2)
% Q2: [NT,Nmu,2,2] -> Q2_{(ab)}
Q2s = 0.5 .* (Q2 + permute(Q2, [1 2 4 3]));
end

function Q3s = symmetrize_Q3_local(Q3)
% Q3: [NT,Nmu,2,2,2] -> Q3_{(abc)}
ordered_abcs = [
    1 1 1
    1 1 2
    1 2 2
    2 2 2
];
Q3s = zeros(size(Q3));
for iord = 1:size(ordered_abcs,1)
    abc = ordered_abcs(iord,:);
    plist = unique(perms(abc), 'rows');
    nperm = size(plist,1);

    acc = zeros(size(Q3,1), size(Q3,2));
    for ip = 1:nperm
        a = plist(ip,1);
        b = plist(ip,2);
        c = plist(ip,3);
        acc = acc + Q3(:,:,a,b,c);
    end
    acc = acc ./ nperm;

    for ip = 1:nperm
        a = plist(ip,1);
        b = plist(ip,2);
        c = plist(ip,3);
        Q3s(:,:,a,b,c) = acc;
    end
end
end

function Q2_vec = flatten_Q2_local(Q2)
% Q2: [NT,Nmu,2,2] -> [NT,Nmu,4] = [xx,xy,yx,yy]
NT  = size(Q2,1);
Nmu = size(Q2,2);
Q2_vec = zeros(NT,Nmu,4);
Q2_vec(:,:,1) = Q2(:,:,1,1);
Q2_vec(:,:,2) = Q2(:,:,1,2);
Q2_vec(:,:,3) = Q2(:,:,2,1);
Q2_vec(:,:,4) = Q2(:,:,2,2);
end

function Q3_vec = flatten_Q3_local(Q3)
% Q3: [NT,Nmu,2,2,2] -> [NT,Nmu,8]
% order = [xxx,xxy,xyx,xyy,yxx,yxy,yyx,yyy]
NT  = size(Q3,1);
Nmu = size(Q3,2);
Q3_vec = zeros(NT,Nmu,8);
idx = 1;
for a = 1:2
    for b = 1:2
        for c = 1:2
            Q3_vec(:,:,idx) = Q3(:,:,a,b,c);
            idx = idx + 1;
        end
    end
end
end
