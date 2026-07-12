function chi_mu = ForthHall(Ham, klist, mu_list, optsParallel, opts)
% ForthHall
%   Loop over klist, call ForthHall_ik on each k, accumulate χ^{(4,1)}(μ)
%
% Inputs:
%   Ham, klist, mu_list, parameters, optsParallel.ncore, opts (T, eps...)
%
% Output:
%   chi_mu : 5 x length(mu_list)

arguments
    Ham
    klist double
    mu_list double
    optsParallel.ncore = 4
    opts.T = 50
    opts.eps = 1e-4
    opts.dk = 1e-6
end

nk = size(klist, 1);
nmu = length(mu_list);
eps = opts.eps;
T = opts.T;
dk = opts.dk;
chi_mu = zeros( nmu,12);
% -------- Volume factor and prefactor ----------
volume = dot(cross(Ham.Rm(1,:),Ham.Rm(2,:)), Ham.Rm(3,:));

const_factor = (constants.charge_C^3) ...
    * (volume / Nk_effect) ...
    / (constants.hbar_eV_s^2);
%% --- Parallel setup  -------------------------------------
use_parallel = (optsParallel.ncore > 1);
if use_parallel
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= optsParallel.ncore
        if ~isempty(pool); delete(pool); end
        pool = parpool(optsParallel.ncore);
    end
end

%% --- Loop k ------------------------------------------------
if use_parallel
    parfor ki = 1:nk
        chi_mu_loc = ForthHall_ik_2D_sym(Ham, klist(ki,:), mu_list,dk, T, eps);
        chi_mu = chi_mu + chi_mu_loc;
    end
else
    for ki = 1:nk
        chi_mu_loc = ForthHall_ik_2D_sym(Ham, klist(ki,:), mu_list,dk,T, eps);
        chi_mu = chi_mu + chi_mu_loc;
    end
end
chi_mu = chi_mu * const_factor;
end
