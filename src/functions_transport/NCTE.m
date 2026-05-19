function alpha_mu = NCTE(Ham, tensor_index, klist, mu_list, optionsParallel, options,optionsAdapt)
% PARAMETERS
% ----------
% Ham    : TBkit HK/HR object
% tensor_index : [1×3] or [] for all 27 components
% klist  : [Nk × 3] k-point list
% mu_list: [1 × Nμ] chemical potentials
%
% optionsParallel.ncore : number of cores
%
% options.T    : temperature (scalar)
% options.eps  : degeneracy threshold
% options.batch_size : (reserved for future multi-batch integration)
%
% RETURNS
% -------
% alpha_mu :
%     If tensor_index is single (a,b,c): [1 × Nμ]
%     If tensor_index = []:             [27 × Nμ]
%
% ---------------------------------------------------------------
arguments
    Ham TBkit
    tensor_index  double = [];
    klist double = [];
    mu_list double = [];
    optionsParallel.ncore = 4
    options.T = 50 % Kelvin
    options.eps = 1e-4
    options.batch_size = 1e6  % 默认批次大小 100*100*100
    options.Coeffs = [1 1 -1];
    options.eta = 0.025;
    optionsAdapt.AdapteEnable (1,1) logical = false
    optionsAdapt.min_keep (1,1) double {mustBePositive} = 32
    optionsAdapt.mesh_shape double = []
    optionsAdapt.refine_factor (1,1) double {mustBePositive} = 1
    optionsAdapt.Rm = [];
end

% optionscell = namedargs2cell(options);
%% prepare dH_dk
switch class(Ham)
    case "HK"
    case "HR"
        Ham = Ham.tjmti_gen();
end

clear fft;
% ------------- Basic info -----------------------
Nmu = numel(mu_list);
Nk  = size(klist,1);
Nk_effect = Nk*(optionsAdapt.refine_factor^3);
% Output shape:
if isempty(tensor_index)
    alpha_mu = zeros(27, Nmu);
else
    alpha_mu = zeros(1, Nmu);
end

if isvector(klist) && length(klist) ==3
    optionsAdapt.mesh_shape = klist;
    kcube_bulk = krange .* [ ...
    -0.5  -0.5  -0.5   % v0
     1     0     0     % v0 + a1
     0     1     0     % v0 + a2
     0     0     1 ];  % v0 + a3

Nk = klist;

%% ========================= 生成 k 点 ================================
klist = kmeshgen( ...
    H_kp_n.Rm, ...
    kcube_bulk, ...
    "Nk1", Nk1, ...
    "Nk2", Nk1, ...
    "Nk3", Nk1, ...
    "full_edge", false);
end

if isempty(optionsAdapt.Rm) 
    optionsAdapt.Rm = Ham.Rm;
end

T   = options.T;
eps = options.eps;

% -------- Volume factor and prefactor ----------
volume = dot(cross(Ham.Rm(1,:),Ham.Rm(2,:)), Ham.Rm(3,:));

const_factor = (constants.charge_C / T / constants.hbar_eV_s) ...
    * (volume / Nk_effect) ...
    / constants.hbar_eV_s;
%%
% ------------------------------------------------
%   MAIN LOOP OVER k-POINTS
% ------------------------------------------------
kfun = @(kpt) NCTE_k(Ham, tensor_index, kpt, mu_list, T, eps, options.Coeffs, options.eta);

options_par = namedargs2cell(optionsParallel);
options_adp = namedargs2cell(optionsAdapt);

alpha_mu = kloop_accumulate(klist, kfun, alpha_mu, options_par{:}, options_adp{:});

%% -------- Multiply physical prefactor ----------
alpha_mu = alpha_mu * const_factor;
end
