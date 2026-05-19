function accumulator = kloop_accumulate_final( ...
        klist, kfun, accumulator, optionsParallel, optionsAdapt)

% ========================================================================
%   kloop_accumulate_final
%   - 内存最小：不存 A(k)
%   - coarse: 只计算 weight(k)，立即丢弃 A(k)
%   - selection: 挑选重要 k 点
%   - refine: 构造 refined_klist
%   - final accumulate: 仅计算 refined 区域
% ========================================================================

arguments
    klist double
    kfun (1,1) function_handle
    accumulator
    optionsParallel.ncore (1,1) double = 1

    optionsAdapt.AdapteEnable (1,1) logical = true
    optionsAdapt.min_keep (1,1) double = 32
    optionsAdapt.coverage_target (1,1) double = 0.9999
    optionsAdapt.adapt_scaling_factor (1,1) double = 1.0
    optionsAdapt.weight_method (1,:) char {mustBeMember(optionsAdapt.weight_method, ...
        {'maxabs','fro','sumabs','absmean'})} = 'maxabs'

    optionsAdapt.mesh_shape double = []
    optionsAdapt.refine_factor (1,1) double = 1
    optionsAdapt.Rm double = []
end

Nk = size(klist,1);
fprintf("=============================================================\n");
fprintf("[kloop] Nk = %d\n", Nk);

if Nk == 0
    return;
end

%% =======================================================================
% 并行池初始化
% ========================================================================
use_parallel = (optionsParallel.ncore > 1);

if use_parallel
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= optionsParallel.ncore
        if ~isempty(pool), delete(pool); end
        parpool('Processes', optionsParallel.ncore);
    end
end


%% =======================================================================
%  (1) coarse：只计算 weight(k)，不存 A(k)
% ========================================================================
fprintf("[kloop] Coarse evaluating weights (%s)...\n", optionsAdapt.weight_method);
weight = zeros(Nk,1);

t0 = tic;

if use_parallel
    parfor k = 1:Nk
        Ak = kfun(klist(k,:));
        weight(k) = compute_weight(Ak, optionsAdapt.weight_method);
    end
else
    for k = 1:Nk
        Ak = kfun(klist(k,:));
        weight(k) = compute_weight(Ak, optionsAdapt.weight_method);
    end
end

t_coarse = toc(t0);
fprintf("[kloop] coarse time = %.3f s\n", t_coarse);

if ~optionsAdapt.AdapteEnable
    fprintf("[kloop] No adaptation. Direct full accumulation...\n");

    if use_parallel
        parfor k = 1:Nk
            accumulator = accumulator + kfun(klist(k,:));
        end
    else
        for k = 1:Nk
            accumulator = accumulator + kfun(klist(k,:));
        end
    end
    fprintf("[kloop] Full accumulate finished.\n");
    return;
end


%% =======================================================================
%  (2) scale + log + contribution selection
% ========================================================================
asf = optionsAdapt.adapt_scaling_factor;
weight = weight * (asf / max(weight));
weight = log(1 + weight);

total_w = sum(weight);

[sorted_w, ord] = sort(weight, 'descend');
cum_w = cumsum(sorted_w);

idx_cut = find(cum_w >= optionsAdapt.coverage_target * total_w, 1,'first');
idx_cut = max(idx_cut, optionsAdapt.min_keep);

selected_idx = ord(1:idx_cut);
coverage_real = sum(weight(selected_idx))/total_w*100;

fprintf("[kloop][adapt] selected = %d (%.3f%%)\n", ...
    numel(selected_idx), numel(selected_idx)/Nk*100);
fprintf("[kloop][adapt] coverage = %.6f%%\n", coverage_real);

%% =======================================================================
%  (3) refine：构建 refined_klist
% ========================================================================
refine_factor = optionsAdapt.refine_factor;

if refine_factor == 1
    refined_klist = klist(selected_idx,:);
else
    refined_klist = refined_klist_gen( ...
        klist, selected_idx, ...
        optionsAdapt.mesh_shape, refine_factor, optionsAdapt.Rm);
end

N_ref = size(refined_klist,1);
fprintf("[kloop] refined k-count = %d\n", N_ref);


%% =======================================================================
%  (4) final accumulation：仅 refined 区域
% ========================================================================
fprintf("[kloop] Final accumulate...\n");
t1 = tic;

if use_parallel
    parfor k = 1:N_ref
        accumulator = accumulator + kfun(refined_klist(k,:));
    end
else
    for k = 1:N_ref
        accumulator = accumulator + kfun(refined_klist(k,:));
    end
end

fprintf("[kloop] Final accumulate time = %.3f s\n", toc(t1));
fprintf("=============================================================\n");

end


%% =======================================================================
%  Weight function (No A_list stored)
%% =======================================================================
function w = compute_weight(Ak, method)
switch method
    case 'maxabs'
        w = max(abs(Ak(:)));
    case 'fro'
        w = norm(Ak,'fro');
    case 'sumabs'
        w = sum(abs(Ak(:)));
    case 'absmean'
        w = mean(abs(Ak(:)));
    otherwise
        error("Unknown weight method");
end
end


