
function H_hr_out = symmetrization(H_hr, OperObj, opt, options)
% SYMMETRIZATION   Symmetrize the real-space Hamiltonian using group theory.
%
%   Implements:
%     H̃^{ij}[R′] = (1/|G|) ∑_{g∈G} ∑_{l,m}
%                      D_il(g) H^{lm}[S_g^{-1}(R′ - T_{ij}^{ml})] D_mj(g^{-1})
%
%     with T_{ij}^{ml} = S_g(t_m - t_l) - (t_j - t_i), integer lattice vector

arguments
    H_hr
    OperObj
    opt (1,1) string {mustBeMember(opt,["numerical","symbolic"])} = "numerical"
    options.generator = 1
    options.center = [0,0,0];
    options.Ugen = false;
    options.Accuracy = 1e-5;
    options.Thershold = 1e-4;
    options.associate_orbL  = [];
    options.ncore = 1;
end
if options.ncore == 1
    ParallelFlag = false;

else
    ParallelFlag = true;

    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= options.ncore
        if ~isempty(pool)
            delete(pool);
        end
        pool = parpool(options.ncore);
    end

end
tol = options.Accuracy;

%% Step 0: Optionally generate representation matrices U(g)
if options.Ugen
    try
        BasisFunction = BasisFunc(H_hr);
        OperObj = OperObj.Ugen(BasisFunction, 'Rm', H_hr.Rm, 'center', options.center);
    catch
        warning("Basis generation failed, continue with original OperObj.");
    end
end

nSym = numel(OperObj);

%% Ensure Hamiltonian is in 'list' form
if ~strcmp(H_hr.Type, 'list')
    H_hr = H_hr.rewrite();
end

%% Symbolic mode not implemented
if opt ~= "numerical"
    warning("Symbolic mode not implemented.");
    H_hr_out = H_hr;
    return;
end

%% Extract Hamiltonian data
if options.Thershold == -1
    HnumL      = H_hr.HnumL;
    vectorList = H_hr.vectorL; % [R(3), l, m]
else
    mask       = abs(H_hr.HnumL) > options.Thershold;
    HnumL      = H_hr.HnumL(mask);
    vectorList = H_hr.vectorL(mask,:); % [R(3), l, m]
end

orbL       = H_hr.orbL;
Nbands     = H_hr.WAN_NUM;

%% Precompute fast hash map for row lookup
keyFunc = @(row) sprintf('%d_%d_%d_%d_%d', row);  % encode [Rx,Ry,Rz,l,m] to string
indexMap = containers.Map('KeyType','char','ValueType','int32');
for idx = 1:size(vectorList,1)
    indexMap(keyFunc(vectorList(idx,:))) = idx;
end


%% === Main loop over symmetry sets ===
for n = 1:nSym
    fprintf('*** Applying symmetry set (%d/%d) ***\n', n, nSym);
    SymSet = OperObj(n).generate_group();
    nOps   = numel(SymSet);

    HnumL_sym = zeros(size(HnumL,1), nOps);
for g = 1:nOps
    gOp = SymSet(g);

    % Identity → copy directly
    if gOp.identity == gOp
        HnumL_sym(:, g) = HnumL;
        continue
    end

    % Extract symmetry data
    Sg = round(inv(gOp.Rf));
    Dg = gOp.U;
    Dg_inv = inv(Dg);
    nBands = size(Dg, 1);

    [Tij_ml, ilmj_index, U_factor] = Tij2lm(orbL, Sg, Dg, Dg_inv, nBands, tol, options.associate_orbL);

    newRows = [];
    newVals = [];

    % === Loop over Hamiltonian elements H^{lm}[R] ===
    for in = 1:numel(HnumL)
        R = vectorList(in, 1:3);
        l = vectorList(in, 4);
        m = vectorList(in, 5);

        maskIJ = (ilmj_index(:, 2) == l) & (ilmj_index(:, 3) == m);
        if ~any(maskIJ), continue; end

        ij_pairs = ilmj_index(maskIJ, [1, 4]);
        U_vals = U_factor(maskIJ);
        Tlist = Tij_ml(maskIJ, :);

        % --- 向量化 Rprime 和 val ---
        % Rprime_all = R + zeros(size(Tlist,1),1);       % 初始化
        Rprime_all = repmat(R, [size(Tlist, 1), 1]) * Sg + Tlist; % 向量化批量计算

        val_all = HnumL(in) * U_vals; % 批量计算

        for k = 1:numel(val_all)
            if gOp.conjugate, val_all(k) = conj(val_all(k)); end
            if gOp.antisymmetry, val_all(k) = -val_all(k); end
        end

        % --- 累加或新增行 ---
        for k = 1:size(ij_pairs, 1)
            i = ij_pairs(k, 1);
            j = ij_pairs(k, 2);

            rowCandidate = [Rprime_all(k, :), i, j];
            key = keyFunc(rowCandidate);

            if isKey(indexMap, key)
                idx = indexMap(key);
                HnumL_sym(idx, g) = HnumL_sym(idx, g) + val_all(k);
            else
                newRows = [newRows; rowCandidate];
                newVals = [newVals; val_all(k)];
            end

        end

    end

    % Append new rows
    if ~isempty(newRows)
        nNew = size(newRows, 1);
        vectorList = [vectorList; newRows];
        HnumL = [HnumL; zeros(nNew, 1)];
        HnumL_sym = [HnumL_sym; zeros(nNew, nOps)];

        for nr = 1:nNew
            key = keyFunc(newRows(nr, :));
            idx = size(vectorList, 1) - nNew + nr;
            indexMap(key) = idx;
            HnumL_sym(idx, g) = HnumL_sym(idx, g) + newVals(nr);
        end

    end

end
   % 5. Average over group
    HnumL = sum(HnumL_sym,2) / nOps;
end
%% Build output
H_hr_out = H_hr;
H_hr_out.HnumL   = HnumL;
H_hr_out.vectorL = vectorList;
end

function [valid_T, valid_ilmjL, valid_UinvU] = Tij2lm(orbL, Rf, U, invU, Nbands, eta, associate_orbL)
%TIJ2LM  Compute valid (i,l,m,j) index quadruples and corresponding translations
%         according to group action formula.
%
%   Input:
%     orbL   : orbital positions (Norb x 3)
%     Rf     : inverse rotation matrix S_g^{-1} in lattice basis (3x3)
%     U      : representation matrix D(g)   (Nbands x Nbands)
%     invU   : representation matrix D(g^-1) (Nbands x Nbands)
%     Nbands : number of bands (size of U)
%     eta    : numerical tolerance for integer lattice vectors
%     associate_orbL : (optional) alternative orbital list (unused here)
%
%   Output:
%     valid_T      : integer lattice translations T_{ij}^{ml} (K x 3)
%     valid_ilmjL  : index quadruples [i, l, m, j] (K x 4)
%     valid_UinvU  : prefactor D_il(g) * D_mj(g^-1) (K x 1)
%
%   Optimized version: uses sparse index extraction, avoids loops, minimizes allocations.

arguments
    orbL
    Rf
    U
    invU
    Nbands
    eta = 1e-6;
    associate_orbL = [];
end

%% Step 1: Find all nonzero elements of U and invU
% This already encodes allowed (i,l) and (m,j) pairs.
[Possible_i, Possible_l, U_vals]     = find(U);
[Possible_m, Possible_j, invU_vals]  = find(invU);

% Number of nonzero entries
nU    = numel(U_vals);
nInvU = numel(invU_vals);

%% Step 2: Build all combinations of (i,l) with (m,j)
% Cartesian product via kron/reshape (avoids loops)
[grid_i, grid_m] = ndgrid(1:nU, 1:nInvU);
grid_i = grid_i(:);
grid_m = grid_m(:);

ilmj = [ Possible_i(grid_i), ...
    Possible_l(grid_i), ...
    Possible_m(grid_m), ...
    Possible_j(grid_m) ];   % (K x 4)

% Prefactors: D_il(g) * D_mj(g^-1)
U_factor = U_vals(grid_i) .* invU_vals(grid_m);

%% Step 3: Compute translation vectors T_{ij}^{ml}
% T = (t_m - t_l)*Rf + (t_i - t_j)
T_raw = (orbL(ilmj(:,3),:) - orbL(ilmj(:,2),:)) * Rf + ...
    orbL(ilmj(:,1),:) - orbL(ilmj(:,4),:);

% Keep only integer-lattice ones (within tolerance)
isLattice = all(abs(T_raw - round(T_raw)) < eta, 2);

valid_ilmjL = ilmj(isLattice,:);
valid_T     = round(T_raw(isLattice,:));
valid_UinvU = U_factor(isLattice);

end




% function [valid_T, valid_ilmjL,valid_UinvU] = Tij2lm(orbL, Rf,U,invU,Nbands,eta,associate_orbL)
% arguments
%     orbL
%     Rf
%     U
%     invU
%     Nbands
%     eta = 1e-6;
%     associate_orbL = [];
% end
% % Optimized TIJ2LM with vectorization
% % if isempty(associate_orbL)
% %     orblist2 = H_hr.orbL;
% % else
% %     orblist1= H_hr.orbL;
% %     orblist2 = H_hr.associate_orbL;
% % end
% Norb = size(orbL, 1);
%
% [I,J] = ndgrid(1:Nbands,1:Nbands);
% I = I(:); J = J(:);
%
% % Possible_ilL = [I,J];
% % Possible_mjL = [I,J];
%
%
% NonZeroUIndex = find(U~=0);
% [Possible_iL,Possible_lL] = ind2sub([Nbands,Nbands],NonZeroUIndex);
%
%
% NonZeroInvUIndex = find(invU~=0);
% [Possible_mL,Possible_jL] = ind2sub([Nbands,Nbands], NonZeroInvUIndex);
%
% Possible_ilL = [Possible_iL,Possible_lL];
%
% Possible_mjL = [Possible_mL,Possible_jL];
% %  Possible_mjL =U(sub2ind([Nbands,Nbands], I,J)) ~=0 ;
% %  % invU(sub2ind([Nbands,Nbands], Possible_ilmjL(:,3),Possible_ilmjL(:,4))) ~=0 ;
% %
% Possible_ilmjL = gridhorzcat(Possible_ilL,Possible_mjL);
%
%
% % Possible_ilmjL = Possible_ilmjL(HasValue,:);
% % i l m j
% % 1 2 3 4
% %
% %T_ij_ml = (m - l)-(j -i) ; 3 2   + 1
% T_ij_ml = double( (orbL(Possible_ilmjL(:,3), :) - orbL(Possible_ilmjL(:,2), :) )* Rf + ...
%     orbL(Possible_ilmjL(:,1), :) - orbL(Possible_ilmjL(:,4), :) );
% lattice_mask = all(abs(T_ij_ml - round(T_ij_ml)) < eta, 2);
% valid_ilmjL = Possible_ilmjL(lattice_mask, :);
% valid_T     = round(T_ij_ml(lattice_mask, :));
%
% valid_U  = U(sub2ind([Nbands,Nbands],valid_ilmjL(:,1),valid_ilmjL(:,2)));
% valid_invU  = invU(sub2ind([Nbands,Nbands],valid_ilmjL(:,3),valid_ilmjL(:,4)));
%
% valid_UinvU = valid_U.*valid_invU;
%
% end
% %
%
%

