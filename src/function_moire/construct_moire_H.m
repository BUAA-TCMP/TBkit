%% construct_moire_H.m
% Construct the Moiré Hamiltonian from a tight-binding model from HR class.
%
% This function builds the Moiré Hamiltonian Hm given a parent Hamiltonian HR,
% a Moiré potential matrix Vmijs, and mBZ/BZ ratio R.
% The resulting Hamiltonian includes the central block and six surrounding
% translated copies (±Gm shifts), along with coupling terms that connect them, if restore_C3z = True.
%
% Usage:
%   Hm = construct_moire_H(HR, Vmijs, R, restore_C3z)
%
% Inputs:
%   HR           - HR object from TBkit, containing the parent TB Hamiltonian
%   Vmijs        - Moiré coupling potential matrix of size (N × N), typically
%                  block-diagonal, e.g. blkdiag(v1*I, v2*I).
%   R            - Scaling factor between mBZ and BZ (default: 0.25).
%   restore_C3z  - Logical flag, if true, restore C3z rotational symmetry by
%                  coupling between neighboring satellite blocks (default: true).
%
% Outputs:
%   Hm           - Function handle @(kx,ky,kz) returning the Moiré Hamiltonian
%                  matrix of dimension N*(1+6), i.e. center + 6 neighbors.
%
% Notes:
%   - The function builds symbolic expressions in k-space, then converts them
%     into a MATLAB function handle using matlabFunction.
%   - Coupling terms are inserted between central and satellite Hamiltonians,
%     and optionally between satellites to enforce C3 symmetry.
%   - Hermiticity is enforced explicitly at the end.
%  The theoretical formulation of the construction is in the supplementary
%  materials of the paper.

function Hm = construct_moire_H(HR, Vmijs, R, restore_C3z)
arguments
    HR
    Vmijs
    R = 0.25;                 % Ratio mBZ/BZ
    restore_C3z logical = true;
end

% Extract reciprocal lattice vectors from parent HR
Gk = HR.Gk;

% Construct Moiré reciprocal lattice vectors
Gm = R * Gk; 
Gm(3,:) = -Gm(1,:) + Gm(2,:);   % Third vector for hexagonal symmetry
N = HR.WAN_NUM;                 % Number of Wannier orbitals per unit cell

% Define symbolic momentum variables
syms k_x k_y k_z real;

% Parent Hamiltonian as function of k
Hfun = matlabFunction(HR.sym, 'Vars', [k_x, k_y, k_z]);

% Extend Gm to include both +Gm and -Gm shifts
Gm = [Gm; -Gm];

% Construct block-diagonal Hamiltonian: center + 6 neighbors
Hblocks = cell(7,1);
Hblocks{1} = Hfun(k_x, k_y, k_z);  % central block
for i = 1:6
    k_shift = [k_x, k_y, k_z] - Gm(i,:);
    Hblocks{i+1} = Hfun(k_shift(1), k_shift(2), k_shift(3));
end
Hm = blkdiag(Hblocks{:}); 

% Add couplings between central block and each satellite
for i = 1:6
    row_idx = 1:N;
    col_idx = (N*i + 1):(N*(i+1));
    
    % Insert Hermitian coupling
    Hm(row_idx, col_idx) = Vmijs;
    Hm(col_idx, row_idx) = Vmijs';
end

% Optionally restore C3z symmetry by adding couplings between satellites
for i = 1:6
    j = mod(i,6) + 1;  % next neighbor index (1..6 cyclic)
    idx_i = (N*i + 1):(N*(i+1));
    idx_j = (N*j + 1):(N*(j+1));
    
    if restore_C3z
        Hm(idx_i, idx_j) = Vmijs;
        Hm(idx_j, idx_i) = Vmijs';
    end
end

% Enforce Hermiticity explicitly
Hm = (Hm + Hm')/2;

% Convert symbolic Hamiltonian into function handle @(kx,ky,kz)
Hm = matlabFunction(Hm, 'Vars', [k_x, k_y, k_z]);
end
