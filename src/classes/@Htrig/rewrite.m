function H_htrig2 = rewrite(H_htrig,mode)
%REWRITE Convert between different representations of a Htrig object.
%
%   H_htrig2 = REWRITE(H_htrig) converts the representation of the Htrig
%   object H_htrig to its alternative form. The function automatically
%   chooses the target representation based on the current Type:
%       'sincos' -> 'exp'
%       'exp'    -> 'sincos'
%       'mat'    -> 'list'
%       'list'   -> 'mat'
%
%   H_htrig2 = REWRITE(H_htrig, mode) converts H_htrig into the specified
%   representation mode. Valid modes are:
%       'sincos' - Sine-cosine symbolic expression (NOT IMPLEMENTED)
%       'exp'    - Exponential symbolic expression (e.g., exp(±i·k·a))
%       'mat'    - Dense matrix (3D array) representation
%       'list'   - Sparse list representation
%
%   Input:
%       H_htrig - An object of class Htrig
%       mode    - Target representation type (optional)
%
%   Output:
%       H_htrig2 - A new Htrig object in the desired representation
%
%   Example:
%       H2 = rewrite(H1, 'mat');   % Convert to dense matrix form
%
%   See also: Htrig, setup, NumOrCoe, split_sym_eq

arguments
    H_htrig Htrig;
    mode {mustBeMember(mode,{'sincos','exp','mat','list',''})}= '';
end

if nargin < 2 || strcmp(mode , '')
    if strcmp(H_htrig.Type,'sincos')
        mode = 'exp';
    elseif strcmp(H_htrig.Type,'exp')
        mode = 'sincos';
    elseif strcmp(H_htrig.Type,'mat')
        mode = 'list';
    elseif strcmp(H_htrig.Type,'list')
        mode = 'mat';
    end
end

H_htrig2 = H_htrig;
VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);

switch mode
    case 'exp'
        % Convert sin/cos form to exp(i*k·a) form
        HsymL_trig_tmp  = combine(rewrite(H_htrig.HsymL_trig,'exp'));
        H_htrig2.HsymL_trig = sym([]);
        H_htrig2.HsymL_trig_bk = [exp(1i*VarUsing) exp(-1i*VarUsing)];
        H_htrig2.HcoeL = sym([]);
        H_htrig2.Type = 'exp';
        count = 0;
        for i = 1:H_htrig.Kinds
            [coeff_trig, symvar_list_trig, H_htrig2] = split_sym_eq(H_htrig2, HsymL_trig_tmp(i));
            for j = 1:numel(coeff_trig)
                count = count + 1;
                k_cell{count} = symvar_list_trig(j);
                mat_cell{count} = H_htrig.HcoeL(:,:,i);
                Var_cell{count} = coeff_trig(j);
            end
        end
        H_htrig2 = H_htrig2.setup(Var_cell, k_cell, mat_cell);

    case 'sincos'
        % Placeholder for future sincos implementation

    case 'mat'
        % Convert list form to 3D matrix form
        [~,~,H_htrig2] = H_htrig2.NumOrCoe();
        WANNUM = H_htrig2.Basis_num;
        if H_htrig2.num
            [vectorList,~,icL] = unique(H_htrig2.HsymL_numL(:,1:H_htrig.Dim),'rows');
            KINDS = size(vectorList,1);
            HnumLtmp = zeros(WANNUM,WANNUM,KINDS);
            sizemesh = [WANNUM,WANNUM,KINDS];
        else
            [vectorList,~,icL] = unique(H_htrig2.HsymL_coeL(:,1:H_htrig.Dim),'rows');
            KINDS = size(vectorList,1);
            HcoeLtmp = sym(zeros(WANNUM,WANNUM,KINDS));
            sizemesh = [WANNUM,WANNUM,KINDS];
        end

        if H_htrig2.num
            iL = double(H_htrig2.HsymL_numL(:,H_htrig.Dim+1));
            jL = double(H_htrig2.HsymL_numL(:,H_htrig.Dim+2));
            indL = sub2ind(sizemesh,iL,jL,icL);
            HnumLtmp(indL) = H_htrig2.HnumL;
            H_htrig2.HnumL = HnumLtmp;
            H_htrig2.HsymL_numL = vectorList;
        else
            iL = double(H_htrig2.HsymL_coeL(:,H_htrig.Dim+1));
            jL = double(H_htrig2.HsymL_coeL(:,H_htrig.Dim+2));
            indL = sub2ind(sizemesh,iL,jL,icL);
            HcoeLtmp(indL) = H_htrig2.HcoeL;
            H_htrig2.HcoeL = HcoeLtmp;
            H_htrig2.HsymL_coeL = vectorList;
        end
        H_htrig2.Type = 'mat';

    case 'list'
        % Convert 3D matrix form to list form
        if strcmp(H_htrig2.Type,'sincos') || strcmp(H_htrig2.Type,'exp')
            [num_label,coe_label,H_htrig2] = H_htrig2.NumOrCoe();
        else
            [num_label,coe_label,H_htrig2] = H_htrig2.NumOrCoe();
            if H_htrig2.num
                if isvector(H_htrig2.HnumL)
                    warning('May not need to rewrite. Do nothing');
                    return;
                end
                NRPTS_ = numel(H_htrig2.HnumL);
                sizeHcoeL = size(H_htrig2.HnumL);
                HnumLtmp = reshape(H_htrig2.HnumL,[NRPTS_,1]);
                [iL,jL,kL] = ind2sub(sizeHcoeL,1:NRPTS_);
                vectorList = [H_htrig2.HsymL_numL(kL,:),iL.',jL.'];
                H_htrig2.HnumL = HnumLtmp;
                H_htrig2.HsymL_numL = vectorList;
            elseif coe_label && ~num_label
                if isvector(H_htrig2.HcoeL)
                    warning('May not need to rewrite. Do nothing');
                    return;
                end
                NRPTS_ = numel(H_htrig2.HcoeL);
                sizeHcoeL = size(H_htrig2.HcoeL);
                HcoeLtmp = reshape(H_htrig2.HcoeL,[NRPTS_,1]);
                [iL,jL,kL] = ind2sub(sizeHcoeL,1:NRPTS_);
                vectorList = [H_htrig2.HsymL_coeL(kL,:),iL.',jL.'];
                H_htrig2.HcoeL = HcoeLtmp;
                H_htrig2.HnumL = zeros(size(HcoeLtmp));
                H_htrig2.HsymL_coeL = vectorList;
            end
        end
        H_htrig2.Type = 'list';

    otherwise
        % Unrecognized mode, do nothing
end
end
