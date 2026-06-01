function H_hr = horzcat(A,B)
% HORZCAT Horizontal concatenation of HR objects
%
%   H_hr = HORZCAT(A,B) combines two HR objects horizontally,
%   creating a block-diagonal Hamiltonian.
%
%   INPUT ARGUMENTS:
%       A - First HR object or matrix
%       B - Second HR object or matrix
%
%   OUTPUT ARGUMENTS:
%       H_hr - Combined HR object
%
%   NOTES:
%       - Handles HR-HR, HR-matrix, and matrix-HR combinations
%       - Preserves both numerical and symbolic components
%
%   SEE ALSO:
%       HR, blkdiag
%

if isa(A,'HR') && isa(B,'HR')
    H_hr1 = A;
    H_hr2 = B;
    H_hr = A;
    H_hr.vectorL = unique([H_hr1.vectorL;H_hr2.vectorL],'rows');
    H_hr.orbL = [H_hr1.orbL;H_hr2.orbL];
    H_hr.orbL_cart = [H_hr1.orbL_cart;H_hr2.orbL_cart];
    H_hr.quantumL = [H_hr1.quantumL;H_hr2.quantumL];
    H_hr.elementL = [H_hr1.elementL;H_hr2.elementL];
    H_hr.orb_symL = [H_hr1.orb_symL;H_hr2.orb_symL];
    H_hr.HnumL = zeros(H_hr1.WAN_NUM+H_hr2.WAN_NUM,H_hr1.WAN_NUM+H_hr2.WAN_NUM,H_hr.NRPTS);
    H_hr.HcoeL = sym(H_hr.HnumL);
    for i = 1:H_hr.NRPTS
        vector = H_hr.vectorL(i,:);
        [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
        [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
        if seq1 ~= 0 && seq2 ~=0
            if H_hr1.num
                amp = blkdiag(H_hr1.HnumL(:,:,seq1) ,...
                    H_hr2.HnumL(:,:,seq2));
                H_hr = H_hr.set_hop_mat(amp,vector,'set');
            end
            if H_hr1.coe
                amp_sym = blkdiag(H_hr1.HcoeL(:,:,seq1) ,...
                    H_hr2.HcoeL(:,:,seq2));
                H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
            end
        end
    end
elseif isa(A,'HR') && ~isa(B,'HR')
    H_hr1 = A;
    H_hr = H_hr1;
    if isa(B,'sym')
        for i = 1:H_hr.NRPTS
            H_hr.HcoeL(:,:,i) = blkdiag(H_hr.HcoeL(:,:,i),B);
        end
    elseif isa(B,'numeric')
        for i = 1:H_hr.NRPTS
            H_hr.HnumL(:,:,i) = blkdiag(H_hr.HnumL(:,:,i),B);
        end
    else
    end
elseif ~isa(A,'HR') && isa(B,'HR')
    H_hr1 = B;
    H_hr = H_hr1;
    if isa(A,'sym')
        for i = 1:H_hr.NRPTS
            H_hr.HcoeL(:,:,i) = blkdiag( A , H_hr.HcoeL(:,:,i));
        end
    elseif isa(A,'numeric')
        for i = 1:H_hr.NRPTS
            H_hr.HnumL(:,:,i) = blkdiag( A , H_hr.HnumL(:,:,i));
        end
    else
    end
end
end
