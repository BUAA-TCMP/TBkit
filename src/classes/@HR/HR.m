classdef HR < TBkit & matlab.mixin.CustomDisplay

    % HR Tight-binding Hamiltonian class with real-space representation
    %
    %   The HR class represents tight-binding Hamiltonians in real space,
    %   supporting both numerical and symbolic coefficients, multiple storage
    %   formats, and various operations including symmetry applications,
    %   k-space conversions, and nanostructure generation.
    %
    %   PROPERTIES:
    %       vectorL - R-vector list
    %       HnumL - Numerical Hamiltonian
    %       HcoeL - Symbolic Hamiltonian
    %       Type - Storage type ('mat','list','sparse')
    %       overlap - Include overlap matrix (logical)
    %
    %   METHODS:
    %       Construction:
    %           from_POSCAR_SE, from_wannier90, from_Hstruct
    %       Manipulation:
    %           set_hop, add_empty_one, reseq, rewrite
    %       Conversion:
    %           HR2HK, HR2Htrig, HR2Hckt
    %       Symmetry:
    %           symmetrize, dualize, hermitize
    %       Visualization:
    %           show, printout
    %       Nanostructures:
    %           Hnanowire_gen, slab, surf
    %
    %   SEE ALSO:
    %       HK, Htrig, Hckt
    %

    properties
        vectorL ;
        HnumL   ;
        HcoeL   ;
    end
    properties
        Duality_vector_dist;
    end
    properties
        Type;
        overlap logical= false;
        num logical= false;
        coe logical= true;
        % soc logical= false;
        AvectorL;
        BvectorL;
        CvectorL;
        vectorhopping = false;
    end
    properties (Transient,Hidden = true)
        R_vector_dist;
    end
    properties(Dependent = true)
        NRPTS    ;
        WAN_NUM  ;
        Line_000 ;
        homecell ;
    end
    properties(Dependent = true,Hidden = true)
    end
    properties (Hidden = true)
        nn_store_smart   ;
        nn_sparse_n      ;
        Atom_store_smart ;
        Rnn_map          ;
        vectorL_map      ;
        % try mex
        mex_initialized = false;
        mex_handle = [];
    end
    methods (Access = protected)
        propgrp = getPropertyGroups(~)
    end
    %% Construction
    methods
        function H_hr = HR(WAN_NUM,vectorL,options,propArgs)
            % HR H_hr = HR(WAN_NUM,vectorL) construct a empty TB obj
            %   H_hr = HR() construct a 4-orbs empty TB obj
            %   H_hr = HR(WAN_NUM) construct a WAN_NUM-orbs empty TB obj
            %   H_hr = HR(WAN_NUM,vectorL) construct a WAN_NUM-orbs with vectorL H_R empty TB obj
            %   H_hr = HR(WAN_NUM,vectorL,HnumL,HcoeL) construct a TB obj
            %   with full information
            %   H_hr = HR(WAN_NUM,vectorL,HnumL,HcoeL,Type) indicate the
            %   Type of this TB obj.
            %   See also FROM_HSTRUCT, FROM_HSPARSE, FROM_HDF5, FROM_WANNIER90.

            % ----------- nargin ----------
            arguments
                WAN_NUM double{mustBeInteger} =4;
                vectorL = ([]);
                options.HnumL double=[]   ;
                options.HcoeL sym=sym([]) ;
                options.Type char = 'mat' ;
                options.overlap logical = false;
                options.sym = true;
                propArgs.?TBkit;
            end
            propArgsCell = namedargs2cell(propArgs);
            H_hr = H_hr@TBkit(propArgsCell{:});
            H_hr.Basis_num = WAN_NUM;
            % -------------check---------------
            Type = options.Type;
            % vectorL
            if strcmp(Type,'list')
                if nargin < 2
                    vectorL = ([zeros(1,H_hr.Dim),WAN_NUM,WAN_NUM]);
                end
            else
                if size(vectorL,2) == 1
                    switch vectorL
                        case 1
                            tmp_vectorL = zeros(1,H_hr.Dim);
                        case 3
                            tmp_vectorL = [-1,0,0;0 0 0;1 0 0 ];
                        case 5
                            tmp_vectorL = [1,0,0;0 1 0;0 0 0;0 -1 0;-1 0 0];
                        case 7
                            tmp_vectorL = [0 0 1;1,0,0;0 1 0;0 0 0;0 -1 0;-1 0 0;0 0 -1];
                        case 9
                            tmp_vectorL = [1,0,0;0 1 0;0 0 0;0 -1 0;-1 0 0;...
                                1 1 0;1 -1 0;-1 1 0;-1 -1 0;];
                        case 11
                            tmp_vectorL = [0,0,0];
                        case 27
                            tmp_vectorL = [0,0,0];
                    end
                    vectorL = (tmp_vectorL);
                end
            end
            if isempty(vectorL)
                vectorL = zeros(1,H_hr.Dim);
            end
            % HnumL
            if isempty(options.HnumL)
                [NRPTS,~]  = size(vectorL);
                if strcmp(Type,'sparse')
                    HnumL{NRPTS} = sparse(WAN_NUM,WAN_NUM);
                    for i =1:NRPTS-1
                        HnumL{i} = sparse(WAN_NUM,WAN_NUM);
                    end
                elseif strcmp(Type,'list')
                    HnumL = zeros(1,1);
                else
                    HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
                end
            else
                HnumL = options.HnumL;
            end
            % HcoeL
            if options.sym
                if isempty(options.HcoeL)
                    if strcmp(Type,'sparse')
                        HcoeL = sym([]);
                    else
                        HcoeL = sym(HnumL);
                    end
                else
                    HcoeL = options.HcoeL;
                end
            else
                HcoeL = options.HcoeL;
            end

            %
            % H_hr.NRPTS   =  NRPTS; % the total number of H(Rn)
            % H_hr.WAN_NUM =  WAN_NUM ; % the num of wannier like orbs
            H_hr.vectorL = vectorL; % the Rvector list
            H_hr.HnumL = HnumL  ; % Hnum_list
            H_hr.HcoeL = HcoeL  ; % Hcoe_list
            H_hr.Type  = Type   ;
            %             if strcmp(Type,'sparse')
            %                 H_hr.orbL  = sparse(WAN_NUM,3);
            %             else
            %                 H_hr.orbL  = zeros(WAN_NUM,3);
            %             end
            H_hr.orbL  = zeros(WAN_NUM,H_hr.Dim);
            H_hr.overlap = options.overlap;
            if options.overlap
                %H_hr = HR(WAN_NUM,vectorL,options,propArgs)
            end
            %H_hr.Duality_vector_dist = containers.Map('KeyType','double','ValueType','double');
        end
    end
    methods (Static)
        H_hr = from_Hstruct(Hstruct)
        H_hr = from_Hsparse(Hsparse)
        H_hr = from_POSCAR_SE(POSCAR_file,options)
        H_hr = from_hdf5(filename)
        H_hr = from_wannier90(filename,Type,options)
    end
    methods(Static,Hidden,Access= protected)
        [dataArray,NRPT_list,NRPTS,NUM_WAN]=hrdat_read(filename)
    end
    %% setup
    methods
        H_hr = add_empty_one(H_hr,vector)
        H_hr = expand_empty_one(H_hr,orbOne,QuantumOne,elementOne)
        H_hr = set_hop(H_hr,amp,hi,hj,vector_list,mode)
        H_hr = set_hop_mat(H_hr,amp,vector,mode)
        H_hr = set_hop_single(H_hr,amp,hi,hj,vector,mode)
    end
    %% get
    % ----------------  get property method --------------------
    methods
        function NRPTS = get.NRPTS(H_hr)
            NRPTS = size(H_hr.vectorL,1);
        end
        function WAN_NUM = get.WAN_NUM(H_hr)
            if strcmp(H_hr.Type,'sparse')
                WAN_NUM = size(H_hr.HnumL{1},1);
            elseif strcmp(H_hr.Type,'list')
                try
                    [Sparse_vector,~,~] = unique(H_hr.vectorL(:,H_hr.Dim+1:H_hr.Dim+2),'rows');
                    WAN_NUM = double(max(Sparse_vector(:)));
                catch
                    WAN_NUM = 0;
                end
            else
                %WAN_NUM = max(size(H_hr.HnumL,1),size(H_hr.HcoeL,1));
                if H_hr.num && ~H_hr.coe
                    WAN_NUM = size(H_hr.HnumL,1);
                elseif ~H_hr.num && H_hr.coe
                    WAN_NUM = size(H_hr.HcoeL,1);
                elseif H_hr.num && H_hr.coe
                    WAN_NUM = max(size(H_hr.HnumL,1),size(H_hr.HcoeL,1));
                else
                    WAN_NUM = H_hr.Basis_num;
                end
            end
        end
        function Line_000 = get.Line_000(H_hr)
            vector = [0,0,0];
            [~,Line_000]=ismember((vector),H_hr.vectorL(:,1:H_hr.Dim),'rows');
        end
        function homecell = get.homecell(H_hr)
            if strcmp(H_hr.Type,'list')
                homecell = H_hr.list([0,0,0]);
            else
                if H_hr.coe
                    homecell = H_hr.HcoeL(:,:,H_hr.Line_000);
                else
                    homecell = H_hr.HnumL(:,:,H_hr.Line_000);
                end
            end
        end
    end
    methods
        Type = type(H_hr)
        vectorSeq = Getvector(H_hr,vector)
        H_hr = autohermi(H_hr,mode,options)
    end
    %% overload
    methods
        H_hr = plus(A,B)
        H_hr = minus(A,B)
        H_hr = uminus(H_hr)
        H_hr = times(A,B)
        H_hr = mrdivide(A,B)
        H_hr = premtimes(A,B)
        H_hr = mtimes(A,B)
        H_hr = power(H_hr,b)
        H_hr = mpower(A,B)
        varargout = eq(H_hr1,H_hr2)
        logicalal_num = ne(H_hr1,H_hr2)
        H_hr = ctranspose(H_hr)
        H_hr = transpose(H_hr)
        H_hr = horzcat(A,B)
        H_hr = vertcat(A,B)
        C = gt(B,A)
        C = lt(A,B)
        C = le(A,B)
        H_hr = kron(A,B)
        H_hr = conj(H_hr)
        H_hr = full(H_hr)
        H_hr = sparse(H_hr)
        [H_hr,EQL] = subs(H_hr,varargin)
        H_hr = simplify(H_hr,Accuracy)
        H_hr = filter(H_hr,Accuracy)
        Hsym = sym(H_hr,options)
        H_hr = sum(H_hr_list)
        [H_hr,Sublist,Unique_term] = unique(H_hr,seed,checklist,Accuracy)
    end
    %% reshape 
    methods
        H_hr = reseq(H_hr,wan_list,nrpt_list,nrpt_list_S)
        H_hr = reseq_spin_basis(H_hr,old2new)        
        H_hr = cut_orb(H_hr,rm_list,options)
        H_hr = clean(H_hr,WANNUM)
        H_hr = project(H_hr,BASIS_MAT)
        H_hr = charalize(H_hr)
        H_hr = ForceToMat(H_hr)
        H_hr = ForceTosparse(H_hr)
        H_hr = ForceTolist(H_hr)
        H_hr = ForceToType(H_hr,Type)
        H_hr = OpenBoundary(H_hr,OBC_list)
        H_hr = translation(H_hr,translation_vector,options)
    end

    methods
        H_hr = enlarge(H_hr,dir,amp)
        H_hr = addorb(H_hr,orblist,options)
        H_hr = add_orb(H_hr,hop_struct,orbOne,QuantumOne,elementOne)
        H_hr = deltarule(H_hr,level_cut,mode,options)
        H_hr = alpharule(H_hr,level_cut,mode,options)
        H_hr = shift_Fermi_energy(H_hr, Efermi)
        dH_dk_xyz = dH_dk(H_hr, kpoint)
        [H_soc_sym, lambda_syms] = SOC_on_site_gen(H_hr)
    end
    %% expand
    methods
        H_hr = Hnanowire_gen(H_hr,Nslab,np,vacuum_mode,options);
        H_hr = cut_piece(H_hr,repeatnum,fin_dir,glue_edges,vacuum_mode)
        H_hr = supercell_hr(H_hr,Ns,options)
        H_hr = unfold_hr(H_hr,Ns,options)
        [sc_orb,sc_vec,sc_elementL,sc_quantumL] = supercell_orb(H_hr,Ns,Accuracy)
        [pc_orb,pc_orbL_full,pc_elementL,pc_quantumL,orb_id_L,pc_orb_id_L,pc_orb_selectL] = unfold_orb(H_hr,Ns,Accuracy,orb_id_L)
        H_hr = descritize(H_hr,Nslab,options)
        orbital_out = nanowire_orb(H_hr,fin_dir,vacuum_mode,options)
        H_hr = supercell(H_hr,Ns,filename,Rm,sites,Atom_name,Atom_num,findir)
        [H_hr_out,H_hr_pi_plus,H_hr_pi_minus] = realmap(H_hr)
    end
    %%
    methods
        [H_hr_R,OperObj] = rotate(H_hr,OperObj)
        H_hr = symmetrization(H_hr,OperObj,opt)
    end
    %%
    methods (Static,Hidden,Access= protected)
        Hnum_list_wire_iz =  Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab,Type)
        Poly_priciplayer_mat = Poly_priciplayer_mat_gen(principle_layer)
        cutlist= unique_label2cutlist(unique_label,NRPTS)
    end
    methods
        H_hr= hrz_gen(H_hr,kpoints_f,fin_dir,mode)
        [H00_H11_cell_list_1,H00_H11_cell_list_2] = H00_H11_cell_list_gen(H_hr,fin_dir,principle_layer)
        varargout = Green_prepare(H_hrz,principle_layer,fin_dir)
        [klist_cart,klist_frac,gap_list,fig] = findnodes2(H_hr,kzf_list,Noccupy,tolerance)
        H_hr = Subsall(H_hr,mode)
        H_hr = input_orb_struct(H_hr,filename,mode,options)
    end
    methods
        H_hr = GenfromOrth(H_hr,seed_r,seed_i,Accuracy,options)
        hrdat = Gen_hr(H_hr,filename,mode)
        [H_hr_bk] = POSCAR_gen(H_hr,filename,options)
        H_hr = tbbox_in_gen(H_hr,options)
        H_hr = pythtb_gen(H_hr,filename,options)
    end
    methods
        HcktObj = HR2Hckt(H_hr,options,options_homecell,options_para)
        [H_hr_forHckt,maxOnsite] = HRforHckt(H_hr,options)
        H_htrig = HR2Htrig(H_hr,options)
        H_hk = HR2HK(H_hr,kpoints_frac,options)
    end
    methods
        varargout = EIGENCAR_gen(H_hr,options)
        EIGENCARout = EIGENCAR_gen_sparse(H_hr,fermi,norb_enforce,klist_s_tmp)
    end
    methods
        function obj = init_mex(obj)
            if isempty(obj.mex_handle) || ~obj.mex_initialized
                % 准备初始化数据
                init_data = struct(...
                    'HnumL', obj.HnumL, ...
                    'vectorList_R', double(obj.vectorL(:,1:obj.Dim)) * obj.Rm);

                % 调用 MEX 初始化
                obj.mex_handle = mex_hamiltonian('init', init_data);
            end
        end
        [W,D,dH_dk_R,vec] = fft_wrapper(H_hr, klist)
        [W, D, dH_dk_R,vec] = fft(H_hr,klist_cart, rotate_cart)
        [W,D,dH_dk_R,dH_dk_dk_R] = fft_2(H_hr, klist_cart)
        %     obj.init_mex();
        %     [W, D, dH_dk_R] = mex_hamiltonian_calc(obj.mex_handle, klist);
        % end
    end
    methods
        H_hr = rewrite(H_hr,options)
        H_hr = rewind(H_hr)
        H_hr = init(H_hr,options)
        H_hr = applyOper(H_hr,SymOper,options)
        H_hr = symmetrize(H_hr,SymOper,options)
        H_hr = dualize(H_hr)
        [H_hr,R_vector_dist_] = dualizeR(H_hr,Rf)
        [ml_cell,ij_list] = Tij2lm(H_hr,Rf)
        [H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper)
        H_hr = hermitize(H_hr)
        H_hr_bk = subsOper(H_hr,SymOper)
        [H_hr_R,H_hr] = applyRU(H_hr,SymOper )
        [H_hr_R,H_hr] = applyR(H_hr,R)
        [H_hr,Smat] = Smatgen(H_hr,R,Accuracy)
        Factorlist_parity = factorlist_parity(H_hr)
        H_hr = nn(H_hr,search_range,Accuracy,Rlength_cut,options)
    end
    methods
        [EIGENCAR_3D,klist1,klist2,WEIGHTCAR_3D,WAVECAR_3D] = EIGENCAR_gen_3D(H_hr,kmesh,k3d,options)
        [EIGENCAR_slab,klist_l,kpoints_l,kpoints_name] = slab(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
        [EIGENCAR_slab,WEIGHTCAR_slab,klist1,klist2] = slab_fermiarc(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
        [EIGENCAR,orb_list,WEIGHTCAR,klist_l,kpoints_l,kpoints_name] = EIGENCAR_gen_wire(H_hr,Nslab,fermi,norb_enforce,KPOINTS_wire,vacuum_mode,np)
        [EIGENCAR,orb_list,WAVECAR] = EIGENCAR_gen_disk(H_hr,Nslab,fermi,norb_enforce,kpoints,vacuum_mode,np)
        [DOSCAR_l,DOSCAR_b,DOSCAR_r,w_list,klist_l,kpoints_l,kpoints_name] = surf(H_hr,w_range,fin_dir,KPOINTS_surf,principle_layer,eta,fermi,mode)
        [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2] = fermiarc(H_hr,w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi,mode)
        [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2,E_list] = fermiarc3D(H_hr,w_range,fin_dir,kmesh,kfermiarc,options)
    end
    methods(Static)
        H_hr = from_POSCAR_Script(varargin)
    end
    methods
        A = list(H_hr,vectorL,options)
        H_hr = rmhop(H_hr,vectorL,options)
        [fig,ax] = show(H_hr,mode,options)
        Hout = printout(H_hr,print_list,mode)
        [VarInit,EQL2] =GetInit(H_hr,H_hr2,vectorL)
        H_atom_soc = H_atom_soc(H_hr)
        kloop = kloop_gen(H_hr,input_mat,mode)
        %stucture
        [Rnn,nn_store_smart,Atom_store_smart,Rnn_map] = nn_information(H_hr,silence)
    end
    methods
        H_hr = H_TBSK_gen(H_hr,options)
        H_hr = H_TBSK_gen_sparse(H_hr,options)
    end
    methods
        H_hr_n = Connect2D(H_hr_n,H_hr_unitcell,opt)
    end
    methods(Static,Hidden,Access= protected)
        varargout = TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,options)
        varargout = TBSK_Var_gen_sparse(L_1_L,L_2_L,m_1_L,m_2_L,Rlength_L,l_L,m_L,n_L,options)
        varargout = TBSK_Var_gen_single(L_1,L_2,m_1,m_2,nn_level,l,m,n,options)
        Coeff = TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n,options)
    end
    methods
        A = GetProperty(H_hr,name)
    end
    methods (Static)
        EQ_list = Test_TBSK_Var_gen(testmode)
    end
    methods(Hidden,Access= protected)
        H_hr = nn_sk_smart(H_hr,search_range,Accuracy,Rlength_cut)
        H_hr = nn_sk_sparse(H_hr,search_range,Accuracy,Rlength_cut)
        H_hr = H_TB_gen_SK(H_hr,options)
        H_hr = H_TB_gen_SK_sparse(H_hr,level_cut,para_filename,onsite_mode)
        Atom_smart_t = Atom_smart_t_gen(site1,site2)
        [nn_sparse_temp,Rnn_list] = nn_sparse_t_gen(site1,site2,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
        [Rnn_list,nn_smart_t] = nn_smart_t_gen(Atom_smart_t,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
        hop = hop_gen(hop_pre,nn_level)
        [TBSK_hop,Coff] = TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orbsym1,orbsym2)
        [TBSK_hop,Coff] =  TBSK_hop_gen_sparse(site1,site2,Rij_cart,Rlength,nn_level)
    end
end
