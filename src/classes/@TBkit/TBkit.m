classdef TBkit < matlab.mixin.CustomDisplay
        %TBkit Tight-binding toolkit for electronic structure calculations
    %
    %   This class provides a comprehensive toolkit for tight-binding model
    %   calculations, including band structure computation, topological
    %   invariants, and various visualization tools.
    %
    %   Properties:
    %       Dim          - Dimension of the system (default=3)
    %       Basis_num    - Number of basis orbitals
    %       Rm           - Real-space lattice vectors
    %       Hermitian    - Flag for Hermitian Hamiltonian (default=true)
    %       VarsSeqLcart - Cartesian k-space symbolic variables
    %       VarsSeqLfrac - Fractional k-space symbolic variables
    %       Hsym         - Symbolic Hamiltonian
    %       Basis        - Basis orbital information
    %       orbL         - Orbital list
    %       elementL     - Element list for each orbital
    %       quantumL     - Quantum numbers for orbitals
    %       klist_cart   - Cartesian k-points list
    %       klist_frac   - Fractional k-points list
    %
    %   Dependent Properties:
    %       Gk           - Reciprocal lattice vectors
    %       symvar_list  - List of symbolic variables in Hamiltonian
    %       Nbands       - Number of bands
    %       Hfun         - Function handle for Hamiltonian
    %
    %   Methods:
    %       TBkit            - Constructor
    %       TBkit_init       - Initialize default properties
    %       TBkitCopy        - Copy properties from another TBkit object
    %       timtj_gen        - Generate tij-tj terms
    %       tjmti_gen        - Generate tj-ti terms
    %       input_Rm         - Input real-space lattice vectors
    %       bandplot         - Plot band structure
    %       Chern            - Calculate Chern number
    %       WilsonLoop       - Compute Wilson loop spectrum
    %       WindingNumber    - Calculate winding number
    %
    %   Static Methods:
    %       POSCAR_read      - Read VASP POSCAR file
    %       BerryConnection  - Calculate Berry connection
    %       EIGENSOLVE       - Solve eigenvalue problem
    properties
        Dim = 3;
        Basis_num;
        Rm = [];
        Hermitian = true;
    end
    properties
        VarsSeqLcart = [];
        VarsSeqLfrac = [];
        Hsym;
    end
    properties (GetAccess = protected, Hidden = true)
    end
    properties (GetAccess = protected)
    end
    properties (Transient, Hidden = true)
    end
    properties(Dependent = true)
        Gk;
        symvar_list;
        Nbands;
    end
    properties
        Hfun;
        Basis;
        orbL = [];
        orbL_cart = [];
        elementL = [];
        quantumL = [];
        orb_symL = sym([]);
        sgn = 1;
        Atom_name = [];
        Atom_num = [];
        sites;
        symmetry_operation = [];
        klist_cart;
        klist_l;
        klist_frac;
        kpoints_frac;
        kpoints_l;
        kpoints_name;
        Rnn;
        nn_store;
        timtj;
        tjmti;
        Sparse_vector = [];
        N_Sparse_vector = [];
        CutList = [];
    end
    methods (Access = protected)
        propgrp = getPropertyGroups(~)
    end
    methods
        function TBkitobj =TBkit(propArgs)
            arguments
                propArgs.?TBkit;
            end
            Fieldnames = fieldnames(propArgs);
            if isempty(Fieldnames)
            else
                for iFieldname = 1:numel(Fieldnames)
                    TBkitobj.(Fieldnames{iFieldname}) = propArgs.(Fieldnames{iFieldname});
                end
            end
            TBkitobj = TBkitobj.TBkit_init();
        end
        function TBkitobj = TBkit_init(TBkitobj)
            if isempty(TBkitobj.Rm)
                TBkitobj.Rm = eye(TBkitobj.Dim);
            end
            if isempty(TBkitobj.VarsSeqLcart)
                syms k_x k_y k_z k_w real;
                TBkitobj.VarsSeqLcart = [k_x k_y k_z k_w];
            end
            if isempty(TBkitobj.VarsSeqLfrac)
                syms k_1 k_2 k_3 k_4 real;
                TBkitobj.VarsSeqLfrac = [k_1 k_2 k_3 k_4];
            end
        end
        function TBkitobj = TBkitCopy(TBkitobj,TBkitobj_in)
            CopyItem = [...
                "Rm","orbL","Dim","elementL","quantumL","orb_symL","sgn","Atom_name","Atom_num",...
                "sites","symmetry_operation","klist_cart","klist_l","klist_frac","kpoints_l","kpoints_name",...
                "Rnn","nn_store"...
                ];
            for i = 1:numel(CopyItem)
                TBkitobj.(CopyItem(i)) = TBkitobj_in.(CopyItem(i));
            end
        end
    end
    %% get
    methods
        %function Hfun = get.Hfun(TBkitobj)
        %    Hfun = matlabFunction(TBkitobj.Hsym,'Vars',TBkitobj.VarsSeqLcart(1:TBkitobj.Dim));
        %end
        function Gk = get.Gk(TBkitobj)
            Gk = (eye(length(TBkitobj.Rm))*2*pi/(TBkitobj.Rm)).';
        end
        function symvar_list = get.symvar_list(TBkitobj)
            symvar_list = symvar(TBkitobj.HcoeL);
        end
        function Nbands = get.Nbands(TBkitobj)
            if isprop(TBkitobj,"WAN_NUM")
                Nbands = TBkitobj.WAN_NUM;
            else
                Nbands = TBkitobj.Basis_num;
            end
        end
    end
    methods
        TBkitobj = timtj_gen(TBkitobj,mode)
        TBkitobj = tjmti_gen(TBkitobj,mode)
        TBkitobj = SliceGen(TBkitobj)
    end
    methods
        TBkitobj = input_Rm(TBkitobj,Rm)
        TBkitobj = input_orb_struct(TBkitobj,filename,mode,options)
        TBkitobj = nn(TBkitobj,search_range,Accuracy,Rlength_cut,options)
        [klist_l,kpoints_l,kpoints_name] = kpath_information(TBkitobj)
        TBkitobj = kpathgen3D(TBkitobj,KPOINTS_name,nodes,kpoints_name_tmp)
        [Rm,sites,Atom_name,Atom_num] = POSCAR_gen(TBkitobj,filename,Rm,sites,Atom_name,Atom_num)
        [num_label,coe_label] = NumOrCoe(TBkitobj)
        varargout = klist_show(TBkitobj,options)
        varargout = PARCHG_gen(varargin)
    end
    methods(Static)
        [kloop1_frac,kloop2_frac,kloop1_cart,kloop2_cart,klist_l,kstart_frac,kstart_cart] = kloop2D(Rm,options)
        [klist_cart,klist_frac] = kloop1D(kpoint_frac,Orientation,radius,opt)
    end
    methods
        fid = kpath_card_gen(TBkitobj,mode,filename)
        fid = write_pj(TBkitobj,mode,filename)
    end
    methods(Static)
        [orbital_out,Rm_s_fin] = AddVacuumLayer(orbital_init,POSCAR_file,fin_dir_list,options)
        POSCAR = POSCAR_cell_read(filename,formatSpec)
        [Rm,sites,Atom_name,Atom_num,elements]=POSCAR_read(filename,mode,options)
        tmpsites = MakeUpSites(sites,spintype)
        [klist_cart,klist_frac,klist_l,kpoints_l,kpoints_frac] = kpathgen(kpoints,nodes,Gk,Gk_,options)
    end
    methods(Static)
        DOSCAR = DOSCAR_gen(GREENCAR,mode)
        GREENCAR = GREENCAR_gen(w_list,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode,mu_max)
        [Gwl,Gwb,Gwr] = Tmatrix_iter(H00,H01,w,eta,mu_max,infinity_small)
        [Gwl,Gwb,Gwr] = GW_iter(H00,H01,w,eta,mu_max,infinity_small)
        Tmatrix = Tmatrix_gen(H00,H01,w,eta)
        [Green_00,Green_00s1,Green_00s2]= Tmatrix2Green00(Tmatrix,H00,H01,w,eta,n)
    end
    methods(Static)
        [Hue,surf_level,hing_level] = orbone2hsv(orb_one,discrimination,center,orientation)
        HSVCAR = HSVCAR_gen(orb_list,mode,discrimination,center,orientation)
        OberserveValue = Observecar_gen(WAVECAR,Oper)
        [COLORCAR,WEIGHTCAR] = COLORCAR_gen(WAVECAR,HSVCAR,signlist)
        [COLOR_one,SIGN_one] = COLOR_one_gen(WF,HSVCAR,signlist)
    end
    methods(Static)
        str = relaceSpaceByComma(str)
        str = mat2str_python(mat)
        SIGN_one = SIGN_one_gen(WF,signlist)
        Symble = SymbolicVarible(seeds,superscript,subscript,level)
    end
    methods
        varargout = bandplot(TBkitobj,varargin)
    end
    methods(Static)
    end
    methods
        Chern = ChernNumber(TBkitobj, tensor_index, optionsK, options)
        [BCCAR, Grid] = BerryCurvature_2D(TBkitobj, optionsK, options)
        [BFCAR,BF_WAVECAR,klist_l,WAVELOOPCAR] = WilsonLoop(TBkitobj,optionsK,options)
        [BFCAR,WEIGHTCAR,klist_l] = ProjectedWilsonLoop(TBkitobj,options)
        [nested_BFCAR,nested_BF_ALL,klist_l] = nested_WilsonLoop(TBkitobj,optionsK,options,optionsNested)
    end
    methods
        [pmulist,klist_l] = nestedWilsonLoop(TBkitobj,optionsK,options,optionsNested)
        [WannierCenterCAR,WannierCenterWAVECAR,klist_l,WAVELOOPCAR] = WannierCenter(TBkitobj,optionsK,options)
    end
    methods
        [WindingNumber,WL] = WindingNumber(Ham_obj,GammaOper,kloop,options)
        [HL,pHL] = Topo1DpreHL(TBkitobj,klist,options)
        WAVECAR_loop = Topo1DpreWAVECAR(TBkitobj,klist,options)
        [BP,WAVECAR_loop] = BP_1D(TBkitobj,klist,options,optionsWAVE)
    end
    methods(Static)
        WAVECAR = cleanWAVECAR(WAVECAR,EIGENCAR,V,Accuracy)
        WAVEFuncL = smoothDegen2(WAVEFuncL)
        WAVEFuncL = smoothDegen(WAVEFuncL,V)
        DegenPair = checkDegen(EIGEN,Accuracy)
        [WAVECAR_loop_] = modify_WAVECAR(WAVECAR_loop,BF_WAVECAR)
        [WCC,WCCvec,HWan] = WannierCenter1D(WAVECAR_loop)
        [BF,BF_WAVE,Wan] = wancenter_1D(WAVECAR_loop,mode)
        F = BerryConnection(W1,W2)
        F = Berryphase_2D()
    end
    methods(Static)
        F = BerryPhaseLine_fun(Fun,kloopr,options)
        F = BerryPhaseLine_definition_sym_i(Hsym,options)
        [BFCAR,BF_WAVECAR,klist_l] = WilsonLoop_fun(Hfun,optionsK,options)
        BC = BC_definition(Hsym,para1,para2,epsilon,options)
        Bc = BC_kubo_sym(Hsym,para1,para2,epsilon,options)
        A_1 = BerryConnection_definition(Eigenvector_sym,para1)
    end
    methods(Static)
        HoutL = HCAR_gen(Hfun,klist_cart,Norb)
        [EIGENCAR,WAVECAR,HoutL] = EIGENSOLVE(Hfun,klist_cart,Norb,opt)
        EIGENCAR = arrangeEIGENCAR(EIGENCAR,REFCAR,method,opt)
    end
    methods(Static)
        pagenew = page_mtimes_matrix(page,mat)
        pagenew = matrix_mtimes_page(mat,page)
        pagenew = matrixtimespage(mat,page)
        [c,t] = coeffsAll(p,varargin)
        Equation_list = isolateAll(Equation_list,Symvar_list)
        [factor_list_1,factor_list_2] = factorAll(SymListRow)
        SymListRow = cleanVar(SymListRow,Accuracy)
        Pmat = Pvector(A)
        Pmat = Pplane(A,mode)
        Gknew = CartisianMat(Gk,dir_seq,kstart)
        M  = P2M(P)
        Eigenvetor=NomalizeEigenvector(Eigenvetor)
    end
    methods(Static)
        [H_sym_Gamma,H_sym_Gamma_L,H_latex_Gamma] = GammaDecomposition(H_sym)
        [CoeForPauli] = pauliDecompositionNumerial(H_double)
        [H_sym_pauli,H_sym_pauli_L,H_latex_pauli] = pauliDecomposition(H_sym)
    end
    methods(Static)
        Value = loss_func(parameters,extra_parm,options)
        Varbayes = VarBayes(Varlist,VarGuess,VarWidth,VarFix)
        ValueTotal = EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR,extra_parm,options)
        options_extra = FitOptionHelper(EIGENCAR_DFT,algorithm,options)
        IMG = EIGENCAR2IMG(EIGENCAR,dE,Erange,kselect)
        [DEGENCAR, NODEINSECT]= Degeneracy_EIGENCAR(EIGENCAR,highK,dE)
    end
    methods(Static,Hidden)
        [ChirdrenCell,Type] = fixedchildren(SymVar,mode)
        [Asort,Usort] = sorteig(U,A)
        SymVar = shelling(SymVar)
        R = Sph2Cart(S)
    end
    methods (Static,Hidden,Access= protected)
        C = isclose(A,B)
        C = allclose(A,B)
        TrueOrFalse = LatticeVectorTest(Vector,Accuracy)
        angle = name_angle(theta, Latex)
        symPage = subspage(symPage,lhs,rhs)
        [nn_sparse_temp,Rnn_list] = nn_sparse_gen(orb1,orb2,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut,onsite)
        labelcut_list = labelcut_list_gen(Atom_num)
        l_num = orb2l(input)
        m_num = orb_sym2m(input)
        orb_sym = Ymlsym(l,m,orb_symbk)
        orbsym_n = subs_xyz(orbsym,Rlmn)
        out = delta_orb(orb1,orb2)
        out = delta_orb_sym(orb_sym1,orb_sym2)
        R_struct = Rm2abc(Rm)
        flag = strcontain(str,containlist)
        rc_plus = plusrc(rc)
        To_red_sc=to_red_sc(red_vec_orig ,Ns)
        To_red_pc=to_red_pc(red_vec_sc ,Ns)
        [orb_one_incell,translation_vector]=translation_orb(orb_one)
        [list_obj_unique,sorted_label,cut_slice] = cut_tools(list_obj)
    end
end
