classdef Oper < group
    %%OPER Symmetry group operation for real and Hilbert space actions.
    % This class represents symmetry operations (continuous/discrete groups)
    % including spatial transformations and Hilbert space unitary/antiunitary operations.
    %
    % Properties:
    %   R           - Spatial rotation matrix (double/sym)
    %   t           - Translation vector (double)
    %   conjugate   - Antiunitary conjugation flag (logical)
    %   antisymmetry- Hamiltonian sign flip flag (logical)
    %   Rf          - Local rotation matrix (double/sym)
    %   tf          - Local translation vector (double)
    %
    % Inherits from: group
    %% Note:
    % This class is partially inspired by the python code 'Qsymm'
    % <https://github.com/quantum-tinkerer/qsymm>
    % Dániel Varjas, Tómas Ö Rosdahl, and Anton R Akhmerov
    % Qsymm: algorithmic symmetry finding and symmetric Hamiltonian generation
    % New J. Phys. 20 093026 (2018)
    properties
        e
        R = []               % Spatial rotation matrix
        t = [0 0 0]          % Translation vector
        conjugate = false    % Complex conjugation flag
        antisymmetry = false % Antisymmetry operation flag
        Rf = []              % Local frame rotation matrix
        tf = [0 0 0]         % Local translation
        TBkitObj             % ConnectTBkitObj
        Rm
    end

    properties(Dependent=true)
        unitary {mustBeMember(unitary, [-1, 1])}
    end

    properties (GetAccess = protected, Hidden = true)
        continuous = false   % Continuous group flag
        strict_eq = false    % Strict equality check flag
    end
    %% Constuction
    methods
        function SymOper = Oper(R, U, t, options, optionsbasis)
            %%OPER Constructor for symmetry operator
            % Constructs symmetry operation with spatial, Hilbert space, and symmetry properties
            %
            % Inputs:
            %   R         - Rotation matrix (default: identity)
            %   U         - Unitary matrix (default: nan)
            %   t         - Translation vector (default: zeros)
            %   options   - Key-value pairs for additional parameters:
            %       * Rlocal          - Local rotation matrix
            %       * conjugate       - Antiunitary flag
            %       * antisymmetry    - Antisymmetry flag
            %       * strict_eq       - Strict equality checking
            arguments
                % (,)size class {Functions} = default
                R  =[1 0 0;0 1 0;0 0 1];
                U  = nan;
                t  {mustBeVector} = zeros(1,length(R));
                options.R = [];
                options.U = [];
                options.t = [];
                options.Rlocal = [];
                options.conjugate logical = false;
                options.unitary {mustBeMember(options.unitary, [-1, 1])} = 1 
                options.antisymmetry logical = false;
                options.strict_eq logical = false;
                options.Rm = [];
                optionsbasis.TBkitObj  = [];   
            end
            %
            if isempty(options.U)
                SymOper_U = U;
            else
                SymOper_U = options.U;
            end
            SymOper = SymOper@group(SymOper_U);
            if isempty(options.R)
                SymOper.R = R;
            else
                SymOper.R = options.R;
            end
            SymOper.R = integer(SymOper.R);
            if isempty(options.t)
                SymOper.t = t;
            else
                SymOper.t = options.t;
            end
            SymOper.t = (SymOper.t);

            SymOper.Rf = options.Rlocal;

            if options.conjugate
                SymOper.conjugate = options.conjugate;
            elseif options.unitary == -1
                SymOper.conjugate = true;
            end

            SymOper.antisymmetry = options.antisymmetry;
            SymOper.strict_eq = options.strict_eq;
            if isnan(SymOper.U)
                if ~isempty(optionsbasis.TBkitObj)
                    BasisFunction = BasisFunc(optionsbasis.TBkitObj);
                    SymOper  = SymOper.attachRm(optionsbasis.TBkitObj.Rm);
                    SymOper.U = BasisFunction.rotation('Oper',SymOper,'Rm',optionsbasis.TBkitObj.Rm);
                end
            end
            if ~isempty(options.Rm)
                SymOper  = SymOper.attachRm(options.Rm);
            end
            %SymOper.BasisHandle = options.BasisHandle ;
        end
    end
    %% get
    methods
        function e = get.e(OperObj)
            e = OperObj.isIdentity();
        end

        function unitary_value = get.unitary(OperObj)
            if OperObj.conjugate
                unitary_value = -1;
            else
                unitary_value =  1;
            end
        end

        function OperObj = set.unitary(OperObj, input)
            arguments
                OperObj Oper
                input {mustBeMember(input, [-1, 1])}
            end
            if input == 1
                OperObj.conjugate = false;
            else
                OperObj.conjugate = true;
            end
        end
    end
    %% 
    methods(Static)
        SymOper = identity(dim, shape,propArgs)
        SymOper = time_reversal(realspace_dim, U, spin,propArgs)
        SymOper = particle_hole(realspace_dim, U,propArgs)
        SymOper = chiral(realspace_dim, U,propArgs)
        SymOper = inversion(realspace_dim, U,quantumL,propArgs)
        SymOper = rotation(angle, axis, inversion, U, spin,options,propArgs)
        SymOper = spaceRotation(angle, axis,t, inversion, U, spin,options,propArgs)
        SymOper = C3z(realspace_dim, inversion, U, spin)
        SymOper = C4z(realspace_dim, inversion, U, spin)
        SymOper = C6z(realspace_dim, inversion, U, spin)
        SymOper = mirror(axis, U, spin,options,propArgs)
        group = square(tr, ph, generators, spin,options,propArgs)
        group = cubic(tr, ph, generators, spin,options,propArgs)
        group = hexagonal_2D(tr, ph, generators, spin,options, propArgs)
        group = hexagonal(tr, ph, generators, spin,options,propArgs)
    end
    methods
        [TrueOrNot,result] = commute(SymOper1,SymOper2)
    end
    methods
        SymOper = Ugen(SymOper,Basis,options)
    end
    methods
        [SYMCAR,OperObj] = Character_gen(OperObj,TBkitobj,klist,options)
        EIGENCAR_SYM = EIGEN(OperObj,WAVECAR,klist,orbL,options)
    end
    methods(Static)
        EIGENCAR_Kone = EIGEN_Kone(WAVECAR_one,D,V)
        EIGENCAR_Kone = EIGEN_Kone2(WAVECAR_one,U)
        EIGEN = EIGEN_one(WAVEFUNC,U)
    end
    methods(Static)
        U = Umat(SymOper1,SymOper2,options)
        P = similarity(SymOper1,SymOper2,options)
    end
    methods
        str = disp(SymOper,options)
        [SymMat,SymR] = sym(SymOper)
        SymOper = E(SymOper)
    end
    methods
        [basic_eq,U_eq]=eq(SymOper1,SymOper2)
        result = lt(SymOper1,SymOper2)
        SymOper_out = times(SymOper1,SymOper2)
        SymOper = inv(SymOper)
    end
    methods
        OperObj = attachRm(OperObj,Rm)
    end
    methods
        SymOper_Str = string(SymOper)
        SymOper_latex= latex(SymOper)
        SymOper_pretty= repr_pretty(SymOper,Cycle)
        SymOper_latex= repr_latex(SymOper)
        name= pretty(SymOper,options)
    end
    methods
        out = symmetry_from_permutation()
    end
    methods (Access= protected)
    end
    methods (Static)
        S_mat = spin_matrices(s, include_0)
        J_mat = spin_matrices_from_orb(quantumL,include_0,strict)
        L_mat = L_matrices(d, l)
        U = spin_rotation(n, s, roundint)
        phi = braket(phi1,Oper,phi2)
    end
    methods(Static)
        angle = name_angle(theta, Latex)
        str = mat2latex(mat,accuracy)
        str = num2latex(num,accuracy)
        C = isclose(A,B)
        C = allclose(A,B)
        n = round_axis(n)
        M = tensordot_naive(A,B,sizeLast)
        RotationMat = nThetad2Rotation(thetad,n)
        RotationMat = nTheta2Rotation(theta,n)
        [n,theta]= Rotation2nTheta(R,Rm)
        [n,thetad]= Rotation2nThetad(R,Rm)
        Rf = Rc2Rf(Rc,Rm)
        Rft = Rct2Rft(Rct,Rm)
        Rc = Rf2Rc(Rf,Rm)
        Rct = Rft2Rct(Rft,Rm)
        [prop,Coeff] = prop_to_id(A)
        mustBeOfClass(input,className)
        mustBeEqualSize(a,b)
        mustBeSize(a,b)
        mustBeDims(input,numDims)
        mustBeHalfInteger(A)
        [Asort,Usort] = sorteig(U,A)
    end
    methods(Static)
        eul = Rotation2eul( R, Rm )
        R = axang2rotm( axang )
        eul = axang2eul(axang)
        quat = rotm2quat( R )
    end
end
