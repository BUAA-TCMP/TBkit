function H_hr = applyOper(H_hr, SymOper, options)
%APPLYOPER Apply symmetry operations to Hamiltonian object
%   This function enforces symmetry constraints on tight-binding Hamiltonian
%   through specified symmetry operations while maintaining Hermiticity
%
%   Inputs:
%       H_hr      - Initial HR object representing Hamiltonian
%       SymOper   - Array of symmetry operation objects (Oper class)
%       options   - Configuration parameters (name-value pairs):
%         * generator : Generate group elements from generators (default=false)
%         * HcoeL     : Handle symbolic coefficients (default=false)
%         * Accuracy  : Tolerance for numerical simplification (default=1e-6)
%         * fast      : Enable fast mode for large systems (default=false)
%         * center    : Rotation center [x,y,z] (default=[0,0,0])
%         * Ugen      : Generate symmetry-adapted basis (default=false)
%
%   Output:
%       H_hr      - Symmetrized Hamiltonian with reduced degrees of freedom
%
%   Features:
%       - Supports both numerical and symbolic Hamiltonian manipulation
%       - Automatic Hermiticity enforcement
%       - Progress tracking for large symmetry groups
%       - Recursive symmetry application for generator mode

    % Input validation and default parameter setup
    arguments
        H_hr HR;
        SymOper Oper = Oper();
        options.generator = false;
        options.HcoeL = false;
        options.Accuracy = 1e-6;
        options.fast = false;
        options.center = [0,0,0];
        options.Ugen = false;
    end
    
    options2 = options;
    nSymOper = length(SymOper);
    
    % Generate symmetry-adapted basis if requested
    if options.Ugen
        try
            BasisFunction = BasisFunc(H_hr);
            SymOper = SymOper.Ugen(BasisFunction, 'Rm', H_hr.Rm, 'center', options.center);
        catch
            % Gracefully handle basis generation failures
        end
    end
    
    % Fast processing mode for large systems
    if options.fast
        % Initialize Hamiltonian data structures if empty
        if isempty(H_hr.AvectorL) && isempty(H_hr.BvectorL)
            H_hr = H_hr.init('fast', true);
            H_hr = H_hr.hermitize();
        end
        
        % Generator mode: recursively apply group elements
        if options.generator
            options2.generator = false;
            optionsCell = namedargs2cell(options2);
            for i = 1:nSymOper
                fprintf('******** apply (%d/%d)symmetry ********\n', i, nSymOper);
                disp(SymOper(i));
                SymOper_tmp = SymOper(i).generate_group();
                H_hr = applyOper(H_hr, SymOper_tmp, optionsCell{:});
                fprintf('----------   SymVarNum: %d   ----------\n', rank(H_hr.CvectorL));
            end
            H_hr = H_hr.hermitize;
            H_hr = H_hr.simplify(options.Accuracy);
        else
            % Direct symmetry application with progress tracking
            H_hr_R = H_hr;
            nSymOper = length(SymOper);
            pb = CmdLineProgressBar('Applying Symmetry ...');
            for j = 1:nSymOper
                [H_hr_R(j), H_hr] = applyRU(H_hr, SymOper(j));
                pb.print(j, nSymOper);
            end
            pb.delete();
            H_hr = sum(H_hr_R);
            H_hr = H_hr.hermitize;
            H_hr = H_hr.simplify(options.Accuracy);
        end
        return;
    end
    
    % Standard initialization sequence
    if ~H_hr.coe && ~H_hr.num
        H_hr = H_hr.init();
        H_hr = H_hr.hermitize();
    end
    
    % Format conversion if needed
    if ~strcmp(H_hr.Type, 'list')
        H_hr = H_hr.rewrite();
    end
    
    % Single symmetry operation handler
    if length(SymOper) == 1
        % Skip identity operations
        if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R, eye(3))
            return;
        end
        
        nSymOper = length(SymOper);
        fprintf('******** apply (%d/%d)symmetry ********\n', 1, nSymOper);
        disp(SymOper);
        
        % Symbolic coefficient handling
        if options.generator
            SymOper_tmp = SymOper.generate_group();
            nSymOper_tmp = length(SymOper_tmp);
            pb = CmdLineProgressBar('Applying Symmetry ...');
            H_hr_R = H_hr;
            for j = 1:nSymOper_tmp
                pb.print(j, nSymOper_tmp);
                [H_hr_R(j), H_hr] = applyRU(H_hr, SymOper_tmp(j));
            end
            pb.delete();
            H_hr = sum(H_hr_R)/nSymOper_tmp;
            H_hr = H_hr.simplify();
        elseif options.HcoeL
            % Symbolic coefficient manipulation
            [H_hr_R, H_hr] = applyRU(H_hr, SymOper);
            if ~isequal(H_hr_R.HcoeL, H_hr.HcoeL)
                Equationlist_r = (real(H_hr.HcoeL - H_hr_R.HcoeL) == 0);
                Equationlist_i = (imag(H_hr.HcoeL - H_hr_R.HcoeL) == 0);
                Equationlist_r = HR.isolateAll(Equationlist_r);
                Equationlist_i = HR.isolateAll(Equationlist_i);
                HcoeLtmp = H_hr.HcoeL;
                HcoeLtmp_r = subs(real(HcoeLtmp), lhs(Equationlist_r), rhs(Equationlist_r));
                HcoeLtmp_i = subs(imag(HcoeLtmp), lhs(Equationlist_i), rhs(Equationlist_i));
                H_hr.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
            end
            H_hr = H_hr.simplify();
        end
    else
        % Multiple symmetry operations processing
        nSymOper = length(SymOper);
        for i = 1:nSymOper
            fprintf('******** apply (%d/%d)symmetry ********\n', i, nSymOper);
            disp(SymOper(i));
            
            if options.generator
                SymOper_tmp = SymOper(i).generate_group();
                nSymOper_tmp = length(SymOper_tmp);
                pb = CmdLineProgressBar('Applying Symmetry ...');
                H_hr_R = H_hr;
                for j = 1:nSymOper_tmp
                    pb.print(j, nSymOper_tmp);
                    [H_hr_R(j), H_hr] = applyRU(H_hr, SymOper_tmp(j));
                end
                pb.delete();
                H_hr = sum(H_hr_R)/nSymOper_tmp;
                H_hr = H_hr.simplify(options.Accuracy);
            else
                % Recursive symmetry application
                H_hr = H_hr.applyOper(SymOper(i), 'generator', 'false');
            end
            fprintf('----------   SymVarNum: %d   ----------\n', length(H_hr.symvar_list));
        end
    end
end