function symbolic_term = soc_term_gen(l1, l2, m1, m2, spin1, spin2, element_seq)
    % SOC_TERM_GEN Generates the spin-orbit coupling term symbolically.
    %   symbolic_term = SOC_TERM_GEN(l1, l2, m1, m2, spin1, spin2, element_seq)
    %   computes the symbolic expression for the spin-orbit coupling term
    %   based on the provided orbital and spin quantum numbers and element sequence.
    %
    %   Inputs:
    %       l1, l2    - Orbital angular momentum quantum numbers.
    %       m1, m2    - Magnetic quantum numbers.
    %       spin1, spin2 - Spin quantum numbers.
    %       element_seq - Sequence of elements involved.
    %
    %   Outputs:
    %       symbolic_term - Symbolic expression for the spin-orbit coupling term.

    % Initialize the term string
    term_str = 'Esoc_';

    % Return zero if orbital angular momenta are different
    if l1 ~= l2
        symbolic_term = sym(0);
        return;
    end

    % Append orbital type to the term string
    switch l1
        case 0
            term_str = strcat(term_str, 'S_');
        case 1
            term_str = strcat(term_str, 'P_');
        case 2
            term_str = strcat(term_str, 'D_');
        case 3
            term_str = strcat(term_str, 'F_');
        otherwise
            error('Unsupported orbital angular momentum quantum number.');
    end

    % Append element sequence to the term string
    term_str = strcat(term_str, string(element_seq));

    % Generate orbital factors
    [alphaB1, alphaB2] = alpha_yml_gen(l1, m1);
    mb1 = m1;
    mb2 = -m1;

    % Generate orbital factors for L+ operator
    [alphaC1, alphaC2, mc1, mc2] = L_plus(l2, m2);
    orb_factor_L_plus = bra_phiAC_cet(mb1, mb2, mc1, mc2, alphaB1, alphaB2, alphaC1, alphaC2);

    % Generate orbital factors for L- operator
    [alphaC1, alphaC2, mc1, mc2] = L_minus(l2, m2);
    orb_factor_L_minus = bra_phiAC_cet(mb1, mb2, mc1, mc2, alphaB1, alphaB2, alphaC1, alphaC2);

    % Generate orbital factors for Lz operator
    [alphaC1, alphaC2, mc1, mc2] = L_z(l2, m2);
    orb_factor_L_z = bra_phiAC_cet(mb1, mb2, mc1, mc2, alphaB1, alphaB2, alphaC1, alphaC2);

    % Combine orbital factors
    orb_factors = [orb_factor_L_minus, orb_factor_L_z, orb_factor_L_plus];

    % Generate spin factors
    spin_factor_L_plus = delta(spin1, -spin2) * delta(spin1, 1);
    spin_factor_L_z = delta(spin1, spin2) * spin1;
    spin_factor_L_minus = delta(spin1, -spin2) * delta(spin1, -1);

    % Combine spin factors
    spin_factors = [spin_factor_L_minus, spin_factor_L_z, spin_factor_L_plus];

    % Compute the symbolic term
    symbolic_term = dot(orb_factors, spin_factors) * str2sym(term_str);
end


function factor = bra_phiAC_cet(mb1,mb2,mc1,mc2,alphaB1,alphaB2,alphaC1,alphaC2)
    %mode1 = delta(mb1,mc1) && delta(mb2,mc2);
    %mode2 = delta(mb1,mc2) && delta(mb2,mc1);
    factor1 = delta(mb1,mc1)*alphaB1*alphaC1;
    factor2 = delta(mb1,mc2)*alphaB1*alphaC2;
    factor3 = delta(mb2,mc1)*alphaB2*alphaC1;
    factor4 = delta(mb2,mc2)*alphaB2*alphaC2;
    factor = factor1+factor2+factor3+factor4;
    
    % bug?
%     if mode1 || mode2
%         if mode1
%             check1 = (alphaB1*alphaC2-alphaC1*alphaB2) == 0;
%             if check1
%                 factor = alphaB1/alphaC1;
%             else
%                 factor = 0;
%             end
%         elseif mode2
%             check2 = (alphaB1*alphaC1-alphaC2*alphaB2) == 0;
%             if check2
%                 factor = alphaB1/alphaC2;
%             else
%                 factor = 0;
%             end
%         else
%             factor = 0;
%         end
%     else
%         factor = 0;
%     end
end

function [alpha1,alpha2] = alpha_yml_gen(l,m)
    switch l
        case 0
            switch m
                case 0
                    alpha1 = sym(1/2);
                    alpha2 = sym(1/2);
                otherwise
                    warning('check your input POSCAR');
            end
        case 1
            switch m
                case 0
                    alpha1 = sym(1/2);
                    alpha2 = sym(1/2);
                case -1
                    alpha1 = sym(1i*2^(-1/2));
                    alpha2 = sym(1i*2^(-1/2));
                case 1
                    alpha1 = sym(-1*2^(-1/2));
                    alpha2 = sym(1*2^(-1/2));
                otherwise
                    warning('check your input POSCAR');
            end
        case 2
            switch m
                case 0
                    alpha1 = sym(1/2);
                    alpha2 = sym(1/2);
                case -1
                    alpha1 = sym(1i*2^(-1/2));
                    alpha2 = sym(1i*2^(-1/2));
                case 1
                    alpha1 = sym(-1*2^(-1/2));
                    alpha2 = sym(1*2^(-1/2)) ;
                case -2
                    alpha1 = sym(1i*2^(-1/2));
                    alpha2 = sym(-1i*2^(-1/2));
                case 2
                    alpha1 = sym(1*2^(-1/2));
                    alpha2 = sym(1*2^(-1/2)) ;
                otherwise
                    warning('check your input POSCAR');
            end
        case 3
            switch m
                case 0
                    alpha1 = sym(1);
                    alpha2 = sym(0);
                case -1
                    alpha1 = sym(1i*2^(-1/2));
                    alpha2 = sym(1i*2^(-1/2));
                case 1
                    alpha1 = sym(-1*2^(-1/2));
                    alpha2 = sym(1*2^(-1/2)) ;
                case -2
                    alpha1 = sym(1i*2^(-1/2)) ;
                    alpha2 = sym(-1i*2^(-1/2));
                case 2
                    alpha1 = sym(1*2^(-1/2)) ;
                    alpha2 = sym(1*2^(-1/2)) ;
                case -3
                    alpha1 = sym(1i*2^(-1/2));
                    alpha2 = sym(1i*2^(-1/2));
                case 3
                    alpha1 = sym(-1*2^(-1/2));
                    alpha2 = sym(1*2^(-1/2)) ;
                otherwise
                    warning('check your input POSCAR');
            end
        case -1 %wait for 
            switch m
                case 1
                    orb_sym = str2sym('1+x');
                case 2
                    orb_sym = str2sym('1-x');
                otherwise
                    warning('check your input POSCAR');
            end 
        case -2 %wait for 
            switch m
                case 1
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                case 2
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                case 3
                    orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');                 
                otherwise
                    warning('check your input POSCAR');
            end
        case -3 %wait for 
            switch m
                case 1
                    orb_sym = str2sym('1+x+y+z');
                case 2
                    orb_sym = str2sym('1+x-y-z');
                case 3
                    orb_sym = str2sym('1-x+y-z'); 
                case 4
                    orb_sym = str2sym('1-x-y+z');
                otherwise
                    warning('check your input POSCAR');
            end
        case -4 %wait for 
            switch m
                case 1
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                case 2
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                case 3
                    orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');
                case 4
                    orb_sym = str2sym('2^(-1/2)*z+2^(-1/2)*z^2');
                case 5
                    orb_sym = str2sym('-2^(-1/2)*z+2^(-1/2)*z^2');
                otherwise
                    warning('check your input POSCAR');
            end
        case -5 %wait for 
            switch m
%                 case 1
%                     orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
%                 case 2
%                     orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
%                 case 3
%                     orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
%                 case 4
%                     orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
%                 case 5
%                     orb_sym = str2sym('6^(-1/2)-2^(-1/2)*z+3^(-1)*(z^2)');
%                 case 6
%                     orb_sym = str2sym('6^(-1/2)+2^(-1/2)*z+3^(-1)*(z^2)');
%                 otherwise
%                     warning('check your input POSCAR');
            end
        otherwise
            warning('check your input POSCAR');
    end

end

function [alpha1_prime,alpha2_prime,ma1,ma2] = L_plus(l,m)
    ma1 = m+1;
    ma2 = -m+1;
    %
    [alpha1,alpha2] = alpha_yml_gen(l,m);
    alpha1_prime = alpha1*l_plus_minus_factor(l,m,'+');
    alpha2_prime = alpha2*l_plus_minus_factor(l,-m,'+');
end

function [alpha1_prime,alpha2_prime,ma1,ma2] = L_minus(l,m)
    ma1 = m-1;
    ma2 = -m-1;
    %
    [alpha1,alpha2] = alpha_yml_gen(l,m);
    alpha1_prime = alpha1*l_plus_minus_factor(l,m,'-');
    alpha2_prime = alpha2*l_plus_minus_factor(l,-m,'-');

end

function [alpha1_prime,alpha2_prime,ma1,ma2] = L_z(l,m)
    ma1 = m;
    ma2 = -m;
    %
    [alpha1,alpha2] = alpha_yml_gen(l,m);
    alpha1_prime = alpha1*ma1;
    alpha2_prime = alpha2*ma2;
end

function result = delta(m1,m2)
    if m1 == m2
        result = 1;
    else
        result = 0;        
    end
end

function factor = l_plus_minus_factor(l,m,mode)
switch mode
    case '+'
        factor = sym(((l-m)*(l+m+1))^(1/2));
    case '-'
        factor = sym(((l+m)*(l-m+1))^(1/2));
end


end
