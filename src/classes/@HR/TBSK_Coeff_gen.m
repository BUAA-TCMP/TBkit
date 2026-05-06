function Coeff = TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n,options)
% TBSK_COEFF_GEN Generate Slater-Koster coefficients
%
%   COEFF = TBSK_COEFF_GEN(L_1,L_2,M_1,M_2,L,M,N,OPTIONS) calculates SK coefficients
%
%   Inputs:
%       L_1,L_2 - Angular momentum quantum numbers
%       m_1,m_2 - Magnetic quantum numbers
%       l,m,n - Direction cosines
%       options.sym_mode - Symbolic mode [default: false]
%   Output:
%       Coeff - Slater-Koster coefficients
%
%   Notes:
%       - Handles s, p, and d orbitals
%       - Supports both numeric and symbolic output
%       - Implements standard SK parameterization
arguments
    L_1 double{mustBeInteger};
    L_2 double{mustBeInteger};
    m_1 double{mustBeInteger};
    m_2 double{mustBeInteger};
    l  = 0;
    m  = 0;
    n  = 0;
    options.sym_mode = false;
end
if options.sym_mode
    Coeff = sym(zeros(1,3));
else
    Coeff = zeros(1,3);
end
switch L_1
    case 0
        switch L_2
            case 0
                Coeff(1) = 1;
            case 1
                switch m_2
                    case -1
                        Coeff(1) = m;
                    case 0
                        Coeff(1) = n;
                    case 1
                        Coeff(1) = l;
                end
            case 2
                SQRT3 = sqrt(3) ;
                switch m_2
                    case -2
                        Coeff(1) = SQRT3*l*m;
                    case -1
                        Coeff(1) = SQRT3*m*n;
                    case 0
                        Coeff(1) = (n^2 - (l^2 + m^2)/2);
                    case 1
                        Coeff(1) = SQRT3*n*l;
                    case 2
                        Coeff(1) = SQRT3/2 * (l^2- m^2);
                end
            case 3
        end
    case 1
        switch L_2
            case 1
                if m_1 == m_2
                    switch m_1
                        case -1
                            Coeff(1) = m^2 ; Coeff(2) = 1 - Coeff(1);
                        case 0
                            Coeff(1) = n^2 ; Coeff(2) = 1 - Coeff(1);
                        case 1
                            Coeff(1) = l^2 ; Coeff(2) = 1 - Coeff(1);
                    end
                else
                    switch m_1
                        case -1
                            switch m_2
                                case 0
                                    Coeff(1) = m*n ; Coeff(2) = - Coeff(1) ;
                                case 1
                                    Coeff(1) = m*l ; Coeff(2) = - Coeff(1);
                            end
                        case 0
                            switch m_2
                                case 1
                                    Coeff(1) = n*l ; Coeff(2) = - Coeff(1) ;
                                case -1
                                    Coeff(1) = n*m ; Coeff(2) = - Coeff(1) ;
                            end
                        case 1
                            switch m_2
                                case -1
                                    Coeff(1) = l*m ; Coeff(2) = - Coeff(1) ;
                                case 0
                                    Coeff(1) = l*n ; Coeff(2) = - Coeff(1) ;
                            end
                    end
                end
            case 2
                SQRT3 = sqrt(3) ;
                switch m_1
                    case -1
                        switch m_2
                            case -2
                                lM2 = m^2*l;
                                Coeff(1) = SQRT3*lM2 ; Coeff(2) = l-2*lM2;
                            case -1
                                nM2 = n*m^2 ;
                                Coeff(1) = SQRT3*nM2 ; Coeff(2) = n-2*nM2;
                            case 1
                                LMN = l*m*n ;
                                Coeff(1) = SQRT3*LMN ; Coeff(2) = -2*LMN;
                            case 0
                                mN2 = m*n^2 ;
                                Coeff(1) = mN2-m^3/2-m*l^2/2 ; Coeff(2) = -SQRT3*mN2;
                            case 2
                                M3 = m^3;mL2 = l*m^2;
                                Coeff(1) = SQRT3/2 * (mL2 - M3);Coeff(2) = -m + M3+mL2;
                        end
                    case 0
                        switch m_2
                            case -2
                                LMN = l*m*n ;
                                Coeff(1) = SQRT3*LMN ; Coeff(2) = -2*LMN;
                            case -1
                                mN2 = m*n^2 ;
                                Coeff(1) = SQRT3*mN2 ; Coeff(2) = m-mN2;
                            case 1
                                lN2 = l*n^2 ;
                                Coeff(1) = SQRT3*lN2 ; Coeff(2) = l-2*lN2;
                            case 0
                                nL2 = n*l^2;nM2 = n*m^2;
                                Coeff(1) = n^3 - nL2/2 - nM2/2 ; Coeff(2) = SQRT3*(nM2+nL2);
                            case 2
                                nL2 = n*l^2;nM2 = n*m^2;
                                Coeff(1) = SQRT3/2 * (nL2 - nM2);Coeff(2) = nM2 - nL2;
                        end
                    case 1
                        switch m_2
                            case -2
                                mL2 = l^2*m;
                                Coeff(1) = SQRT3*mL2 ; Coeff(2) = m-2*mL2;
                            case -1
                                LMN = l*m*n ;
                                Coeff(1) = SQRT3*LMN ; Coeff(2) = -2*LMN;
                            case 1
                                nL2 = n*l^2 ;
                                Coeff(1) = SQRT3*nL2 ; Coeff(2) = n-2*nL2;
                            case 0
                                lN2 = l*n^2 ;
                                Coeff(1) = lN2-l^3/2-l*m^2/2 ; Coeff(2) = -SQRT3*lN2;
                            case 2
                                L3 = l^3;lM2 = l*m^2;
                                Coeff(1) = SQRT3/2 * (L3 - lM2);Coeff(2) = l-L3+lM2;
                        end
                end
            case 3
        end
    case 2 % d
        switch L_2
            case 2 % d-d
                switch m_1
                    case {-2,-1,1}
                        switch m_2 % t2g
                            case {-2,-1,1} % t2g-t2g
                                L2 = l^2;M2 = m^2;N2 = n^2;
                                if m_1 == -2 && m_2 == -2 % dxy-dxy
                                    L2M2 = L2*M2;
                                    Coeff(1) = 3*L2M2 ;
                                    Coeff(2) = L2 + M2 - 4 * L2M2;
                                    Coeff(3) = N2 + L2M2 ;
                                elseif m_1 == -1 && m_2 == -1 % dyz-dyz
                                    M2N2 = M2*N2 ;
                                    Coeff(1) = 3*M2N2 ;
                                    Coeff(2) = M2 + N2 - 4 * M2N2;
                                    Coeff(3) = L2 + M2N2;
                                elseif m_1 == 1 && m_2 == 1 % dzx-dzx
                                    N2L2 = N2*L2;
                                    Coeff(1) = 3*N2L2 ;
                                    Coeff(2) = N2 + L2 - 4 * N2L2;
                                    Coeff(3) = M2 + N2L2;
                                elseif (m_1 == -2 && m_2 == -1) || (m_1 == -1 && m_2 == -2) % dxy-dyz
                                    LN = l*n;LNM2 = LN*M2;
                                    Coeff(1) = 3*LNM2 ;
                                    Coeff(2) = LN*(1-4*M2);
                                    Coeff(3) = LNM2 - LN;
                                elseif (m_1 == -2 && m_2 == 1) || (m_1 == 1 && m_2 == -2) % dxy-dzx
                                    MN = m*n;MNL2 = MN*L2;
                                    Coeff(1) = 3*MNL2 ;
                                    Coeff(2) = MN*(1-4*L2);
                                    Coeff(3) = MNL2 - MN;
                                elseif (m_1 == -1 && m_2 == 1) || (m_1 == 1 && m_2 == -1) % dyz-dzx
                                    LM = l*m;LMN2 = LM*N2;
                                    Coeff(1) = 3*LMN2 ;
                                    Coeff(2) = LM*(1-4*N2);
                                    Coeff(3) = LMN2 - LM;
                                else
                                    error('unexpected error!');
                                end
                            case {0,2}
                                if m_1 == -2 && m_2 ==0
                                    SQRT3 =sqrt(3); LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2pM2 = L2+M2;
                                    Coeff(1) = SQRT3/2 *LM*(N2-L2pM2);
                                    Coeff(2) = -2*LM*N2 ;
                                    Coeff(3) = LM*(1+N2)/2;
                                elseif m_1 == -1 && m_2 ==0
                                    SQRT3 =sqrt(3); MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2pM2 = L2+M2;
                                    Coeff(1) = SQRT3/2 *MN*(N2-L2pM2);
                                    Coeff(2) = MN*(L2pM2-N2) ;
                                    Coeff(3) = -MN*(1+L2pM2)/2;
                                elseif m_1 == 1 && m_2 ==0
                                    SQRT3 =sqrt(3); NL =n*l;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2pM2 = L2+M2;
                                    Coeff(1) = SQRT3/2 *NL*(N2-L2pM2);
                                    Coeff(2) = NL*(L2pM2-N2) ;
                                    Coeff(3) = NL*(1+L2pM2)/2;
                                elseif m_1 == -2 && m_2 ==2
                                    LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2mM2 = L2-M2;
                                    Coeff(1) = 1.5*LM*L2mM2;
                                    Coeff(2) = -2*LM*L2mM2 ;
                                    Coeff(3) = LM*L2mM2/2;
                                elseif m_1 == -1 && m_2 ==2
                                    MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2mM2 = L2-M2;
                                    Coeff(1) = 1.5*MN*L2mM2;
                                    Coeff(2) = -2*MN*L2mM2 - MN ;
                                    Coeff(3) = MN*L2mM2/2 + MN;
                                elseif m_1 == 1 && m_2 ==2
                                    NL =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2mM2 = L2-M2;
                                    Coeff(1) = 1.5*NL*L2mM2;
                                    Coeff(2) = -2*NL*L2mM2 + NL ;
                                    Coeff(3) = NL*L2mM2/2 - NL;
                                end
                        end
                    case {0,2}
                        switch m_2
                            case {0,2}
                                L2 = l^2;M2 = m^2;N2 = n^2;
                                if m_1 == 2 && m_2 == 2
                                    L2mM2_2 = (L2-M2)^2;
                                    Coeff(1) = 0.75*L2mM2_2;
                                    Coeff(2) = L2 + M2 - L2mM2_2;
                                    Coeff(3) = N2 + L2mM2_2/4;
                                elseif m_1 == 0 && m_2 == 0
                                    L2pM2 = (L2+M2);
                                    Coeff(1) = (N2 - L2pM2/2)^2;
                                    Coeff(2) = 3*N2*L2pM2;
                                    Coeff(3) = 0.75*L2pM2^2;
                                else
                                    SQRT3 = sqrt(3); L2mM2 = (L2-M2); L2pM2 = (L2+M2);
                                    Coeff(1) = SQRT3/2 * L2mM2*(N2-L2pM2/2 );
                                    Coeff(2) = -N2*L2mM2;
                                    Coeff(3) = (1+N2)*L2mM2/4;
                                end
                            case {-2,-1,1}
                                if m_1 == 0 && m_2 ==-2
                                    SQRT3 =sqrt(3); LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2pM2 = L2+M2;
                                    Coeff(1) = SQRT3/2 *LM*(N2-L2pM2);
                                    Coeff(2) = -2*LM*N2 ;
                                    Coeff(3) = LM*(1+N2)/2;
                                elseif m_1 == 0 && m_2 ==-1
                                    SQRT3 =sqrt(3); MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2pM2 = L2+M2;
                                    Coeff(1) = SQRT3/2 *MN*(N2-L2pM2);
                                    Coeff(2) = MN*(L2pM2-N2) ;
                                    Coeff(3) = -MN*(1+L2pM2)/2;
                                elseif m_1 == 0 && m_2 ==1
                                    SQRT3 =sqrt(3); NL =n*l;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2pM2 = L2+M2;
                                    Coeff(1) = SQRT3/2 *NL*(N2-L2pM2);
                                    Coeff(2) = NL*(L2pM2-N2) ;
                                    Coeff(3) = NL*(1+L2pM2)/2;
                                elseif m_1 == 2 && m_2 ==-2
                                    LM =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2mM2 = L2-M2;
                                    Coeff(1) = 1.5*LM*L2mM2;
                                    Coeff(2) = -2*LM*L2mM2 ;
                                    Coeff(3) = LM*L2mM2/2;
                                elseif m_1 == 2 && m_2 ==-1
                                    MN =m*n;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2mM2 = L2-M2;
                                    Coeff(1) = 1.5*MN*L2mM2;
                                    Coeff(2) = -2*MN*L2mM2 - MN ;
                                    Coeff(3) = MN*L2mM2/2 + MN;
                                elseif m_1 == 2 && m_2 ==1
                                    NL =l*m;L2 = l^2;M2 = m^2;N2 = n^2;
                                    L2mM2 = L2-M2;
                                    Coeff(1) = 1.5*NL*L2mM2;
                                    Coeff(2) = -2*NL*L2mM2 + NL ;
                                    Coeff(3) = NL*L2mM2/2 - NL;
                                end
                        end
                end
            case 3
        end
    case 3
        switch L_2
            case 3
        end
end
end
