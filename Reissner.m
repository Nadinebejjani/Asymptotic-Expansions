function out = Reissner(Plaque,F,PHI,Matrices_Override)

% ---------------- PARAMETRES USER -------------------------------


    % Vecteur fréquence
        NF = length(F) ; 
        
    % Vecteur angle de propagation de l'onde
        NPHI = length(PHI) ;
        
    % COUPLAGE ACOUSTIQUE ?
        acou = 0 ;
        epsilon_Acou = 1e-3 ;
        rho_air = 1.225 ; % kg/m^3
        Temperature = 20 ; % degrés
        c_air = 20.05*sqrt(Temperature + 273.15) ;

    %------------------PLAQUE------------------------------
        NCouches = Plaque.NCouches ;

        e = Plaque.e ;
        EL = Plaque.EL ;
        nuLT = Plaque.nuLT ;
        rho = Plaque.rho ;
        Kappa = Plaque.Kappa ;
        THETA = Plaque.THETA ;
        ET = Plaque.ET ;
        GLT = Plaque.GLT ;
        GNL = Plaque.GNL ;
        GTN = Plaque.GTN ;
        nuTL = Plaque.nuTL ;
    %------------------------------------------------------ 
 

%-------------------- CALCULS DE BASE ------------------------------------

    % Vecteur fréquence
        omega = 2*pi*F ;

    % Calcul des Altitudes (z=0 sur la face inf)    
        z0 = sum(e)/2 ;
        zplus(1) = e(1) - z0 ;
        zmoins(1) = -z0 ;
        if (NCouches>1)
            for n = 2:NCouches
                zplus(n) = sum(e(1:n)) - z0 ;
                zmoins(n) = sum(e(1:n-1)) - z0 ;
            end
        end

% ------------------------------------------------------------------------        
        

% -------------------- PROBLEME AUX VALEURS PROPRES ----------------------

Matrices = {} ;
            
% Matrices Constantes
        
    O = zeros(1,1,NCouches) ;
    
    % Raideur       
        QLL = EL./(1-nuLT.*nuTL) ; QLL = reshape(QLL,[1 1 NCouches]) ;
        QTT = ET./(1-nuLT.*nuTL) ; QTT = reshape(QTT,[1 1 NCouches]) ;
        QLT = nuLT.*ET./(1-nuLT.*nuTL) ; QLT = reshape(QLT,[1 1 NCouches]) ;
        GLT = reshape(GLT,[1 1 NCouches]) ;
        RotMatQ = @(c,s)[c^4 s^4 2*s^2*c^2 4*s^2*c^2 ;...
                        s^4 c^4 2*s^2*c^2 4*s^2*c^2 ;...
                        s^2*c^2 s^2*c^2 c^4+s^4 -4*s^2*c^2 ;...
                        c^3*s -s^3*c s^3*c-c^3*s 2*(s^3*c-c^3*s) ;...
                        s^3*c -c^3*s c^3*s-s^3*c 2*(c^3*s-s^3*c) ;...
                        s^2*c^2 s^2*c^2 -2*s^2*c^2 c^4+s^4-2*s^2*c^2] ;
        SQ44 = 1./reshape(GTN,[1 1 NCouches]) ;
        SQ55 = 1./reshape(GNL,[1 1 NCouches]) ;
        RotMatSQ = @(c,s)[c^2 s^2 ; s^2 c^2 ; -c*s c*s] ;
        
    % Inertie
        M = sum(e.*rho) ;
        I = sum((zplus.^2-zmoins.^2)/2.*rho) ;
        J = sum((zplus.^3-zmoins.^3)/3.*rho) ;
        MM = [M 0 0 I 0 ;...
              0 M 0 0 I ;...
              0 0 M 0 0 ;...
              I 0 0 J 0 ;...
              0 I 0 0 J ] ;
    
    % Matrices Comportement
        QQXY = zeros(3,3,NCouches) ;
        SQ = zeros(2,2,NCouches) ;
        for n = 1:NCouches
            QQ = RotMatQ(cos(THETA(n)),sin(THETA(n)))*[QLL(n) ; QTT(n) ; QLT(n) ; GLT(n)] ;
            QQXY(:,:,n) = [QQ(1) QQ(3) QQ(4) ; QQ(3) QQ(2) QQ(5) ; QQ(4) QQ(5) QQ(6)] ;
            SSQQ = RotMatSQ(cos(THETA(n)),sin(THETA(n)))*[SQ44(n) ; SQ55(n)] ;
            SQ(:,:,n) = [SSQQ(1) SSQQ(3) ; SSQQ(3) SSQQ(2)] ;
        end
            A = sum(repmat(reshape(e,[1 1 NCouches]),[3,3]).*QQXY,3) ;
            B = sum(repmat(reshape((zplus.^2-zmoins.^2)/2,[1 1 NCouches]),[3,3]).*QQXY,3) ;
            D = sum(repmat(reshape((zplus.^3-zmoins.^3)/3,[1 1 NCouches]),[3,3]).*QQXY,3) ;
            % Formule à Arthur
            F = Kappa^2*inv(sum(repmat(reshape(e,[1,1,NCouches]),[2,2]).*SQ/sum(e)^2,3)) ;
            F = [F(2,2) F(1,2) ; F(1,2) F(1,1)] ; % Inversion de V1,V2,d1 et d2 par rapport au poly
            
        % OVERRIDE DES MATRICES (SI MATRICES PRECEDEMENT)
        if exist('Matrices_Override')
            Matrices_Override
            for field = ['A','B','D','F','M','I','J']
                override = '' ;
                if isfield(Matrices_Override,field)
                    eval([field,' = Matrices_Override.',field,';']) ;
                    override = strcat(override,' ',field) ;
                end
            end
            display(['Override Matrices :',override])
        end
            
            Matrices.A = A ;
            Matrices.B = B ;
            Matrices.D = D ;
            Matrices.F = F ;
            Matrices.M = M ;
            Matrices.I = I ;
            Matrices.J = J 
            
            
    % Matrices du problème indépendantes de la fréquence et de l'angle
        % K00 à freq=0
            K00 = zeros(5) ;
            K00(4:5,4:5) = F ;
        % K10
            K10 = zeros(5) ;
            K10(3,4)=1i*F(1,1) ;
            K10(3,5)=1i*F(1,2) ;
            K10(4,3)=-1i*F(1,1) ;
            K10(5,3)=-1i*F(1,2) ;
        % K20
            K20 = zeros(5) ;
            K20(3,4)=1i*F(1,2) ;
            K20(3,5)=1i*F(2,2) ;
            K20(4,3)=-1i*F(1,2) ;
            K20(5,3)=-1i*F(2,2) ;
        % K11
            K11 = [ A(1,1) A(1,3) 0 B(1,1) B(1,3) ;...
                    A(1,3) A(3,3) 0 B(1,3) B(3,3) ;...
                    0 0 F(1,1) 0 0 ;...
                    B(1,1) B(1,3) 0 D(1,1) D(1,3) ;...
                    B(1,3) B(3,3) 0 D(1,3) D(3,3) ] ;
        % K22
            K22 = [ A(3,3) A(2,3) 0 B(3,3) B(2,3) ;...
                    A(2,3) A(2,2) 0 B(2,3) B(2,2) ;...
                    0 0 F(2,2) 0 0 ;...
                    B(3,3) B(2,3) 0 D(3,3) D(2,3) ;...
                    B(2,3) B(2,2) 0 D(2,3) D(2,2) ] ;
        % K12
            K12 = [ 2*A(1,3) A(3,3)+A(1,2) 0 2*B(1,3) B(3,3)+B(1,2) ;...
                    A(3,3)+A(1,2) 2*A(2,3) 0 B(3,3)+B(1,2) 2*B(2,3) ;...
                    0 0 2*F(1,2) 0 0 ;...
                    2*B(1,3) B(3,3)+B(1,2) 0 2*D(1,3) D(3,3)+D(1,2) ;...
                    B(3,3)+B(1,2) 2*B(2,3) 0 D(3,3)+D(1,2) 2*D(2,3) ] ;

    O = zeros(5) ;
    Id = eye(5) ;
    
    wtbr = waitbar(0,'Reissner') ;
    tic ;
    
    U = zeros(NF,NPHI,10,5) ;
    k = zeros(NF,NPHI,10) ;
    eta = zeros(NF,NPHI,10) ;
    Vp = zeros(NF,NPHI,10) ;
    Ve = zeros(NF,NPHI,10,3) ;
              
for p = 1:NPHI
    
    % Matrices indépendantes de la fréquence
        K1 = cos(PHI(p))*K10 + sin(PHI(p))*K20 ;
        K2 = cos(PHI(p))^2*K11 + sin(PHI(p))^2*K22 + cos(PHI(p))*sin(PHI(p))*K12 ;
        invK2 = inv(K2) ;
     
        
  % Résolution du Prob. aux V.P et post-traitement
    for w = 1:NF
        % Problème aux V.P
            K0 = K00-omega(w)^2*MM ;
            M1 = [O Id ; -invK2*K0 -invK2*K1] ;
            M2 = [Id O ; O Id] ;
            [Utemp,ktemp] = eig(M1,M2,'vector') ; 
            U(w,p,:,:) = Utemp(1:5,:).' ;
            k(w,p,:) = ktemp ;
        % Norme, Amortissement , Vp et Ve
            for u = 1:10
                UU = squeeze(U(w,p,u,:)) ;
                % Norme
                Ec = omega(w)^2*UU'*MM*UU ;
                U(w,p,u,:) = U(w,p,u,:)/sqrt(real(Ec)) ; 
                UU = squeeze(U(w,p,u,:)) ;  
                kr = k(w,p,u) ;
                kx = cos(PHI(p))*kr ;
                ky = sin(PHI(p))*kr ;
                DEF = [-1i*kx*UU(1) ; ...
                            -1i*ky*UU(2) ; ...
                            -1i*kx*UU(2)-1i*ky*UU(1) ; ...
                            -1i*kx*UU(4) ; ...
                            -1i*ky*UU(5) ; ...
                            -1i*kx*UU(5) - 1i*ky*UU(4) ; ...
                            UU(4) - 1i*kx*UU(3) ; ... 
                            UU(5) - 1i*ky*UU(3) ] ;
                EFF = [A B zeros(3,2) ; B D zeros(3,2) ; zeros(2,6) F]*DEF ;
                Ed = DEF'*EFF ;
                if acou
                    kz = sqrt(omega(w).^2/c_air^2 - kr^2) ;
                    Ed = DEF'*EFF + 2i*omega(w)^2*rho_air/kz*UU(3)'*UU(3) ;
                end
                eta(w,p,u) = imag(Ed)/real(Ed) ; % amortissement    
                Vp(w,p,u) = omega(w)/real(k(w,p,u)) ; % Vp
                Ve(w,p,u,1:2) = -real([EFF(1) EFF(3) EFF(4) EFF(6) EFF(7) ; EFF(3) EFF(1) EFF(6) EFF(5) EFF(8) ]*conj(-1i*omega(w)*UU)) ; 
                
            end
            
        if(toc>.1)
            waitbar(((p-1)*NF+w)/(NF*NPHI),wtbr) ;
            tic ;
        end
    
    end


        
end
delete(wtbr) ;  
    
    
% Sorties
    out.name = '$Reissner \, Angle$' ;
    out.U = U ;
    out.k = k ;
    out.eta = eta ;
    out.Ve = Ve ;
    out.Vp = Vp ;
    out.F = omega/2/pi ;
    out.Matrices = Matrices ;
    
end