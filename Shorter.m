function out = Shorter(Plaque,F,PHI)

%-------------------- CALCULS DE BASE ------------------------------------     

    %------DONNEES PLAQUE----------------------------------
        NCouches = Plaque.NCouches ;
        NElmts = Plaque.NElmts ; % Elements poutre par couche
        THETA = Plaque.THETA ;
        e = Plaque.e ; % Epaisseurs
        EL = Plaque.EL ;
        ET = Plaque.ET ;
        EN = Plaque.EN ;
        nuLN = Plaque.nuLN ;
        nuTN = Plaque.nuTN ;
        nuLT = Plaque.nuLT ;
        GLT = Plaque.GLT ;
        GNL = Plaque.GNL ;
        GTN = Plaque.GTN ;
%         etaL = Plaque.etaL ;
%         etaT = Plaque.etaT ;
%         etaN = Plaque.etaN ;
%         etaLT = Plaque.etaLT ;
%         etaLN = Plaque.etaLN ;
%         etaTN = Plaque.etaTN ;
%         etaGLT = Plaque.etaGLT ;
%         etaGNL = Plaque.etaGNL ;
%         etaGTN = Plaque.etaGTN ;
        rho = Plaque.rho ;
        
%         nuTL = nuLT.*ET./EL ;
%         nuNL = nuLN.*EN./EL ;
%         nuNT = nuTN.*EN./ET ;  
    %------------------------------------------------------ 
    
    
    % Vecteur fréquence
        NF = length(F) ; 
        omega = 2*pi*F ;
        
    % Vecteur angle de propagation de l'onde
        NPHI = length(PHI) ; 
    % Fonctions de forme de poutre
        N{1} = @(e)1-e ;
        N{2} = @(e)e ;

    % Coordonnées des noeuds (z=0 sur la face inf)    
        NElmtsTotal = sum(NElmts) ;    
        h = sum(e) ;
        Nodes = [] ;
        Nodes(1) = -h/2 ;
        for k = 1:NCouches
            for i = 1:NElmts(k)
                Nodes = [Nodes ; Nodes(end)+e(k)/NElmts(k)] ;
            end
        end
        L = Nodes(2:end)-Nodes(1:end-1) ;
        NNodes = length(Nodes) ;

% ------------------------------------------------------------------------        
        

% -------------------- PROBLEME AUX VALEURS PROPRES ----------------------

% Matrice des raideurs

    O = zeros(1,1,NCouches) ;
    
    % -------- ORTHOTROPE -----------------------     
      % Facteurs de perte
%         EL = EL.*(1+1i*etaL) ;
%         ET = ET.*(1+1i*etaT) ;
%         EN = EN.*(1+1i*etaN) ;
%         GLT = GLT.*(1+1i*etaGLT) ; 
%         GNL = GNL.*(1+1i*etaGNL) ;
%         GTN = GTN.*(1+1i*etaGTN) ;
        
% %         etaL = reshape(etaL,[1 1 NCouches]) ; 
% %         etaT = reshape(etaT,[1 1 NCouches]) ;
% %         etaN = reshape(etaN,[1 1 NCouches]) ;
% %         etaLT = reshape(etaLT,[1 1 NCouches]) ;
% %         etaLN = reshape(etaLN,[1 1 NCouches]) ;
% %         etaTN = reshape(etaTN,[1 1 NCouches]) ;
% %         etaGLT = reshape(etaGLT,[1 1 NCouches]) ;
% %         etaGNL = reshape(etaGNL,[1 1 NCouches]) ;
% %         etaGTN = reshape(etaGTN,[1 1 NCouches]) ;
        
        nuTL = nuLT.*real(ET)./real(EL) ;
        nuNL = nuLN.*real(EN)./real(EL) ;
        nuNT = nuTN.*real(EN)./real(ET) ;  
% %         nuTL = nuLT.*ET./EL ;
% %         nuNL = nuLN.*EN./EL ;
% %         nuNT = nuTN.*EN./ET ;  
        
        d = 1-nuTN.*nuNT-nuNL.*nuLN-nuLT.*nuTL-nuTL.*nuLN.*nuNT-nuLT.*nuTN.*nuNL ;
        d11 = EL.*(1-nuTN.*nuNT)./d ;
        d22 = ET.*(1-nuNL.*nuLN)./d ; 
        d33 = EN.*(1-nuLT.*nuTL)./d ; 
        d12 = real(ET).*(nuLT+nuNT.*nuLN)./d ; 
        d13 = real(EN).*(nuLN+nuTN.*nuLT)./d ; 
        d23 = real(EN).*(nuTN+nuLN.*nuTL)./d ; 
        d21 = d12 ;%EL.*(nuTL+nuNL.*nuTN)./d ; 
        d31 = d13 ;%EL.*(nuNL+nuTL.*nuNT)./d ; 
        d32 = d23 ;%ET.*(nuNT+nuLT.*nuNL)./d ;
            d11 = reshape(d11,[1 1 NCouches]) ;
            d22 = reshape(d22,[1 1 NCouches]) ;
            d33 = reshape(d33,[1 1 NCouches]) ;
            d12 = reshape(d12,[1 1 NCouches]) ;
            d13 = reshape(d13,[1 1 NCouches]) ;
            d23 = reshape(d23,[1 1 NCouches]) ;
            d21 = reshape(d21,[1 1 NCouches]) ;
            d31 = reshape(d31,[1 1 NCouches]) ;
            d32 = reshape(d32,[1 1 NCouches]) ;
            d44 = reshape(GTN,[1 1 NCouches]) ;
            d55 = reshape(GNL,[1 1 NCouches]) ;
            d66 = reshape(GLT,[1 1 NCouches]) ;
        DLTN = [d11 d12 d13 O O O ; ... Raideur
                d21 d22 d23 O O O ; ...
                d31 d32 d33 O O O ; ...
                O O O d44 O O ; ...
                O O O O d55 O ; ...
                O O O O O d66 ] ;
%         dDLTN = [etaL.*d11 etaLT.*d12 etaLN.*d13 O O O ; ... Amortissement
%                  etaLT.*d21 etaT.*d22 etaTN.*d23 O O O ; ...
%                  etaLN.*d31 etaTN.*d32 etaN.*d33 O O O ; ...
%                  O O O etaGTN.*d44 O O ; ...
%                  O O O O etaGNL.*d55 O ; ...
%                  O O O O O etaGLT.*d66 ] ; 
        Os = @(c,s)[c^2, s^2, 0, 0, 0, 2*c*s ;...
                    s^2, c^2, 0, 0, 0, -2*c*s ;...
                    0, 0, 1, 0, 0, 0 ;...
                    0, 0, 0, c, -s, 0 ;...
                    0, 0, 0, s, c, 0 ;...
                    -c*s, c*s, 0, 0, 0, c^2-s^2 ] ;
    % -------------------------------------------
     
    
% Matrices de Masse et Raideur
    Mi = cell(NElmtsTotal,1) ;
    Kr2i = cell(NElmtsTotal,NPHI) ;
    Kr1i = cell(NElmtsTotal,NPHI) ;
    Kr0i = cell(NElmtsTotal,NPHI) ;
    D = zeros(6,6,NCouches,NPHI) ;
    for p = 1:length(PHI)
        for k = 1:NCouches
            D(:,:,k,p) = Os(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p)))*DLTN(:,:,k)*Os(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p))).' ;
%             dD(:,:,k) = Os(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p)))*dDLTN(:,:,k)*Os(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p)))' ;
            for i = 1:NElmts(k)
                
                elmt = i+sum(NElmts(1:k-1)) ;
           
           % Masse
                Mi{elmt} = rho(k)*L(elmt)/6*[2 0 0 1 0 0 ; ...
                                             0 2 0 0 1 0 ; ...
                                             0 0 2 0 0 1 ; ...
                                             1 0 0 2 0 0 ; ...
                                             0 1 0 0 2 0 ; ...
                                             0 0 1 0 0 2 ] ;
           % Raideur                         
                Kr0i{elmt,p} = 1/L(elmt) * [D(5,5,k,p) D(4,5,k,p) 0 (-D(5,5,k,p)) (-D(4,5,k,p)) 0;...
                                            D(4,5,k,p) D(4,4,k,p) 0 (-D(4,5,k,p)) (-D(4,4,k,p)) 0;...
                                            0 0 D(3,3,k,p) 0 0 (-D(3,3,k,p));...
                                            (-D(5,5,k,p)) (-D(4,5,k,p)) 0 D(5,5,k,p) D(4,5,k,p) 0;...
                                            (-D(4,5,k,p)) (-D(4,4,k,p)) 0 D(4,5,k,p) D(4,4,k,p) 0;...
                                            0 0 (-D(3,3,k,p)) 0 0 D(3,3,k,p) ];

                Kr1i{elmt,p} = 1i/2 *  [0 0 (-D(1,3,k,p) + D(5,5,k,p)) 0 0 (D(1,3,k,p) + D(5,5,k,p));...
                                        0 0 (-D(3,6,k,p) + D(4,5,k,p)) 0 0 (D(3,6,k,p) + D(4,5,k,p));...
                                        (D(1,3,k,p) - D(5,5,k,p)) (D(3,6,k,p) - D(4,5,k,p)) 0 (D(1,3,k,p) + D(5,5,k,p)) (D(3,6,k,p) + D(4,5,k,p)) 0;...
                                        0 0 (-D(1,3,k,p) - D(5,5,k,p)) 0 0 (D(1,3,k,p) - D(5,5,k,p));...
                                        0 0 (-D(3,6,k,p) - D(4,5,k,p)) 0 0 (D(3,6,k,p) - D(4,5,k,p));...
                                        (-D(1,3,k,p) - D(5,5,k,p)) (-D(3,6,k,p) - D(4,5,k,p)) 0 (-D(1,3,k,p) + D(5,5,k,p)) (-D(3,6,k,p) + D(4,5,k,p)) 0 ];

                Kr2i{elmt,p} = L(elmt)/6 * [2 * D(1,1,k,p) 2 * D(1,6,k,p) 0 D(1,1,k,p) D(1,6,k,p) 0;...
                                            2 * D(1,6,k,p) 2 * D(6,6,k,p) 0 D(1,6,k,p) D(6,6,k,p) 0;...
                                            0 0 2 * D(5,5,k,p) 0 0 D(5,5,k,p);...
                                            D(1,1,k,p) D(1,6,k,p) 0 2 * D(1,1,k,p) 2 * D(1,6,k,p) 0;...
                                            D(1,6,k,p) D(6,6,k,p) 0 2 * D(1,6,k,p) 2 * D(6,6,k,p) 0;...
                                            0 0 D(5,5,k,p) 0 0 2 * D(5,5,k,p) ];
            end
        end
    end

    
M = zeros(3*NNodes) ;
    for elmt = 1:NElmtsTotal
        M(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6) = M(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6)+ Mi{elmt} ;
    end
    
Kr2 = cell(NPHI,1) ;
Kr1 = cell(NPHI,1) ;
Kr0 = cell(NPHI,1) ;
for p = 1:NPHI
    Kr2{p} = zeros(3*NNodes) ;
    Kr1{p} = zeros(3*NNodes) ;
    Kr0{p} = zeros(3*NNodes) ;
    for elmt = 1:NElmtsTotal
        Kr2{p}(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6) = Kr2{p}(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6)+ Kr2i{elmt,p} ;
        Kr1{p}(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6) = Kr1{p}(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6)+ Kr1i{elmt,p} ;
        Kr0{p}(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6) = Kr0{p}(3*(elmt-1)+1:3*(elmt-1)+6,3*(elmt-1)+1:3*(elmt-1)+6)+ Kr0i{elmt,p} ;
    end
end
    
% Construction du système (A(w)-k.B){q} = 0
    U = zeros(NF,NPHI,6*NNodes,3*NNodes) ;
    k = zeros(NF,NPHI,6*NNodes) ;
    wtbr = waitbar(0,'Calcul des Solutions') ;
    tic ;
    for p = 1:NPHI
        C2 = Kr2{p} ;
        %invC2 = inv(Kr2{p}) ;
        C1 = Kr1{p} ;
        O = zeros(size(C1)) ;
        I = eye(size(C1)) ;
        for i = 1:NF
            C0 = Kr0{p} - omega(i)^2*M ;
            A = [O I ; -C2\C0 -C2\C1] ;
            %A = [O I ; -invC2*C0 -invC2*C1] ;
            B = [I O ; O I] ;
            [Utemp,ktemp] = eig(A,B,'vector') ;
            k(i,p,:) = ktemp ;
            U(i,p,:,:) = Utemp(1:3*NNodes,:).' ;
            if(toc>0.3)
                waitbar(((p-1)*NF+i)/(NF*NPHI),wtbr) ;
                tic ;
            end
        end
    end
    delete(wtbr)
    
    
% Amortissement
        wtbr = waitbar(0,'Calcul de l''Amortissement') ;
        eta = zeros(NF,NPHI,6*NNodes) ;
        for p = 1:NPHI
            for i = 1:NF
                for u = 1:6*NNodes
                    K = Kr0{p} + real(k(i,p,u))*Kr1{p} + real(k(i,p,u))^2*Kr2{p} ;
                    UU = squeeze(U(i,p,u,1:3*NNodes)) ;
                    E = UU'*K*UU ;
                    eta(i,p,u) = imag(E)/real(E) ;
                    % Au passage, on norme par rapport à l'energie cinétique...
                        %norme = abs(real(omega(i)^2*UU'*M*UU)) ; 
                        %U(i,p,u,1:3*NNodes) = U(i,p,u,1:3*NNodes)/sqrt(norme) ;
                        U(i,p,u,1:3*NNodes) = U(i,p,u,1:3*NNodes)/sqrt(real(E)) ;
                end
                if(toc>0.3)
                    waitbar(((p-1)*NF+i)/(NF*NPHI),wtbr) ;
                    tic ;
                end
            end
        end
        delete(wtbr);
        
    % Deformations
        B  = @(z,k,L) 1/L* [1i*k*(z-L) 0 0 -1i*k*z 0 0 ;...
                            0 0 0 0 0 0 ;...
                            0 0 -1 0 0 1 ;...
                            0 -1 0 0 1 0 ;...
                            -1 0 1i*k*(z-L) 1 0 -1i*k*z ;...
                            0 1i*k*(z-L) 0 0 -1i*k*z 0 ] ;
                            
    % Contraintes
        H = @(k,z,L,R) 1/L * [ 1i*k*(z-L)*R(1,1) , 1i*k*(z-L)*R(1,6) , -R(1,3) ,-1i*k*z*R(1,1),-1i*k*z*R(1,6), R(1,3) ;...
                                1i*k*(z-L)*R(1,2) , 1i*k*(z-L)*R(2,6) , -R(2,3) ,-1i*k*z*R(1,2),-1i*k*z*R(2,6), R(2,3) ;...
                                1i*k*(z-L)*R(1,3) , 1i*k*(z-L)*R(3,6) , -R(3,3) ,-1i*k*z*R(1,3),-1i*k*z*R(3,6), R(3,3) ;...
                                -R(4,5) , -R(4,4) , 1i*k*(z-L)*R(4,5) , R(4,5) , R(4,4) ,-1i*k*z*R(4,5);...
                                -R(5,5) , -R(4,5) , 1i*k*(z-L)*R(5,5) , R(5,5) , R(4,5) ,-1i*k*z*R(5,5);...
                                1i*k*(z-L)*R(1,6) , 1i*k*(z-L)*R(6,6) , -R(3,6) ,-1i*k*z*R(1,6),-1i*k*z*R(6,6), R(3,6) ] ;
                            
        
    % SORTIES
        out.name = '$SFEM$' ;
        out.Plaque = Plaque ;
        out.U = U ;
        out.k = k ;
        out.eta = eta ;
        out.Nodes = Nodes ;
        out.L = L ;
        out.D = D ;
        out.DLTN = DLTN ;
        out.B = B ;
        out.H = H ;
        out.Kr0i = Kr0i ;
        out.Kr1i = Kr1i ;
        out.Kr2i = Kr2i ;
        out.F = F ;
        out.Ve = zeros(NF,NPHI,6*NNodes,3) ;
        out.Vp = zeros(NF,NPHI,6*NNodes) ;
    
    
    % Debug
%         imD_reD = imag(DLTN)./real(DLTN)
%         imd_red = imag(d)/real(d)
end