function Out = SortBranches3(Data)
    
    F = Data.F ;
    NF = size(Data.k,1)  ;
    NPHI = size(Data.k,2) ;
    N = size(Data.k,3) ; % Number of branches

    Out = Data ;
    wtbr = waitbar(0,'Tri..') ;
    

    
% Sens de la plus petite à la plus grande fréquence
        % Valeur et dérivée de la valeur les + proches prises en compte
    tic
    for p = 1:NPHI
        for w = 2:NF    
            % Tri sur les donnees :
                %DataTri = Out.k ; % de nombre d'onde
                DataTri = repmat(2*pi*F(:),[1 size(Data.k,2) size(Data.k,3)])./Out.k ; % de vitesse de phase
            if (w==2) % deuxième point, le premier n'a pas besoin d'être trié
                Vitesse = zeros(N,1) ;
            else
                Vitesse = squeeze((DataTri(w-1,p,:)-DataTri(w-2,p,:))/(F(w-1)-F(w-2))) ;
            end
            % Calcul du point théorique suivant
                Ideal = squeeze(DataTri(w-1,p,:)) + Vitesse*(F(w)-F(w-1)) ;
            % Calcul de la matrice des différences
                DiffMat = abs(repmat(squeeze(DataTri(w,p,:)),[1 N])-repmat(Ideal.',[N 1])) ;
            % Calcul des indices de permutation
                ind = zeros(N,1) ;
                for i = 1:N
                    % Recherche du min
                        [~,indice] = min(DiffMat(:)) ;
                    % conversion en subindices 
%                         col = floor((indice-1)/N)+1 ;
%                         row = indice-(col-1)*N ;
                        [row,col] = ind2sub([N N],indice) ;
                    % Sauvegarde
                        ind(col) = row ;
                    % Mise des lignes et colonnes concernées en "oubli"
                        DiffMat(row,:) = Inf ;
                        DiffMat(:,col) = Inf ;
                    %ind(i) = find(abs(Out2.k(w,:,:)-Out2.k(w-1,:,i))==min(abs(Out2.k(w,:,:)-Out2.k(w-1,:,i)))) ;
                    %ind(i) = find(abs(Out.k(w,:,:)-Ideal)==min(abs(Out.k(w,:,:)-Ideal))) ;
                end
            % Permutation
                Out.k(w,p,:) = Out.k(w,p,ind) ;
                Out.U(w,p,:,:) = Out.U(w,p,ind,:) ;
                Out.eta(w,p,:) = Out.eta(w,p,ind) ;
                Out.Vp(w,p,:) = Out.Vp(w,p,ind) ;
                Out.Ve(w,p,:,:) = Out.Ve(w,p,ind,:) ;
            if(toc>0.2)
                wtbr = waitbar(((p-1)*NF + (w-1))/NF/NPHI,wtbr) ;
                tic ;
            end
        end
    end
    delete(wtbr)
end