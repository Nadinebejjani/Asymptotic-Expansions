% Fonction qui trie les ondes les plus propagatives

function Out = TriPropagativ3(Data,GammaMax)
    Out = Data ;
    for w = 1:size(Data.k,1)
        for p = 1:size(Data.k,2)
            % Copie locale des nb d'ondes
                K = Data.k(w,p,:) ;
            % Indices des ondes à conserver
                ind = logical(ones(size(K))) ; % indices des ondes valides
                ind = ind & imag(K)./real(K)<=1e-5 ; % Ondes avec amortissement qui a le bon signe (solutions physiques)
                ind = ind & real(K)>0 ;
                ind = ind & imag(K)./real(K)>=-abs(GammaMax) ; % Ondes avec amortissement inférieur à GammaMax
            % Application à la sortie
                Out.k(w,p,~ind) = -Inf ;
            % Tri par k décroissant
                [~,indSort] = sort(real(Out.k(w,p,:)),'descend') ;
            % Application
                Out.k(w,p,:) = Out.k(w,p,indSort) ;
                Out.U(w,p,:,:) = Out.U(w,p,indSort,:) ;
                if (isfield(Data,'eta'))
                    Out.eta(w,p,:) = Out.eta(w,p,indSort) ;
                end
%             %kk = Data.k(w,p,ind) ;
%             %kk = kk(find(real(kk)>0)) ;
%             [kk ind2] = sort(Data.k(w,p,ind),'descend') ;
%             ind = ind(ind2) ;
%             k(w,p,:) = Data.k(w,p,ind) ;
%             U(w,p,:,:) = Data.U(w,p,ind,:) ;
%                 if (isfield(Data,'eta'))
%                     eta(w,p,:) = Data.eta(w,p,ind) ;
%                 end
            
        end
    end
end
