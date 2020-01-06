%% REINITIALISATION DES VARIABLES ET FIGURES
clc
clear all
close all

%% exemples
 %% Plaque (-30,30,30,-30)
    
        Plaque = [] ;
        Plaque.e = [10 10 10 10]*1e-3 
        Plaque.NCouches = 4
        Plaque.THETA = [-30 30 30 -30]*pi/180 ;
        Plaque.NElmts = [10 10 10 10];
        Plaque.EL = [1.72e+11 1.72e+11 1.72e+11 1.72e+11 ]
        Plaque.ET = [6.89e+09 6.89e+09 6.89e+09 6.89e+09]
        Plaque.EN = [6.89e+09 6.89e+09 6.89e+09  6.89e+09]
        Plaque.GLT = [3.45e+09 3.45e+09 3.45e+09 3.45e+09]
        Plaque.GTN = [2.75e+09 2.75e+09 2.75e+09 2.75e+09 ]
        Plaque.GNL = [3.45e+09 3.45e+09 3.45e+09 3.45e+09 ]
        Plaque.nuLT = [0.25 0.25 0.25 0.25]
        Plaque.nuTN = [0.25 0.25 0.25 0.25 ]
        Plaque.nuLN = [0.25 0.25 0.25 0.25 ]
        Plaque.rho = [2260 2260 2260 2260]
        Plaque.nuTL = Plaque.nuLT.*Plaque.ET./Plaque.EL
        Plaque.nuNL = Plaque.nuLN.*Plaque.EN./Plaque.EL
        Plaque.nuNT = Plaque.nuTN.*Plaque.EN./Plaque.ET
        fmin =50 
        fmax =50000;
        %% Calculs EF
%% CALCUL DES NOMBRES D'ONDE
% Parametres
nF = 100 ;
scaleF = 'lin' ; 
phi0 =0;% % angle de propagation 
nPhi =1; %  Nombre de points sur la diagramme polaire
gammaMax = .1 ; % Critère sur l'amortissement pour le tri
        
    % Premiers calculs
        switch scaleF
            case 'lin'
                F = linspace(fmin,fmax,nF) ;
            case 'log'
                F = logspace(log10(fmin),log10(fmax),nF) ;
        end
        PHI = phi0+(0:nPhi-1)*2*pi/nPhi ;
        
    % CALCUL DES NOMBRES D'ONDE
        out = Shorter(Plaque,F,PHI) ; % Calcul nb ondes
        out = SortBranches3(out) ; % tri des branches
        out = TriPropagativ3(out,gammaMax) ; % tri des nb donde propagtifs
    
    % Résultats
        K = out.k ;  %tenseur dordre trois k(f,phi,mode) = K[f,phi,mode]
        U = out.U ;
    % Calculs
        h = sum(Plaque.e) ;
        starK = h*out.k;
        C = 2*pi*repmat(reshape(out.F,[nF 1 1]),[1 1 size(out.k,3)])./out.k ;
        C2 = (C).^2 ;
 %%
         iPhi = 1;
        ondes = 3; %ondes de compression
        c=real(squeeze(C(:,iPhi,ondes)));
        c2=real(squeeze(C2(:,iPhi,ondes)));
        KK=real(squeeze(K(:,iPhi,ondes)));
        strK = real(squeeze(starK(:,iPhi,ondes)));
%% Ajout valeurs Reissner
         iPhi = 1;
        ondes = 3;
        Plaque.Kappa = pi/sqrt(12) ;
        rei = Reissner(Plaque,F,PHI) ; % Calcul nb ondes
        rei = SortBranches3(rei) ; % tri des branches
        rei = TriPropagativ3(rei,gammaMax) ; % tri des nb donde propagtifs
        starK_Reissner = h*rei.k;
        C_Reissner = 2*pi*repmat(reshape(rei.F,[nF 1 1]),[1 1 size(rei.k,3)])./rei.k ;       
        creiss=real(squeeze(C_Reissner(:,iPhi,ondes)));
          % Ajout valeurs Kirchhoff
        C_Bend_LK =real(sqrt(2*pi*F).*(rei.Matrices.D(1,1)/rei.Matrices.M)^.25) ;
        C_Shear_LK = repmat(sqrt(rei.Matrices.A(3,3)/rei.Matrices.M),size(F)) ;
        C_Comp_LK = repmat(sqrt(rei.Matrices.A(1,1)/rei.Matrices.M),size(F)) ;
  %% Facteur de normalisation
res  = calculus(Plaque,PHI);
Nodes=res.Nodes; 
rhobar=res.rhobar %moyenne de rho
Gbar=res.Gbar %moyenne de G_LN
cs=sqrt(Gbar/rhobar) %facteur de normalisation
%% Compression
n=1000;
v=25;
res  = fct_compression(Plaque,PHI,n,v );
alpha0=res.alpha0
alpha0bar=res.alpha0bar
alpha_comp=res.alpha
alpha2=res.alpha2
alpha4=res.alpha4
%%  RESULTATS COMPRESSION

       XData=strK/(2*pi);%ratio
       YData =c/cs;
       kdata=strK;
        for k=2:v
        alphaprime(k-1)=alpha_comp(k);
        end 
%ordre 0
c2comp_ord0=transpose(alpha0*ones(1,100));
ccomp_ord0=sqrt(c2comp_ord0)/cs;
%ordre 1
c2comp_ord1=c2comp_ord0+alpha2*kdata.^2;
ccomp_ord1=sqrt(c2comp_ord1)/cs;
%ordre 2
c2comp_ord2=c2comp_ord1+alpha4*kdata.^4;
ccomp_ord2=sqrt(c2comp_ord2)/cs;

% Ordre p
c2comp_ordp=alpha0;
t=2;
for k=1:t
c2comp_ordp=c2comp_ordp+real(alphaprime(k))*kdata.^(2*k);
end
ccomp_ordp=sqrt(c2comp_ordp)/cs;%;/gltmoy;

pfem=plot(XData,YData,'k','LineWidth',2);
set(gca,'FontSize',25)
hold on
pord0=plot(XData,ccomp_ord0,'--','Color',[.5 .00 .5],'LineWidth',2);
hold on 
pord1=plot(XData,ccomp_ord1,'--r','LineWidth',2);
hold on
pordp=plot(XData,ccomp_ordp,'--','Color',[.22 .45 .2],'LineWidth',2);
hold off
%xlim([0 1/(2*pi)])
xlim([0 0.1])
ylabel('$\frac{\mbox{Compressional wave velocity}}{\mbox{Normalization factor}}=\frac{c}{c_S}$','Interpreter','latex','FontSize',20);
xlabel(' Thickness-wavelength ratio $(h/\lambda)$','Interpreter','latex','FontSize',20);%lgd={'FEM','p=0','p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','p=10'};
lgd={'FEM','p=0','p=1','p=2'};
legendflex([pfem,pord0,pord1,pordp],lgd,'ref', gcf ,'anchor', {'s','s'},'buffer',[+80 +80], 'ncol',4,'fontsize',25);
 
