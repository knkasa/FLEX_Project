clc
clear
close all
delete  ./output/*.mat

%{
=========================================================================
 Fully self-consistent two-orbital FLEX code in superconducting state.
 This is the main file.  Run this to get the result. All other files are
 sub-routines.  Down the line, you may choose the band of your choice.
 See my dissertation for results.  
   
    Ken Nakatsukasa, Apr. 2018.   

 For bilayer model, set Nelm=2.  
 For general two orbital, set Nelm=3.  
=========================================================================
%}


Nk = [ 16   16 ];       % 2*nk=total_k
Nelm =   2    ;       % # total elements in matrix
gap_symm = 1   ;      % 1=swave,  2=dwave   (see below)
gap_input = [  0.01  0.01  0.01   ];   % initial guess of the gap 
fix_fill = 0 ;        %  1=fix_filling, check mu
filling0 = 1.91  ;    % choose if fix_fill=1,  half_filling=2.0

U = 0.8    ;    % Hubbard U  [eV],  

Tlist =  30 % linspace( 10, 60, 15 )  ;   % [K] this can be a vector
eps = 1e-10 ;     % convergence criteria.  default=1e-5

maxiter = 10000;   % max iteration
Ncut = 3 ;       % cutoff matsubara freq  (3=default)
Power2NumFreq =   0 ;  % fix # of matsubara freq (numwi) 

%--------------------------------------------------------------------------

% profile on

if Nelm==3
   %J1 = 0.0;  U2=0;  J2=0;   % equivalent to bilayer
   J1=1.0;   U2=0.2;    J2=J1;  
   %J1= U/2;  U2=U-2*J1;  J2=J1;   % spin-rotational invariance (by Kontani)
   
   Us = [ U 0 0 J1 ;  0 U2 J2 0 ; 0 J2 U2 0 ; J1 0 0 U ];
   Uc = [ U 0 0 2*U2-J1 ; 0 -U2+2*J1 J2 0 ; 0 J2 -U2+2*J1 0 ; 2*U2-J1 0 0 U ];
end


maxchi_cutoff = 0.999; %0.9999;
Norb = 2  ;     % # of eps(i,i), diagonal elements
AddHFSelfE = 1  ;    % 1=add Hartree, 0=no Hartree
kb = 8.6173303e-5;     % Boltzman constant
Tlist = kb*Tlist;
T = Tlist(1);
Ulist = U;

runInSerial = 1;    %switch between for-loop and parfor-loop
saveSelfE = 1;      %save self-energy
loadSelfE = 1;      %load self-energy at different T or U as the input
flySelfE = 0;       %load self-energy calculated on the fly
if (loadSelfE == 1) && (saveSelfE ~= 1)
    flySelfE = 1;
end

checkFreqSymm = 0;  %check frequency even/odd symmetry
% FFT can produce asymmetry due to small numerical errors

useSymmetry = 1;    % 1 = C_2(inversion) symmetry is used
NCorder = 2;        % correction order in fourier transform
% cpus = feature('numCores'); % for parallel computing
% pool = parpool(cpus)


nktot = 4*Nk(1)*Nk(2);
for ii = 1:2
    K{ii} = (0:1:2*Nk(ii)-1)*pi/Nk(ii);
end
[KX,KY] = meshgrid(K{2},K{1});      %row: y, column: x


%--------- Model band structure -------------------------------------------

    % To reproduce S+- bilayer, 
    % filling0=2.0, tz=-0.6, t=[0.2,-0.02,-0.05]   U=0.8, 
    
    %
   % **** Bilayer model *****  
     t = [  0.2  -0.02   0.05 ]  ;   %t=[eV]        
     ek(:,:,1) = -2*t(1)*(cos(KX)+cos(KY)) - 4*t(2)*cos(KX).*cos(KY) - t(3) ;
    % cuprate band from ARPES
    %t =  [  -0.598   0.0962  -0.1306  -0.0507  0.0939  0.1 ];  
    %ek(:,:,1) = 0.5*t(1)*(cos(KX)+cos(KY)) + t(2)*cos(KX).*cos(KY) ...
    %  + 0.5*t(3)*(cos(2*KX)+cos(2*KY)) + 0.5*t(4)*(cos(2*KX).*cos(KY) ...
    % + cos(2*KY).*cos(KX)) + t(5)*cos(2*KX).*cos(2*KY) + t(6)  ;
       tz = -0.6 ;     % default (negative)
       ek(:,:,2) = tz*ones(size(KX));  
     %}
    
    %{
    % **** General two-orbital *****
     t = [ -0.33   0.384   -0.234  -0.26  ];   mu = 0.0  ;  % paper 
    % t = [ -0.38   0.3 -0.3  -0.2  ];   mu=0.4 ;  % other example
     ek(:,:,1) =   -2*t(1)*cos(kx)-2*t(2)*cos(ky)  - 4*t(3)*cos(kx).*cos(ky) - mu ;
     ek(:,:,2) =   -4*t(4)*sin(kx).*sin(ky)  ; 
     ek(:,:,3) =   -2*t(2)*cos(kx)-2*t(1)*cos(ky)  - 4*t(3)*cos(kx).*cos(ky) - mu ;
     %}
    
    
%--------------------------------------------------------------------------

      
  for k1 = 1:2*Nk(1)
   for k2 = 1:2*Nk(2)    
       if Nelm==2
          mat_eps = [ ek(k1,k2,1)  ek(k1,k2,2) ; ek(k1,k2,2)  ek(k1,k2,1)  ] ;
       elseif Nelm==3
          mat_eps = [ ek(k1,k2,1)  ek(k1,k2,2) ; ek(k1,k2,2)  ek(k1,k2,3)  ] ;
       end
          band(k1,k2,:) = eig(mat_eps) ;
   end
  end
    
% Getting initial filling0. 
% Note: fermi-func is slightly different from summing Gn
if fix_fill==0 
  fk = fermi(band(:),1/T);
  filling0 = 2*sum(fk(:))/nktot  % 2 for two spins
end


muLim = [min(band(:)) max(band(:))];
Wband = muLim(2) - muLim(1);
nuscale = Ncut*Wband/2;  % matsubara freq cutoff
% xifill0 = T*log(2*Norb/filling0 - 1);    

% sig=smearing parameter for delta function
 sig = 0.03  ;
 mixing = 1  ; 
 wgt = 0.2  ;  % mixing weight
 
 %mix_method: 0=mix self-energy;  1=mix Green's function; 
 %            2=mix Green's function &  Anderson acceleration
 if Nelm==2
     mix_method = 2;
 elseif Nelm==3
     mix_method = 1;  
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fileDir = './';
fileDir = './output/';  % make sure to creat folder "output"
filenamestr = ['_Nk=' num2str(Nk(1)) '_' num2str(Nk(2)) '_U=' num2str(U) '.dat'];
fileout = ['out' filenamestr];
filegaps = ['gaps' filenamestr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid = fopen([fileDir,fileout],'w');
% fclose(fid);
diary([fileDir, fileout]);
disp(' ');
fprintf('  # of element = %4.0f  \n', Nelm  );
fprintf('  grid [nky nkx] = [%4d,%4d], convergence critera = %g \n',  Nk(1),Nk(2),eps)
fprintf('  U = %g [eV]\n',U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fidgaps = fopen([fileDir filegaps],'a');



T = [];   % initializing T  (do not remove)
for nt = 1:numel(Tlist)
    for iiu = 1:numel(Ulist)
        
        tic;
        Told = T;
        Uold = U;
        T = Tlist(nt);
        U = Ulist(iiu);
        beta = 1/T;
        wn = pi/beta;
        numwi = round(nuscale/wn);
        if Power2NumFreq == 1
           numwi = 2048;   % fix max_numwi = 1024 (default)
        end
        
        xifill0 = T*log(2*Norb/filling0 - 1);  
        
        fprintf('\n')
        fprintf('**************************************************************\n')
        fprintf('  T = %g [K], U = %g [eV]\n',T/kb,U)
        %fprintf('  xifill0 = %g [eV]\n',xifill0)
        fprintf(' numwi =%6d, Power2NumFreq = %1d\n',numwi,Power2NumFreq)
        
        if nt>1
            WNold = WN;
        end
        
        %Define the freqency grid
        WN = pi*(2*(-numwi:numwi-1)+1)/beta;
        WNU = pi*2*(-numwi:numwi-1)/beta;
        if saveSelfE == 1              
            fileSelfE{nt,iiu} = ['selfE_Nk=' num2str(Nk(1)) '_' num2str(Nk(2)) '_U=' num2str(U) '_T=' num2str(T) '.mat'];
            fprintf('  self-energy %s%s will be saved.\n',fileDir,fileSelfE{nt,iiu})
        end
        if (saveSelfE == 1) && exist([fileDir fileSelfE{nt,iiu}],'file')
            fprintf('  self-energy %s%s exists. Go to next T and U.\n',fileDir,fileSelfE{nt,iiu})
            continue
        end
        
        
        if (loadSelfE == 1) && (nt>1 || iiu>1)
            if numel(Tlist)>1 && numel(Ulist)>1 && saveSelfE~=1
                error('  set saveSelfE = 1!')
            end
            if numel(Tlist) == 1 && numel(Ulist) > 1
                if flySelfE ~= 1
                    load([fileDir fileSelfE{nt,iiu-1}],'S','P','WN','mu','-mat')
                end
            elseif  numel(Tlist) > 1 && numel(Ulist) == 1
                if flySelfE == 1
                    [S P mu] = interpSelfE(WN,1,S,P,WNold,mu);
                else
                    [S P mu] = interpSelfE(WN,0,[fileDir fileSelfE{nt-1,1}]);
                end
            else
                if (Told~=T)
                    [S P mu] = interpSelfE(WN,0,[fileDir fileSelfE{nt-1,iiu}]);
                elseif (Uold~=U)
                    load([fileDir fileSelfE{nt,iiu-1}],'S','P','WN','mu','-mat')
                end
            end
              wn = pi/beta;
            for jj = 1:Nelm
                gap(jj) = max(max( abs(real(P(:,:,numwi+1,jj))) ));
                if abs(gap(jj)) < 0.001
                    for nn = 1:length(WN)
                        P(:,:,nn,jj) = gap_input(jj);
                    end
                end
            end
        else
            
            clear('S','P')
              S(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
              P(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
           
            %Initial guess of self-energies
            if Nelm==2
                for nn=1:2*numwi
                   %S(:,:,nn,1) = -1i*sign(WN(nn));
                    S(:,:,nn,1) = -1i*WN(nn);
                  if gap_symm==1
                      P(:,:,nn,:) = gap_input(1) ;
                  elseif gap_symm==2
                      P(:,:,nn,:) = gap_input(1)*(cos(KX)-cos(KY));
                  end
                end
            elseif Nelm==3
                for nn=1:2*numwi
                    % S(:,:,nn,[1 3]) = -1i*WN(nn);  
                    S(:,:,nn,[1 3]) = -1i*sign(WN(nn));
                   if gap_symm==1
                      P(:,:,nn,1:3) = gap_input(1) ;
                    % P(:,:,nn,1) = gap_input(1)*(cos(KX).*cos(KY) ).^1;
                    %P(:,:,nn,2) = gap_input(1)*(cos(KX).*cos(KY) ).^1;
                    %P(:,:,nn,3) = gap_input(1)*(cos(KX).*cos(KY) ).^1;
                   elseif gap_symm==2
                       P(:,:,nn,1:3) = gap_input(1)*(cos(KX)-cos(kY));
                   end
                end
            end      
              mu = get_mu(S,P,WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
        end
        
        %surf( P(:,:,1,1), 'edgecolor','none' )
        %return
        
        mu0 = get_mu(zeros(size(S)),zeros(size(P)),WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
        filling = get_filling(S,P,WN,ek,mu,xifill0,Norb,beta,useSymmetry);

        fprintf('  filling0= %12.8f, mu0= %12.8f \n',filling0,mu0)        
        fprintf('  filling = %12.8f, mu = %12.8f \n',filling,mu)

        dk = gauss(band-mu,sig);
        Nf0 = sum(dk(:))/nktot;        
        fprintf('  N(E_f,T =%8.6f) = %12.8f (1/eV/per spin/Volume)\n',sig,Nf0)     
        
      
        calculate_im_axis;
        wn = pi/beta;

        if Nelm==3
            maxchi=0;
        end
        
        if (U == 0), maxchi = 0; end
        if (maxchi < 1)
            for jj = 1:Nelm
                gapwn = real(P(:,:,numwi+1,jj));
                gapmax(jj) = max(gapwn(:));
                gapmin(jj) = min(gapwn(:));
            end
            if Nelm == 2
                fprintf(fidgaps,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %g\n', ...
                    T,mu,gapmax(1),gapmin(1),gapmax(2),gapmin(2),eps);
            elseif Nelm == 3
                fprintf(fidgaps,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %g\n', ...
                    T,mu,gapmax(1),gapmin(1),gapmax(2),gapmin(2),gapmax(3),gapmin(3),eps);
            end
        end
        if (U==0), clear('maxchi'); end
      
        
        %fprintf(fidgaps,'\n');
        if saveSelfE == 1
            save([fileDir fileSelfE{nt,iiu}],'S','P','WN','mu','-mat')
        end
        
        vtm = toc;
        
        chis00 = permute(chis0(:,:,numwi+1,:),[1 2 4 3]);
        chic00 = permute(chic0(:,:,numwi+1,:),[1 2 4 3]);
        chisR0 = permute(chisR(:,:,numwi+1,:),[1 2 4 3]);
        Vnm0 = permute(Vnm(:,:,numwi+1,:),[1 2 4 3]);
        S0 = permute(S(:,:,numwi+1,:),[1 2 4 3]);
        P0 = permute(P(:,:,numwi+1,:),[1 2 4 3]);
        % w dependence q = (0 0)
        chis0w1 = reshape(chis0(1,1,:,:),1,[]);
        chic0w1 = reshape(chic0(1,1,:,:),1,[]);
        chisRw1 = reshape(chisR(1,1,:,:),1,[]);
        Vnmw1 = reshape(Vnm(1,1,:,:),1,[]);
        Sw1 = reshape(S(1,1,:,:),1,[]);
        Pw1 = reshape(P(1,1,:,:),1,[]);
        % w dependence q = (pi pi)
        chis0w2 = reshape(chis0(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        chic0w2 = reshape(chic0(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        chisRw2 = reshape(chisR(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        Vnmw2 = reshape(Vnm(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        Sw2 = reshape(S(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        Pw2 = reshape(P(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        % w dependence q = (0 pi)
        chis0w3 = reshape(chis0(Nk(1)+1,1,:,:),1,[]);
        chic0w3 = reshape(chic0(Nk(1)+1,1,:,:),1,[]);
        chisRw3 = reshape(chisR(Nk(1)+1,1,:,:),1,[]);
        Vnmw3 = reshape(Vnm(Nk(1)+1,1,:,:),1,[]);
        Sw3 = reshape(S(Nk(1)+1,1,:,:),1,[]);
        Pw3 = reshape(P(Nk(1)+1,1,:,:),1,[]);
        save([fileDir 'selfE_Nk=' num2str(Nk(1)) '_' num2str(Nk(2))  '_T=' num2str(T) '_U=' num2str(U) '_freqFFT.mat'],...
            'Nk','numwi','T','mu','WN','WNU',...
            'S0','P0','chi*00','chisR0','Vnm0',...
            'Sw*','Pw*','chi*0w*','chisRw*','Vnmw*','ek','-mat')

        fprintf('Done gap calculation. Total Time = %.2f s\n',sum(vtm));
        %save([fileDir fileout(1:end-3) 'mat'],'-mat')
    end  % end U
    
    gapval(nt) = abs(max(max(real(P(:,:,numwi+1,1)) )) ) ; 
    
end  % end T 

delta = real(P(:,:,numwi+1,:));

figure(2)
  surf( K{1}, K{2}, delta(:,:,1), 'edgecolor','none' ) 
  %shading interp
  colorbar
  axis([ 0 2*pi  0 2*pi ])
  xlabel('k_x [0,\pi)'); ylabel('k_y [0,\pi)');
  view(2)


figure(3)
  plot( Tlist/kb, gapval, '-ob' )
  xlabel('T [K]');  ylabel('\phi_{max}(1)');  
  
nnt = 1 ;
for nt = 1:length(Tlist)
    if( gapval(nt) > 0.001 )  % ignore gap<1meV for better fit
        tempgap(nnt) = gapval(nt);
        tp(nnt) = Tlist(nt)/kb  ;
        nnt = nnt + 1 ;
    end
end


if( max(tempgap)>0.001 && length(Tlist)>=2 )
 x0=[  0.04;  3;  30  ]; % initial guess
 %lb=[  0  0  10    ];  ub=[ 1  10  50    ];  
 fac = lsqcurvefit( @BCS_fit, x0, tp, tempgap  )
 temperature=linspace(0,100,500);
 temp_gap = fac(1)*tanh( fac(2)*sqrt( 1-temperature./(fac(3) )) );
figure(3)
   plot( temperature, real(temp_gap), tp, tempgap, 'o' )
end


fclose(fidgaps);
diary off
