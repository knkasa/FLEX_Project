
% Getting chi matrix.  
% G(kx,ky,wn,3) dimension.
% chi_(ijmn) =  G_im G_nj  +-  F_in G_mj    
% Note: f(-tau) -> -f(beta-tau), and  F_ji(tau)=F_ij(-tau)=-F_ij(beta-tau)

Xs0(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:4,1:4) = 0;
Xc0(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:4,1:4) = 0;
 Vn(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:4,1:4) = 0;
 Va(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:4,1:4) = 0;
matXc0(1:4,1:4)=0;
matXc0(1:4,1:4)=0;

%X(1,1)    G11 G11  -+ F11 F11
  Gbeta = -Gn(:,:,1,1);   
  Gbeta(1,1) = Gbeta(1,1) - 1;  
Xs0(:,:,1,1,1) =  Gn(:,:,1,1).*Gbeta - Fn(:,:,1,1).*Fn(:,:,1,1);
Xc0(:,:,1,1,1) =  Gn(:,:,1,1).*Gbeta + Fn(:,:,1,1).*Fn(:,:,1,1);
Xs0(:,:,2:end,1,1) = Gn(:,:,2:end,1).*Gn(:,:,end:-1:2,1) - Fn(:,:,2:end,1).*Fn(:,:,2:end,1);
Xc0(:,:,2:end,1,1) = Gn(:,:,2:end,1).*Gn(:,:,end:-1:2,1) + Fn(:,:,2:end,1).*Fn(:,:,2:end,1);

%X(4,4)    % G22 G22 -+ F22 G22 
  Gbeta = -Gn(:,:,1,3);   
  Gbeta(1,1) = Gbeta(1,1) - 1;  
Xs0(:,:,1,4,4) =  Gn(:,:,1,3).*Gbeta - Fn(:,:,1,3).*Fn(:,:,1,3);
Xc0(:,:,1,4,4) =  Gn(:,:,1,3).*Gbeta + Fn(:,:,1,3).*Fn(:,:,1,3);
Xs0(:,:,2:end,4,4) = Gn(:,:,2:end,3).*Gn(:,:,end:-1:2,3) - Fn(:,:,2:end,3).*Fn(:,:,2:end,3);
Xc0(:,:,2:end,4,4) = Gn(:,:,2:end,3).*Gn(:,:,end:-1:2,3) + Fn(:,:,2:end,3).*Fn(:,:,2:end,3);

%X(1,4)     chi_1122   G12 G21 +- F12 F21
 Gbeta = -Gn(:,:,1,2);  
 Fbeta = -Fn(:,:,1,2);  
Xs0(:,:,1,1,4) = Gn(:,:,1,2).*Gbeta + Fn(:,:,1,2).*Fbeta;
Xc0(:,:,1,1,4) = Gn(:,:,1,2).*Gbeta - Fn(:,:,1,2).*Fbeta;
Xs0(:,:,2:end,1,4) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,2) + Fn(:,:,2:end,2).*Fn(:,:,end:-1:2,2);
Xc0(:,:,2:end,1,4) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,2) - Fn(:,:,2:end,2).*Fn(:,:,end:-1:2,2);

%X(1,3)        chi_1112      G11 G21  -+  F12 F11 
  Gbeta = -Gn(:,:,1,2);     
Xs0(:,:,1,1,2) =  Gn(:,:,1,1).*Gbeta - Fn(:,:,1,2).*Fn(:,:,1,1);
Xc0(:,:,1,1,2) =  Gn(:,:,1,1).*Gbeta + Fn(:,:,1,2).*Fn(:,:,1,1);
Xs0(:,:,2:end,1,2) = Gn(:,:,2:end,1).*Gn(:,:,end:-1:2,2) - Fn(:,:,2:end,2).*Fn(:,:,2:end,1);
Xc0(:,:,2:end,1,2) = Gn(:,:,2:end,1).*Gn(:,:,end:-1:2,2) + Fn(:,:,2:end,2).*Fn(:,:,2:end,1);

%X(1,2)           chi_1121     G12 G11 +- F11 F21
  Gbeta = -Gn(:,:,1,1);    
  Gbeta(1,1) = Gbeta(1,1) - 1;  
  Fbeta = -Fn(:,:,1,2);  
Xs0(:,:,1,1,3) = Gn(:,:,1,2).*Gbeta + Fn(:,:,1,1).*Fbeta;
Xc0(:,:,1,1,3) = Gn(:,:,1,2).*Gbeta - Fn(:,:,1,1).*Fbeta;
Xs0(:,:,2:end,1,3) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,1) + Fn(:,:,2:end,1).*Fn(:,:,end:-1:2,2);
Xc0(:,:,2:end,1,3) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,1) - Fn(:,:,2:end,1).*Fn(:,:,end:-1:2,2);

%X(3,3)           chi_1212    G11 G22  -+ F12 F12
  Gbeta = -Gn(:,:,1,3);     
  Gbeta(1,1) = Gbeta(1,1) - 1; 
Xs0(:,:,1,2,2) = Gn(:,:,1,1).*Gbeta - Fn(:,:,1,2).*Fn(:,:,1,2);
Xc0(:,:,1,2,2) = Gn(:,:,1,1).*Gbeta + Fn(:,:,1,2).*Fn(:,:,1,2);
Xs0(:,:,2:end,2,2) = Gn(:,:,2:end,1).*Gn(:,:,end:-1:2,3) - Fn(:,:,2:end,2).*Fn(:,:,2:end,2);
Xc0(:,:,2:end,2,2) = Gn(:,:,2:end,1).*Gn(:,:,end:-1:2,3) + Fn(:,:,2:end,2).*Fn(:,:,2:end,2);

%X(2,3)   chi_2112     G21 G21 -+ F22 F11
  Gbeta = -Gn(:,:,1,2); 
Xs0(:,:,1,3,2) =  Gn(:,:,1,2).*Gbeta - Fn(:,:,1,3).*Fn(:,:,1,1);
Xc0(:,:,1,3,2) =  Gn(:,:,1,2).*Gbeta + Fn(:,:,1,3).*Fn(:,:,1,1);
Xs0(:,:,2:end,3,2) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,2) - Fn(:,:,2:end,3).*Fn(:,:,2:end,1);
Xc0(:,:,2:end,3,2) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,2) + Fn(:,:,2:end,3).*Fn(:,:,2:end,1);

%X(3,4)    chi_1222    G12 G22 -+ F12 F22
  Gbeta = -Gn(:,:,1,3);   
  Gbeta(1,1) = Gbeta(1,1) - 1; 
Xs0(:,:,1,3,4) =  Gn(:,:,1,2).*Gbeta - Fn(:,:,1,2).*Fn(:,:,1,3);
Xc0(:,:,1,3,4) =  Gn(:,:,1,2).*Gbeta + Fn(:,:,1,2).*Fn(:,:,1,3);
Xs0(:,:,2:end,3,4) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,3) - Fn(:,:,2:end,2).*Fn(:,:,2:end,3);
Xs0(:,:,2:end,3,4) = Gn(:,:,2:end,2).*Gn(:,:,end:-1:2,3) + Fn(:,:,2:end,2).*Fn(:,:,2:end,3);

%X(2,2)   chi_2121   G22 G11  -+ F21 F21
  Gbeta = -Gn(:,:,1,1);   
  Gbeta(1,1) = Gbeta(1,1) - 1; 
  Fbeta = -Fn(:,:,1,2);
Xs0(:,:,1,3,3) = Gn(:,:,1,3).*Gbeta - Fbeta.*Fbeta;
Xc0(:,:,1,3,3) = Gn(:,:,1,3).*Gbeta + Fbeta.*Fbeta;
Xs0(:,:,2:end,3,3) = Gn(:,:,2:end,3).*Gn(:,:,end:-1:2,1) - Fn(:,:,end:-1:2,2).*Fn(:,:,end:-1:2,2);
Xc0(:,:,2:end,3,3) = Gn(:,:,2:end,3).*Gn(:,:,end:-1:2,1) + Fn(:,:,end:-1:2,2).*Fn(:,:,end:-1:2,2);

%X(2,4)    chi_2122    G22 G21 +- F22 F21
  Gbeta = -Gn(:,:,1,2);
  Fbeta = -Fn(:,:,1,2);
Xs0(:,:,1,3,4) = Gn(:,:,1,3).*Gbeta + Fn(:,:,1,3).*Fbeta;
Xc0(:,:,1,3,4) = Gn(:,:,1,3).*Gbeta - Fn(:,:,1,3).*Fbeta;
Xs0(:,:,2:end,3,4) = Gn(:,:,2:end,3).*Gn(:,:,end:-1:2,2) + Fn(:,:,2:end,3).*Fn(:,:,end:-1:2,2);
Xc0(:,:,2:end,3,4) = Gn(:,:,2:end,3).*Gn(:,:,end:-1:2,2) - Fn(:,:,2:end,3).*Fn(:,:,end:-1:2,2);


 %backward Fourier transform susceptibility (r,tau) to (k,w_n)
 %careful about discontinuity.
 % gamma(ij,mn) = -G_im(0)del_nj  +  G_nj(0)del_im
   Xs0(:,:,:,1,1) = fourier_bwd(Xs0(:,:,:,1,1),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,1,1) = fourier_bwd(Xc0(:,:,:,1,1),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   
   Xs0(:,:,:,1,2) = fourier_bwd(Xs0(:,:,:,1,2),Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,1,2) = fourier_bwd(Xc0(:,:,:,1,2),Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
 
   Xs0(:,:,:,1,3) = fourier_bwd(Xs0(:,:,:,1,3),-Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,1,3) = fourier_bwd(Xc0(:,:,:,1,3),-Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
 
   Xs0(:,:,:,1,4) = fourier_bwd(Xs0(:,:,:,1,4),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,1,4) = fourier_bwd(Xc0(:,:,:,1,4),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   
   Xs0(:,:,:,2,2) = fourier_bwd(Xs0(:,:,:,2,2),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,2,2) = fourier_bwd(Xc0(:,:,:,2,2),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);

   Xs0(:,:,:,2,3) = fourier_bwd(Xs0(:,:,:,2,3),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,2,3) = fourier_bwd(Xc0(:,:,:,2,3),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);

   Xs0(:,:,:,2,4) = fourier_bwd(Xs0(:,:,:,2,4),-Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,2,4) = fourier_bwd(Xc0(:,:,:,2,4),-Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);

   Xs0(:,:,:,3,3) = fourier_bwd(Xs0(:,:,:,3,3),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,3,3) = fourier_bwd(Xc0(:,:,:,3,3),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);

   Xs0(:,:,:,3,4) = fourier_bwd(Xs0(:,:,:,3,4),Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,3,4) = fourier_bwd(Xc0(:,:,:,3,4),Gn(1,1,1,2),NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);

   Xs0(:,:,:,4,4) = fourier_bwd(Xs0(:,:,:,4,4),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
   Xc0(:,:,:,4,4) = fourier_bwd(Xc0(:,:,:,4,4),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);

  
% chi matrix are symmetric
Xs0(:,:,:,2,1) = Xs0(:,:,:,1,2);
Xc0(:,:,:,2,1) = Xc0(:,:,:,1,2);
Xs0(:,:,:,3,1) = Xs0(:,:,:,1,3);
Xc0(:,:,:,3,1) = Xc0(:,:,:,1,3);
Xs0(:,:,:,4,1) = Xs0(:,:,:,1,4);
Xc0(:,:,:,4,1) = Xc0(:,:,:,1,4);
Xs0(:,:,:,3,2) = Xs0(:,:,:,2,3);
Xc0(:,:,:,3,2) = Xc0(:,:,:,2,3);
Xs0(:,:,:,4,2) = Xs0(:,:,:,2,4);
Xc0(:,:,:,4,2) = Xc0(:,:,:,2,4);
Xs0(:,:,:,4,3) = Xs0(:,:,:,3,4);
Xc0(:,:,:,4,3) = Xc0(:,:,:,3,4);

%Xs0(:,:,numwi+1,:,:) = Xs0(:,:,numwi,:,:);
%Xc0(:,:,numwi+1,:,:) = Xs0(:,:,numwi,:,:);

% Getting V_eff matrix
%{
for ii = 1:2*Nk(1)
 for jj = 1:2*Nk(2)
  for n = 1:2*numwi
      
      matXs0 = reshape(Xs0(ii,jj,n,:,:),4,4); % make it 4x4 matrix
      matXc0 = reshape(Xc0(ii,jj,n,:,:),4,4);
      
      matXs = ( eye(4)-matXs0*Us )\matXs0;  % RPA form. inv(B)*A order.  Not A*inv(B)
      matXc = ( eye(4)+matXc0*Uc )\matXc0;
      
    % No hartree term.  It can not use fourier transform.  
   Vn(ii,jj,n,:,:) =  1.5*Us*matXs*Us + 0.5*Uc*matXc*Uc ...
       -3/8*Us*(matXs0+matXc0)*Us - 1/8*Uc*(matXs0+matXc0)*Uc  ;
   
   Va(ii,jj,n,:,:) = 1.5*Us*matXs*Us - 0.5*Uc*matXc*Uc ...
       -3/8*Us*(matXs0-matXc0)*Us - 1/8*Uc*(matXs0-matXc0)*Uc  ;
   
  end
 end
end
%}

%
Vn = eval_Vn( Us, Uc, eval_Xs( Us, Xs0, size(Xs0) ), eval_Xc( Uc, Xc0, size(Xc0) ), Xs0, Xc0, size(Xs0)  );
Va = eval_Va( Us, Uc, eval_Xs( Us, Xs0, size(Xs0) ), eval_Xc( Uc, Xc0, size(Xc0) ), Xs0, Xc0, size(Xs0)  );
%}

% fourier transfrom  wn -> tau
for ii = 1:4
 for jj = 1:4
    Vn(:,:,:,ii,jj) = fourier_fwd(Vn(:,:,:,ii,jj),[],beta,Nk,numwi,0,0,useSymmetry,runInSerial);
    Va(:,:,:,ii,jj) = fourier_fwd(Va(:,:,:,ii,jj),[],beta,Nk,numwi,0,0,useSymmetry,runInSerial);
 end
end

Vn = real(Vn);
Va = real(Va);


%Getting self-energies.  Note: it cant be matrix multiplications.
S(:,:,:,1) = Vn(:,:,:,1,1).*Gn(:,:,:,1) + Vn(:,:,:,1,2).*Gn(:,:,:,2) + Vn(:,:,:,2,1).*Gn(:,:,:,2) + Vn(:,:,:,2,2).*Gn(:,:,:,3) ;
S(:,:,:,2) = Vn(:,:,:,1,3).*Gn(:,:,:,1) + Vn(:,:,:,1,4).*Gn(:,:,:,2) + Vn(:,:,:,2,3).*Gn(:,:,:,2) + Vn(:,:,:,2,4).*Gn(:,:,:,3) ;
S(:,:,:,3) = Vn(:,:,:,3,3).*Gn(:,:,:,1) + Vn(:,:,:,3,4).*Gn(:,:,:,2) + Vn(:,:,:,4,3).*Gn(:,:,:,2) + Vn(:,:,:,4,4).*Gn(:,:,:,3) ;

%{
P(:,:,:,1) = Va(:,:,:,1,1).*Fn(:,:,:,1) + Va(:,:,:,1,2).*Fn(:,:,:,2) + Va(:,:,:,3,1).*Fn(:,:,:,2) + Va(:,:,:,3,2).*Fn(:,:,:,3) ;
P(:,:,:,2) = Va(:,:,:,1,3).*Fn(:,:,:,1) + Va(:,:,:,1,4).*Fn(:,:,:,2) + Va(:,:,:,3,3).*Fn(:,:,:,2) + Va(:,:,:,3,4).*Fn(:,:,:,3) ;
P(:,:,:,3) = Va(:,:,:,2,3).*Fn(:,:,:,1) + Va(:,:,:,2,4).*Fn(:,:,:,2) + Va(:,:,:,4,3).*Fn(:,:,:,2) + Va(:,:,:,4,4).*Fn(:,:,:,3) ;
%}

%
Fbeta = -Fn(:,:,1,2);
P(:,:,1,1) = Va(:,:,1,1,1).*Fn(:,:,1,1) + Va(:,:,1,1,3).*Fn(:,:,1,2)   ... 
              - Va(:,:,1,2,1).*Fbeta + Va(:,:,1,2,3).*Fn(:,:,1,3) ;
P(:,:,2:end,1) = Va(:,:,2:end,1,1).*Fn(:,:,2:end,1) + Va(:,:,2:end,1,3).*Fn(:,:,2:end,2) ...
          - Va(:,:,2:end,2,1).*Fn(:,:,end:-1:2,2) + Va(:,:,2:end,2,3).*Fn(:,:,2:end,3) ;

P(:,:,1,2) = Va(:,:,1,1,2).*Fn(:,:,1,1) + Va(:,:,1,1,4).*Fn(:,:,1,2) ...
               - Va(:,:,1,2,3).*Fbeta + Va(:,:,1,2,4).*Fn(:,:,1,3) ;
P(:,:,2:end,2) = Va(:,:,2:end,1,2).*Fn(:,:,2:end,1) + Va(:,:,2:end,1,4).*Fn(:,:,2:end,2) ...
          - Va(:,:,2:end,2,3).*Fn(:,:,end:-1:2,2) + Va(:,:,2:end,2,4).*Fn(:,:,2:end,3) ;
      
P(:,:,1,3) = Va(:,:,1,3,2).*Fn(:,:,1,1) + Va(:,:,1,3,4).*Fn(:,:,1,2) ... 
               - Va(:,:,1,4,2).*Fbeta + Va(:,:,1,4,4).*Fn(:,:,1,3) ;
P(:,:,2:end,3) = Va(:,:,2:end,3,2).*Fn(:,:,2:end,1) + Va(:,:,2:end,3,4).*Fn(:,:,2:end,2) ... 
          - Va(:,:,2:end,4,2).*Fn(:,:,end:-1:2,2) + Va(:,:,2:end,4,4).*Fn(:,:,2:end,3) ;
%}
      
% Fourier transform. Check the discontinuity
 S(:,:,:,1) = fourier_bwd(S(:,:,:,1),-Vn(1,1,1,1,1)-Vn(1,1,1,2,2),NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
 S(:,:,:,2) = fourier_bwd(S(:,:,:,2),-Vn(1,1,1,1,3)-Vn(1,1,1,2,4),NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
 S(:,:,:,3) = fourier_bwd(S(:,:,:,3),-Vn(1,1,1,3,3)-Vn(1,1,1,4,4),NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
 
  for jj = 1:3 
     P(:,:,:,jj) = fourier_bwd(P(:,:,:,jj),0,NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
  end













