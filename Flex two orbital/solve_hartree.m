function  [ S_HF P_HF ] = solve_hartree(Gn,Fn,Us,Uc,pre)  
% Getting hartree for general case. (not bilayer)

 S_HF(1:3) = 0;
 P_HF(1:3) = 0;

 Uself = 1.5*Us - 0.5*Uc ;
  Uphi = 0.5*Us + 0.5*Uc ;
 
  S_HF(1) = Uself(1,1)+Uself(2,2) + pre*Uself(1,1)*sum(sum(sum(Gn(:,:,:,1)))) ...
            + pre*Uself(2,2)*sum(sum(sum(Gn(:,:,:,3)))) ...
            + 2*pre*Uself(1,2)*sum(sum(sum(Gn(:,:,:,2))))  ;

  S_HF(3) = Uself(3,3)+Uself(4,4) +  pre*Uself(3,3)*sum(sum(sum(Gn(:,:,:,1))))  ...
        + pre*Uself(4,4)*sum(sum(sum(Gn(:,:,:,3)))) ...
        + 2*pre*Uself(3,4)*sum(sum(sum(Gn(:,:,:,2))))  ;
 
  S_HF(2) = Uself(1,3)+Uself(2,4) + pre*Uself(1,3)*sum(sum(sum(Gn(:,:,:,1))))  ...
          + pre*Uself(2,4)*sum(sum(sum(Gn(:,:,:,3)))) ...
        + pre*Uself(1,4)*sum(sum(sum(Gn(:,:,:,2))))  ...  
        + pre*Uself(2,3)*sum(sum(sum(Gn(:,:,:,2))))   ;
           
S_avg   = 0.5*(S_HF(1)+S_HF(3))  ;
S_HF(1) =  S_HF(1) - S_avg;
S_HF(3) = S_HF(3) - S_avg;

P_HF(1) = pre*Uphi(1,1)*sum(sum(sum(Fn(:,:,:,1)))) + pre*Uphi(2,3)*sum(sum(sum(Fn(:,:,:,3)))) ...
    + pre*Uphi(1,3)*sum(sum(sum(Fn(:,:,:,2)))) + pre*Uphi(2,1)*sum(sum(sum(Fn(:,:,:,2)))) ; 

P_HF(3) = pre*Uphi(3,2)*sum(sum(sum(Fn(:,:,:,1)))) + pre*Uphi(4,4)*sum(sum(sum(Fn(:,:,:,3)))) ...
    + pre*Uphi(3,4)*sum(sum(sum(Fn(:,:,:,2)))) + pre*Uphi(4,2)*sum(sum(sum(Fn(:,:,:,2)))) ; 

P_HF(2) = pre*Uphi(1,2)*sum(sum(sum(Fn(:,:,:,1)))) + pre*Uphi(2,4)*sum(sum(sum(Fn(:,:,:,3)))) ...
    + pre*Uphi(1,4)*sum(sum(sum(Fn(:,:,:,2)))) + pre*Uphi(2,2)*sum(sum(sum(Fn(:,:,:,2)))) ; 





