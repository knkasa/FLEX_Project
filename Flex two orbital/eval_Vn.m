function V = eval_Vn( Us, Uc, Xs, Xc, Xs0, Xc0, dim  )

% evaluatin Vn

V = zeros( dim(1), dim(2),dim(3), 4, 4 );

for i = 1:4
 for j = i:4
     V(:,:,:,i,j) = 0 ; 
   for m=1:4
   for n=1:4
       V(:,:,:,i,j) = V(:,:,:,i,j) + 1.5*Us(i,m)*Xs(:,:,:,m,n)*Us(n,j)  ...
          + 0.5*Uc(i,m)*Xc(:,:,:,m,n)*Uc(n,j)    ...
          - 3/8*Us(i,m)*( Xs0(:,:,:,m,n)+Xc0(:,:,:,m,n) )*Us(n,j)  ...
           - 1/8*Uc(i,m)*(  Xs0(:,:,:,m,n)+Xc0(:,:,:,m,n) )*Uc(n,j)  ;
   end
   end
 end
end

V(:,:,:,2,1) = V(:,:,:,1,2);
V(:,:,:,3,1) = V(:,:,:,1,3);
V(:,:,:,3,2) = V(:,:,:,2,3);
V(:,:,:,4,1) = V(:,:,:,1,4);
V(:,:,:,4,2) = V(:,:,:,2,4);
V(:,:,:,4,3) = V(:,:,:,3,4);

