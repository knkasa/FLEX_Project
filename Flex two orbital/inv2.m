function A = inv2(A)

a11 = A(:,1); a21 = A(:,2); a12 = A(:,3); a22 = A(:,4);
invd = 1./(a11.*a22 - a21.*a12);
A(:,1) = a22.*invd;
A(:,2) = -a21.*invd;
A(:,3) = -a12.*invd;
A(:,4) = a11.*invd;
