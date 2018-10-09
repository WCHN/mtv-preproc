function X = A(Y,dat)
Y(~isfinite(Y)) = 0;
Y               = fftn(Y);    

X = cell(1,dat.N);
for n=1:dat.N
    T = dat.mat\dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);
     
    tmp  = real(ifftn(Y.*dat.A(n).S));  
    
    X{n} = samp0(tmp,y); 
    
    vx1  = sqrt(sum(dat.mat(1:3,1:3).^2));
    vx0  = sqrt(sum(dat.A(n).mat(1:3,1:3).^2));
    scl  = prod(vx1./vx0);
    X{n} = scl*X{n};
end
%==========================================================================  