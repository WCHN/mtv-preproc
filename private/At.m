function Y = At(X,tau,dat)  
Y = single(0);   
for n=1:dat.N      
    T = dat.mat\dat.A(n).mat;
    y = apply_affine(T,dat.A(n).dm);
    
    tmp = spm_diffeo('push',X{n},y,dat.dm);     
    tmp(~isfinite(tmp)) = 0;
       
    tmp = real(ifftn(fftn(tmp).*dat.A(n).S));  
     
    Y = Y + tau(n).*tmp;           
end
%==========================================================================