function Nii_H = approx_hessian(Nii_H,Nii_x,dat)
% Compute approximation to the diagonal of the Hessian 
%
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(dat);

for c=1:C        
    A1 = A(ones(dat(c).dm,'single'),dat(c));
    N  = numel(A1);
    for n=1:N
        X           = get_nii(Nii_x{c}(n));
        msk         = get_msk(X);
        clear X
        A1{n}(~msk) = NaN;
%         figure(111);imagesc3d(reshape(msk,dat(c).A.dm))
    end
    H        = At(A1,dat(c));   
%     msk      = get_msk(H);
%     figure(111);imagesc3d(H)
%     H        = max(H(:))*ones(size(H),'single'); % Infinity norm
%     H(~msk)  = 0;
    Nii_H(c) = put_nii(Nii_H(c),H);
end    
%==========================================================================