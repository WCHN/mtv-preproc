lat_ref = [172 221 168];
vs_ref  = [1 1 1];
lat_sub = [28 221 168];
vs_sub  = [6 1 1];

A = [6 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0];
J = [6 0 0; 0 1 0; 0 0 1];
J = single(reshape(J, [1 1 1 3 3]));

y = spm_warps('affine', A, lat_sub);

ref = randn(lat_ref, 'single');
sub = zeros(lat_sub, 'single');
sub(14,111,84) = 1;

kernel = blur_fun(lat_ref, eye(3), [6 1 1]);

%%
tic, ref2sub = spm_diffeo('pull', ref, y); toc
tic, sub2ref = spm_diffeo('push', sub, y, lat_ref); toc

%%

tic, ref2sub = pushpull('pull', ref, y, J); toc
tic, sub2ref = pushpull('push', sub, y, J, lat_ref); toc

%%

tic, ref2sub = spm_diffeo('pull',single(real(ifftn(fftn(ref).*kernel))),single(y)); toc
tic, sub2ref = single(real(ifftn(fftn(spm_diffeo('push',sub,single(y),lat_ref)).*kernel))); toc

%%

u = randn(lat_ref, 'single');
v = randn(lat_sub, 'single');
Au = pushpull('pull', u, y, J);
Atv = pushpull('push', v, y, J, lat_ref);
v_dot_Au = v(:)' * Au(:);
Atv_dot_v = Atv(:)' * u(:);
v_dot_Au - Atv_dot_v