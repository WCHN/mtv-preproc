function compile_mex(fname_mex)
f = fileparts(which('compile_mex'));
mex(fullfile(f,fname_mex),'-outdir',f);
%==========================================================================