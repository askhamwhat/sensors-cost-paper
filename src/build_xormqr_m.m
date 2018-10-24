fortheaderdir = fullfile(matlabroot,'extern','examples','refbook');
fortfile = fullfile(matlabroot,'extern','examples','refbook','fort.c');

mex('-v','-largeArrayDims',['-I' fortheaderdir],'xormqr_m.c',fortfile,'-lmwlapack')
