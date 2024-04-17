function recon = reconstruction_MBIR(rawdata,param)

    P=param.P;
    P{2}=param.P{2}*max(abs(col(rawdata{param.shot})));
    recon{param.shot} = gather(regularizedReconstruction(param.F{param.shot},rawdata{param.shot},P{:},'maxit',20,...
        'verbose_flag', 0,'tol',1e-5));
    printf('complete')
        
end
