function []=SPEAD_1D_save(dat14,dat15,dat2020a,dat1010,...
    dat1020,dat1022,dat1030,dat2010,dat2020,dat2030,dat70,dat80,dat24)
%Save output from SPEAD_1D.
%
%Usage : 
%
% [dat14,dat15,dat2020a,dat1010,dat1020,dat1022,dat1030,...
%  dat2010,dat2020,dat2030,dat70,dat80,dat24]=SPEAD_1D();
% 
% SPEAD_1D_save( dat14,dat15,dat2020a,dat1010,dat1020,dat1022,dat1030,...
%                       dat2010,dat2020,dat2030,dat70,dat80,dat24)

    fil14=SPEAD_1D_save_fig14(dat14)
end

function [fil]=SPEAD_1D_save_one(fil,nam,dat)
    tmp1=struct();
    for i=1:length(nam)
        tmp1.(nam{i})=dat(i);
    end
    save(fil,'-struct','tmp1');
end

function [fil]=SPEAD_1D_save_fig14(dat)
    fil=fullfile(tempdir(),"SPEAD_1D_fig14.mat");
    nam={'itemp','iparz0','PAR2D','imld','iKZ','fignum','mypackages'};
    SPEAD_1D_save_one(fil,nam,dat);
end
