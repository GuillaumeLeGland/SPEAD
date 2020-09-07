function [] = mysavecell2text(filenameout,cellarray)
%====================================================================
% Use: mysavecell2text(filenameout,cellarray) 
%--------------------------------------------------------------------
% <http://stackoverflow.com/questions/26125280/matlab-saving-cell-array-to-text-file> 
%Example: 
% array = cell(10,1);
% for i=1:10
%     array{i} = ['someText ' num2str(i)];
% end
% fileID = fopen('celldata.dat','w');
% [nrows,ncols] = size(array);
% for row = 1:nrows
%     temp_str = array{row,:};
%     fprintf(fileID ,'%s\n', temp_str);
% end
% fclose(fileID);
%--------------------------------------------------------------------
%====================================================================
%.................................................................... 
fileID = fopen(filenameout,'w');
%.................................................................... 
[nrows,ncols] = size(cellarray);
%.................................................................... 
for irow = 1:nrows 
    cell_str = cellarray{irow,:};
    fprintf(fileID ,'%s\n', cell_str);
end
%.................................................................... 
fclose(fileID);
%.................................................................... 
%====================================================================
%********************************************************************
return

