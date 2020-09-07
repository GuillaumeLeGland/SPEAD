function [irow,jcol] = myrawiseindexing(counter,irow,jcol,msize,keyRawise)

%========================================================================
%........................................................................
if strcmp(keyRawise,'not') %Col wise indexing
    if mod(counter,msize) == 1 %colwise indexing.
	jcol = jcol + 1;
	irow = 0;
    end 
    irow = irow + 1; 
%........................................................................
elseif strcmp(keyRawise,'yes') %Row wise indexing
    if mod(counter,msize) == 1 %rowwise indexing.
	irow = irow + 1;
	jcol = 0;
    end 
    jcol = jcol + 1; 
end
%........................................................................
%========================================================================
return


