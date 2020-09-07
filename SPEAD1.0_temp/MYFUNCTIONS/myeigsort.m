function [vec,val] = myeigsort(A,k)
    
%****************************************************************************
% Use: [vec,val] = myeigsort(A,ndims)
%****************************************************************************
%============================================================================
%............................................................................
%COMPUTE USING EIG + SORT FUNCTIONS:
[VECTOR_nxn,LAMBDA_nxn] = eig(full(A));
[Lambda_diag] = diag(LAMBDA_nxn);
[Lambda_sort,Indexe_sort] = sort(Lambda_diag);
[Lambda_ndim] = Lambda_sort(end-k+1:end);
[Column_ndim] = Indexe_sort(end-k+1:end);
[Vector_eigsort] = VECTOR_nxn(:,Column_ndim);
[Lambda_eigsort] = diag(Lambda_ndim);
%............................................................................
%COMPUTE USING EIGS FUNCTION:
N = size(A,1);
opts.v0 = rand(N,1); %Starting vector for eigs (see QUESTION/ANSWER below)
[Vector_eigs,Lambda_eigs] = eigs(A,k,'LR',opts);
%............................................................................
%TESTING:
Vector_eig = [Vector_eigsort,Vector_eigs]; %Should be the same.
Lambda_eig = [Lambda_eigsort,Lambda_eigs]; %Should be the same.
%............................................................................
%OUTPUTS:
vec = Vector_eigsort;
val = Lambda_eigsort;
%............................................................................
%============================================================================
%****************************************************************************
return

%============================================================================
%----------------------------------------------------------------------------
% <https://www.mathworks.com/matlabcentral/answers/224516-why-matlab-function-eigs-has-different-results-for-the-same-input-data>
%----------------------------------------------------------------------------
%QUESTION: Why Matlab function eigs has different results for the same input data?
%Latest activity Commented on by Steven Lord on 10 Apr 2017
%Asked by Honggui Li on 19 Jun 2015
%Accepted Answer by John D'Errico
%----------------------------------------------------------------------------
%ANSWER: This is because eigs uses a random start. So it need not generate
%exactly the same result every time it is called. Worse, since
%eigenvectors are not unique(!), they may vary in the sign of those
%vectors generated. That is, you can multiply all elements of an
%eigenvector by -1 (actually, any constant will suffice) and still have
%an equally valid eigenvector for that eigenvalue. As far as getting a
%stable answer, you can set the seed for the random generator that eigs
%will use to some fixed value. That will cause eigs to start from the
%same point every time, so the result will be stable (by some arbitrary
%definition of stability.) Alternatively, you can pass in an options
%structure to eigs, with a starting vector already supplied. 
%Thus "opts.v0" must be an Nx1 vector.
%----------------------------------------------------------------------------
% k = 5;
% opts.v0 = rand(npoints,1);
% [V1,D1] = eigs(A,k,'LM',opts);
% [V2,D2] = eigs(A,k,'LM',opts);
%----------------------------------------------------------------------------
%============================================================================

