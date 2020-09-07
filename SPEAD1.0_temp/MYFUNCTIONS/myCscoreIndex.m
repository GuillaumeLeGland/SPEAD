function [] = myCscoreIndex()

%------------------------------------------------------------------------
% Cij :: C-score for the two species sp1, sp2 in the given set of islands
% Sij :: The number of co-occurrences of sp1, sp2
% ri :: Number of islands in which sp1 has 1
% rj :: Number of islands in which sp2 has 1
% The checkerboard score (c-score) for the colonisation pattern is then
% calculated as the mean number of checkerboard units per species-pair in
% the community: For M species, there are P = M(M-1)/2 species-pairs, so C-score is calculated:
% The C-score is sensitive to the proportion of islands that are
% occupied, thereby confounding comparisons between matrices or 
% sets of species pairs within them. An extension of the C-score 
% therefore standardizes by the number of islands each species-pair occupies using:
%------------------------------------------------------------------------

for i = 1:nspecies

    ri = find()

    for j = 1:nspecies 

	Cij = (ri - Sij) * (rj - Sij)

    end
end

%************************************************************************
return

