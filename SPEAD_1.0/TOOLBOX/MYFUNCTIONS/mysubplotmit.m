function handle = mysubplot(ny,nx,it,bor,pos)
%
% function handle = mysubplot(ny,nx,it,bor,pos)

  if nargin < 4, bor = 0.05; end
  if nargin < 5, pos = []; end
  if (length(bor) == 1) 
    bor(2) = bor(1);
  end
  if (length(bor) == 2)
    bor = [bor(1) bor(2) 1-2*bor(1) 1-2*bor(2)];
  end
  if (length(pos) < 1) pos(1)=0; end
  if (length(pos) < 2) pos(2)=0; end
  if (length(pos) < 3) pos(3)=1; end
  if (length(pos) < 4) pos(4)=1; end
  x0=pos(1);
  y0=pos(2);
  ww=pos(3);
  hh=pos(4);
  i = mod(it-1,nx)+1;
  j = (it-i)/nx+1;
  ps=[x0+ww*(i-1+bor(1))/nx y0+hh*(ny-j+bor(2))/ny ww*bor(3)./nx hh*bor(4)./ny];
  handle = subplot('Position',ps);
end

