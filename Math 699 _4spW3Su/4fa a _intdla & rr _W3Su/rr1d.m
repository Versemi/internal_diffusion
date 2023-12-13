function [cnt,imax] = rr1d(Nbugs,Ngrid)
%RR1D   Rotor-router walk in one dimension.
%   [CNT,IMAX] = RR2D(NBUGS) performs a rotor-router walk for NBUGS bugs and
%   returns a vector CNT that counts the number of times each bin was
%   visited.  The bins are ordered {-1,0,1,2,...}, with the -1 and 0 bins
%   the "cups" where the bug exits.  CNT has NBUGS rows.
%
%   The optional return argument IMAX gives the largest excursion from the
%   starting point for each bug.
%
%   References:
%
%   M. Kleber, "Goldbug Variations," The Mathematical Intelligencer,
%   Vol. 27, 55-63 (2005).
%
%   See also RR2D.

if nargin < 1, Nbugs = 100; end
if nargin < 2, Ngrid = 10; end

% Convert bin counts to arrow directions (+1/-1).
toarw = @(i) 1 - 2*mod(i,2);

% Bin counts, including the bins at -1 and 0 at the start.
cnt = zeros(Nbugs,Ngrid+2);
% Maximum excursion for each bug.
imax = ones(Nbugs,1);

for b = 1:Nbugs
  i = 1;  % bugs start at 1
  while i > 0
    if toarw(cnt(b,i+2)) == 1
      cnt(b,i+2) = cnt(b,i+2) + 1;
      i = i-2;
    elseif toarw(cnt(b,i+2)) == -1
      cnt(b,i+2) = cnt(b,i+2) + 1;
      i = i+1;
    end
    if i > imax(b), imax(b) = i; end
  end
  % If i <= 0, then the bug falls into the cups at -1 (cnt(1))and 0 (cnt(2)).
  cnt(b,i+2) = cnt(b,i+2) + 1;
  % Copy bin counts as starting point for the next bug.
  if b < Nbugs, cnt(b+1,:) = cnt(b,:); end
end

figure(1)
plot(imax,'k.','MarkerSize',10)
hold on
plot(cummax(imax),'r-','LineWidth',2)
hold off
xlabel('bug')
ylabel('max distance')


%[P,d] = hist(imax,unique(imax))

[Pm,dm] = hist(cummax(imax),unique(cummax(imax)))
