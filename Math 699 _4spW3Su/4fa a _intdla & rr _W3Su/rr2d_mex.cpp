#include <cmath>
#ifdef RR2D_DEBUG
#include <iostream>
#include <cassert>
#endif
#include "mex.h"

extern void _main();

mwIndex ind(const mwIndex i, const mwIndex j, const mwSize N)
{
#ifdef RR2D_DEBUG
  assert(i >= 0 && i < N && j >= 0 && j < N);
#endif
  return (i + j*N);
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mwSize Nbugs = (mwSize) mxGetScalar(prhs[0]);
  double *grid0 = mxGetPr(prhs[1]);
  const mwSize N = mxGetN(prhs[1]);
  const mwSize Ngrid = (N-1)/2;		// size of grid quadrant
  plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
  double *grid = mxGetPr(plhs[0]);
  for (mwIndex i = 0; i < N*N; ++i) grid[i] = grid0[i];
  const int X0 = Ngrid, Y0 = Ngrid;	// center of grid

  for (mwIndex i = 1; i <= Nbugs; ++i)
    {
      if (i % 1000 == 0)
	{ mexPrintf("bug %8d\n",i); mexEvalString("drawnow;"); }

#ifdef RR2D_DEBUG
      std::cout << "\nbug " << i << std::endl;
#endif
      mwIndex X = X0, Y = Y0;
      while (true)
	{
#ifdef RR2D_DEBUG
	  std::cout << "X=" << X << " Y=" << Y << " ind=" << ind(X,Y,N) << " ";
#endif
	  if (!grid[ind(X,Y,N)])
	    {
	      // Grid cell is empty.
	      // Occupy grid by setting grid to 1.
	      grid[ind(X,Y,N)] = 1;
#ifdef RR2D_DEBUG
	      std::cout << std::endl;
#endif
	      break;
	    }
	  else
	    {
	      // Grid cell has an arrow.
	      // Rotate arrow ccw then move in that direction.
	      ++grid[ind(X,Y,N)];
	      int direc = ((int)(grid[ind(X,Y,N)]-1) % 4) + 1;
#ifdef RR2D_DEBUG
	      std::cout << "grid=" << grid[ind(X,Y,N)] << " direc=" << direc;
	      std::cout << std::endl;
#endif
	      switch (direc)
		{
		case 1: ++X; break; // E
		case 2: ++Y; break; // N
		case 3: --X; break; // W
		case 4: --Y; break; // S
		}
	    }
	}
    }
}
