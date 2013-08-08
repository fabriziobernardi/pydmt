from numpy import array,empty
import sys

def column(matrix, i):
   return [row[i] for row in matrix]


def gaussj(A,v):

   v=column(v, 0)
   N = len(v)

   # Gaussian elimination
   for m in range(N):

     # Divide by the diagonal element
     div = A[m,m]
     A[m,:] /= div
     v[m] /= div

     # Now subtract from the lower rows
     for i in range(m+1,N):
         mult = A[i,m]
         A[i,:] -= mult*A[m,:]
         v[i] -= mult*v[m]

   # Backsubstitution

   #define length of x
   x = empty(N+1,float)
   for m in range(N-1,-1,-1):
       x[m] = v[m]
       for i in range(m+1,N):
           x[m] -= A[m,i]*x[i]

   return x
