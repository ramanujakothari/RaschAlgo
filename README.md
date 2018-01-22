# RaschAlgo
The code is basically a Fortran Module that uses Fortran GSL to first evaluate Wigner 3J symbols and then it stores them on an array. Then there is a separate method for the symbols retrieval. The method is based on Rasch Algorithm (J. Rasch & A.C.H. Yu EFFICIENT STORAGE SCHEME FOR PRECALCULATED WIGNER 3J, 6J AND GAUNT COEFFICIENTS, SIAM Journal of Scientific Computing 25, 4). It is very helpful if your work needs evaluation of same wigner symbols again and again. This scheme exploits the Regge symmetry satisfied by Wigner symbols. Just to appreciate the efficiency let us take the following example: let the Wigner symbol be

l1  l2  l3
m1  m2  m3

Assuming that l1 and l2 can go from l_min to l_max, then imposing triangular condition and m1+m2+m3 = 0, means evaluation of  1,973,989 symbols for l_min = 2 and l_max = 16. Thus number can be significantly reduced to just 205,647 which is about 10% of the original number if we use Rasch Algorithm. For any further queries and if you find any ways of improving the code please email me at quantummechanicskothari@gmail.com
