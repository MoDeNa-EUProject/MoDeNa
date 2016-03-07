

                     EXPLANATION FILE OF MODULE LU
                     =============================


  LU Decomposition
  ----------------

  Suppose we are able to write the matrix A as a product of two matrices,

                 L . U = A           (2.3.1)

  where L is a lower triangular (has elements only on the diagonal ans below)
  and U is upper triangular ((has elements only on the diagonal ans above).

  For the case of a 4 x 4 matrix A, for example, equation (2.3.1) would look
  like this:

       |a11 0   0   0  | |b11 b12 b13 b14|   |A11 A12 A13 A14|
       |a21 a22 0   0  | |0   b22 b23 b24| = |A21 A22 A23 A24|  (2.3.2)
       |a31 a32 a33 0  | |0   0   b33 b34|   |A31 A32 A33 A34|
       |a41 a42 a43 a44| |0   0   0   b44|   |A41 A42 A43 A11|  
 	 
  We can use a decomposition such as (2.3.1) to solve the linear set
   
           A . x = (L . U) . x = L . (U . x) = b                (2.3.3)
            
  by first solving for the vector y such that
   
           L . y = b                                            (2.3.4) 
 
  and then solving 

           U . x = y                                            (2.3.5)
   
    What is the advantage of breaking up one linear set into two successive 
  ones? The advantage is that the solution of a triangular set of equations is 
  quite trivial, as we have already seen before. Thus, equation (2.3.4) can be
  solved by forward substitution as follows, 

           y  = b  / a 
            1    1    11
		              i-1       
           y  = (1/a ) [ b - Sum a   Y  ]    i = 2,3,...,N      (2.3.6)
	    i       11    i  j=1  ij  j

  while (2.3.5) can then be solved by backsubstitution: 

           x  = y  / b
            N    N    NN
                               N
           x  = (1/b ) [ y -  Sum  b   x ]   i=N-1,N-2,...,1    (2.3.7)
            i       ii    i  j=i+1  ij  j

    Equations (2.3.6) and (2.3.7) total (for each right-hand side b) N2 exe- 
  cutions of an inner loop containing one multiply and one add. If we have N
  right-hand sides which are the unit column vectors (which is the case when we
  are inverting a matrix), then taking into account the leading zeros reduces
  the total execution count of (2.3.6) from N^3/2 to N^3/6, while (2.3.7) is
  unchanged at N^3/2.
   
    Notice that, once we have the LU decomposition of A we can solve with as 
  many right-hand sides as we then care to, one at a time. This is a distinct 
  advantage over the methods seen before.
   
  Performing the LU Decomposition
  -------------------------------
   
    How then can we solve for Land U, given A? First, we write out the i,j th 
  component of equation (2.3.1) or (2.3.2). That component always is a sum
  beginning with

            a   b   + ... = A             with notations of (2.3.2)
             i1  1j          ij

  The number of terms in the sum depends, however, on whether i or j is the 
  smaller number. We have, in fact, the three cases, 

          i < j:     a  b  + a  b  + ... + a   b   = A      (2.3.8)
                      i1 1j   i2 2j         ii  ij    ij

          i = j:     a  b  + a  b  + ... + a   b   = A      (2.3.9)
                      i1 1j   i2 2j         ii  jj    ij

          i > j:     a  b  + a  b  + ... + a   b   = A     (2.3.10)
                      i1 1j   i2 2j         ij  jj    ij

    
    Equations (2.3.8)-(2.3.10) total N^2 equations for the N^2 + N unknown
  a's and b's (the diagonal being represented twice). Since the number of 
  unknowns is greater than the number of equations, we are invited to specify 
  N of the unknowns arbitrarily and then try to solve for the others. ln fact, 
  as we shall see, it is always possible to take

            a  = 1       i = 1,...,N                       (2.3.11)
             ii

    A surprising procedure, now, is Crout's algorithm, which quite trivially 
  solves the set of N^2 + N equations (2.3.8)-(2.3.11) for all the a's and b's
  by just arranging the equations in a certain order! That order is as follows: 

       • Set a  = 1,  i = 1,...,N                  (equation 2.3.11).
              ii

       • For each j = 1,2,3, ... , N do these two procedures: First, for i = 
         1,2, ... .i, use (2.3.8), (2.3.9), and (2.3.11) to solve for b  , 
         namely                                                        ij
                              i-1
                   b   = A  - Sum a   b
                    ij    ij  k=1  ik  kj                  (2.3.12)

    (When i = 1 in 2.3.12 the summation term is taken to mean zero.) Second,
    for i = j + l,j + 2, ... ,N use (2.3.10) to solve for a  , namely 
		j-1 	)                                ij
                                        j-1
                   a   = (1/b)  ( aij - Sum a   b  )       (2.3.13) 
                    ij       jj         k=1  ik  kj
 
    Be sure to do both procedures before going on to the next j.
     
    If you work through a few iterations of the above procedure, you will see 
  that the a's and b's that occur on the right-hand side of equations (2.3.12) 
  and (2.3.13) are already determined by the time they are needed. You will 
  also see that every aij is used only once and never again. This means that 
  the corresponding aij or bij can be stored in the location that the A used 
  to occupy: the decomposition is "in place." (The diagonal unity elements aii 
  (equation 2.3.11) are not stored at aIl). ln brief, Crout's method fills in
  the combined matrix of a's and b's,

                  | b11 b12 b13 b14 |
                  | a21 b22 b23 b24 |                     (2.3.14)
                  | a31 a32 b33 b34 |
                  | a41 a42 a43 b44 |
   
  by columns from left to right, and within each column from top to bottom. 
 
    What about pivoting? Pivoting (i.e. selection of a salubrious pivot ele- 
  ment for the division in equation 2.3.13) is absolutely essential for the sta-
  bility of Crout's method. Only partial pivoting (interchange of rows) can be
  implemented efficiently. However this is enough to make the method stable.
  This means, incidentally, that we don't actually decompose the matrix A into
  LU form, but rather we decompose a rowwise permutation of A. (If we keep track
  of what that permutation is, this decomposition is just as use ful as the 
  original one would have been).
   
    Pivoting is slightly subtle in Crout's algorithm. The key point to notice 
  is that equation (2.3.12) in the case of i = J (its final application) is
  exactly the same as equation (2.3.13) except for the division in the latter
  equation; in both cases the upper limit of the sum is k = j - 1 (= i - 1).
  This means that we don't have to commit ourselves as to whether the diagonal
  element bjj is the one which happens to fall on the diagonal in the first
  instance, or whether one of the (undivided) aij's below it in the column,
  i = j + 1, ... ,N, is to be "promoted" to become the diagonal b. This can be
  decided after all the candidates in the column are in hand. As you should be
  able to guess by now, we will choose the largest one as the diagonal b (pivot
  element), then do all the divisions by that element en masse. This is Crout's
  method with partial pivoting. Our implementation has one additional subtlety:
  it initially finds the largest element in each row, and subsequently (when it
  is looking for the maximal pivot element) scales the comparison as if we had
  initially scaled all the equations to make their maximum coefficient equal to
  unity; this is the implicit pivoting mentioned before.

    The LU decomposition in routine LUDCMP requires about N^3/3 executions of
  the inner loops (each with one multiply and one add). This is thus the ope-
  ration count for solving one (or a few) right-hand sides, and is a factor of 3
  better than the Gauss-Jordan routine GAUSSJ which was given in §2.1, and a
  factor of 1.5 better than a Gauss-Jordan routine which does not compute the
  inverse matrix. For inverting a matrix, the total count (inc1uding the forward
  and backsubstitution as discussed following equation 2.3.7 above) is (1/3 +
  1/6 + 1/2) N^3 = N^3, the same as GAUSSJ. 

    To summarize, this is the preferred way to solve the linear set of equations 
  A·x = b: 

    CALL LUDCMP(A,N,NP,INDX,D) 
    CALL LUBKSB(A,N,NP,INDX,B) 

    The answer x will be returned in B. Your original matrix A will have been
  destroyed. 

    If you subsequently want to solve a set of equations with the same A but 
  a different right-hand side b, you repeat only 

    CALL LUBKSB(A,N,NP,INDX,B) 

  not, of course, with the original matrix A, but with A and INDX as were 
  already returned from LUDCMP.

From [BIBLI 08].
--------------------------------------------------
End of file lu.txt
