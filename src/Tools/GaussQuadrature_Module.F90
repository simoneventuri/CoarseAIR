Module GaussQuadrature_Module

  use Parameters_Module           ,only:  rkp, Zero, Half, One, Two, Three, Four, Pi

  implicit none

  contains 

  Pure Subroutine GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)
        implicit real*8 (A-H,O-Z)
  !     b is temp, t is nodes, w is weights.
  !
  !           THIS SET OF ROUTINES COMPUTES THE NODES X(I) AND WEIGHTS
  !        C(I) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED
  !        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE
  !
  !                 INTEGRAL (FROM A TO B)  F(X) W(X) DX
  !
  !                              N
  !        BY                   SUM C  F(X )
  !                             I=1  I    I
  !
  !        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT
  !        FUNCTIONS (LISTED BELOW), AND F(X) IS THE
  !        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY
  !        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT
  !        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.
  !
  !           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF
  !        ORTHOGONAL POLYNOMIALS.  THE NODES X(I) ARE JUST THE ZEROES
  !        OF THE PROPER N-TH DEGREE POLYNOMIAL.
  !
  !     INPUT PARAMETERS
  !
  !        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF
  !                 QUADRATURE RULE
  !
  !        KIND = 1=  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)
  !        KIND = 2=  CHEBYSHEV QUADRATURE OF THE FIRST KIND
  !                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)
  !        KIND = 3=  CHEBYSHEV QUADRATURE OF THE SECOND KIND
  !                   W(X) = SQRT(1 - X*X) ON (-1, 1)
  !        KIND = 4=  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON
  !                   (-INFINITY, +INFINITY)
  !        KIND = 5=  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**
  !                   BETA ON (-1, 1), ALPHA, BETA  >  -1.
  !                   NOTE= KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.
  !        KIND = 6=  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*
  !                   X**ALPHA ON (0, +INFINITY), ALPHA  >  -1
  !
  !        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE
  !        ALPHA    DOUBLE PRECISION PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-
  !                 LAGUERRE QUADRATURE (OTHERWISE USE 0.).
  !        BETA     DOUBLE PRECISION PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE--
  !                 (OTHERWISE USE 0.).
  !        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-
  !                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A
  !                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO
  !                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED
  !                 ENDPOINTS (1 OR 2).
  !        ENDPTS   DOUBLE PRECISION ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF
  !                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.
  !        B        DOUBLE PRECISION SCRATCH ARRAY OF LENGTH N
  !
  !     OUTPUT PARAMETERS (BOTH ARRAYS OF LENGTH N)
  !
  !        T        WILL CONTAIN THE DESIRED NODES X(1),,,X(N)
  !        W        WILL CONTAIN THE DESIRED WEIGHTS C(1),,,C(N)
  !
  !   SubroutineS REQUIRED
  !
  !        GBSLVE, CLASS, AND GBTQL2 ARE PROVIDED. UNDERFLOW MAY SOMETIMES
  !        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE
  !        TURNED OFF AS THEY ARE ON THIS MACHINE.
  !
  !     ACCURACY
  !
  !        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,
  !        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP
  !        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,
  !        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT
  !        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR
  !        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME
  !        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,
  !        COMPLETELY HARMLESS.
  !
  !     METHOD
  !
  !           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION
  !        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE
  !        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE
  !        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH
  !        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF
  !        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,
  !        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A
  !        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.
  !        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF
  !        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.
  !
  !     REFERENCES
  !
  !        1.  GOLUB, G. H., AND WELSCH, J. H.,  CALCULATION OF GAUSSIAN
  !            QUADRATURE RULES,  MATHEMATICS OF COMPUTATION 23 (APRIL,
  !            1969), PP. 221-230.
  !        2.  GOLUB, G. H.,  SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,
  !            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).
  !        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE-
  !            HALL, ENGLEWOOD CLIFFS, N.J., 1966.
  !
  !     ..................................................................
  !
        integer                 ,intent(in)    :: KIND
        integer                 ,intent(in)    :: N
        real(rkp)               ,intent(in)    :: ALPHA
        real(rkp)               ,intent(in)    :: BETA
        integer                 ,intent(in)    :: KPTS
        real(rkp) ,dimension(2) ,intent(in)    :: ENDPTS
        real(rkp) ,dimension(N) ,intent(inout) :: B
        real(rkp) ,dimension(N) ,intent(out)   :: T
        real(rkp) ,dimension(N) ,intent(out)   :: W

        real(rkp)                              :: MUZERO
        integer                                :: K1, I, IERR

  !
  !$$$      if ( kind /= 132)stop
        CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)
  !$$$      if ( kind /= 132)stop
  !
  !           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.
  !           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY
  !           B THE OFF-DIAGONAL ELEMENTS.
  !           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2
  !           SUBMATRIX.
  !
        IF (KPTS == 0 )  GO TO 100
        IF (KPTS == 2)  GO TO  50
  !
  !           IF KPTS=1, ONLY T(N) MUST BE CHANGED
  !
        T(N) =GBSLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)
        GO TO 100
  !
  !           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED
  !
     50 GAM =GBSLVE(ENDPTS(1), N, T, B)
        T1 = ((ENDPTS(1) - ENDPTS(2))/(GBSLVE(ENDPTS(2), N, T, B) - GAM))
        B(N-1) =  SQRT(T1)
        T(N) = ENDPTS(1) + GAM*T1
  !
  !           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
  !           AND THUS THE VALUE OF B(N) IS ARBITRARY.
  !           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL
  !           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.
  !           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING
  !
    100 W(1) = One
        DO 105 I = 2, N
    105    W(I) = Zero
  !
        CALL GBTQL2 (N, T, B, W, IERR)
  !$$$      if ( kind /= 6 ) then
        DO 110 I = 1, N
    110    W(I) = MUZERO * W(I) * W(I)
  !$$$      else
  !$$$       write(6,*)('returning sqrt of weights ')
  !$$$      end if
  End Subroutine


  real(rkp) Pure Function GBSLVE(SHIFT, N, A, B)
        implicit real*8 (A-H,O-Z)
  !
  !       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE
  !       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION
  !
  !             (JN - SHIFT*IDENTITY) * DELTA  = EN,
  !
  !       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN
  !       THE N-TH POSITION.
  !
  !       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL
  !       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION
  !       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER
  !       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.
        real(rkp)               ,intent(in)  :: SHIFT
        integer                 ,intent(in)  :: N
        real(rkp) ,dimension(N) ,intent(in)  :: A
        real(rkp) ,dimension(N) ,intent(in)  :: B

        integer                              :: I, NM1

        ALPHA = A(1) - SHIFT
        NM1 = N - 1
        DO 10 I = 2, NM1
     10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA
        GBSLVE = One  /ALPHA
  End Function


  Pure Subroutine CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)
        implicit real*8 (A-H,O-Z)
  !
  !           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE
  !        RECURRENCE RELATION
  !
  !             B P (X) = (X - A ) P   (X) - B   P   (X)
  !              J J            J   J-1       J-1 J-2
  !
  !        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,
  !        AND THE ZERO-TH MOMENT
  !
  !             MUZERO = INTEGRAL W(X) DX
  !
  !        OF THE GIVEN POLYNOMIAL   WEIGHT FUNCTION W(X).  SINCE THE
  !        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS
  !        GUARANTEED TO BE SYMMETRIC.
  !
  !           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND
  !        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR
  !        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS
  !        REQUIRE THE GAMMA FUNCTION.
  !
  !     ..................................................................
        integer                 ,intent(in)  :: KIND
        integer                 ,intent(in)  :: N
        real(rkp)               ,intent(in)  :: ALPHA
        real(rkp)               ,intent(in)  :: BETA
        real(rkp) ,dimension(N) ,intent(out) :: A
        real(rkp) ,dimension(N) ,intent(out) :: B
        real(rkp)               ,intent(out) :: MUZERO

        integer                              :: K1, I, NM1

  !$$$      if ( n /= 132)stop
  !
        NM1 = N - 1
        GO TO (10, 20, 30, 40, 50, 60), KIND
  !
  !              KIND = 1=  LEGENDRE POLYNOMIALS P(X)
  !              ON (-1, +1), W(X) = 1.
  !
     10 MUZERO = Two
        DO 11 I = 1, NM1
           A(I) = Zero
           ABI = I
     11    B(I) = ABI/ SQRT(4*ABI*ABI - One  )
        A(N) = Zero
        RETURN
  !
  !              KIND = 2=  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)
  !              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)
  !
     20 MUZERO = PI
        DO 21 I = 1, NM1
           A(I) = Zero
     21    B(I) = Half
        B(1) =  SQRT(Half)
        A(N) = Zero
        RETURN
  !
  !              KIND = 3=  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)
  !              ON (-1, +1), W(X) = SQRT(1 - X*X)
  !
     30 MUZERO = PI/Two
        DO 31 I = 1, NM1
           A(I) = Zero
     31    B(I) = Half
        A(N) = Zero
        RETURN
  !
  !              KIND = 4=  HERMITE POLYNOMIALS H(X) ON (-INFINITY,
  !              +INFINITY), W(X) = EXP(-X**2)
  !
     40 MUZERO =  SQRT(PI)
        DO 41 I = 1, NM1
           A(I) = Zero
     41    B(I) =  SQRT(I/Two)
        A(N) = Zero
        RETURN
  !
  !              KIND = 5=  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON
  !              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND
  !              BETA GREATER THAN -1
  !
     50 AB = ALPHA + BETA
        ABI = Two   + AB
        if ( abi > 60d0 ) then
         x1=gammln(ALPHA + One)
         x2=gammln(BETA  + One)
         x3=gammln(abi)
         x4=log(2d0)
         muzero=exp((ab+One)*x4+x1+x2-x3)
        else
        MUZERO = Two   ** (AB + One  ) * APPROXGAMMA(ALPHA + One  ) * &
        APPROXGAMMA(BETA + One  ) / APPROXGAMMA(ABI)
        end if
        A(1) = (BETA - ALPHA)/ABI
        B(1) =  SQRT(Four  *(One   + ALPHA)*(One   + BETA)/((ABI + One  )*ABI*ABI))
        A2B2 = BETA*BETA - ALPHA*ALPHA
        DO 51 I = 2, NM1
           ABI = Two  *I + AB
           A(I) = A2B2/((ABI - Two  )*ABI)
     51    B(I) =  SQRT (Four  *I*(I + ALPHA)*(I + BETA)*(I + AB)/   ((ABI*ABI - 1)*ABI*ABI))
        ABI = Two  *N + AB
        A(N) = A2B2/((ABI - Two  )*ABI)
        RETURN
  !
  !              KIND = 6=  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON
  !              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER
  !              THAN -1.
  !
  !$$$   60 MUZERO = APPROXGAMMA(ALPHA + One  )
     60 continue
  !$$$      if ( nm1 /= 132)stop
        if ( ALPHA + One > 60d0 ) then
         !write(6,*)('in class for generalized laguerre polynomials')
         !write(6,*)('replace muzero with unity to avoid overflow')
         muzero=One
        else
         MUZERO = APPROXGAMMA(ALPHA + One  )
        end if
  !$$$      if ( nm1 /= 132)stop
        DO 61 I = 1, NM1
           A(I) = Two  *I - One   + ALPHA
     61    B(I) =  SQRT(I*(I + ALPHA))
        A(N) = Two  *N - 1 + ALPHA
  End Subroutine


  Pure Subroutine GBTQL2(N, D, E, Z, IERR)
        implicit real*8 (A-H,O-Z)
  !
  !     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
  !     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
  !     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
  !
  !     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE
  !     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL
  !     METHOD, AND IS ADAPTED FROM THE EISPAK ROUTINE IMTQL2
  !
  !     ON INPUT=
  !
  !        N IS THE ORDER OF THE MATRIX;
  !
  !        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
  !
  !        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
  !          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;
  !
  !        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.
  !
  !      ON OUTPUT=
  !
  !        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
  !          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
  !          UNORDERED FOR INDICES 1, 2, ..., IERR-1;
  !
  !        E HAS BEEN DESTROYED;
  !
  !        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
  !          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
  !          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
  !          EIGENVALUES;
  !
  !        IERR IS SET TO
  !
  !        IERR IS SET TO
  !          ZERO       FOR NORMAL RETURN,
  !          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
  !                     DETERMINED AFTER 30 ITERATIONS.
  !
  !     ------------------------------------------------------------------
  !
        integer                 ,intent(in)     :: N
        real(rkp) ,dimension(N) ,intent(inout)  :: D
        real(rkp) ,dimension(N) ,intent(inout)  :: E
        real(rkp) ,dimension(N) ,intent(inout)  :: Z
        integer                 ,intent(out)    :: IERR

        INTEGER                                 :: I, J, K, L, M, II, MML
        real(rkp)                               :: MACHEP
  !
  !     ========== MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
  !                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
  !                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
  !                ON S360 ==========
         MACHEP=1.0e-14_rkp
         xold=machep
     87  continue
         try=One+machep
         if ( try /= One ) then
          xold=machep
          machep=machep*Half
          go to 87
         end if
         machep=xold
  !$$$       write(6,*)('machep '),machep
  !
        IERR = 0
        IF (N  ==  1) GO TO 1001
  !
        E(N) = Zero
        DO 240 L = 1, N
           J = 0
  !     ========== LOOK FOR SMALL SUB-DIAGONAL ELEMENT ==========
    105    DO 110 M = L, N
              IF (M  ==  N) GO TO 120
              IF ( ABS(E(M))  <=  MACHEP * ( ABS(D(M)) +  ABS(D(M+1)))) GO TO 120
    110    CONTINUE
  !
    120    P = D(L)
           IF (M  ==  L) GO TO 240
           IF (J  ==  30) GO TO 1000
           J = J + 1
  !     ========== FORM SHIFT ==========
           G = (D(L+1) - P) / (Two   * E(L))
           R =  SQRT(G*G+One )
           G = D(M) - P + E(L) / (G +  SIGN(R, G))
           S = One
           C = One
           P = Zero
           MML = M - L
  !     ========== FOR I=M-1 STEP -1 UNTIL L DO -- ==========
           DO 200 II = 1, MML
              I = M - II
              F = S * E(I)
              B = C * E(I)
              IF ( ABS(F)  <   ABS(G)) GO TO 150
              C = G / F
              R =  SQRT(C*C+One )
              E(I+1) = F * R
              S = One   / R
              C = C * S
              GO TO 160
    150       S = F / G
              R =  SQRT(S*S+One )
              E(I+1) = G * R
              C = One   / R
              S = S * C
    160       G = D(I+1) - P
              R = (D(I) - G) * S + Two   * C * B
              P = S * R
              D(I+1) = G + P
              G = C * R - B
  !     ========== FORM FIRST COMPONENT OF VECTOR ==========
              F = Z(I+1)
              Z(I+1) = S * Z(I) + C * F
              Z(I) = C * Z(I) - S * F
  !
    200    CONTINUE
  !
           D(L) = D(L) - P
           E(L) = G
           E(M) = Zero
           GO TO 105
    240 CONTINUE
  !     ========== ORDER EIGENVALUES AND EIGENVECTORS ==========
        DO 300 II = 2, N
           I = II - 1
           K = I
           P = D(I)
  !
           DO 260 J = II, N
              IF (D(J) >= P) GO TO 260
              K = J
              P = D(J)
    260    CONTINUE
  !
           IF (K  ==  I) GO TO 300
           D(K) = D(I)
           D(I) = P
  !
           P = Z(I)
           Z(I) = Z(K)
           Z(K) = P
  !
    300 CONTINUE
  !
        GO TO 1001
  !     ========== SET ERROR -- NO CONVERGENCE TO AN
  !                EIGENVALUE AFTER 30 ITERATIONS ==========
   1000 IERR = L
   1001 RETURN
  !     ========== LAST CARD OF GBTQL2 ==========
  End Subroutine


  real(rkp) Pure Function  APPROXGAMMA(ZZ)
        implicit real*8 (A-H,O-Z)
!     MODIFIED 4/4/88 BY DWS TO BE ACCURATE FOR ZZ GT 3.
  !  THIS IS A PROCEDURE THAT EVALUATES GAMMA(Z) FOR
  !     0 LT Z LE 3 TO 16 SIGNIFICANT FIGURES
  !    IT IS BASED ON A CHEBYSHEV-TYPE POLYNOMIAL
  !   APPROXIMATION GIVEN IN H. WERNER AND R. COLLINGE, MATH. COMPUT.
  !    15 (1961), PP. 195-97.
  !     according to this reference, these parameters are for the
  !     approximation of gamma (2+x) with 0 le x le 1
  !     an odd choice, but what can I say.
  !   APPROXIMATIONS TO THE GAMMA FUNCTION, ACCURATE UP TO 18 SIGNIFICANT
  !   DIGITS, MAY BE FOUND IN THE PAPER QUOTED ABOVE
  !
        real(rkp) ,intent(in) :: ZZ
        real(rkp) ,dimension(18) , parameter :: A =[                    One,  .4227843350984678_rkp,  .4118403304263672_rkp, &
                                                      .0815769192502609_rkp,  .0742490106800904_rkp, -.0002669810333484_rkp, &
                                                      .0111540360240344_rkp, -.0028525821446197_rkp,  .0021036287024598_rkp, &
                                                     -.0009184843690991_rkp,  .0004874227944768_rkp, -.0002347204018919_rkp, &
                                                      .0001115339519666_rkp, -.0000478747983834_rkp,  .0000175102727179_rkp, &
                                                     -.0000049203750904_rkp,  .0000009199156407_rkp, -.0000000839940496_rkp]
        integer               :: K, K1

        PREFAC=One
        Z=ZZ
        if ( z > Three ) then
     44 CONTINUE
        if ( Z < Three)GO TO 45
        Z=Z-One
        PREFAC=PREFAC*Z
        GO TO 44
        else if ( z < Two ) then
     46  continue
          if ( z >= Two)go to 45
          prefac=prefac/z
          z=z+One
         go to 46
        end if
     45 CONTINUE
  !
        if ( Z <= One  ) GO TO 10
        if ( Z <= Two  ) GO TO 20
        T=Z-Two
        GO TO 30
  10    T=Z
        GO TO 30
  20    T=Z-One
  30    P=A(18)
        DO 40 K1=1,17
        K=18-K1
        P=T*P+A(K)
  40    CONTINUE
  !
        P=P*PREFAC
        if ( Z > Two  ) GO TO 50
        if ( Z > One  ) GO TO 60
        APPROXGAMMA=P/(Z*(Z+One  ))
        RETURN
  60    APPROXGAMMA=P/Z
        RETURN
  50    APPROXGAMMA=P
  End Function


  ! log gamma to about 1 part in 1d10 from numerical recipes, 2nd ed, page 207.
  real(rkp) Pure Function gammln(xx)
    implicit real*8 (a-h,o-z)

    real(rkp) ,intent(in) :: xx
    integer               :: j

    !save !!! WARNING: SIMONE REMOVED THE SAVE AND MADE IT PURE; SEEMS WORKING FINE WITHOUT
    real(rkp) ,dimension(6) ,parameter :: cof = [76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0,-1.231739572450155d0,0.1208650973866179d-2, -0.5395239384953d-5]
    real(rkp)               ,parameter :: stp = 2.5066282746310005d0

    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+Half)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
     y=y+One
     ser=ser+cof(j)/y
    end do
    gammln=tmp+log(stp*ser/x)

  End Function


End Module 