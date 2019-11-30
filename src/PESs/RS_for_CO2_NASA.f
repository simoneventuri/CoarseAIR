      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)                       04460004
C                                                                       04460005
      INTEGER N,NM,IERR,MATZ                                            04460006
      REAL*8 A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)                         04460007
C                                                                       04460008
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 04460009
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     04460010
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             04460011
C     OF A REAL SYMMETRIC MATRIX.                                       04460012
C                                                                       04460013
C     ON INPUT:                                                         04460014
C                                                                       04460015
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    04460016
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            04460017
C        DIMENSION STATEMENT;                                           04460018
C                                                                       04460019
C        N  IS THE ORDER OF THE MATRIX  A;                              04460020
C                                                                       04460021
C        A  CONTAINS THE REAL SYMMETRIC MATRIX;                         04460022
C                                                                       04460023
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              04460024
C        ONLY EIGENVALUES ARE DESIRED;  OTHERWISE IT IS SET TO          04460025
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    04460026
C                                                                       04460027
C     ON OUTPUT:                                                        04460028
C                                                                       04460029
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER;                04460030
C                                                                       04460031
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO;              04460032
C                                                                       04460033
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN            04460034
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE           04460035
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO;            04460036
C                                                                       04460037
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   04460038
C                                                                       04460039
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        04460040
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         04460041
C                                                                       04460042
C     ------------------------------------------------------------------04460043
C                                                                       04460044
      IF (N .LE. NM) GO TO 10                                           04460045
      IERR = 10 * N                                                     04460046
      GO TO 50                                                          04460047
C                                                                       04460048
   10 IF (MATZ .NE. 0) GO TO 20                                         04460049
C     :::::::::: FIND EIGENVALUES ONLY ::::::::::                       04460050
      CALL  TRED1(NM,N,A,W,FV1,FV2)                                     04460051
      CALL  TQLRAT(N,W,FV2,IERR)                                        04460052
      GO TO 50                                                          04460053
C     :::::::::: FIND BOTH EIGENVALUES AND EIGENVECTORS ::::::::::      04460054
   20 CALL  TRED2(NM,N,A,W,FV1,Z)                                       04460055
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                     04460056
   50 RETURN                                                            04460057
C     :::::::::: LAST CARD OF RS ::::::::::                             04460058
      END                                                               04460059
C                                                                       77440001
C     ------------------------------------------------------------------77440002
C                                                                       77440003
      SUBROUTINE TRED1(NM,N,A,D,E,E2)                                   77440004
C                                                                       77440005
      INTEGER I,J,K,L,N,II,NM,JP1                                       77440006
      REAL*8 A(NM,N),D(N),E(N),E2(N)                                    77440007
      REAL*8 F,G,H,SCALE                                                77440008
      REAL*8 DSQRT,DABS,DSIGN                                           77440009
C                                                                       77440010
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,    77440011
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   77440012
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   77440013
C                                                                       77440014
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX                   77440015
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING                           77440016
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            77440017
C                                                                       77440018
C     ON INPUT:                                                         77440019
C                                                                       77440020
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         77440021
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          77440022
C          DIMENSION STATEMENT;                                         77440023
C                                                                       77440024
C        N IS THE ORDER OF THE MATRIX;                                  77440025
C                                                                       77440026
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          77440027
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               77440028
C                                                                       77440029
C     ON OUTPUT:                                                        77440030
C                                                                       77440031
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             77440032
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER         77440033
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED;        77440034
C                                                                       77440035
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;    77440036
C                                                                       77440037
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         77440038
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;      77440039
C                                                                       77440040
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    77440041
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        77440042
C                                                                       77440043
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        77440044
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         77440045
C                                                                       77440046
C     ------------------------------------------------------------------77440047
C                                                                       77440048
      DO 100 I = 1, N                                                   77440049
  100 D(I) = A(I,I)                                                     77440050
C     :::::::::: FOR I=N STEP -1 UNTIL 1 DO -- ::::::::::               77440051
      DO 300 II = 1, N                                                  77440052
         I = N + 1 - II                                                 77440053
         L = I - 1                                                      77440054
         H = 0.0D0                                                      77440055
         SCALE = 0.0D0                                                  77440056
         IF (L .LT. 1) GO TO 130                                        77440057
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::       77440058
         DO 120 K = 1, L                                                77440059
  120    SCALE = SCALE + DABS(A(I,K))                                   77440060
C                                                                       77440061
         IF (SCALE .NE. 0.0D0) GO TO 140                                77440062
  130    E(I) = 0.0D0                                                   77440063
         E2(I) = 0.0D0                                                  77440064
         GO TO 290                                                      77440065
C                                                                       77440066
  140    DO 150 K = 1, L                                                77440067
            A(I,K) = A(I,K) / SCALE                                     77440068
            H = H + A(I,K) * A(I,K)                                     77440069
  150    CONTINUE                                                       77440070
C                                                                       77440071
         E2(I) = SCALE * SCALE * H                                      77440072
         F = A(I,L)                                                     77440073
         G = -DSIGN(DSQRT(H),F)                                         77440074
         E(I) = SCALE * G                                               77440075
         H = H - F * G                                                  77440076
         A(I,L) = F - G                                                 77440077
         IF (L .EQ. 1) GO TO 270                                        77440078
         F = 0.0D0                                                      77440079
C                                                                       77440080
         DO 240 J = 1, L                                                77440081
            G = 0.0D0                                                   77440082
C     :::::::::: FORM ELEMENT OF A*U ::::::::::                         77440083
            DO 180 K = 1, J                                             77440084
  180       G = G + A(J,K) * A(I,K)                                     77440085
C                                                                       77440086
            JP1 = J + 1                                                 77440087
            IF (L .LT. JP1) GO TO 220                                   77440088
C                                                                       77440089
            DO 200 K = JP1, L                                           77440090
  200       G = G + A(K,J) * A(I,K)                                     77440091
C     :::::::::: FORM ELEMENT OF P ::::::::::                           77440092
  220       E(J) = G / H                                                77440093
            F = F + E(J) * A(I,J)                                       77440094
  240    CONTINUE                                                       77440095
C                                                                       77440096
         H = F / (H + H)                                                77440097
C     :::::::::: FORM REDUCED A ::::::::::                              77440098
         DO 260 J = 1, L                                                77440099
            F = A(I,J)                                                  77440100
            G = E(J) - H * F                                            77440101
            E(J) = G                                                    77440102
C                                                                       77440103
            DO 260 K = 1, J                                             77440104
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)                  77440105
  260    CONTINUE                                                       77440106
C                                                                       77440107
  270    DO 280 K = 1, L                                                77440108
  280    A(I,K) = SCALE * A(I,K)                                        77440109
C                                                                       77440110
  290    H = D(I)                                                       77440111
         D(I) = A(I,I)                                                  77440112
         A(I,I) = H                                                     77440113
  300 CONTINUE                                                          77440114
C                                                                       77440115
      RETURN                                                            77440116
C     :::::::::: LAST CARD OF TRED1 ::::::::::                          77440117
      END                                                               77440118
C                                                                       35440001
C     ------------------------------------------------------------------35440002
C                                                                       35440003
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                    35440004
      implicit real*8 (a-h,o-z)                                         8d3s93
C                                                                       35440005
      INTEGER I,J,L,M,N,II,L1,MML,IERR                                  35440006
      REAL*8 D(N),E2(N)                                                 35440007
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                     35440008
      REAL*8 DSQRT,DABS,DSIGN                                           35440009
C                                                                       35440010
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,   35440011
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                35440012
C                                                                       35440013
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              35440014
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     35440015
C                                                                       35440016
C     ON INPUT:                                                         35440017
C                                                                       35440018
C        N IS THE ORDER OF THE MATRIX;                                  35440019
C                                                                       35440020
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          35440021
C                                                                       35440022
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     35440023
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. 35440024
C                                                                       35440025
C      ON OUTPUT:                                                       35440026
C                                                                       35440027
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          35440028
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          35440029
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            35440030
C          THE SMALLEST EIGENVALUES;                                    35440031
C                                                                       35440032
C        E2 HAS BEEN DESTROYED;                                         35440033
C                                                                       35440034
C        IERR IS SET TO                                                 35440035
C          ZERO       FOR NORMAL RETURN,                                35440036
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               35440037
C                     DETERMINED AFTER 30 ITERATIONS.                   35440038
C                                                                       35440039
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        35440040
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         35440041
C                                                                       35440042
C     ------------------------------------------------------------------35440043
C                                                                       35440044
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     35440045
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   35440046
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        35440047
C                ON S360 ::::::::::                                     35440048
c$$$      DATA MACHEP/Z3410000000000000/                                    35440049
      factm=0.5d0
      machep=1d-12
  333 continue
       if(machep+1.d0.eq.1d0)go to 353
       machep=machep*factm
       go to 333
  353 continue
      machep=machep*2d0
C                                                                       35440050
      IERR = 0                                                          35440051
      IF (N .EQ. 1) GO TO 1001                                          35440052
C                                                                       35440053
      DO 100 I = 2, N                                                   35440054
  100 E2(I-1) = E2(I)                                                   35440055
C                                                                       35440056
      F = 0.0D0                                                         35440057
      B = 0.0D0                                                         35440058
      E2(N) = 0.0D0                                                     35440059
C                                                                       35440060
      DO 290 L = 1, N                                                   35440061
         J = 0                                                          35440062
         H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))                       35440063
         IF (B .GT. H) GO TO 105                                        35440064
         B = H                                                          35440065
         C = B * B                                                      35440066
C     :::::::::: LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT :::::::::: 35440067
  105    DO 110 M = L, N                                                35440068
            IF (E2(M) .LE. C) GO TO 120                                 35440069
C     :::::::::: E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              35440070
C                THROUGH THE BOTTOM OF THE LOOP ::::::::::              35440071
  110    CONTINUE                                                       35440072
C                                                                       35440073
  120    IF (M .EQ. L) GO TO 210                                        35440074
  130    IF (J .EQ. 30) GO TO 1000                                      35440075
         J = J + 1                                                      35440076
C     :::::::::: FORM SHIFT ::::::::::                                  35440077
         L1 = L + 1                                                     35440078
         S = DSQRT(E2(L))                                               35440079
         G = D(L)                                                       35440080
         P = (D(L1) - G) / (2.0D0 * S)                                  35440081
         R = DSQRT(P*P+1.0D0)                                           35440082
         D(L) = S / (P + DSIGN(R,P))                                    35440083
         H = G - D(L)                                                   35440084
C                                                                       35440085
         DO 140 I = L1, N                                               35440086
  140    D(I) = D(I) - H                                                35440087
C                                                                       35440088
         F = F + H                                                      35440089
C     :::::::::: RATIONAL QL TRANSFORMATION ::::::::::                  35440090
         G = D(M)                                                       35440091
         IF (G .EQ. 0.0D0) G = B                                        35440092
         H = G                                                          35440093
         S = 0.0D0                                                      35440094
         MML = M - L                                                    35440095
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::             35440096
         DO 200 II = 1, MML                                             35440097
            I = M - II                                                  35440098
            P = G * H                                                   35440099
            R = P + E2(I)                                               35440100
            E2(I+1) = S * R                                             35440101
            S = E2(I) / R                                               35440102
            D(I+1) = H + S * (H + D(I))                                 35440103
            G = D(I) - E2(I) / G                                        35440104
            IF (G .EQ. 0.0D0) G = B                                     35440105
            H = G * P / R                                               35440106
  200    CONTINUE                                                       35440107
C                                                                       35440108
         E2(L) = S * G                                                  35440109
         D(L) = H                                                       35440110
C     :::::::::: GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST :::::::::: 35440111
         IF (H .EQ. 0.0D0) GO TO 210                                    35440112
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210                      35440113
         E2(L) = H * E2(L)                                              35440114
         IF (E2(L) .NE. 0.0D0) GO TO 130                                35440115
  210    P = D(L) + F                                                   35440116
C     :::::::::: ORDER EIGENVALUES ::::::::::                           35440117
         IF (L .EQ. 1) GO TO 250                                        35440118
C     :::::::::: FOR I=L STEP -1 UNTIL 2 DO -- ::::::::::               35440119
         DO 230 II = 2, L                                               35440120
            I = L + 2 - II                                              35440121
            IF (P .GE. D(I-1)) GO TO 270                                35440122
            D(I) = D(I-1)                                               35440123
  230    CONTINUE                                                       35440124
C                                                                       35440125
  250    I = 1                                                          35440126
  270    D(I) = P                                                       35440127
  290 CONTINUE                                                          35440128
C                                                                       35440129
      GO TO 1001                                                        35440130
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN                      35440131
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::              35440132
 1000 IERR = L                                                          35440133
 1001 RETURN                                                            35440134
C     :::::::::: LAST CARD OF TQLRAT ::::::::::                         35440135
      END                                                               35440136
C                                                                       78440001
C     ------------------------------------------------------------------78440002
C                                                                       78440003
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    78440004
C                                                                       78440005
      INTEGER I,J,K,L,N,II,NM,JP1                                       78440006
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)                                  78440007
      REAL*8 F,G,H,HH,SCALE                                             78440008
      REAL*8 DSQRT,DABS,DSIGN                                           78440009
C                                                                       78440010
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    78440011
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   78440012
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   78440013
C                                                                       78440014
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              78440015
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               78440016
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            78440017
C                                                                       78440018
C     ON INPUT:                                                         78440019
C                                                                       78440020
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         78440021
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          78440022
C          DIMENSION STATEMENT;                                         78440023
C                                                                       78440024
C        N IS THE ORDER OF THE MATRIX;                                  78440025
C                                                                       78440026
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          78440027
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               78440028
C                                                                       78440029
C     ON OUTPUT:                                                        78440030
C                                                                       78440031
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;    78440032
C                                                                       78440033
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         78440034
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;      78440035
C                                                                       78440036
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                78440037
C          PRODUCED IN THE REDUCTION;                                   78440038
C                                                                       78440039
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            78440040
C                                                                       78440041
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        78440042
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         78440043
C                                                                       78440044
C     ------------------------------------------------------------------78440045
C                                                                       78440046
      DO 100 I = 1, N                                                   78440047
C                                                                       78440048
         DO 100 J = 1, I                                                78440049
            Z(I,J) = A(I,J)                                             78440050
  100 CONTINUE                                                          78440051
C                                                                       78440052
      IF (N .EQ. 1) GO TO 320                                           78440053
C     :::::::::: FOR I=N STEP -1 UNTIL 2 DO -- ::::::::::               78440054
      DO 300 II = 2, N                                                  78440055
         I = N + 2 - II                                                 78440056
         L = I - 1                                                      78440057
         H = 0.0D0                                                      78440058
         SCALE = 0.0D0                                                  78440059
         IF (L .LT. 2) GO TO 130                                        78440060
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::       78440061
         DO 120 K = 1, L                                                78440062
  120    SCALE = SCALE + DABS(Z(I,K))                                   78440063
C                                                                       78440064
         IF (SCALE .NE. 0.0D0) GO TO 140                                78440065
  130    E(I) = Z(I,L)                                                  78440066
         GO TO 290                                                      78440067
C                                                                       78440068
  140    DO 150 K = 1, L                                                78440069
            Z(I,K) = Z(I,K) / SCALE                                     78440070
            H = H + Z(I,K) * Z(I,K)                                     78440071
  150    CONTINUE                                                       78440072
C                                                                       78440073
         F = Z(I,L)                                                     78440074
         G = -DSIGN(DSQRT(H),F)                                         78440075
         E(I) = SCALE * G                                               78440076
         H = H - F * G                                                  78440077
         Z(I,L) = F - G                                                 78440078
         F = 0.0D0                                                      78440079
C                                                                       78440080
         DO 240 J = 1, L                                                78440081
            Z(J,I) = Z(I,J) / H                                         78440082
            G = 0.0D0                                                   78440083
C     :::::::::: FORM ELEMENT OF A*U ::::::::::                         78440084
            DO 180 K = 1, J                                             78440085
  180       G = G + Z(J,K) * Z(I,K)                                     78440086
C                                                                       78440087
            JP1 = J + 1                                                 78440088
            IF (L .LT. JP1) GO TO 220                                   78440089
C                                                                       78440090
            DO 200 K = JP1, L                                           78440091
  200       G = G + Z(K,J) * Z(I,K)                                     78440092
C     :::::::::: FORM ELEMENT OF P ::::::::::                           78440093
  220       E(J) = G / H                                                78440094
            F = F + E(J) * Z(I,J)                                       78440095
  240    CONTINUE                                                       78440096
C                                                                       78440097
         HH = F / (H + H)                                               78440098
C     :::::::::: FORM REDUCED A ::::::::::                              78440099
         DO 260 J = 1, L                                                78440100
            F = Z(I,J)                                                  78440101
            G = E(J) - HH * F                                           78440102
            E(J) = G                                                    78440103
C                                                                       78440104
            DO 260 K = 1, J                                             78440105
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)                  78440106
  260    CONTINUE                                                       78440107
C                                                                       78440108
  290    D(I) = H                                                       78440109
  300 CONTINUE                                                          78440110
C                                                                       78440111
  320 D(1) = 0.0D0                                                      78440112
      E(1) = 0.0D0                                                      78440113
C     :::::::::: ACCUMULATION OF TRANSFORMATION MATRICES ::::::::::     78440114
      DO 500 I = 1, N                                                   78440115
         L = I - 1                                                      78440116
         IF (D(I) .EQ. 0.0D0) GO TO 380                                 78440117
C                                                                       78440118
         DO 360 J = 1, L                                                78440119
            G = 0.0D0                                                   78440120
C                                                                       78440121
            DO 340 K = 1, L                                             78440122
  340       G = G + Z(I,K) * Z(K,J)                                     78440123
C                                                                       78440124
            DO 360 K = 1, L                                             78440125
               Z(K,J) = Z(K,J) - G * Z(K,I)                             78440126
  360    CONTINUE                                                       78440127
C                                                                       78440128
  380    D(I) = Z(I,I)                                                  78440129
         Z(I,I) = 1.0D0                                                 78440130
         IF (L .LT. 1) GO TO 500                                        78440131
C                                                                       78440132
         DO 400 J = 1, L                                                78440133
            Z(I,J) = 0.0D0                                              78440134
            Z(J,I) = 0.0D0                                              78440135
  400    CONTINUE                                                       78440136
C                                                                       78440137
  500 CONTINUE                                                          78440138
C                                                                       78440139
      RETURN                                                            78440140
C     :::::::::: LAST CARD OF TRED2 ::::::::::                          78440141
      END                                                               78440142
C                                                                       90220001
C     ------------------------------------------------------------------90220002
C                                                                       90220003
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  90220004
C                                                                       90220005
      implicit real*8 (a-h,o-z)                                         8d3s93
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR                             90220006
      REAL*8 D(N),E(N),Z(NM,N)                                          90220007
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                     90220008
      REAL*8 DSQRT,DABS,DSIGN                                           90220009
C                                                                       90220010
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     90220011
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     90220012
C     WILKINSON.                                                        90220013
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   90220014
C                                                                       90220015
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            90220016
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               90220017
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              90220018
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  90220019
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  90220020
C                                                                       90220021
C     ON INPUT:                                                         90220022
C                                                                       90220023
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         90220024
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          90220025
C          DIMENSION STATEMENT;                                         90220026
C                                                                       90220027
C        N IS THE ORDER OF THE MATRIX;                                  90220028
C                                                                       90220029
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          90220030
C                                                                       90220031
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        90220032
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;               90220033
C                                                                       90220034
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           90220035
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      90220036
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        90220037
C          THE IDENTITY MATRIX.                                         90220038
C                                                                       90220039
C      ON OUTPUT:                                                       90220040
C                                                                       90220041
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          90220042
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          90220043
C          UNORDERED FOR INDICES 1,2,...,IERR-1;                        90220044
C                                                                       90220045
C        E HAS BEEN DESTROYED;                                          90220046
C                                                                       90220047
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           90220048
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     90220049
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       90220050
C          EIGENVALUES;                                                 90220051
C                                                                       90220052
C        IERR IS SET TO                                                 90220053
C          ZERO       FOR NORMAL RETURN,                                90220054
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               90220055
C                     DETERMINED AFTER 30 ITERATIONS.                   90220056
C                                                                       90220057
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        90220058
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         90220059
C                                                                       90220060
C     ------------------------------------------------------------------90220061
C                                                                       90220062
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     90220063
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   90220064
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        90220065
C                ON S360 ::::::::::                                     90220066
c      DATA MACHEP/Z3410000000000000/                                    90220067
      factm=0.5d0
      machep=1d-12
  333 continue
       if(machep+1.d0.eq.1d0)go to 353
       machep=machep*factm
       go to 333
  353 continue
      machep=machep*2d0
C                                                                       90220068
      IERR = 0                                                          90220069
      IF (N .EQ. 1) GO TO 1001                                          90220070
C                                                                       90220071
      DO 100 I = 2, N                                                   90220072
  100 E(I-1) = E(I)                                                     90220073
C                                                                       90220074
      F = 0.0D0                                                         90220075
      B = 0.0D0                                                         90220076
      E(N) = 0.0D0                                                      90220077
C                                                                       90220078
      DO 240 L = 1, N                                                   90220079
         J = 0                                                          90220080
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))                         90220081
         IF (B .LT. H) B = H                                            90220082
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::         90220083
         DO 110 M = L, N                                                90220084
            IF (DABS(E(M)) .LE. B) GO TO 120                            90220085
C     :::::::::: E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               90220086
C                THROUGH THE BOTTOM OF THE LOOP ::::::::::              90220087
  110    CONTINUE                                                       90220088
C                                                                       90220089
  120    IF (M .EQ. L) GO TO 220                                        90220090
  130    IF (J .EQ. 30) GO TO 1000                                      90220091
         J = J + 1                                                      90220092
C     :::::::::: FORM SHIFT ::::::::::                                  90220093
         L1 = L + 1                                                     90220094
         G = D(L)                                                       90220095
         P = (D(L1) - G) / (2.0D0 * E(L))                               90220096
C         R = DSQRT(P*P+1.0D0)                                           90220097
         IF(ABS(P).LT.1.D18)THEN
          R=DSQRT(P*P+1.0D0)
         ELSE
          R=ABS(P)
         END IF
         D(L) = E(L) / (P + DSIGN(R,P))                                 90220098
         H = G - D(L)                                                   90220099
C                                                                       90220100
         DO 140 I = L1, N                                               90220101
  140    D(I) = D(I) - H                                                90220102
C                                                                       90220103
         F = F + H                                                      90220104
C     :::::::::: QL TRANSFORMATION ::::::::::                           90220105
         P = D(M)                                                       90220106
         C = 1.0D0                                                      90220107
         S = 0.0D0                                                      90220108
         MML = M - L                                                    90220109
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::             90220110
         DO 200 II = 1, MML                                             90220111
            I = M - II                                                  90220112
            G = C * E(I)                                                90220113
            H = C * P                                                   90220114
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150                      90220115
            C = E(I) / P                                                90220116
            R = DSQRT(C*C+1.0D0)                                        90220117
            E(I+1) = S * P * R                                          90220118
            S = C / R                                                   90220119
            C = 1.0D0 / R                                               90220120
            GO TO 160                                                   90220121
  150       C = P / E(I)                                                90220122
            R = DSQRT(C*C+1.0D0)                                        90220123
            E(I+1) = S * E(I) * R                                       90220124
            S = 1.0D0 / R                                               90220125
            C = C * S                                                   90220126
  160       P = C * D(I) - S * G                                        90220127
            D(I+1) = H + S * (C * G + S * D(I))                         90220128
C     :::::::::: FORM VECTOR ::::::::::                                 90220129
            DO 180 K = 1, N                                             90220130
               H = Z(K,I+1)                                             90220131
               Z(K,I+1) = S * Z(K,I) + C * H                            90220132
               Z(K,I) = C * Z(K,I) - S * H                              90220133
  180       CONTINUE                                                    90220134
C                                                                       90220135
  200    CONTINUE                                                       90220136
C                                                                       90220137
         E(L) = S * P                                                   90220138
         D(L) = C * P                                                   90220139
         IF (DABS(E(L)) .GT. B) GO TO 130                               90220140
  220    D(L) = D(L) + F                                                90220141
  240 CONTINUE                                                          90220142
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::          90220143
      DO 300 II = 2, N                                                  90220144
         I = II - 1                                                     90220145
         K = I                                                          90220146
         P = D(I)                                                       90220147
C                                                                       90220148
         DO 260 J = II, N                                               90220149
            IF (D(J) .GE. P) GO TO 260                                  90220150
            K = J                                                       90220151
            P = D(J)                                                    90220152
  260    CONTINUE                                                       90220153
C                                                                       90220154
         IF (K .EQ. I) GO TO 300                                        90220155
         D(K) = D(I)                                                    90220156
         D(I) = P                                                       90220157
C                                                                       90220158
         DO 280 J = 1, N                                                90220159
            P = Z(J,I)                                                  90220160
            Z(J,I) = Z(J,K)                                             90220161
            Z(J,K) = P                                                  90220162
  280    CONTINUE                                                       90220163
C                                                                       90220164
  300 CONTINUE                                                          90220165
C                                                                       90220166
      GO TO 1001                                                        90220167
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN                      90220168
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::              90220169
 1000 IERR = L                                                          90220170
 1001 RETURN                                                            90220171
C     :::::::::: LAST CARD OF TQL2 ::::::::::                           90220172
      END                                                               90220173
