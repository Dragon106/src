C***********************************************************************
      SUBROUTINE D01BCF(ITYPE,AA,BB,CC,DD,NPNTS,WEIGHT,ABSCIS,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9C REVISED. IER-370 (JUN 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14A REVISED. IER-677 (DEC 1989).
C     MARK 14B REVISED. IER-840 (MAR 1990).
C     SUBROUTINE FOR THE DETERMINATION OF GAUSSIAN QUADRATURE RULES
C     **************************************************************
C
C     INPUT PARAMETERS
C
C     ITYPE  INTEGER WHICH SPECIFIES THE RULE TYPE CHOSEN
C               WEIGHT W(X)          INTERVAL     RESTRICTIONS
C     0            1                    A,B          B.GT.A
C     1    (B-X)**C*(X-A)**D            A,B     B.GT.A,C,D.GT.-1
C     2   ABS(X-0.5*(A+B))**C           A,B     C.GT.-1,B.GT.A
C     3  ABS(X-A)**C*EXP(-B*X)          A,INF   C.GT.-1,B.GT.0
C     3  ABS(X-A)**C*EXP(-B*X)       -INF,A     C.GT.-1,B.LT.0
C     4 ABS(X-A)**C*EXP(-B*(X-A)**2) -INF,INF   C.GT.-1,B.GT.0
C     5  ABS(X-A)**C/ABS(X+B)**D      A,INF A.GT.-B,C.GT.-1,D.GT.C+1
C     5  ABS(X-A)**C/ABS(X+B)**D     -INF,A A.LT.-B,C.GT.-1,D.GT.C+1
C     ABS(ITYPE) MUST BE LESS THAN 6. IF ITYPE IS GIVEN LESS THAN
C     ZERO THEN THE ADJUSTED WEIGHTS ARE CALCULATED. IF NPNTS IS
C     ODD AND ITYPE EQUALS -2 OR -4 AND C IS NOT ZERO, THERE MAY BE
C     PROBLEMS.
C
C     AA     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     BB     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     CC     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     DD     REAL PARAMETER USED TO SPECIFY RULE TYPE.  SEE ITYPE.
C
C     IFAIL  NAG FAILURE PARAMETER.  SEE NAG DOCUMENTATION.
C
C     NPNTS  INTEGER THAT DETERMINES DIMENSION OF WEIGHT AND ABSCIS
C
C     OUTPUT PARAMETERS
C
C     WEIGHT  REAL ARRAY OF DIMENSION NPNTS WHICH CONTAINS
C     RULE WEIGHTS
C
C     ABSCIS  REAL ARRAY OF DIMENSION NPNTS WHICH CONTAINS
C     RULE ABSCISSAE
C
C     IFAIL INTEGER NAG FAILURE PARAMETER
C       IFAIL=0 FOR NORMAL EXIT
C       IFAIL=1 FOR FAILURE IN NAG ROUTINE F02AVF
C       IFAIL=2 FOR PARAMETER NPNTS OR ITYPE OUT OF RANGE
C       IFAIL=3 FOR PARAMETER AA OR BB OR CC OR DD OUT OF
C               ALLOWED RANGE
C       IFAIL=4 FOR OVERFLOW IN CALCULATION OF WEIGHTS
C       IFAIL=5 FOR UNDERFLOW IN CALCULATION OF WEIGHTS
C       IFAIL=6 FOR ITYPE=-2 OR -4, NPNTS ODD, C NOT ZERO
C
C     **************************************************************
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01BCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AA, BB, CC, DD
      INTEGER           IFAIL, ITYPE, NPNTS
C     .. Array Arguments ..
      DOUBLE PRECISION  ABSCIS(NPNTS), WEIGHT(NPNTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ABSPNC, B, BN, C, CN, CNO, D, FACN, FN, FOUR,
     *                  GAMMA, GAMMAB, GAMMB, HALF, ONE, PNA, PNB, PNC,
     *                  PONORM, PSQRD, REALMX, SMALL, SQRTCN, STORE,
     *                  TWNAPB, TWO, WTSUM, Y, ZERO
      INTEGER           IERROR, ISUB, J, MITYPE, N, NBUG, NFAC, NHALF
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S14AAF, X02AJF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          S14AAF, X02AJF, X02ALF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02AVF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MOD, LOG, EXP, DBLE, SQRT, INT
C     .. Data statements ..
      DATA              ZERO, ONE, TWO, FOUR/0.0D0, 1.0D0, 2.0D0, 4.0D0/
      DATA              HALF/0.5D0/
C     .. Executable Statements ..
C
C        INITIALISATION AND PARAMETER CHECKING
C
      SMALL = X02AJF()
      IF (NPNTS.LE.0) GO TO 780
      DO 20 J = 1, NPNTS
         ABSCIS(J) = ZERO
         WEIGHT(J) = ZERO
   20 CONTINUE
      MITYPE = ABS(ITYPE) + 1
      IF (MITYPE.GT.6) GO TO 780
      A = AA
      B = BB
      C = CC
      D = DD
      GO TO (40,60,100,120,140,160) MITYPE
   40 C = ZERO
      D = ZERO
   60 IF (C.LE.-ONE .OR. D.LE.-ONE) GO TO 800
      IF (B.LE.A) GO TO 800
      PONORM = (HALF*(B-A))**(C+D+ONE)
      IF (ITYPE.LT.0) PONORM = PONORM/(HALF*(B-A))**(C+D)
   80 IERROR = 1
      GAMMA = S14AAF(C+ONE,IERROR)
      IF (IERROR.GT.0) GO TO 800
      IERROR = 1
      GAMMB = S14AAF(D+ONE,IERROR)
      IF (IERROR.GT.0) GO TO 800
      IERROR = 1
      GAMMAB = S14AAF(C+D+TWO,IERROR)
      IF (IERROR.GT.0) GO TO 800
      PONORM = PONORM*TWO**(C+D+ONE)*GAMMA*GAMMB/GAMMAB
      ABSCIS(1) = (D-C)/(C+D+TWO)
      GO TO 180
  100 IF (C.LE.-ONE .OR. B.LE.A) GO TO 800
      PONORM = TWO*(HALF*(B-A))**(C+ONE)/(C+ONE)
      IF (ITYPE.LT.0) PONORM = PONORM/(HALF*(B-A))**C
      GO TO 180
  120 IF (C.LE.-ONE .OR. B.EQ.ZERO) GO TO 800
      IERROR = 1
      PONORM = S14AAF(C+ONE,IERROR)*EXP(-B*A)/ABS(B)**(C+ONE)
      IF (ITYPE.LT.0) PONORM = PONORM/EXP(-B*A)*ABS(B)**C
      IF (IERROR.GT.0) GO TO 800
      ABSCIS(1) = C + ONE
      GO TO 180
  140 IF (C.LE.-ONE .OR. B.LE.ZERO) GO TO 800
      IERROR = 1
      PONORM = S14AAF((C+ONE)/TWO,IERROR)/B**((C+ONE)/TWO)
      IF (ITYPE.LT.0) PONORM = PONORM*B**(C/TWO)
      IF (IERROR.GT.0) GO TO 800
      GO TO 180
  160 IF (A+B.EQ.ZERO) GO TO 800
      IF (C.LE.-ONE .OR. D.LE.C+ONE) GO TO 800
      D = D - C - TWO
      PONORM = ONE/(TWO**(C+D+ONE))/(ABS(A+B)**(D+ONE))
      IF (ITYPE.LT.0) PONORM = PONORM*(TWO**(C+D+TWO))*(ABS(A+B)
     *                         **(D+TWO))
      GO TO 80
C
C       COMPUTE DIAGONAL AND OFF-DIAGONAL OF SYMMETRIC TRI-DIAGONAL
C         MATRIX WHICH HAS ABSCISSAE AS EIGENVALUES
C
  180 IF (NPNTS.EQ.1) GO TO 320
      DO 300 N = 2, NPNTS
         FN = N - 1
         GO TO (200,200,220,240,260,200) MITYPE
  200    TWNAPB = FN + FN + C + D
         ABSCIS(N) = (D+C)*(D-C)/(TWNAPB*(TWNAPB+TWO))
         CN = FOUR*(FN+C)*(FN+D)*FN/(TWNAPB**2*(TWNAPB+ONE))
         IF (N.GT.2) CN = CN*((C+D+FN)/(TWNAPB-ONE))
         GO TO 280
  220    ABSCIS(N) = ZERO
         CN = (FN+C*MOD(FN,TWO))**2/((FN+FN+C)**2-ONE)
         GO TO 280
  240    ABSCIS(N) = C + FN + FN + ONE
         CN = FN*(C+FN)
         GO TO 280
  260    ABSCIS(N) = ZERO
         CN = (FN+C*MOD(FN,TWO))/TWO
  280    WEIGHT(N) = SQRT(CN)
  300 CONTINUE
C
C        USE NAG ROUTINE TO FIND EIGENVALUES WHICH ARE ABSCISSAE
C
  320 IERROR = 1
      CALL F02AVF(NPNTS,X02AJF(),ABSCIS,WEIGHT,IERROR)
      IF (IERROR.GT.0) GO TO 760
C
C        LOOP TO DETERMINE WEIGHTS
C             EVALUATE EACH ORTHONORMAL POLYNOMIAL OF DEGREE
C         LESS THAN NPNTS AT ABSCIS(J) AND SUM SQUARES OF
C         RESULTS TO DETERMINE WEIGHT(J)
      IERROR = 0
      REALMX = X02ALF()
      DO 700 J = 1, NPNTS
C
C        INITIALISE INNER LOOP AND SCALE WEIGHT(J) AND ABSCIS(J)
C        DIVIDE EXPONENTIAL TERMS INTO FACTORS THAT DON'T UNDERFLOW
C
         WEIGHT(J) = ZERO
         Y = ABSCIS(J)
         PNA = ZERO
         CNO = ZERO
         NFAC = 0
         PNB = ONE/SQRT(PONORM)
         GO TO (340,340,360,400,420,440) MITYPE
  340    ABSCIS(J) = Y*(HALF*(B-A)) + (HALF*(A+B))
         IF (ITYPE.GT.0) GO TO 460
         PNB = PNB*(ONE-Y)**(C*HALF)*(ONE+Y)**(D*HALF)
         GO TO 460
  360    ABSCIS(J) = Y*(HALF*(B-A)) + (HALF*(A+B))
         IF (ITYPE.GT.0 .OR. C.EQ.ZERO) GO TO 460
         IF (Y.EQ.ZERO .AND. C.GT.ZERO) GO TO 660
         IF (C.GT.ZERO) GO TO 380
         IF (PONORM.GE.ONE) GO TO 380
         IF (ABS(Y).LE.(ONE/(REALMX*PONORM))**(-ONE/C)) GO TO 680
  380    PNB = PNB*ABS(Y)**(C*HALF)
         GO TO 460
  400    ABSCIS(J) = Y/B + A
         IF (ITYPE.GT.0) GO TO 460
         PNB = PNB*Y**(C*HALF)
         NFAC = INT(Y/LOG(HALF*REALMX)) + 1
         FACN = EXP(-HALF*Y/DBLE(NFAC))
         GO TO 460
  420    ABSCIS(J) = Y/SQRT(B) + A
         IF (ITYPE.GT.0) GO TO 460
         NFAC = INT(Y*Y/LOG(HALF*REALMX)) + 1
         FACN = EXP(-HALF*Y*Y/DBLE(NFAC))
         IF (C.EQ.ZERO) GO TO 460
         IF (Y.EQ.ZERO .AND. C.GT.ZERO) GO TO 660
         IF (Y.EQ.ZERO .AND. C.LT.ZERO) GO TO 680
         PNB = PNB*ABS(Y)**(C*HALF)
         GO TO 460
  440    ABSCIS(J) = TWO*(A+B)/(Y+ONE) - B
         IF (ITYPE.GT.0) GO TO 460
         PNB = PNB*(ONE-Y)**(C*HALF)*(ONE+Y)**(HALF*(D+TWO))
  460    WTSUM = PNB*PNB
         IF (NPNTS.EQ.1) GO TO 640
C
C          LOOP TO EVALUATE ORTHONORMAL POLYNOMIALS USING THREE
C           TERM RECURRENCE RELATION.
C
         DO 620 N = 2, NPNTS
            FN = N - 1
            GO TO (480,480,500,520,540,480) MITYPE
  480       TWNAPB = FN + FN + C + D
            BN = (D-C)/TWNAPB
            IF (N.GT.2) BN = BN*(C+D)/(TWNAPB-TWO)
            CN = FOUR*FN*(C+FN)*(D+FN)/(TWNAPB**2*(TWNAPB+ONE))
            IF (N.GT.2) CN = CN*((C+D+FN)/(TWNAPB-ONE))
            GO TO 560
  500       BN = ZERO
            CN = (FN+C*MOD(FN,TWO))**2/((FN+FN+C)**2-ONE)
            GO TO 560
  520       BN = C + FN + FN - ONE
            CN = FN*(FN+C)
            GO TO 560
  540       BN = ZERO
            CN = (FN+C*MOD(FN,TWO))/TWO
  560       SQRTCN = SQRT(CN)
            PNC = ((Y-BN)*PNB-CNO*PNA)/SQRTCN
            CNO = SQRTCN
            ABSPNC = ABS(PNC)
            IF (ABSPNC.LE.ONE) GO TO 580
            IF (ABSPNC.LE.REALMX/ABSPNC) GO TO 580
            IF (ITYPE.GT.0) GO TO 680
            IF (NFAC.LE.0) GO TO 680
            PNB = PNB*FACN
            PNC = PNC*FACN
            WTSUM = WTSUM*FACN*FACN
            NFAC = NFAC - 1
  580       PSQRD = PNC*PNC
            IF (WTSUM.LE.REALMX-PSQRD) GO TO 600
            IF (ITYPE.GT.0) GO TO 680
            IF (NFAC.LE.0) GO TO 680
            PNB = PNB*FACN
            PNC = PNC*FACN
            WTSUM = WTSUM*FACN*FACN
            PSQRD = PSQRD*FACN*FACN
            NFAC = NFAC - 1
  600       WTSUM = WTSUM + PSQRD
            PNA = PNB
            PNB = PNC
  620    CONTINUE
C
C          END LOOP FOR POLYNOMIAL EVALUATION
C
C Richard Brankin - NAG, Oxford - 26th July 1989
C replaced the following line ....
C
C  640    IF (NFAC.GT.0) WTSUM = WTSUM*FACN**(2*NFAC)
C
C so as not to get needless underflow to zero when powering up FACN
C for 0.0 < FACN << 1.0. The error was brought to light in a VAX
C double precision implementation when a user tried to compute modified
C Laguerre weights (ITYPE = -3) for more than 25 abscissae (N > 25).
C As a result, before the assignment in the above line
C   WTSUM = O(1.0e+38), FACN = O(1.0e-10), NFAC = 2
C WTSUM was assigned a value of 0.0 since O(1.0e-10)**4 underflows
C although WTSUM should have been assigned O(1.0e+2). This correction
C also applies for other values of ITYPE.
C
  640    IF (NFAC.GT.0) THEN
            DO 650 NBUG = 1, 2*NFAC
               WTSUM = WTSUM*FACN
  650       CONTINUE
         END IF
C
C End of correction
C
         IF (WTSUM.EQ.ZERO) GO TO 660
         WEIGHT(J) = ONE/WTSUM
         GO TO 700
  660    IERROR = 4
         WEIGHT(J) = REALMX
         GO TO 700
  680    IERROR = 5
  700 CONTINUE
C
C        END LOOP FOR WEIGHTS
C
C        REVERSE RATIONAL OR LAGUERRE POINTS
C
      IF ((MITYPE.NE.6 .OR. A+B.LT.ZERO)
     *    .AND. (MITYPE.NE.4 .OR. B.GT.ZERO)) GO TO 740
      NHALF = NPNTS/2
      IF (NHALF.LE.1) GO TO 740
      DO 720 J = 1, NHALF
         ISUB = NPNTS + 1 - J
         STORE = ABSCIS(J)
         ABSCIS(J) = ABSCIS(ISUB)
         ABSCIS(ISUB) = STORE
         STORE = WEIGHT(J)
         WEIGHT(J) = WEIGHT(ISUB)
         WEIGHT(ISUB) = STORE
  720 CONTINUE
C
C        ASSIGNMENT OF IFAIL PARAMETER
C
  740 IF ((ITYPE.EQ.-2 .OR. ITYPE.EQ.-4) .AND. MOD(NPNTS,2)
     *    .EQ.1 .AND. C.NE.ZERO) IERROR = 6
      GO TO 820
  760 IERROR = 1
      GO TO 820
  780 IERROR = 2
      GO TO 820
  800 IERROR = 3
  820 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F02AVF(N,ACHEPS,D,E,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 9 REVISED. IER-326 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     TQL1
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A TRIDIAGONAL
C     MATRIX,
C     T, GIVEN WITH ITS DIAGONAL ELEMENTS IN THE ARRAY D(N) AND
C     ITS SUBDIAGONAL ELEMENTS IN THE LAST N - 1 STORES OF THE
C     ARRAY E(N), USING QL TRANSFORMATIONS. THE EIGENVALUES ARE
C     OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE ARRAY D IN
C     ASCENDING ORDER. THE SUBROUTINE WILL FAIL IF ALL
C     EIGENVALUES TAKE MORE THAN 30*N ITERATIONS.
C     1ST APRIL 1972
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AVF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACHEPS
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), E(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, F, G, H, P, R, S
      INTEGER           I, I1, II, ISAVE, J, L, M, M1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I = 2, N
         E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      J = 30*N
      DO 340 L = 1, N
         H = ACHEPS*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C        LOOK FOR SMALL SUB DIAGONAL ELEMENT
         DO 60 M = L, N
            IF (ABS(E(M)).LE.B) GO TO 80
   60    CONTINUE
   80    IF (M.EQ.L) GO TO 260
  100    IF (J.LE.0) GO TO 360
         J = J - 1
C        FORM SHIFT
         G = D(L)
         H = D(L+1) - G
         IF (ABS(H).GE.ABS(E(L))) GO TO 120
         P = H*0.5D0/E(L)
         R = SQRT(P*P+1.0D0)
         H = P + R
         IF (P.LT.0.0D0) H = P - R
         D(L) = E(L)/H
         GO TO 140
  120    P = 2.0D0*E(L)/H
         R = SQRT(P*P+1.0D0)
         D(L) = E(L)*P/(1.0D0+R)
  140    H = G - D(L)
         I1 = L + 1
         IF (I1.GT.N) GO TO 180
         DO 160 I = I1, N
            D(I) = D(I) - H
  160    CONTINUE
  180    F = F + H
C        QL TRANSFORMATION
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         M1 = M - 1
         DO 240 II = L, M1
            I = M1 - II + L
            G = C*E(I)
            H = C*P
            IF (ABS(P).LT.ABS(E(I))) GO TO 200
            C = E(I)/P
            R = SQRT(C*C+1.0D0)
            E(I+1) = S*P*R
            S = C/R
            C = 1.0D0/R
            GO TO 220
  200       C = P/E(I)
            R = SQRT(C*C+1.0D0)
            E(I+1) = S*E(I)*R
            S = 1.0D0/R
            C = C/R
  220       P = C*D(I) - S*G
            D(I+1) = H + S*(C*G+S*D(I))
  240    CONTINUE
         E(L) = S*P
         D(L) = C*P
         IF (ABS(E(L)).GT.B) GO TO 100
  260    P = D(L) + F
C        ORDER EIGENVALUE
         IF (L.EQ.1) GO TO 300
         DO 280 II = 2, L
            I = L - II + 2
            IF (P.GE.D(I-1)) GO TO 320
            D(I) = D(I-1)
  280    CONTINUE
  300    I = 1
  320    D(I) = P
  340 CONTINUE
      IFAIL = 0
      RETURN
  360 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*1             REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
      DOUBLE PRECISION FUNCTION S14AAF(X,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 7C REVISED IER-184 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     GAMMA FUNCTION
C
C     **************************************************************
C
C     TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
C     ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING  CDD ,
C     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
C     DIGITS REPRESENTED BY THE MACHINE
C     DELETE THE ILLEGAL DUMMY STATEMENTS OF THE FORM
C     * EXPANSION (NNNN) *
C
C     ALSO INSERT APPROPRIATE DATA STATEMENTS TO DEFINE CONSTANTS
C     WHICH DEPEND ON THE RANGE OF NUMBERS REPRESENTED BY THE
C     MACHINE, RATHER THAN THE PRECISION (SUITABLE STATEMENTS FOR
C     SOME MACHINES ARE CONTAINED IN COMMENTS BEGINNING CRD WHERE
C     D IS A DIGIT WHICH SIMPLY DISTINGUISHES A GROUP OF MACHINES).
C     DELETE THE ILLEGAL DUMMY DATA STATEMENTS WITH VALUES WRITTEN
C     *VALUE*
C
C     **************************************************************
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S14AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 G, GBIG, T, XBIG, XMINV, XSMALL,
     *                                 Y
      INTEGER                          I, M
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SIGN, DBLE
C     .. Data statements ..
C08   DATA XSMALL/1.0D-8/
C09   DATA XSMALL/3.0D-9/
C12   DATA XSMALL/1.0D-12/
      DATA XSMALL/3.0D-15/
C17   DATA XSMALL/1.0D-17/
C19   DATA XSMALL/1.7D-18/
C
      DATA XBIG,GBIG,XMINV/ 1.70D+2,4.3D+304,2.23D-308 /
C     XBIG = LARGEST X SUCH THAT  GAMMA(X) .LT. MAXREAL
C                            AND  1.0/GAMMA(X+1.0) .GT. MINREAL
C             (ROUNDED DOWN TO AN INTEGER)
C     GBIG = GAMMA(XBIG)
C     XMINV = MAX(1.0/MAXREAL,MINREAL)  (ROUNDED UP)
C     FOR IEEE SINGLE PRECISION
CR0   DATA XBIG,GBIG,XMINV /33.0E0,2.6E+35,1.2E-38/
C     FOR IBM 360/370 AND SIMILAR MACHINES
CR1   DATA XBIG,GBIG,XMINV /57.0D0,7.1D+74,1.4D-76/
C     FOR DEC-10, HONEYWELL, UNIVAC 1100 (S.P.)
CR2   DATA XBIG,GBIG,XMINV /34.0D0,8.7D+36,5.9D-39/
C     FOR ICL 1900
CR3   DATA XBIG,GBIG,XMINV /58.0D0,4.0D+76,1.8D-77/
C     FOR CDC 7600/CYBER
CR4   DATA XBIG,GBIG,XMINV /164.0D0,2.0D+291,3.2D-294/
C     FOR UNIVAC 1100 (D.P.)
CR5   DATA XBIG,GBIG,XMINV /171.0D0,7.3D+306,1.2D-308/
C     FOR IEEE DOUBLE PRECISION
CR7   DATA XBIG,GBIG,XMINV /170.0D0,4.3D+304,2.3D-308/
C     .. Executable Statements ..
C
C     ERROR 1 AND 2 TEST
      T = ABS(X)
      IF (T.GT.XBIG) GO TO 160
C     SMALL RANGE TEST
      IF (T.LE.XSMALL) GO TO 140
C     MAIN RANGE REDUCTION
      M = X
      IF (X.LT.0.0D0) GO TO 80
      T = X - DBLE(M)
      M = M - 1
      G = 1.0D0
      IF (M) 20, 120, 40
   20 G = G/X
      GO TO 120
   40 DO 60 I = 1, M
         G = (X-DBLE(I))*G
   60 CONTINUE
      GO TO 120
   80 T = X - DBLE(M-1)
C     ERROR 4 TEST
      IF (T.EQ.1.0D0) GO TO 220
      M = 1 - M
      G = X
      DO 100 I = 1, M
         G = (DBLE(I)+X)*G
  100 CONTINUE
      G = 1.0D0/G
  120 T = 2.0D0*T - 1.0D0
C
C      * EXPANSION (0026) *
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 08E.09
C08   Y = ((((((((((((+1.88278283D-6)*T-5.48272091D-6)
C08  *    *T+1.03144033D-5)*T-3.13088821D-5)*T+1.01593694D-4)
C08  *    *T-2.98340924D-4)*T+9.15547391D-4)*T-2.42216251D-3)
C08  *    *T+9.04037536D-3)*T-1.34119055D-2)*T+1.03703361D-1)
C08  *    *T+1.61692007D-2)*T + 8.86226925D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 09E.10
C09   Y = (((((((((((((-6.463247484D-7)*T+1.882782826D-6)
C09  *    *T-3.382165478D-6)*T+1.031440334D-5)*T-3.393457634D-5)
C09  *    *T+1.015936944D-4)*T-2.967655076D-4)*T+9.155473906D-4)
C09  *    *T-2.422622002D-3)*T+9.040375355D-3)*T-1.341184808D-2)
C09  *    *T+1.037033609D-1)*T+1.616919866D-2)*T + 8.862269255D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 12E.13
C12   Y = (((((((((((((((-7.613347676160D-8)*T+2.218377726362D-7)
C12  *    *T-3.608242105549D-7)*T+1.106350622249D-6)
C12  *    *T-3.810416284805D-6)*T+1.138199762073D-5)
C12  *    *T-3.360744031186D-5)*T+1.008657892262D-4)
C12  *    *T-2.968993359366D-4)*T+9.158021574033D-4)
C12  *    *T-2.422593898516D-3)*T+9.040332894085D-3)
C12  *    *T-1.341185067782D-2)*T+1.037033635205D-1)
C12  *    *T+1.616919872669D-2)*T + 8.862269254520D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 15E.16
C15   Y = (((((((((((((((-1.243191705600000D-10
C15  *    *T+3.622882508800000D-10)*T-4.030909644800000D-10)
C15  *    *T+1.265236705280000D-9)*T-5.419466096640000D-9)
C15  *    *T+1.613133578240000D-8)*T-4.620920340480000D-8)
C15  *    *T+1.387603440435200D-7)*T-4.179652784537600D-7)
C15  *    *T+1.253148247777280D-6)*T-3.754930502328320D-6)
C15  *    *T+1.125234962812416D-5)*T-3.363759801664768D-5)
C15  *    *T+1.009281733953869D-4)*T-2.968901194293069D-4)
C15  *    *T+9.157859942174304D-4)*T-2.422595384546340D-3
C15   Y = ((((Y*T+9.040334940477911D-3)*T-1.341185057058971D-2)
C15  *    *T+1.037033634220705D-1)*T+1.616919872444243D-2)*T +
C15  *     8.862269254527580D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 17E.18
      Y = (((((((((((((((-1.46381209600000000D-11
     *    *T+4.26560716800000000D-11)*T-4.01499750400000000D-11)
     *    *T+1.27679856640000000D-10)*T-6.13513953280000000D-10)
     *    *T+1.82243164160000000D-9)*T-5.11961333760000000D-9)
     *    *T+1.53835215257600000D-8)*T-4.64774927155200000D-8)
     *    *T+1.39383522590720000D-7)*T-4.17808776355840000D-7)
     *    *T+1.25281466396672000D-6)*T-3.75499034136576000D-6)
     *    *T+1.12524642975590400D-5)*T-3.36375833240268800D-5)
     *    *T+1.00928148823365120D-4)*T-2.96890121633200000D-4
      Y = ((((((Y*T+9.15785997288933120D-4)*T-2.42259538436268176D-3)
     *    *T+9.04033494028101968D-3)*T-1.34118505705967765D-2)
     *    *T+1.03703363422075456D-1)*T+1.61691987244425092D-2)*T +
     *     8.86226925452758013D-1
C
C     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 19E.20
C19   Y = (((((((((((((((+6.7108864000000000000D-13
C19  *    *T-1.6777216000000000000D-12)*T+6.7108864000000000000D-13)
C19  *    *T-4.1523609600000000000D-12)*T+2.4998051840000000000D-11)
C19  *    *T-6.8985815040000000000D-11)*T+1.8595971072000000000D-10)
C19  *    *T-5.6763875328000000000D-10)*T+1.7255563264000000000D-9)
C19  *    *T-5.1663077376000000000D-9)*T+1.5481318277120000000D-8)
C19  *    *T-4.6445740523520000000D-8)*T+1.3931958370304000000D-7)
C19  *    *T-4.1782339907584000000D-7)*T+1.2528422549504000000D-6)
C19  *    *T-3.7549858152857600000D-6)*T+1.1252456510305280000D-5
C19   Y = (((((((((Y*T-3.3637584239226880000D-5)
C19  *    *T+1.0092815021080832000D-4)
C19  *    *T-2.9689012151880000000D-4)*T+9.1578599714350784000D-4)
C19  *    *T-2.4225953843706897600D-3)*T+9.0403349402888779200D-3)
C19  *    *T-1.3411850570596516480D-2)*T+1.0370336342207529018D-1)
C19  *    *T+1.6169198724442506740D-2)*T + 8.8622692545275801366D-1
C
      S14AAF = Y*G
      IFAIL = 0
      GO TO 240
C
C     ERROR 3 TEST
  140 IF (T.LT.XMINV) GO TO 200
      S14AAF = 1.0D0/X
      IFAIL = 0
      GO TO 240
C
C     ERROR EXITS
  160 IF (X.LT.0.0D0) GO TO 180
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      S14AAF = GBIG
      GO TO 240
C
  180 IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      S14AAF = 0.0D0
      GO TO 240
C
  200 IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
      T = X
      IF (X.EQ.0.0D0) T = 1.0D0
      S14AAF = SIGN(1.0D0/XMINV,T)
      GO TO 240
C
  220 IFAIL = P01ABF(IFAIL,4,SRNAME,0,P01REC)
      S14AAF = GBIG
C
  240 RETURN
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /Z'3CA0000000000001' /
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02ALF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1 - B**(-P)) * B**EMAX  (THE LARGEST POSITIVE MODEL
C     NUMBER)
C
      DOUBLE PRECISION X02CON
      DATA X02CON /Z'7FEFFFFFFFFFFFFF' /
C     .. Executable Statements ..
      X02ALF = X02CON
      RETURN
      END
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
C***********************************************************************
