       subroutine permutation(nperm,nrow,ncol,ngeno,model,
     +                        genotRate,data1,Gstat)
       
       implicit none
       integer nperm,nrow,ncol,ngeno,genotRate,model
       integer data1(nrow,ncol)
       double precision Gstat(ncol-1,nperm),GstatAux(ncol-1)
       integer i,j,k
       integer p(nrow),aux(nrow)
       
       do i=1,nrow
        aux(i)=data1(i,1) 
       end do    
 
c initiate the n permutations (nperm)
       do i=1,nperm
c permutation of case-control label
        call rperm(p,nrow) 

        do j=1,nrow
         data1(j,1)=aux(p(j)) 
        end do
c computing G statistic
        call WGassociation(nrow,ncol,ngeno,model,genotRate,
     +          data1,GstatAux)
        do k=1,ncol-1
         Gstat(k,i)=GstatAux(k)
        end do 

       end do

       end subroutine





       subroutine WGassociation(nrow,ncol,ngeno,model,genotRate,
     +                         data1,Gstat)

       implicit none
       integer nrow,ncol,ngeno,genotRate,ipred,i,j,k
       integer data1(nrow,ncol),nn(2,ngeno),model
       integer caco(nrow), snp(nrow)
       double precision Gstat(ncol-1)
       real x(nrow), y(nrow),ymin,ymax,control, r2
       integer nomis
    
       do k=1,nrow
        caco(k)=data1(k,1)
       end do 

       do i=2,ncol
        if (model.ne.5) then
         do k=1,nrow
          snp(k)=data1(k,i)
         end do 
         call table(nrow,ngeno,caco,snp,model,genotRate,
     +            nn,ipred)
         if (ipred.eq.1) then
          call G(ngeno,nn,Gstat(i-1))
         else
          Gstat(i-1)=-2         
         end if

c additive
        else if (model.eq.5) then
         nomis=0
         do j=1,nrow
           if (data1(j,i).ne.0) then
            x(nomis+1)=real(data1(j,1))
            y(nomis+1)=real(data1(j,i))
            nomis=nomis+1
           end if  
          end do

         control=(real(nomis)/real(nrow))*100.0
         call USMNMX(y,nomis,1,ymin,ymax)

         if (control.lt.genotRate) then
          Gstat(i-1)=-2
         else if (ymin.eq.ymax) then
          Gstat(i-1)=-1
         else
          
c VM
          Gstat(i-1)=(nomis-1)*r2(x,y,nomis)
         end if
        end if  

       end do 
      
       end subroutine

    



      subroutine table(n,geno,vec1,vec2,model,genotRate,nn,ipred)
      implicit none

c
c   IN:
c   n -- length of vectors (number of individuals)
c   geno -- 2 or 3 indicating the number of genotypes (to collapse dominant, recessive, etc.)
c   vec1 -- 1 y 0
c   vec2 -- 1,2 or 1,2,3
c   genotRate -- percentage of desired genotyping 
c
c   OUT:
c   nn -- 2x3 1:case  2:control  
c   ipred -- controls wheter the genotyping is less than genotRate
c
   
      integer n,i,j,k,geno,genotRate,ipred,nOK
      integer vec1(n),vec2(n)
      integer nn(2,geno),model
      real control

      do j=1,2
       do k=1,geno
        nn(j,k)=0
       end do
      end do          


c  Codominant
      if (model.eq.1) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if (vec2(i).eq.1) then
            nn(1,1)=nn(1,1)+1
         else if (vec2(i).eq.2) then
            nn(1,2)=nn(1,2)+1
         else if (vec2(i).eq.3) then
            nn(1,3)=nn(1,3)+1 
         end if
        else
         if (vec2(i).eq.1) then
            nn(2,1)=nn(2,1)+1
         else if (vec2(i).eq.2) then
            nn(2,2)=nn(2,2)+1
         else if (vec2(i).eq.3) then
            nn(2,3)=nn(2,3)+1 
         end if
        end if        
       end do   

c  Dominant
      else if (model.eq.2) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if (vec2(i).eq.1) then
            nn(1,1)=nn(1,1)+1
         else if ((vec2(i).eq.2).or.(vec2(i).eq.3)) then
            nn(1,2)=nn(1,2)+1
         end if
        else
         if (vec2(i).eq.1) then
            nn(2,1)=nn(2,1)+1
         else if ((vec2(i).eq.2).or.(vec2(i).eq.3)) then
            nn(2,2)=nn(2,2)+1
         end if
        end if        
       end do

c  Recessive
      else if (model.eq.3) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if ((vec2(i).eq.1).or.(vec2(i).eq.2)) then
            nn(1,1)=nn(1,1)+1
         else if (vec2(i).eq.3) then
            nn(1,2)=nn(1,2)+1
         end if
        else
         if ((vec2(i).eq.1).or.(vec2(i).eq.2)) then
            nn(2,1)=nn(2,1)+1
         else if (vec2(i).eq.3) then
            nn(2,2)=nn(2,2)+1
         end if
        end if        
       end do

c  Overdominant
      else if (model.eq.4) then
       do i=1,n
        if (vec1(i).eq.1) then  
         if (vec2(i).eq.2) then
            nn(1,1)=nn(1,1)+1
         else if ((vec2(i).eq.1).or.(vec2(i).eq.3)) then
            nn(1,2)=nn(1,2)+1
         end if
        else
         if (vec2(i).eq.2) then
            nn(2,1)=nn(2,1)+1
         else if ((vec2(i).eq.1).or.(vec2(i).eq.3)) then
            nn(2,2)=nn(2,2)+1
         end if
        end if        
       end do

      end if



      nOK=0
      do j=1,2
       do k=1,geno
        nOK=nOK+nn(j,k)
       end do
      end do          
   
      ipred=1
      control=(real(nOK)/real(n))*100.0
      if (control.lt.genotRate) then
        ipred=0
      end if
  
      end subroutine table 


      subroutine G(geno,tt,Gstat)
      implicit none
      integer geno
      integer tt(2,geno),i,j
      double precision r(2),s(geno),N,Gstat,Gij,obs,esp

      do i=1,2
       if (geno.eq.3) then 
        r(i)=tt(i,1)+tt(i,2)+tt(i,3)
       else if (geno.eq.2) then
        r(i)=tt(i,1)+tt(i,2)
       end if 
      end do
      N=r(1)+r(2)

      do i=1,geno
       s(i)=tt(1,i)+tt(2,i)
      end do


     
      if ((geno.eq.3).and.(((s(1).eq.0).and.(s(2).eq.0)).or.
     +   ((s(1).eq.0).and.(s(3).eq.0)).or.
     +   ((s(2).eq.0).and.(s(3).eq.0)))) then

         Gstat=-1.0

      else if ((geno.eq.2).and.((s(1).eq.0).or.(s(2).eq.0).or.
     +   (s(3).eq.0))) then
    
         Gstat=-1.0

      else
       Gstat=0.0
       do i=1,2
        do j=1,geno
         if (tt(i,j).gt.0) then
          esp=(r(i)*s(j))/N
          obs=tt(i,j)
          Gij=obs*log(obs/esp)
          Gstat=Gstat+Gij  
         end if
         end do
       end do   
       Gstat=2.0*Gstat  
      end if  

      end subroutine G 




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      REAL FUNCTION r2(x,y,n)

      INTEGER n
      REAL x(n),y(n)

      INTEGER i
      REAL sxx,sxy,syy,sx,sy, nn
      nn=real(n)

      sx=0.
      sy=0.
      sxx=0.
      syy=0.
      sxy=0.
      do 12 i=1,n
        sx=sx+x(i)
        sy=sy+y(i)
        sxx=sxx+x(i)*x(i)
        syy=syy+y(i)*y(i)
        sxy=sxy+x(i)*y(i)
12    continue

      r2=((sxy-sx*sy/nn)**2)/(sxx-sx*sx/nn)/(syy-sy*sy/nn)

      return
      END




      SUBROUTINE RPERM(P,N)

c---Changed RAND by RAND2 to allow g77 (JRG January 2007)

C===Generate a random permutation, P, of the first N integers.
C   (equivalent to sampling WITHOUT REPLACEMENT).
C   Adaptation of Knuth Volume 2, Algorithm 3.4.2P.
      INTEGER N,P(N), K,J,I,IPJ,ITEMP,M
      REAL U(100)
      DO 1 I=1,N
1     P(I)=I
C---Generate up to 100 U(0,1) numbers at a time.
      DO 3 I=1,N,100
        M=MIN(N-I+1,100)
        CALL RAND2(U,M)
        DO 2 J=1,M
          IPJ=I+J-1
          K=INT(U(J)*(N-IPJ+1))+IPJ
          ITEMP=P(IPJ)
          P(IPJ)=P(K)
2       P(K)=ITEMP
3     CONTINUE
      RETURN
      END
      SUBROUTINE IRAND(S,N,LOW,HI)
C===Generate a random integer sequence: S(1),S(2), ... ,S(N)
C   such that each element is in the closed interval <LOW,HI> and
C   sampled WITH REPLACEMENT.                            HDK, JUNE 1971.
      INTEGER N,S(N),LOW,HI,IX,I
      REAL U(1)
      DOUBLE PRECISION X
      DO 1 I=1,N
        CALL RAND2(U,1)
C---Use DP arithmetic to effect a more precise transformation.
        X=DBLE((HI+1)-LOW)*U(1) + DBLE(LOW)
        IX=X
        IF(X.LT.0 .AND. IX.NE.X) IX=X-1.D0
        S(I)=IX
1     CONTINUE
      RETURN
      END
      BLOCK DATA
C=======================================================================
C  Portable pseudo-random integer generator, especially for
C  microcomputers with a limitation of 16 bit integers. Translated from
C  original Pascal version(1) to Fortran 77 by H. D. Knoble, PSU.
C
C   The supporting paper is:
C   (1) B. Wichmann & D. Hill, "Building a Random-Number Generator",
C             BYTE, March, 1987, 12(3):127-128.
C
C   Also see the following related works:
C   (2) Wichmann, B.A. and Hill, I.D., "An Efficient and Portable",
C             Pseudo-random Number Generator", Applied Statistics,
C             Volume 31, 1982, pages 188-190.
C   (3) Haas, Alexander, "The Multiple Prime Random Number Generator",
C             ACM Transactions on Mathematical Software; December,
C             1987, 13(4):368-381.
C   (4) L'Ecuyer, Pierre, "Efficient and Portable Combined Random Number
C             Generators", Communications of the ACM; June, 1988,
C             31(6):742-749,774.
C
C Use...
C      CALL RAND2(U,N)
C          To generate a sequence, U, of N Uniform(0,1) numbers.
C          Cycle length is ((30269-1)*(30307-1)*(30323-1))/4  or
C          6953607871644  > 6.95E+12.
C
C     To access the SEED vector in the calling program use statements:
C     INTEGER SEED(3)
C     COMMON/RANDOM/SEED
C
C  The common variable SEED is the array of three current seeds.
      INTEGER SEED(3)
      COMMON/RANDOM/SEED
      DATA SEED(1),SEED(2),SEED(3)/1,10000,3000/
      END
C=======================================================================
      SUBROUTINE RAND2(U,N)
      INTEGER N,X,Y,Z
      REAL U(N),V
      COMMON/RANDOM/X,Y,Z
      IF(N.LE.0) RETURN
      DO 1 I=1,N
        X=171*MOD(X,177)-2*(X/177)
        IF(X.LT.0) X=X+30269
        Y=172*MOD(Y,176)-35*(Y/176)
        IF(Y.LT.0) Y=Y+30307
        Z=170*MOD(Z,178)-63*(Z/178)
        IF(Z.LT.0) Z=Z+30323
        V=X/30269.0 + Y/30308.0 + Z/30323.0
1       U(I)=V-INT(V)
      RETURN
      END


CUSMNMX
C   IMSL ROUTINE NAME   - USMNMX                                        USMN0010
C                                                                       USMN0020
C-----------------------------------------------------------------------USMN0030
C                                                                       USMN0040
C   COMPUTER            - CDCFT5/SINGLE                                 USMN0050
C                                                                       USMN0060
C   LATEST REVISION     - JANUARY 1, 1978                               USMN0070
C                                                                       USMN0080
C   PURPOSE             - DETERMINATION OF THE MINIMUM AND MAXIMUM      USMN0090
C                           VALUES OF A VECTOR                          USMN0100
C                                                                       USMN0110
C   USAGE               - CALL USMNMX (X,N,INC,XMIN,XMAX)               USMN0120
C                                                                       USMN0130
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N FROM WHICH MINIMUM,  USMN0140
C                           MAXIMUM VALUES ARE TO BE TAKEN.             USMN0150
C                N      - LENGTH OF THE INPUT VECTOR X. (INPUT)         USMN0160
C                INC    - DISPLACEMENT BETWEEN CONSECUTIVE VALUES OF X  USMN0170
C                           TO BE CONSIDERED.                           USMN0180
C                XMIN   - OUTPUT SCALAR CONTAINING MINIMUM VALUE OF X.  USMN0190
C                XMAX   - OUTPUT SCALAR CONTAINING MAXIMUM VALUE OF X.  USMN0200
C                                                                       USMN0210
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         USMN0220
C                       - SINGLE/H36,H48,H60                            USMN0230
C                                                                       USMN0240
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 USMN0250
C                                                                       USMN0260
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           USMN0270
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      USMN0280
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  USMN0290
C                                                                       USMN0300
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       USMN0310
C                                                                       USMN0320
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN USMN0330
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    USMN0340
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        USMN0350
C                                                                       USMN0360
C-----------------------------------------------------------------------USMN0370
C                                                                       USMN0380
      SUBROUTINE USMNMX (X,N,INC,XMIN,XMAX)                             USMN0390
C                                                                       USMN0400
      DIMENSION          X(N)                                           USMN0410
C                                  FIRST EXECUTABLE STATEMENT           USMN0420
      XMIN = X(1)                                                       USMN0430
      XMAX = X(1)                                                       USMN0440
      DO 10 I=1,N,INC                                                   USMN0450
         IF (X(I) .GE. XMIN) GO TO 5                                    USMN0460
         XMIN = X(I)                                                    USMN0470
         GO TO 10                                                       USMN0480
   5     IF (X(I) .GT. XMAX) XMAX = X(I)                                USMN0490
  10  CONTINUE                                                          USMN0500
      RETURN                                                            USMN0510
      END                                                               USMN0520


