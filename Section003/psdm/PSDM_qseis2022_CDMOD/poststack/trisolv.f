C ******************************************************************
C  
C  SUBROUTINE TRISOLV(P,ZA,ZB,ZAA,ZBB,ZBNX1,ZBNXN,ZBBNX1,ZBBNXN,ZE,NX)
C                                                          
C    IMPLICIT COMPLEX TRIDIAGONAL SOLVER:        
C                                                          
C            P  -   NX  (WAVEFIELD VECTOR / INPUT & OUTPUT)
C           ZA  -   NX                  
C           ZB  -   NX                                      
C          ZAA  -   NX                                 
C          ZBB  -   NX                                  
C        ZBNX1  -   SCALAR                                        
C        ZBNXN  -   SCALAR                                         
C       ZBBNX1  -   SCALAR                                   
C       ZBBNXN  -   SCALAR                                   
C           ZE  -   NX  (WORK ARRAY)
C
C
C    SOLVES THE FOLLOWING FIRST-ORDER LINEAR DIFFERENCE EQUATION:
C
C                                                           (K+1)
C    ___                                             ___    __  __ 
C   |                                                   |  |      |
C   | ZBNX1	ZA(2)    0      ...                         |  | P(1) |
C   |                                                   |  |      |
C   | ZA(1)   ZB(2)   ZA(3)    0                        |  | P(2) |
C   |                                                   |  |  .   |
C   |   0     ZA(2)   ZB(3)   ZA(4)                     |  |  .   |
C   |                                                   |  |  .   |
C   |                          ...                      |  |      |
C   |                                                   |  |      |
C   |                                 ZA(NX-1)     0    |  |  .   |
C   |                                                   |  |  .   |
C   |                                 ZB(NX-1)   ZA(NX) |  |P(NX-1)
C   |                                                   |  |      |
C   |                            0    ZA(NX-1)   ZBNXN  |  | P(NX)|
C   |___                                             ___|  |__  __|
C
C
C                                                              (K)
C      _______                                        ____    __  __ 
C     |                                                   |  |      |
C     |ZBBNX1   ZAA(2)    0       ...                     |  | P(1) |
C     |                                                   |  |      |
C     |ZAA(1)   ZBB(2)   ZAA(3)    0                      |  | P(2) |
C     |                                                   |  |  .   |
C     |   0     ZAA(2)   ZBB(3)  ZAA(4)                   |  |  .   |
C     |                                                   |  |  .   |
C     |                          ...                      |  |      |
C  =  |                                                   |  |      |
C     |                                ZAA(NX-1)     0    |  |  .   |
C     |                                                   |  |  .   |
C     |                                ZBB(NX-1)  ZAA(NX) |  |P(NX-1)
C     |                                                   |  |      |
C     |                            0   ZAA(NX-1)  ZBBNXN  |  | P(NX)|
C     |___                                           _____|  |__  __|
C
C
C    NOTE THAT WITH THIS INDEXING SCHEME THE INDICES LINE UP VERTICALLY
C    IN THE MATRIX COLUMNS. THIS IS BECAUSE THE PHYSICAL CONSTANTS
C    WHICH DEFINE THE MATRIX ENTRIES DEPEND ON SPATIAL LOCATION.
C
C     
      SUBROUTINE TRISOLV( P, ZA, ZB, ZAA, ZBB,
     1			    ZBNX1, ZBNXN, ZBBNX1, ZBBNXN,
     2           	    ZE, NX )
C  
      IMPLICIT NONE
      INTEGER  NX, IX, NX1, IXM1, IXP1
      COMPLEX      P(  NX)                                  
      COMPLEX     ZA(  NX)                                  
      COMPLEX     ZB(  NX)                                  
      COMPLEX    ZAA(  NX)                                  
      COMPLEX    ZBB(  NX)                                  
      COMPLEX     ZE(  NX)                                  
      COMPLEX  DEN, ZF, ZBNX1, ZBNXN, ZBBNX1, ZBBNXN              
C                                                           
C                                                           
C     (A)  FORWARD                                          
C                                                           
C     IX=1                                                  
C                                                           
      ZE(1) = -ZA(2)/ZBNX1                                  
      ZF    = (ZBBNX1*P(1)+ZAA(2)*P(2))/ZBNX1               
C       
C     IX=2                                                  
C   
      DEN=1./(ZB(2)+ZA(1)*ZE(1))                            
      ZE(2)=-ZA(3)*DEN                                      
      P(1)=((ZAA(1)*P(1)+ZBB(2)*P(2)                        
     *    + ZAA(3)*P(3))-ZA(1)*ZF)*DEN                      
C                                                           
      NX1 = NX-1
      DO IX=3,NX1                                       
          IXM1=IX-1                                             
          IXP1=IX+1                                             
          DEN=1./(ZB(IX)+ZA(IXM1)*ZE(IXM1))                     
          ZE(IX)=-ZA(IXP1)*DEN                                  
          P(IXM1)=((ZAA(IXM1)*P(IXM1)+ZBB(IX)*P(IX)             
     *    +   ZAA(IXP1)*P(IXP1))-ZA(IXM1)*P(IXM1-1))*DEN      
      ENDDO
C                                                           
C     (B)  LAST ELEMENT      IX=NX                          
C                                                           
      P(NX)=( ZBBNXN*P(NX)+ZAA(NX1)*P(NX1)                  
     *   -ZA(NX1)*P(NX1-1) )/(ZBNXN+ZA(NX1)*ZE(NX1))        
C
C     (C)  BACKWARD                                         
C                                                           
      DO IX=NX1,2,-1                                    
          IXM1=IX-1                                             
          IXP1=IX+1                                             
          P(IX)=ZE(IX)*P(IXP1)+P(IXM1)                          
      ENDDO
C                                                           
C     IX=1                                                  
C                                                           
      P(1)=ZE(1)*P(2)+ZF                                    
C                                                           
C                                                           
      RETURN                                                
      END                                                   
