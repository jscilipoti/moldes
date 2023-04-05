	SUBROUTINE PAUSA
	LOGICAL SIGUE
	WRITE(6,1010)
	READ(5,890)   SIGUE
 890	FORMAT(A)
1010	FORMAT(' ',70X,'<RET>',$)
	RETURN
	END
integer function count_groups(comp)
!------------------------------------------------------------
!   Cuenta la cantidad de grupos en el vector comp
!------------------------------------------------------------
use CONSTANTES
implicit none
integer,intent(in)::comp(DiffStructGroups)
integer::i
!SENTENCIAS
    i=1
    do while(comp(i)/=0)
        i=i+1
    enddo
    count_groups = i-1

endfunction count_groups 


SUBROUTINE GAUSL(ND,NCOL,N,NS,A)
!--------------------------------------------------------------------
!C  SUBROUTINE GAUSL SOLVES N LINEAR ALGEBRAIC EQUATIONS BY GAUSS        
!C  ELIMINATION WITH ROW PIVOTING                                        
!C  TO SOLVE THE PROBLEM QX=U, WHERE Q IS A NXN MATRIX AND U IS NXNS,    
!C  ONE PLACES Q IN THE FIRST N COLUMNS OF A AND U IS PLACED IN THE      
!C  FOLLOWING NS COLUMNS.                                                
!C  THE PROGRAM RETURNS X=Q**(-1)*U AT THE PREVIOUS POSITION OF U.       
!C  *                                                                    
!C  ND IS THE ROW DIMENSION AND NCOL IS THE COLUMN DIMENSION OF A.       
!C  BOTH MUST BE TRANSFERRED TO THE SUBROUTINE.                          
!---------------------------------------------------------------------    
      implicit real*8 (A-H,O-Z)                               
      DIMENSION A(ND,NCOL)                                              
      N1=N+1                                                            
      NT=N+NS                                                           
      IF (N .EQ. 1) GO TO 50                                            
!C      START ELIMINATION                                                
!C                                                                       
!C                                                                       
      DO 10 I=2,N                                                       
        IP=I-1                                                          
        I1=IP                                                           
        X=DABS(A(I1,I1))                                                
        DO 11 J=I,N                                                     
            IF (DABS(A(J,I1)) .LT. X) GO TO 11                          
            X=DABS(A(J,I1))                                             
            IP=J                                                        
   11   CONTINUE                                                        
        IF (IP .EQ. I1) GO TO 13                                        
!C                                                                       
!C     ROW INTERCHANGE                                                   
!C                                                                       
        DO 12 J=I1,NT                                                   
            X=A(I1,J)                                                   
            A(I1,J)=A(IP,J)                                             
   12       A(IP,J)=X                                                   
   13   DO 10 J=I,N                                                     
            X=A(J,I1)/A(I1,I1)                                          
            DO 10 K=I,NT                                                
   10 A(J,K)=A(J,K) - X*A(I1,K)                                         
!C                                                                       
!C      ELIMINATION FINISHED, NOW BACKSUBSTITUTION                       
!C                                                                       
   50 DO 20 IP=1,N                                                      
        I=N1-IP                                                         
        DO 20 K=N1,NT                                                   
            A(I,K) = A(I,K)/A(I,I)                                      
            IF (I .EQ. 1) GO TO 20                                      
            I1=I-1                                                      
            DO 25 J=1,I1                                                
   25           A(J,K) = A(J,K) - A(I,K)*A(J,I)                         
   20 CONTINUE                                                          
      RETURN                                                            
endsubroutine
SUBROUTINE LIMP (IDEV)
IF (IDEV.EQ.6) THEN
   WRITE (IDEV,*) CHAR(27),'[1;1H',CHAR(27),'[J'
END IF
RETURN
END