real*8 function NEWTON (x,KeqExp,Zt0,compuestos,T,P,Ncomponentes)
           
    !x:          conversión
    !Keq:        Contante de equilibrio de reacción
    !Z:          Composición inicial global del sistema
    !compuestos: componentes del sistema
    !T:          Temperatura del sistema
    !P:          Presión del sistema
    !N:          Número de componentes
    
    use constantes
    
    implicit none
    
    real*8,intent(inout)::x
    real*8,intent(in)::KeqExp
    real*8,intent(in)::Zt0(NCOM)
    integer,intent(in)::compuestos(NCOM,DiffStructGroups,2)
    real*8,intent(in)::T,P
    integer,intent(in)::Ncomponentes
!Variables internas
    real*8::act(NCOM)
    real*8::Ztf(NCOM)
    real*8::epsilon,f_value,per,FDeriv_Old,  F_ValueOld
      INTEGER :: NumIterations, N

      real*8:: XOLD,XNEW
      CHARACTER (1) :: A
      real*8,external::F,fd

      ! Get termination values Epsilon and NumIterations and 
      ! initial approximation

      NumIterations=1000
      Epsilon=1E-3
      PRINT *
      ! Initialize F_Value and N; print headings and initial values

      F_Value = F(x,KeqExp,Zt0,compuestos,T,P,Ncomponentes)
      N = 0
      PRINT *, "  N      X(N)       F(X(N))"
      PRINT *, "============================="
      PRINT 10, 0, X, F_Value, 0.0
  10  FORMAT(1X, I3, F11.5, 2E14.5)
    if(isnan(F_Value))then !salgo de la función
        Newton = F_Value
        return
    endif
! Iterate using Newton's method while ABS(F_Value) is greater
! than or equal to Epsilon and N has not reached NumIterations

      !DO
              PER = X*1E-2
        ! If a termination condition met, stop generating approximations

        ! Otherwise continue with the following
        N = N + 1
        FDeriv_Old = FD(X,PER,KeqExp,Zt0,compuestos,T,P,Ncomponentes)

    if(isnan(FDeriv_Old))then !salgo de la función
        Newton = FDeriv_Old
        return
    endif        
        
        ! Terminate if the derivative is 0 at some approximation
        IF (FDeriv_Old == 0) THEN
            PRINT *, "Newton's method fails -- derivative = 0"
        END IF
        XOLD=X
        ! Generate a new approximation
        X = XOLD - (F_Value / FDeriv_Old)
        
        if(x>0.99)then
            x=0.8
        elseif(x<0.001)then
            x=0.2
        endif
        
        F_Value = F(X,KeqExp,Zt0,compuestos,T,P,Ncomponentes)
    if(isnan(F_Value))then !salgo de la función
        Newton = F_Value
        return
    endif        
        F_ValueOld = F(XOLD,KeqExp,Zt0,compuestos,T,P,Ncomponentes)
    if(isnan(F_ValueOld))then !salgo de la función
        Newton = F_ValueOld
        return
    endif        
        PRINT 10, N, X, F_Value, FDeriv_Old

      !END DO
!C
!C-----Método de la secante
!C
      DO
        IF ((ABS(X-XOLD) < Epsilon) .OR. (N > NumIterations)) EXIT
        
        N = N + 1

        
        XNEW = X - ((X-XOLD)/(F_Value-F_ValueOld))*F_Value
        
        XOLD = X
        if(xnew>0.99)then
            x=0.8
        elseif(xnew<0.001)then
            x=0.2
        else
            X = XNEW
        endif
        F_Value = F(X,KeqExp,Zt0,compuestos,T,P,Ncomponentes)
    if(isnan(F_Value))then !salgo de la función
        Newton = F_Value
        return
    endif          
        F_ValueOld = F(XOLD,KeqExp,Zt0,compuestos,T,P,Ncomponentes)
    if(isnan(F_ValueOld))then !salgo de la función
        Newton = F_ValueOld
        return
    endif          
        PRINT 10, N, X, F_Value, F_ValueOld
        
      ENDDO
      
      
      
      newton = f_value
 

    endfunction NEWTON
    
    
    
    
real*8 FUNCTION FD (X,Per,KeqExp,Zt0,compuestos,T,P,Ncomponentes)
!*****************************************************
!     Derivada de la Función Objetivo
!*****************************************************      
    use constantes      
implicit none
    real*8,intent(inout)::x
    real*8,intent(in)::KeqExp
    real*8,intent(in)::Zt0(NCOM)
    integer,intent(in)::compuestos(NCOM,DiffStructGroups,2)
    real*8,intent(in)::T,P,per
    integer,intent(in)::Ncomponentes
      real*8::XP
      real*8,external::f
      
      XP=X+Per
             
      FD=(F(XP,KeqExp,Zt0,compuestos,T,P,Ncomponentes)-F(X,KeqExp,Zt0,compuestos,T,P,Ncomponentes))/Per
      
      END FUNCTION FD