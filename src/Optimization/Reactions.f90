real*8 function F(x,KeqExp,Zt0,compuestos,T,P,N)

    !x:          conversión
    !Keq:        Contante de equilibrio de reacción
    !Z:          Composición inicial global del sistema
    !compuestos: componentes del sistema
    !T:          Temperatura del sistema
    !P:          Presión del sistema
    !N:          Número de componentes
    
    use constantes
    
    implicit none
    
    real*8,intent(in)::x,KeqExp
    real*8,intent(in)::Zt0(NCOM)
    integer,intent(in)::compuestos(NCOM,DiffStructGroups,2)
    real*8,intent(in)::T,P
    integer,intent(in)::N
    !Variables internas
    real*8::act(NCOM)
    real*8::Ztf(NCOM)
    
    Ztf(5) = Zt0(5)
    Ztf(4) = Zt0(4) + Zt0(1)*x
    Ztf(3) = Zt0(3) + Zt0(1)*x
    Ztf(2) = Zt0(2) - Zt0(1)*x
    Ztf(1) = Zt0(1) - Zt0(1)*x
    
    
    call llecalas(T,P,compuestos,N,Ztf,act)
    
    F = (KeqExp - (  (act(4)*act(3))  / (act(1)*act(2)) ))**2
    
    
endfunction F