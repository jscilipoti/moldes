! 
!  Esta subrutina calcula la temperatura de ebullición de un dado número de componentenes.
!
!  Esto se logra resolviendo el siguiente función:

!    f(T) = ln(φ_L) - ln(φ_V) = 0
!
!  utilizando el método de Newton amortiguado. En este, no se permiten cambios mayores a cierto porcentaje del
!  valor actual de la temperatura de ebullición. Los errores en la temperatura de ebullición predichas por tu
!  correlación y el modelo no deberían diferir en más de 50%.
!  Naturalmente, el vector de composición está evaluado sólo para componentes puros. Recomiendo llamar a esta 
!  subrutina con NC = 1.
subroutine boilingPoint (NC, NG, NGA, Tc, Pc, Tsat,  dc, nu, q, Tstr, gstr, gdot, gddot, kstr, kdot, alpha, Tb,fallo)


  implicit none
  integer, parameter          :: Liq = 2, maxIt = 100, min_iter_next_volume = 5, n_deriv = 1, NCM = 15, NGM = 15, &
                                 temp_deriv = 1, Vap = 1
  
! deltaMax  = porcentaje de cambio máximo en T
! tolerance = paso mínimo en cada iteración. No se pretenden presiciones mayores a 0.01 K
! eps       = valor máximo para el cual 1 + eps = 1
  real(8), parameter          :: deltaMax = .05D0, eps = epsilon(1.D0), tolerance = 0.01D0
                             
                             
  integer                     :: error, guess, i, iter, j, k, NC, NG, NGA, NST, nu(NCM,NGM), phase_type
  real(8)                     :: deltaT, deltaNewton, dIsoFug_dT, isoFug, Told
  real(8), dimension(2)       :: dP_dT, dP_dV, v, Z
  real(8), dimension(NCM)     :: dc, DLPHIP, omega, Pc, Psat, Tb, Tc, Tsat, x
  real(8), dimension(2,NCM)   :: lnPhi, dLnPhi_dT
  real(8), dimension(NCM,NCM) :: dLPhi
  
  real(8), dimension(NGM)     :: gstr, gdot, gddot, q, Tstr
  real(8), dimension(NGM,NGM) :: alpha, kdot, kstr
 
  
! Variables dentro de COMMONs:
  integer, dimension(NCM,NGM) :: ny
  real(8)                     :: dPV, dPDT, R, ZZ
  real(8), dimension(NCM)     :: d, dcrit, dt, dPDn, HA, HB, Pcrit, Tcrit
  real(8), dimension(NGM)     :: epx, gs, g1,g2, qArea, Ts, Tspl
  real(8), dimension(NGM,NGM) :: a, dadt, akij, alfa, xkij
  
  COMMON /GROUP1/                GS, G1, G2, TS, TSPL, XKIJ, AKIJ, EPX
  COMMON /GROUP2/                Qarea, A, DADT, ALFA, R, NY
  COMMON /COORD/                 ZZ  
  COMMON /MOL/                   DCrit, D, DT, HA, HB
  COMMON /CRIT/                  TCrit, Pcrit
  common /ScndDer/               dPV, dPDT, dPDn  
  logical::fallo

  NST = 0 !no hay asociación en esta versión
  
  R = 82.05D0 !atm·cm3/mol K
  ZZ = 10D0   !Número de coordinación
  psat(:) = 1.
  
! Asignación de los vectores COMMON:  
  do j = 1, NG
    
    gs(j) = gstr(j)
    g1(j) = gdot(j)
    g2(j) = gddot(j)
    Ts(j) = Tstr(j)
    qArea(j) = q(j)
    Tspl(j) = 1D3
    epx(j) = 0D0
    do k = j + 1, NG
      
      xkij(k,j) = kstr(k,j)
      xkij(j,k) = kstr(k,j)
      akij(k,j) = kdot(k,j)
      akij(j,k) = kdot(k,j)
      alfa(k,j) = alpha(k,j)
      alpha(j,k) = alpha(j,k)
      
    enddo
    ny(:NC,j) = nu(:NC,j)
    
  enddo
  do i = 1, NC
    
    Tcrit(i) = Tc(i)
    Pcrit(i) = Pc(i)
    dcrit(i) = dc(i)
    
  enddo
  
!-------------------------------------------------------------------------------------

! Lazo de cálculo para cada componente:
!
  do i = 1, NC
    
    x(:NC) = eps
    x(i) = 1D0
    Tb(i) = Tsat(i)
    deltaT = 10
    guess = 0
    if (1 < Pc(i)) then !Psat == 1, normal boiling point
      do iter = 1, maxIt
        
        if (iter > min_iter_next_volume .OR. dabs(deltaT) < 1D0) then

  !       Estimación de los volúmenes de las fases que requiere la GCA para reproducir P y T.
  !       Esto se hace sólo luego de cierto número de iteraciones (min_iter_next_volume) o si el paso
  !       en temperatura es < 1 K.
          v(Vap) = Z(Vap)*R*Tb(i)/Psat(i)
          v(Vap) = v(Vap) + dP_dT(Vap)/dP_dV(Vap) *(Tb(i) - Told)
          Z(Vap) = Psat(i)*v(Vap)/R/Tb(i)             
          v(Liq) = Z(Liq)*R*Tb(i)/Psat(i)
          v(Liq) = v(Liq) + dP_dT(Liq)/dP_dV(Liq) *(Tb(i) - Told)
          Z(Liq) = Psat(i)*v(Liq)/R/Tb(i)
          guess = 2
          
        endif  
        call PARAGC (Tb(i), NC, NG, NST, temp_deriv) !evaluación de parámetros dependientes de la temperatura: g(T), d(T) y Delta(T) (si hay asoc.)
        
        phase_type = Liq - 1
        call GCEOS (NC, NG, NST, n_deriv, temp_deriv, Tb(i), Psat(i), x, lnPhi(Liq,:NCM), dLPhi, dLnPhi_dT(Liq,:NCM), &
                    DLPHIP, Z(Liq), guess, phase_type, error,fallo)      
        dP_dV(Liq) = dPV
        dP_dT(Liq) = dPdT
        
        phase_type = Vap - 2
        call GCEOS (NC, NG, NST, n_deriv, temp_deriv, Tb(i), Psat(i), x, lnPhi(Vap,:NCM), dLPhi, dLnPhi_dT(Vap,:NCM), &
                    DLPHIP, Z(Vap), guess, phase_type, error,fallo)
        dP_dV(Vap) = dPV
        dP_dT(Vap) = dPdT      
        
        isoFug = lnPhi(Liq,i) - lnPhi(Vap,i)              !F(T)
        dIsoFug_dT = dLnPhi_dT(Liq,i) - dLnPhi_dT(Vap,i)  !dF/dT
        if(isnan(isofug))return
        deltaNewton = -isoFug/dIsoFug_dT                   !Paso de Newton original
        deltaT = min(dabs(deltaNewton), dabs(deltaMax*Tb(i))) * deltaNewton/dabs(deltaNewton) !elección del mínimo entre el original y el máx permitido
        Told = Tb(i)
        if(isnan(Tb(1)))return
        Tb(i) = Tb(i) + deltaT
        
        if (dabs(deltaT) < tolerance) exit
        
      enddo
      
    endif
    
  enddo
    
  
  
  

  return

endsubroutine boilingPoint