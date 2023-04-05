!MODULE PROP_PDB
!contains
!!==========================================================================
!Subroutine Search_Isomers(nsol)
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!  use constantes
!  use StructuresDesign
!  implicit none
!!Variables de ENTRADA
!  integer,intent(in)::nsol
!!Variables INTERNAS
!  type(FinalStructure),pointer::recorre
!  type(isomers),pointer::recorreIsomer
!  integer::nisom,NI,k,i,n
!  integer,dimension(8)::IDENT
!  integer,dimension(DiffStructGroups,2)::UNIF
!  real*8::aux,VLIQ,BP
!  character(len=20)::nameaux
!  character(len=35)::ANAME,formula(8),form,anamet(8)
!!Sentencias----------------------------------------------------------------
!    write(6,*) ' Checking solvent structures with data base'
!	recorre=>FMSs
!	call ordenar(recorre%Formula)
!	DO i=1,nsol
!		call checkdatos (recorre%Formula,nisom,ident,anamet,formula)
!		recorre%NumberIsomers = nisom
!		if(nisom > 0)then
!		    allocate(recorre%nextisomer)
!		    recorreIsomer => recorre%nextisomer
!			do n=1,nisom
!				READ(7,"(I4,A35,A20,A35,20I3,7D13.6)",REC=IDENT(N)) NI,ANAME,nameaux,form,&
!                    (UNIF(K,1),UNIF(K,2),K=1,10),aux,aux,aux,aux,aux,BP,VLIQ
!                recorreIsomer%index = NI
!                recorreIsomer%name = aname
!                recorreIsomer%Formula(:,:) = UNIF(:,:)
!                recorreIsomer%BoilingPoint = BP
!                recorreIsomer%LiquidMolarVolume = VLIQ
!                recorreIsomer%FormChem = form
!                
!                allocate(recorreIsomer%next)
!                recorreIsomer => recorreIsomer%next
!!                CALL RELVOLPRO (RELV,TB(J,N))
!!                VRELSOLV(J,N)=RELV 
!			enddo
!			nullify(recorreIsomer)
!		endif
!		if(.not.associated(recorre%next))exit
!		recorre => recorre%next
!		call ordenar(recorre%Formula)
!    enddo
!endsubroutine Search_Isomers
!
!logical function check_inter_param (ipareq,imprim,kgrup,iprin)
!!c-----------------------------------------------------------------------------
!!c     Esta subrutina controla la existencia de parametros de interaccion entre
!!c     los subgrupos cuyas identificaciones  estan  almacenadas  en  el  vector
!!c     kgrup.  Toma un elemento de  kgrup  y  controla   que  ese  grupo  tenga
!!c     interaccion con los grupos almacenados desde  el  lugar  iprin  a   iult
!!c     del vector kgrup.  Continua asi tomando cada elemento de kgrup desde  el
!!c     primero hasta iult-1.  Iprin puede valer 2 o si  ya  fueron  controlados
!!c     los n primeros elementos, iprin puede valer desde n+1 a iult.
!!c
!!c     Variables de entrada:
!!c     --------------------
!!c                   ipareq: modelo de parametros de equilibrio:
!!c                           1:liquido-liquido
!!c                           2:liquido-vapor
!!c                           3:liquido-vapor a dilucion infinita.
!!c                   imprim: 0: sin salida.
!!c                           1: salida por pantalla.
!!c                           2: salida por pantalla y por archivo (unidad 2).
!!c                    kgrup: vector de  numeros  de  identificacion  de  grupos
!!c                           principales a ser analizado.  Su dimension es nsel.
!!c                    iprin: elemento a partir del cual se desea analizar kgrup.
!!c
!!c     Variables de salida:
!!c     -------------------
!!c                     prob: true:  todos   los   grupos  de  kgrup  no  tienen
!!c                                  parametros de interaccion
!!c                           false: todos los grupos de kgrup tienen parametros
!!c                                  de interaccion.
!!c
!!c     La dimension de los vectores queda fijada por los siguientes parametros:
!!c
!      use CONSTANTES
!!c
!!c-----------------------------------------------------------------------------
!      implicit real*8 (a-h,o-z)
!      external nom_gru1,Leer_In
!      dimension kgrup(nint)
!      character*8 fs1,fs2,nom_gru1
!      logical prob
!      integer inicio
!!c     common ipant
!!c
!      inicio = iprin
!      prob = .false.
!      ind=iprin
!      do while(kgrup(ind)/=0)
!        if (ind.eq.inicio) then
!            inicio = inicio + 1
!        end if
!        j=inicio
!        do while (kgrup(j)/=0)
!            k1 = kgrup(ind)
!            k2 = kgrup(j)
!            fs1 = nom_gru1 (k1,ipareq)
!            fs2 = nom_gru1 (k2,ipareq)
!            apar = Leer_In (k1,k2,ipareq)
!            if (dabs(apar-9000.).le.0.1) then
!               if (imprim.gt.0) then
!                  write (6,10)  fs1,fs2
!                  if (imprim.eq.2) then
!                     write (2,10)  fs1,fs2
!                  end if
!               end if
!               prob = .true.
!            end if
!            j=j+1
!        enddo
!        ind=ind+1
!      enddo
!      
!      check_inter_param = prob
!      
!!---- formatos
! 10   format (' ',/,' ** Interaction parameters not available for: ',2(a8,2x))
!      
!      
!endfunction check_inter_param
!
!
!!==========================================================================
!      SUBROUTINE CHECK_INTRCN
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------   
!      use SubGrupos
!      use CONSTANTES
!      IMPLICIT real*8 (A-H,O-Z)
!!INTERFACES
!interface
!    subroutine ab_ban1(mod)
!    use input_data, only:Model_S
!    integer, intent(out),optional::mod
!    endsubroutine ab_ban1
!endinterface
!      INTEGER GRUPOS(NMG)      
!      CHARACTER*8 MGV(NMG)
!!Sentencias----------------------------------------------------------------        
!      CALL AB_BAN1(model)
!      CALL TABLA_PARAMETROS(IPAREQ)
!      ngrmax = max_sub_int (ipareq)
!      DO I=1,NGRMAX
!        irec =I+(ipareq-1)*70
!		grupos(i) = i
!		READ(13,24,rec=irec)MGV(I)
!      ENDDO
!      CALL Groups_Present ()
!      CALL PAUSA
!      CALL CI_BAN1
!!C-----FORMATOS
!  24  FORMAT (a8)      
!      RETURN
!      endsubroutine
!!==========================================================================
!      SUBROUTINE CHECK_GRUPOS
!      use SubGrupos
!      use CONSTANTES
!      IMPLICIT real*8 (A-H,O-Z)
!!INTERFACES
!interface
!    subroutine ab_ban1(mod)
!    use input_data, only:Model_S
!    integer, intent(out),optional::mod
!    endsubroutine ab_ban1
!endinterface
!      INTEGER GRUPOS(NMG)      
!      CHARACTER*8 MGV(NMG)
!!Sentencias----------------------------------------------------------------          
!      CALL AB_BAN1(model)
!      CALL TABLA_PARAMETROS(IPAREQ)
!      ngrmax = Size_LSubGroups()
!      DO I=1,NGRMAX
!        irec =I+(ipareq-1)*150
!		grupos(i) = i
!		READ(14,24,rec=irec)MGV(I)
!      ENDDO
!      CALL Groups_Present ()
!      CALL PAUSA
!      CALL CI_BAN1
!!C-----FORMATOS
!  24  FORMAT (12x,a8)      
!      RETURN
!      endsubroutine
!!subroutine Search_Isomers(FMSs)
!!  use StructuresDesign
!!  implicit none
!!!Variables de ENTRADA
!!  type(FinalStructure),pointer,intent(in)::FMSs
!!!Variables INTERNAS
!!  type(FinalStructure),pointer::recorre
!!!SENTENCIAS
!!    recorre=>FMSs      
!!    do while (associate)
!!      
!!endsubroutine Search_Isomers
!
!!==========================================================================      
!subroutine load_pure_comp_prop (puntero)
!!--------------------------------------------------------------------------	
!!   description
!!--------------------------------------------------------------------------	
!    use CONSTANTES
!    use Input
!    implicit none
!!Variable de ENTRADA/SALIDA
!    type(Compound),pointer,intent(inout)::puntero
!!Variables INTERNAS
!    integer::nrcar,nisom
!	integer::ident(8) !compt(10,2),UNIF(10,2),
!	character*35::nombre(8),formula(8) 
!	character*37::nomco3,nomco4
!!SENTENCIAS
!
!	    nrcar=0
!        if (ipareq.eq.2) then
!		    call checkdatos (puntero%Formula,nisom,ident,nombre,formula)! checkea en prop.pdb
!		    
!		    if (nisom /= 0) call elegir_CXR (nisom,ident,nombre,formula,nrcar) !Si existe algún isómero, pregunta al usuario cuál elegir
!        end if
!        
!!-------Carga de propiedades de componente puro
!	    if (nrcar.ne.0) then !si existe en la base de datos
!		    call pvomegaybp (nrcar,puntero)
!
!	    else !carga manual si NO existe en la base de datos 
!            nomco4 = 'component                            '
!!		    if (InputProblem01%mop.eq.1) then ! Lectura de la temperatura de ebullicion del CAR
!!			    nomco4 = 'component to be recovered:           '
!!		    else
!!			    nomco3 = 'more  volatile  component:           '
!!			    nomco4 = 'less  volatile  component:           '
!!		    end if
!3010		write (6,"(1x,/,' Give the boiling point (K) of the ',a37,/,64x,'> ',$)") nomco4
!		    read (5,*,err=3010) puntero%BoilingPoint
!
!		    if (InputProblem01%mop.eq.2.and.nrcar.eq.0) then ! Lectura de las constantes de Antoine
!			    call antoine (nomco4,'A1','A2','A3',puntero)
!		    end if
!        end if
!
!
!
!endsubroutine load_pure_comp_prop
!
!      subroutine antoine (nomco4,c1,c2,c3,puntero)
!!c-----------------------------------------------------------------------
!!c     Esta subrutina carga las constantes de Antoine para un componente
!!c-----------------------------------------------------------------------
!      use Input
!      implicit real*8 (a-h,o-z)
!      type(Compound),pointer,intent(inout)::puntero
!      character*1 iopt
!      character*2 c1,c2,c3
!      character*37 nomco4
!!c
!      idev = 6
!	idevr = 5
!!c     call limp (idev)
!  70  write (idev,50) nomco4
!      write (idev,110) c1,c2,c3
!  80  write (idev,60) c1
!      read (idevr,*,err=80) puntero%a(1)
!  90  write (idev,60) c2
!      read (idevr,*,err=90) puntero%a(2)
! 100  write (idev,60) c3
!      read (idevr,*,err=100) puntero%a(3)
!!c     call limp (idev)
! 170  write (idev,180) (puntero%a(i),i=1,3) 
!      write (idev,970)
!      read (idevr,980,err=170) icomp
!      if (icomp.eq.1) then
! 130     write (idev,120) c1,c2,c3
!         read (idevr,510,err=130) iopt
!         i1 = index ('123',iopt)
!         if (i1.eq.0) then
!!c           call limp (idev)
!            go to 130
!         end if
!!c        call limp (idev)
!         write (idev,110) c1,c2,c3
!         if (iopt.eq.'1') then
! 140        write (idev,60) c1
!            read (idevr,*,err=140) puntero%a(1)
!            go to 170
!         else if (iopt.eq.'2') then
! 150        write (idev,60) c2
!            read (idevr,*,err=150) puntero%a(2)
!            go to 170
!         else if (iopt.eq.'3') then
! 160        write (idev,60) c3
!            read (idevr,*,err=160) puntero%a(3)
!            go to 170
!         end if        
!      end if
!!c
!!c---formatos
!  50  format (' Give the Antoine constants for the ',a37,/)
! 110  format (11x,'ln(Pv/mmHg) = ',a2,' - ',a2,'/ (T/K + ',a2,')',/)
! 180  format (1x,/,1x,'ln(Pv/mmHg) = ',f9.4,' - ',f9.2,'/ (T/K + ',f9.2,')')
!  60  format (1x,33x,a2,'  > ',$)
! 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
! 980  format (i2)
! 120  format(' ',/,' Which constant do you want to change?',//,34x,&
!            a2,': 1',/,34x,a2,': 2',/,34x,a2,': 3',//,42x,'> ',$)
! 510  format (a)
!      return 
!      endsubroutine antoine
!
!
!
!!==========================================================================      
!	subroutine checkdatos (compt,nisomt,identt,nombret,formulat)
!!--------------------------------------------------------------------------	
!!	Esta subrutina checkea o confirma la existencia en el
!!	banco de datos "Prop.pdb" de uno o mas compuestos cuya 
!!	estructura de grupos según UNIFAC sea igual a la que 
!!	contiene el vector de entrada "compt".
!!	La variable de salida "nisomt" indica el número de
!!	diferentes isómeros que satisfacen tal condición.
!!	Sus números identificatorios y nombres se almacenan en 
!!	los vectores "identt" y "nombret" respectivamente.
!!--------------------------------------------------------------------------	
!	parameter(NA=150)
!	integer compt(10,2),UNIF(10,2),identt(8)
!	INTEGER I,J,AUX,AUXB
!	character*35 nombret(8),formulat(8) 
! 	CHARACTER*35 ANAME,FORM2
!	CHARACTER*20 FORM1
!    CHARACTER*8 fst(NA)
!    logical isom
!	COMMON/NOM/fst
!!Sentencias----------------------------------------------------------------    
!    nisomt = 0
!	DO i=22,1097
!	    
!		READ(7,100,REC=i,ERR=120) NI,ANAME,FORM1,FORM2,&
!                                   (UNIF(j,1),UNIF(j,2),j=1,10)
! 100 		FORMAT(I4,A35,A20,A35,20I3)
!
!        call ordenar (UNIF)
!		if (ni.eq.0) cycle
!		call comparar (UNIF,compt,1,isom)
!        if(isom)then
!		    nisomt = nisomt + 1
!			identt(nisomt) = ni
!			nombret(nisomt) = aname
!			formulat(nisomt) = form2
!		end if	
!    enddo
!	GOTO 500
!
!120 WRITE(6,*) '*ERROR* EN LA LECTURA DE "PROP.PDB"'
!
!500	RETURN 
!	endsubroutine
!!==========================================================================
!	subroutine elegir_CXR (nisom,ident,nombre,formula,nr)
!	integer ident(8)
!	character*35 nombre(8),formula(8)
!!Sentencias----------------------------------------------------------------	
!	write (6,9)
!  9	format (/,&
!        '  The data bank "PROP.PDB" contains information about the',/,&
!     	'  following isomers for the UNIFAC structure you have entered:'&
!     	,//,&
!     	'  R.N.	Name				   Formula')
!    do 10 i=1,nisom
!	    write (6,7) ident(i),nombre(i),formula(i)
!   7	format (x,i4,4x,a35,a35)
! 10 continue
!100	write (6,8)
!  8	format (/,&
!     	'  If you wish, enter a Reference Number and its properties',/,&
!     	'  from the data bank will be taken into account.',/,&
!     	'  or choose 0 to work with the MOLDES group contribution',/,&
!     	'  estimations for physical properties.',/,&
!     	60x,'>')
!	read (5,*) nr
!	if (nr.ne.0.and.nr.ne.ident(1).and.nr.ne.ident(2)&
!     	.and.nr.ne.ident(3).and.nr.ne.ident(4).and.nr.ne.ident(5)&
!     	.and.nr.ne.ident(6).and.nr.ne.ident(7).and.nr.ne.ident(8)) &
!     	goto 100
!    RETURN 
!	endsubroutine
!!==========================================================================
!	subroutine pvomegaybp (nr,puntero)
!!--------------------------------------------------------------------------	
!!	Esta subrutina utiliza Tb, Tc, Pc y W de una sustancia determinada
!!	(de la base de datos PROP.PDB) para calcular las constantes A1 y
!!	A2 que sean útiles para estimar la Pv de la sustancia a temperaturas
!!	cercanas a su punto de ebullición según la ecuación:
!!	Ln(Pv)=A1-A2/T 
!!--------------------------------------------------------------------------
!    use Input
!    use SubGrupos
!    implicit real*8 (a-h,o-z)
!    type(Compound),pointer,intent(inout)::puntero
!	integer UNIF(20)
! 	CHARACTER*35 ANAME,FORM2
!	CHARACTER*20 FORM1
!!Sentencias----------------------------------------------------------------		
!
!	READ(7,100,REC=NR,ERR=20) NI,ANAME,FORM1,FORM2,&
!                                   (UNIF(K),K=1,20),TC,PC,VC,ZC,TF,TB,&
!                                   VLIQ,W
!    puntero%Name = ANAME
!    puntero%TC = TC
!    puntero%PC = PC
!    puntero%VC = VC
!    puntero%BoilingPoint = TB
!    puntero%vliq = VLIQ
!	pmt=0
!	i=1
!	do while(puntero%Formula(i,1)/=0)
!		PMT = PMT + Obtain_MW(puntero%Formula(i,1))*puntero%Formula(i,2) 
!        i=i+1
!    enddo 
!    puntero%MW = pmt    
!    puntero%Dens = 0.001*PMT/VLIQ
!100	FORMAT(I4,A35,A20,A35,20I3,8D13.6)
!!c	Pv = Pc * 10** ((7/3)*(1+W)*(T-Tc)/T) !Esto era interpolando entre el punto crítico y Tr=0.7
!!c	Ahora, interpolando entre Tb y Tr=0.7:
!    puntero%a(2) = dlog(Pc*10**(-1-W)/101325)/((1/Tb) - (1/(0.7*Tc)))
!	puntero%a(1) = alog(760.) + puntero%a(2)/Tb
!    GOTO 500
!
! 10	continue
!	GOTO 500
! 20	WRITE(6,*) '*ERROR* EN LA LECTURA DE "PROP.PDB"'
!
!500	RETURN 
!	endsubroutine
!ENDMODULE PROP_PDB