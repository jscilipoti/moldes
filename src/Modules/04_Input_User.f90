module input_data
contains
!==========================================================================
Subroutine Search_Isomers(nsol)
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
  use constantes
  use StructuresDesign
  implicit none
!Variables de ENTRADA
  integer,intent(in)::nsol
!Variables INTERNAS
  type(FinalStructure),pointer::recorre
  type(isomers),pointer::recorreIsomer
  integer::nisom,NI,k,i,n
  integer,dimension(8)::IDENT
  integer,dimension(DiffStructGroups,2)::UNIF
  real*8::aux,VLIQ,BP
  character(len=20)::nameaux
  character(len=35)::ANAME,formula(8),form,anamet(8)
!Sentencias----------------------------------------------------------------
    write(6,*) ' Checking solvent structures with data base'
	recorre=>FMSs
	call ordenar(recorre%Formula)
	DO i=1,nsol
		call checkdatos (recorre%Formula,nisom,ident,anamet,formula)
		recorre%NumberIsomers = nisom
		if(nisom > 0)then
		    allocate(recorre%nextisomer)
		    recorreIsomer => recorre%nextisomer
			do n=1,nisom
				READ(7,"(I4,A35,A20,A35,20I3,7D13.6)",REC=IDENT(N)) NI,ANAME,nameaux,form,&
                    (UNIF(K,1),UNIF(K,2),K=1,10),aux,aux,aux,aux,aux,BP,VLIQ
                recorreIsomer%index = NI
                recorreIsomer%name = aname
                recorreIsomer%Formula(:,:) = UNIF(:,:)
                recorreIsomer%BoilingPoint = BP
                recorreIsomer%LiquidMolarVolume = VLIQ
                recorreIsomer%FormChem = form
                
                allocate(recorreIsomer%next)
                recorreIsomer => recorreIsomer%next
!                CALL RELVOLPRO (RELV,TB(J,N))
!                VRELSOLV(J,N)=RELV 
			enddo
			nullify(recorreIsomer)
		endif
		if(.not.associated(recorre%next))exit
		recorre => recorre%next
		call ordenar(recorre%Formula)
    enddo
endsubroutine Search_Isomers

logical function check_inter_param (ipareq,imprim,kgrup,iprin)
!c-----------------------------------------------------------------------------
!c     Esta subrutina controla la existencia de parametros de interaccion entre
!c     los subgrupos cuyas identificaciones  estan  almacenadas  en  el  vector
!c     kgrup.  Toma un elemento de  kgrup  y  controla   que  ese  grupo  tenga
!c     interaccion con los grupos almacenados desde  el  lugar  iprin  a   iult
!c     del vector kgrup.  Continua asi tomando cada elemento de kgrup desde  el
!c     primero hasta iult-1.  Iprin puede valer 2 o si  ya  fueron  controlados
!c     los n primeros elementos, iprin puede valer desde n+1 a iult.
!c
!c     Variables de entrada:
!c     --------------------
!c                   ipareq: modelo de parametros de equilibrio:
!c                           1:liquido-liquido
!c                           2:liquido-vapor
!c                           3:liquido-vapor a dilucion infinita.
!c                   imprim: 0: sin salida.
!c                           1: salida por pantalla.
!c                           2: salida por pantalla y por archivo (unidad 2).
!c                    kgrup: vector de  numeros  de  identificacion  de  grupos
!c                           principales a ser analizado.  Su dimension es nsel.
!c                    iprin: elemento a partir del cual se desea analizar kgrup.
!c
!c     Variables de salida:
!c     -------------------
!c                     prob: true:  todos   los   grupos  de  kgrup  no  tienen
!c                                  parametros de interaccion
!c                           false: todos los grupos de kgrup tienen parametros
!c                                  de interaccion.
!c
!c     La dimension de los vectores queda fijada por los siguientes parametros:
!c
      use CONSTANTES
!c
!c-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      external nom_gru1,Leer_In
      dimension kgrup(nint)
      character*8 fs1,fs2,nom_gru1
      logical prob
      integer inicio
!c     common ipant
!c
      inicio = iprin
      prob = .false.
      ind=iprin
      do while(kgrup(ind)/=0)
        if (ind.eq.inicio) then
            inicio = inicio + 1
        end if
        j=inicio
        do while (kgrup(j)/=0)
            k1 = kgrup(ind)
            k2 = kgrup(j)
            fs1 = nom_gru1 (k1,ipareq)
            fs2 = nom_gru1 (k2,ipareq)
            apar = Leer_In (k1,k2,ipareq)
            if (dabs(apar-9000.).le.0.1) then
               if (imprim.gt.0) then
                  write (6,10)  fs1,fs2
                  if (imprim.eq.2) then
                     write (2,10)  fs1,fs2
                  end if
               end if
               prob = .true.
            end if
            j=j+1
        enddo
        ind=ind+1
      enddo
      
      check_inter_param = prob
      
!---- formatos
 10   format (' ',/,' ** Interaction parameters not available for: ',2(a8,2x))
      
      
endfunction check_inter_param


!==========================================================================
      SUBROUTINE CHECK_INTRCN
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------   
      use SubGrupos
      use CONSTANTES
      IMPLICIT real*8 (A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
  !  use input_data, only:Model_S
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
      INTEGER GRUPOS(NMG)      
      CHARACTER*8 MGV(NMG)
!Sentencias----------------------------------------------------------------        
      CALL AB_BAN1(model)
      CALL TABLA_PARAMETROS(IPAREQ)
      ngrmax = max_sub_int (ipareq)
      DO I=1,NGRMAX
        irec =I+(ipareq-1)*70
		grupos(i) = i
		READ(13,24,rec=irec)MGV(I)
      ENDDO
      CALL Groups_Present ()
      CALL PAUSA
      CALL CI_BAN1
!C-----FORMATOS
  24  FORMAT (a8)      
      RETURN
      endsubroutine
!==========================================================================
      SUBROUTINE CHECK_GRUPOS
      use SubGrupos
      use CONSTANTES
      IMPLICIT real*8 (A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
  !  use input_data, only:Model_S
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
      INTEGER GRUPOS(NMG)      
      CHARACTER*8 MGV(NMG)
!Sentencias----------------------------------------------------------------          
      CALL AB_BAN1(model)
      CALL TABLA_PARAMETROS(IPAREQ)
      ngrmax = Size_LSubGroups()
      DO I=1,NGRMAX
        irec =I+(ipareq-1)*150
		grupos(i) = i
		READ(14,24,rec=irec)MGV(I)
      ENDDO
      CALL Groups_Present ()
      CALL PAUSA
      CALL CI_BAN1
!C-----FORMATOS
  24  FORMAT (12x,a8)      
      RETURN
      endsubroutine
!subroutine Search_Isomers(FMSs)
!  use StructuresDesign
!  implicit none
!!Variables de ENTRADA
!  type(FinalStructure),pointer,intent(in)::FMSs
!!Variables INTERNAS
!  type(FinalStructure),pointer::recorre
!!SENTENCIAS
!    recorre=>FMSs      
!    do while (associate)
!      
!endsubroutine Search_Isomers

!==========================================================================      
subroutine load_pure_comp_prop (puntero)
!--------------------------------------------------------------------------	
!   description
!--------------------------------------------------------------------------	
    use CONSTANTES
    use Input
    implicit none
!Variable de ENTRADA/SALIDA
    type(Compound),pointer,intent(inout)::puntero
!Variables INTERNAS
    integer::nrcar,nisom
	integer::ident(8) !compt(10,2),UNIF(10,2),
	character*35::nombre(8),formula(8) 
	character*37::nomco3,nomco4
!SENTENCIAS

	    nrcar=0
        if (ipareq.eq.2) then
		    call checkdatos (puntero%Formula,nisom,ident,nombre,formula)! checkea en prop.pdb
		    
		    if (nisom /= 0) call elegir_CXR (nisom,ident,nombre,formula,nrcar) !Si existe alg�n is�mero, pregunta al usuario cu�l elegir
        end if
        
!-------Carga de propiedades de componente puro
	    if (nrcar.ne.0) then !si existe en la base de datos
		    call pvomegaybp (nrcar,puntero)

	    else !carga manual si NO existe en la base de datos 
            nomco4 = 'component                            '
!		    if (InputProblem01%mop.eq.1) then ! Lectura de la temperatura de ebullicion del CAR
!			    nomco4 = 'component to be recovered:           '
!		    else
!			    nomco3 = 'more  volatile  component:           '
!			    nomco4 = 'less  volatile  component:           '
!		    end if
3010		write (6,"(1x,/,' Give the boiling point (K) of the ',a37,/,64x,'> ',$)") nomco4
		    read (5,*,err=3010) puntero%BoilingPoint

		    if (InputProblem01%mop.eq.2.and.nrcar.eq.0) then ! Lectura de las constantes de Antoine
			    call antoine (nomco4,'A1','A2','A3',puntero)
		    end if
        end if



endsubroutine load_pure_comp_prop

      subroutine antoine (nomco4,c1,c2,c3,puntero)
!c-----------------------------------------------------------------------
!c     Esta subrutina carga las constantes de Antoine para un componente
!c-----------------------------------------------------------------------
      use Input
      implicit real*8 (a-h,o-z)
      type(Compound),pointer,intent(inout)::puntero
      character*1 iopt
      character*2 c1,c2,c3
      character*37 nomco4
!c
      idev = 6
	idevr = 5
!c     call limp (idev)
  70  write (idev,50) nomco4
      write (idev,110) c1,c2,c3
  80  write (idev,60) c1
      read (idevr,*,err=80) puntero%a(1)
  90  write (idev,60) c2
      read (idevr,*,err=90) puntero%a(2)
 100  write (idev,60) c3
      read (idevr,*,err=100) puntero%a(3)
!c     call limp (idev)
 170  write (idev,180) (puntero%a(i),i=1,3) 
      write (idev,970)
      read (idevr,980,err=170) icomp
      if (icomp.eq.1) then
 130     write (idev,120) c1,c2,c3
         read (idevr,510,err=130) iopt
         i1 = index ('123',iopt)
         if (i1.eq.0) then
!c           call limp (idev)
            go to 130
         end if
!c        call limp (idev)
         write (idev,110) c1,c2,c3
         if (iopt.eq.'1') then
 140        write (idev,60) c1
            read (idevr,*,err=140) puntero%a(1)
            go to 170
         else if (iopt.eq.'2') then
 150        write (idev,60) c2
            read (idevr,*,err=150) puntero%a(2)
            go to 170
         else if (iopt.eq.'3') then
 160        write (idev,60) c3
            read (idevr,*,err=160) puntero%a(3)
            go to 170
         end if        
      end if
!c
!c---formatos
  50  format (' Give the Antoine constants for the ',a37,/)
 110  format (11x,'ln(Pv/mmHg) = ',a2,' - ',a2,'/ (T/K + ',a2,')',/)
 180  format (1x,/,1x,'ln(Pv/mmHg) = ',f9.4,' - ',f9.2,'/ (T/K + ',f9.2,')')
  60  format (1x,33x,a2,'  > ',$)
 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 980  format (i2)
 120  format(' ',/,' Which constant do you want to change?',//,34x,&
            a2,': 1',/,34x,a2,': 2',/,34x,a2,': 3',//,42x,'> ',$)
 510  format (a)
      return 
      endsubroutine antoine



!==========================================================================      
	subroutine checkdatos (compt,nisomt,identt,nombret,formulat)
!--------------------------------------------------------------------------	
!	Esta subrutina checkea o confirma la existencia en el
!	banco de datos "Prop.pdb" de uno o mas compuestos cuya 
!	estructura de grupos seg�n UNIFAC sea igual a la que 
!	contiene el vector de entrada "compt".
!	La variable de salida "nisomt" indica el n�mero de
!	diferentes is�meros que satisfacen tal condici�n.
!	Sus n�meros identificatorios y nombres se almacenan en 
!	los vectores "identt" y "nombret" respectivamente.
!--------------------------------------------------------------------------	
	parameter(NA=150)
	integer compt(10,2),UNIF(10,2),identt(8)
	INTEGER I,J,AUX,AUXB
	character*35 nombret(8),formulat(8) 
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
    CHARACTER*8 fst(NA)
    logical isom
	COMMON/NOM/fst
!Sentencias----------------------------------------------------------------    
    nisomt = 0
	DO i=22,1097
	    
		READ(7,100,REC=i,ERR=120) NI,ANAME,FORM1,FORM2,&
                                   (UNIF(j,1),UNIF(j,2),j=1,10)
 100 		FORMAT(I4,A35,A20,A35,20I3)

        call ordenar (UNIF)
		if (ni.eq.0) cycle
		call comparar (UNIF,compt,1,isom)
        if(isom)then
		    nisomt = nisomt + 1
			identt(nisomt) = ni
			nombret(nisomt) = aname
			formulat(nisomt) = form2
		end if	
    enddo
	GOTO 500

120 WRITE(6,*) '*ERROR* EN LA LECTURA DE "PROP.PDB"'

500	RETURN 
	endsubroutine
!==========================================================================
	subroutine elegir_CXR (nisom,ident,nombre,formula,nr)
	integer ident(8)
	character*35 nombre(8),formula(8)
!Sentencias----------------------------------------------------------------	
	write (6,9)
  9	format (/,&
        '  The data bank "PROP.PDB" contains information about the',/,&
     	'  following isomers for the UNIFAC structure you have entered:'&
     	,//,&
     	'  R.N.	Name				   Formula')
    do 10 i=1,nisom
	    write (6,7) ident(i),nombre(i),formula(i)
   7	format (x,i4,4x,a35,a35)
 10 continue
100	write (6,8)
  8	format (/,&
     	'  If you wish, enter a Reference Number and its properties',/,&
     	'  from the data bank will be taken into account.',/,&
     	'  or choose 0 to work with the MOLDES group contribution',/,&
     	'  estimations for physical properties.',/,&
     	60x,'>')
	read (5,*) nr
	if (nr.ne.0.and.nr.ne.ident(1).and.nr.ne.ident(2)&
     	.and.nr.ne.ident(3).and.nr.ne.ident(4).and.nr.ne.ident(5)&
     	.and.nr.ne.ident(6).and.nr.ne.ident(7).and.nr.ne.ident(8)) &
     	goto 100
    RETURN 
	endsubroutine
!==========================================================================
	subroutine pvomegaybp (nr,puntero)
!--------------------------------------------------------------------------	
!	Esta subrutina utiliza Tb, Tc, Pc y W de una sustancia determinada
!	(de la base de datos PROP.PDB) para calcular las constantes A1 y
!	A2 que sean �tiles para estimar la Pv de la sustancia a temperaturas
!	cercanas a su punto de ebullici�n seg�n la ecuaci�n:
!	Ln(Pv)=A1-A2/T 
!--------------------------------------------------------------------------
    use Input
    use SubGrupos
    implicit real*8 (a-h,o-z)
    type(Compound),pointer,intent(inout)::puntero
	integer UNIF(20)
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
!Sentencias----------------------------------------------------------------		

	READ(7,100,REC=NR,ERR=20) NI,ANAME,FORM1,FORM2,&
                                   (UNIF(K),K=1,20),TC,PC,VC,ZC,TF,TB,&
                                   VLIQ,W
    puntero%Name = ANAME
    puntero%TC = TC
    puntero%PC = PC
    puntero%VC = VC
    puntero%BoilingPoint = TB
    puntero%vliq = VLIQ
	pmt=0
	i=1
	do while(puntero%Formula(i,1)/=0)
		PMT = PMT + Obtain_MW(puntero%Formula(i,1))*puntero%Formula(i,2) 
        i=i+1
    enddo 
    puntero%MW = pmt    
    puntero%Dens = 0.001*PMT/VLIQ
100	FORMAT(I4,A35,A20,A35,20I3,8D13.6)
!c	Pv = Pc * 10** ((7/3)*(1+W)*(T-Tc)/T) !Esto era interpolando entre el punto cr�tico y Tr=0.7
!c	Ahora, interpolando entre Tb y Tr=0.7:
    puntero%a(2) = dlog(Pc*10**(-1-W)/101325)/((1/Tb) - (1/(0.7*Tc)))
	puntero%a(1) = alog(760.) + puntero%a(2)/Tb
    GOTO 500

 10	continue
	GOTO 500
 20	WRITE(6,*) '*ERROR* EN LA LECTURA DE "PROP.PDB"'

500	RETURN 
	endsubroutine
   ! ENDMODULE PROP_PDB
    
    
    
    
    
    
!module input_data
!contains

subroutine enter_problem(nfunc)
!-------------------------------------------------
!   Descripci�n
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
!    parameter()
    implicit none
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
!Variables de ENTRADA
    integer,intent(in)::nfunc
!Variables INTERNAS
    integer::j,model !contadores
    integer::nopt
    character::ifile
    integer::nsolut
    
!...Sentencias
!---Apertura del banco de datos y selecci�n
    call ab_ban1(model)	


    if(nfunc==1)then
    elseif(nfunc==2)then !EvaSol
        
        call input_method(nopt) !1-keuboard; 2-file

        if (nopt==1) then !from keyboar
            write (6,"(///' ','How many component to be recovered do you want evaluate?',//,60x,'> ',$)") 
            read(5,*) NSOLUT
        else !from file
112         write (6,*) "Give the file name (maximun 35 characters):"
            read (5,"(a)",err=112) ifile
            open(unit=20,file=ifile,status='old',access='direct', form='formatted',RECL=62)
            Do j=1, 100
                !read  (20,"(20i3)",rec=j,err=113) (idrs(j,i),nyrs(j,i),i=1,10)
                !write (6,"(20i3)") (idrs(j,i),nyrs(j,i),i=1,10)
            enddo
113         nsolut = j-1   
        endif
        
        
        
    endif



endsubroutine enter_problem

subroutine input_method(nopt)
!-------------------------------------------------------------
!   Asks the user the input method  
!   - Variables de salida
!       nopt: Integer. 1-keyboard, 2-input file
!-------------------------------------------------------------
    implicit none
    integer,intent(out)::nopt

!SENTENCES    
    nopt = 0
    do while (nopt<1.or.nopt>2)
        write (6,"(////' ','Component read from:',&
                              //,11x,' keyboard:           :  1',&
                               /,11x,' input file          :  2',&
     	                      //,60x,'> ',$)")
        read (5,*) nopt
    enddo      

endsubroutine

subroutine leer_comp (ipareq,nomcom,comp)
!c-----------------------------------------------------------------------------
!c     Esta subrutina lee por pantalla un componente con estructura UNIFAC.
!c
!c     Variables de entrada:
!c     --------------------
!c                   ipareq: modelo de parametros de equilibrio:
!c                           1: liquido-liquido
!c                           2: liquido-vapor
!c                           3: liquido-vapor a dilucion infinita.
!c                   ngrmax: cantidad de grupos en la tabla de par�metros 
!                            ipareq.
!c                   nomcom: nombre del componente. Variable character *37.
!c       
!c     Variables de salida:
!c     -------------------
!c                     comp: vector de numeros de identificacion y n�meros de 
!                            repetici�n de subgrupos del componente 
!                            seleccionado. 
!c
!c     Observacion:
!c     -----------
!c                 Los   componentes  no  podran  tener mas de ncomp subgrupos
!c                 diferentes.
!c
!c     Las dimensiones de los  vectores  quedan  fijadas  por  los  siguientes
!c     parametros:
!c
    use constantes
    use SubGrupos

!c
!c-----------------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      integer,intent(out)::comp(DiffStructGroups,2)
      character*8 fs(NMG)
      character*37 nomcom
!c
!c---Escritura de los grupos seleccionables para el componente 
      idev = 6
	  idevr = 5

!-----Averigua cantidad de elementos en la lista LSubGroups
      ngrmax = Size_LSubGroups()

 960  call Groups_Present ()
!c
!c---Escritura del ejemplo
      if (ipareq.eq.1) then
		write (idev,770)
      else if(ipareq.eq.2) then
		write (idev,1400)
      else
		write (idev,1410)
      end if
!c
!c---Lectura del componente
      comp(:,:)=0
      write (idev,100) nomcom ! 100  format (' ',a37,2x,'(maximum 10 different subgroups)')
      read (idevr,580,err=960) (comp(i,1),comp(i,2),i=1,DiffStructGroups)
      msol = 0
      do 320 i=1,DiffStructGroups+1 
		if (comp(i,1).eq.0) goto 330
		if ((comp(i,1).gt.ngrmax).or.(comp(i,1).lt.0)) then
			write (idev,930)
			read (idevr,510) return
			go to 960
		else
			msol = msol + 1
		end if
 300	if (comp(i,2).le.0) then
3020        write (idev,920) Obtain_SubGroup_Name(comp(i,1))
            read (idevr,*,err=3020) comp(i,2)
		    goto 300
         end if
 320  continue
 330  continue
      if (msol.gt.DiffStructGroups) then
         write (idev,110)
         read (idevr,510) return
         go to 960
      end if
!c
!c---- formatos
 100  format (' ',a37,2x,'(maximum 10 different subgroups)')
 110  format (' ','*Error*: You chose a component with more than ','10 different subgroups.',//,' Please <ret>.',$)
 770  format (' ',/,' ','Example:   Propanoic Acid  (CH3)1(CH2)1','(COOH)1',/,' ','1,1,2,1,23,1')
1400  format (' ',/,' ','Example:   Propanoic Acid  (CH3)1(CH2)1','(COOH)1',/,' ','1,1,2,1,43,1')
1410  format (' ',/,' ','Example:   Ethanol (CH3)1(CH2)1(OH)1',/,' ','1,1,2,1,15,1')
 930  format (' ','You chose a nonexistent group.  Please <ret>.',$)
 510  format (a)
 580  format (22i3)
 920  format (' ',/,' How many ',a8,' groups are there in your',' component?  > ',$)
      return

endsubroutine leer_comp
      
subroutine ingresar_solventes (numsol,mst,nmst,nicomp)
    use Input
    use SubGrupos
    use CONSTANTES
    !use prop_pdb,only:check_inter_param
	parameter (NMAX=60000,Nscm=10)
	integer::comp(DiffStructGroups,2)
	dimension ids(nsel),nys(nsel),MST(NMAX,NSCM,2),nicomp(nmax)
    integer grupo1(NMG),nmst(nmax),maingr(NMG),grupos(nsel),kgrup(nsel)
    CHARACTER*8 FS(NMG)
    character*37 nomco
    character*1 sigue
	logical prob
	common/int/ifin2,kgrup,maingr
    COMMON/gru/GRUPO1,NGRMAX
	COMMON/NOM/FS
!SENTENCIAS
    idev = 6
	idevr = 5
    write(idev,*)'  How many solvent structures do you want ','to evaluate?'
    write(idev,694)
	read(idevr,*)numsol
	nomco = 'SOLVENT:                     '
    nums=numsol
	i1=0
    do 11 i=1,nums
		write(idev,*)'Solvent structure number ',i
		write(idev,541)
 1960	call leer_comp (ipareq,nomco,comp)
!c
!c---Ingreso del Solvente al VGP
		call asignar_grupos_principales (ids,msol,grupos)
		call armar_grupos_distintos (grupos,msol,kgrup,ifin2,ifinS)
!c
!c---Control de parametros de interaccion para el Solvente
		PROB=.FALSE.
		if (ifinS.gt.ifin2) then
		    write (idev,1910)
		    if (model /= 3)then
			    prob = check_inter_param (ipareq,1,kgrup,ifin2+1)
		    else
		        prob = .False.
		    endif
		end if
		if (prob) then
		    numsol=numsol-1
     		write (idev,1000)
		    WRITE (6,695)
		    READ (5,696) SIGUE
			goto 11
		else
!c
!c---------Confirmacion del Solvente
3030			write (idev,1210) nomco,(fs(ids(l)),nys(l),l=1,msol)
			write (idev,970)
			read (idevr,980,err=3030) icompa
			if (icomp.eq.1) then
			   go to 1960
			end if
		i1=i1+1
		end if
		nmst(i1)=msol
		do 12 j=1,msol
			mst(i1,j,1) = ids(j)
			mst(i1,j,2) = nys(j)
  12		continue
		write (idev,693)
		write (idev,694) 
		read (idevr,*) nicomp(i1)
  11  continue	
!c	formatos
 541  format (//,' Give the group composition of the solvent: ',&
     		' id1,ny1,id2,ny2,etc.',/,&
             10x,'id1,id2: subgroup identification number',/,&
             10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
 693	format (' ','Give the structure type of the solvent you have ',&
     		'entered:',//,' ',     &
             10x,'-Aromatic                           : 1',//&
             10x,'-Single substance groups            : 2',//&
             10x,'-Aliphatic                          : 3',//&
             10x,'-Cyclic                             : 4',///&
     	      51x,' > ',$)
 694  FORMAT (' ',//,70X,'>',$)
 695  FORMAT (' ',//,70X,'ret>',$)
 696	format (a)
1910  format (' ','*** Checking interaction parameters ***',/)
1000  format (' ',/,' ','* No interaction parameters for some of the ',&
             'subgroups',/,'   already entered.')
 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 980  format (i2)
1210  format (' ',a37,10(a8,i2))
      RETURN
      endsubroutine ingresar_solventes 
      
subroutine ingresar_componentes (NC,prob)
    use Input
    use SubGrupos
    use CONSTANTES
  !  use prop_pdb, only: check_inter_param
	parameter (Nscm=10)
      implicit real*8 (A-H,O-Z)
      integer::comp(DiffStructGroups,2)
      COMMON/US/MS(NCOM,NSCM,2),nms(NCOM)
	common/int/ifin2,kgrup,maingr
      COMMON/gru/GRUPO1,NGRMAX
	COMMON/NOM/FS
	dimension ids(nsel),nys(nsel)
      integer grupo1(NMG),maingr(NMG),grupos(nsel),kgrup(nsel)
      CHARACTER*8 FS(NMG)
      character*37 nomco
      character*1 sigue
	logical prob
	ifin2=0
	if(nc.eq.-2) then
	nc=2
	goto 31
	end if
      write(6,*)'  Give number of components in the mixture> '
      write(6,694)
	read(5,*)NC
 31	nomco = 'COMPONENT:                   '
      nums=NC
	i1=0
      ngrmax = Size_LSubGroups()
	do 10 i=1,ngrmax
		grupo1(i) = i
 		maingr(i) = mainsg (i,ipareq)
 10	continue
      do 11 i=1,nums
		write(6,*)'Component structure number ',i
		write(6,541)
 1960		call leer_comp (ipareq,nomco,comp)

!---Ingreso del Componente al VGP
		call asignar_grupos_principales (ids,msol,grupos)
		call armar_grupos_distintos (grupos,msol,kgrup,ifin2,ifinS)

!---Control de parametros de interaccion para el Componente
		PROB=.FALSE.
		if (ifinS.gt.ifin2) then
		    write (6,1910)
		    if (model /= 3)then
			    prob = check_inter_param (ipareq,1,kgrup,ifin2+1)
			else
			    prob = .False.
			endif
		end if
		if (prob) then
		    NC=NC-1
     		write (6,1000)
		    WRITE (6,695)
		    READ (5,696) SIGUE
			goto 11
		else

!---------Confirmacion del Componente
3030			write (6,1210) nomco,(fs(ids(l)),nys(l),l=1,msol)
			write (6,970)
			read (5,980,err=3030) icomp
			if (icomp.eq.1) then
			   go to 1960
			end if
		i1=i1+1
		end if
		nms(i1)=msol
		do 12 j=1,msol
			ms(i1,j,1) = ids(j)
			ms(i1,j,2) = nys(j)
  12		continue
	ifin2=ifins
  11	continue	
!	formatos
 541  format (//,' Give the group composition of the component: ',&
      		' id1,ny1,id2,ny2,etc.',/,&
              10x,'id1,id2: subgroup identification number',/,&
              10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
 694  FORMAT (' ',//,70X,'>',$)
 695  FORMAT (' ',//,70X,'ret>',$)
 696	format (a)
1910  format (' ','*** Checking interaction parameters ***',/)
1000  format (' ',/,' ','* No interaction parameters for some of the ',&
             'subgroups',/,'   already entered.')
 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 980  format (i2)
1210  format (' ',a37,10(a8,i2))
      RETURN
endsubroutine

SUBROUTINE Model_S(m)
    use Input
    COMMON/AS/ASOC
    LOGICAL ASOC
    INTEGER M
      
    IDEV=6
    ASOC=.FALSE.
    write (idev,"(//,' Choose the model:',//,&
             15x,'UNIFAC  : 1',/,&
             15x,'A-UNIFAC: 2',/,&
             15x,'GC-EOS  : 3')")  
    M=0
10  do while (M<1 .or. M>3)
        write (idev,"(1x,/,60x,'> ',$)") 
        read (5,"(I1)",err=10) M
    enddo
    	
    if(M==2)ASOC=.TRUE.
    

endsubroutine
      

subroutine rec_inf (nomb,MolecularDesign,ptr_FMSs,arch)
!c-----------------------------------------------------------------------
!c   Esta subrutina recupera la informacion de un archivo de datos '.MDI.'
!c   
!c   Datos de entrada
!c   ----------------
!c           nombre: nombre del archivo de datos a recuperar.
!c-----------------------------------------------------------------------
!c
    use GRPS
    use Input
    use StructuresDesign
    use PropertiesData
   
    implicit none
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
!Variables de ENTRADA
    character*35,intent(in)::nomb    
    logical,intent(in)::MolecularDesign
!Variables de SALIDA
    type(FinalStructure),pointer,intent(out)::ptr_FMSs 
    logical,intent(out)::arch
!Variables internas
    integer::i,j,k,lenfile,ierror,qLimits,SolventsNumber,nousari
    integer,dimension(NCOM,DiffStructGroups,2)::MS
    real*8::sclli,select,solv,pmma,dtc2s,tazeo,x1azeo,top !sacar estas variables cuando se actualice el c�digo para mop==2
    character*1::sigue
    logical::nousarl
    
    type(Compound),pointer::ptr_pcp  
    type(CalculatedProperties),pointer::ListProperties
    type(CalculatedPerformance),pointer::ListPerformance      
!Commons
    common/NGX/NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1 
    integer::NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1   
    COMMON/PUNSUB/npunt,ngrup,ncant
    integer::ncant 
    integer,dimension(NMG)::NPUNT,NGRUP
    COMMON/PUNGRU/npint,nintt,numint
    integer::numint
    integer,dimension(nint)::npint,nintt    

!---Apertura de archivo
    arch=.True.
    lenfile = index(nomb,' ')-1
    open(unit=1,file=nomb(1:lenfile)//'.mdi',STATUS='OLD',FORM='FORMATTED',IOSTAT=IERROR)
    if(ierror/=0)then
        write(6,*) '*****',nomb(1:lenfile)//'.mdi',': file not found *****'
        write(6,"(' ',//,70X,'<RET>',$)")
        read(5,"(A)") SIGUE
        arch = .FALSE.
        return
    endif    
    
!---    
    read (1,"(A70)") InputProblem01%ProblemTitle !l�nea 1
    read (1,"(3I3)") ipareq,InputProblem01%mop,npepis !l�nea 2
    call ab_ban1(model)  !Apertura del banco de datos BADAUN 
    call SubGroups_Characterisation()  
    
!---Datos de la mezcla  
    allocate(InputProblem01%MixtureInput)
    nullify(InputProblem01%MixtureInput%Solutes)    
    nullify(InputProblem01%MixtureInput%PCR)    
   !Solutos
    read(1,"(i2)") InputProblem01%MixtureInput%SolutesNumber !l�nea 3
    do i=1,InputProblem01%MixtureInput%SolutesNumber
        nullify(ptr_pcp)
        allocate(ptr_pcp)
        read (1,"(20I3)",err=100) (ptr_pcp%formula(J,1),ptr_pcp%formula(J,2),J=1,DiffStructGroups) !l�nea 4
100     read (1,"(F9.2)",err=101) ptr_pcp%BoilingPoint !l�nea 5
101     call incorporate_compound(ptr_pcp,InputProblem01%MixtureInput%Solutes)
        !Carga de subgrupos de solutos en NPUNT y NPINT
        k=1
        do while(ptr_pcp%Formula(k,1)/=0)
            call CR_PUNTF(ptr_pcp%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
            k=k+1
        enddo
   
    enddo   
     
   !Ppal comp en el refinado
    nullify(ptr_pcp)
    allocate(ptr_pcp)
    read (1,"(20I3)",err=102) (ptr_pcp%formula(J,1),ptr_pcp%formula(J,2),J=1,DiffStructGroups) ! l�nea 6
102 call incorporate_compound(ptr_pcp,InputProblem01%MixtureInput%PCR)
    !Carga de subgrupos de solutos en NPUNT y NPINT
    k=1
    do while(ptr_pcp%Formula(k,1)/=0)
        call CR_PUNTF(ptr_pcp%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
        k=k+1
    enddo
    ms(2,:,:)=ptr_pcp%formula
   !Operation temperature
    read(1,*)InputProblem01%T
    
!---Datos para el dise�o molecular    
    if(MolecularDesign)then
        read(1,"(20I3)")MDV,MSV,MSV1,MJK,InputProblem01%kilout,MDV2 !l�nea 7
        !		MDV = Cantidad de subgrupos "intermedios" seleccionados por el usuario.
        ! 		MSV	= Cantidad de subgrupos "terminales" seleccionados por el usuario.
        !		MSV1 = MSV que no comparten un mismo grupo con subgrupos de valencia dual.
        !		MJK = Ten�a relaci�n con la modificaci�n de las propiedades de combinaci�n.
        !		KILOUT = Type of output
        !       MDV2 = cantidad de subgrupos intermedios alif�ticos seleccionados para ifam=1
        MDV1=MDV-MDV2
        !       MDV1 = cantidad de subgrupos arom�ticos seleccionados por el usuario    
        
   !    Carca de grupos
   !    Grupos intermedios
        NGDV(:)=0
        if(MDV > 0)then
	        read(1,"(11I3)") NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2 !l�nea 8
	    	read(1,"(20I3)") (NGDV(I),I=1,MDV) !l�nea 9 NGDV = Vector de tama�o MDV que contiene los identificadores de subgr. interm. selecc.
            k=1
            do while(NGDV(k)/=0)
                call CR_PUNTF(NGDV(k),NPUNT,NGRUP,NPINT,NINTT)
                k=k+1
            enddo     
        endif
   !    Grupos terminales
        NGSV(:)=0
        if(MSV > 0)then
	        read(1,"(3I3)") NGK1,NGK1J2,NGM1 !l�nea 10
	    	read(1,"(20I3)") (NGSV(I),I=1,MSV) !l�nea 11 NGSV = Vector de tama�o MSV que contiene los identificadores de subgr. termin. selecc.
            k=1
            do while(NGSV(k)/=0)
                call CR_PUNTF(NGSV(k),NPUNT,NGRUP,NPINT,NINTT)
                k=k+1
            enddo          
        endif  
   !    Grupos terminales que no comparten un mismo grupo ppal con subgrupos de valencia dual.    
        NGSV1(:)=0
        if(MSV1 > 0)then
            read(1,"(20I3)") (NGSV1(I),I=1,MSV1) !NGSV1 = Vector de tama�o MSV1 que contiene los n� de subgr. term. selecc.                           
            k=1
            do while(NGSV1(k)/=0)
                call CR_PUNTF(NGSV1(k),NPUNT,NGRUP,NPINT,NINTT)
                k=k+1
            enddo          
        endif
        
    !---Boundaries
        if(InputProblem01%mop == 1)then !liquid-liquid extraction
            read(1,*) qLimits
            do i=1,qLimits
                read(1,*)lugar(i,1),lugar(i,2)
                if(lugar(i,2)==0) then
                    read(1,*)limits(lugar(i,1))%RunValue
                elseif(lugar(i,2)==1)then
                    read(1,*)limits(lugar(i,1))%LowerBound    
                elseif(lugar(i,2)==2)then        
                    read(1,*)limits(lugar(i,1))%UpperBound   
                endif        
            enddo        
            !variables que saqu� SCLLI,SELECT,SOLV,SLSUP1,SLSUPL,DIST,PMMA
	        read (1,*)	DENS2,IFAM,NSOL,IS,NALAR,(IALAR(I),I=1,NALAR)
                                           			
        
        else !if InputProblem01%mop == 2
	        read(1,*)	SCLLI,SELECT,SOLV,PMMA,DTC2S,BP1,DENS2,TAZEO,&
        				x1azeo, A1,A2,A3,B1,B2,B3,IFAM,NSOL,IS,NALAR,&
        				(IALAR(I),I=1,NALAR)
	    endif   
        !
        !	Para MOP = 1
        ! 
        !		SCLLI : Minimum selectivity of an intermediate structure
        !		SELECT: Minimum selectivity of a final structure (wt)       
        !		SOLV  : Minimum solvent power required (wt)                 
        !		SLSUP1: Solvent loss of an intermediate structure (wt %)      
        !		SLSUPL: Solvent loss of a final structure (wt %)              
        !		DIST  : Minimum distribution coefficient required (wt %)      
        !		PMMA  : Maximum molecular weight of the final structure (wt)
        !		Top	  : Operation temperature (K)
        !		BP1	  : Boiling point (K) of the component to be recovered		
        !		DENS2 : Density of raffinate (gr/ml)
        !		IFAM  : Family of solvents to be generated:
        !					Aromatic solvents						: 1
        !					Single substance groups					: 2
        !					No-aromatics with up to 8 groups 
        !					in the intermediate structure			: 3
        !					Mixed aromatic-aliphatic structures
        !					in which the aromatic  part  is the
        !					structure   (ACH)5(ACCH2)1			or
        !								(ACH)4(ACCL)1(ACCH2)1   or
        !								(ACH)4(ACNH2)1(ACCH2)1  or
        !								(ACH)4(ACNO2)1(ACCH2)1		: 4			
        !		Nsol  : Maximum number of solvents to be listed                 
        !		Is	  : Number of solvents per page
        !		IALAR : Vector que (cuando ifam=4) puede contener los n�meros
        !				1 al 4 indicando que los siguientes subgrupos tienen
        !				par�metros de interacci�n con todos los subgrupos del 
        !				CAR y el CPR, y tambi�n con ACH y ACCH2.
        !					ACH y ACCH2	: 1
        !					ACCl		: 2
        !					ACNH2		: 3
        !					ACNO2		: 4
        !		NALAR : Cantidad de elementos en IALAR (max.4)
        !	Para MOP = 2
        ! 
        !       SCLLI : Minimum volatility obtained with an intermediate str.
        !		SELECT: Minimum volatility obtained with a final structucture
        !		SOLV  : Minimum solvent power required                  
        !		Pmma  : Maximum molecular weight of the final structure 
        !         DTC2S : Minimun boiling temperature difference between  
        !				solvent and less volatile component (K)
        !     	Top   : 1.0
        !         BP1	  :	Boiling point (K) of the less volatile component           
        !		TAZEO :	Azeotropic temperature (K)
        !		x1azeo: Azeotropic composition (molar fraction of the less 
        !				volatile component)
        !		A1,A2,A3: Antoine constants for the less volatile component (CAR)
        !		B1,B2,B3: Antoine constants for the more volatile component	(CPR)
        !	    (When nrcar or nrcpr .ne. 0 the program use two coef. calculated from W,Tb,Tc, and Pc)
        !		DENS2,IFAM,NSOL,IS,NALAR e IALAR idem MOP = 1
        !
!        if(MJK==1)then
!            read(1,"(20I3)")NGM,(NGM1V(I),I=1,NGM)
!            do 5 I=1,NGM
!                read(1,"(5I3)")MM(NPUNT(NGM1V(I))),MJ(NPUNT(NGM1V(I))),&
!                             MK(NPUNT(NGM1V(I))),MI(NPUNT(NGM1V(I))),MH(NPUNT(NGM1V(I)))
!5           continue
!        endif  
    else !not MolecularDesign
        !Solventes
        read(1,"(i2)") SolventsNumber
        do i=1,SolventsNumber
            nullify(ptr_FMSs)
            allocate(ptr_FMSs)
            read (1,"(20I3)") (ptr_FMSs%formula(J,1),ptr_FMSs%formula(J,2),J=1,DiffStructGroups)
            read (1,"(I3)") ptr_FMSs%tsolv
!!Carga de subgrupos del SOLVENTE en NPUNT y NPINT     
!            k=1
!            do while (ptr_FMSs%Formula(k,1)/=0)
!                call CR_PUNTF(ptr_FMSs%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
!                k=k+1
!            enddo
!            call Evaluate_Pure(.True.,ptr_FMSs%Formula,InputProblem01%T,ptr_FMSs%tsolv,nousari,nousarl,ListProperties)
!            ms(3,:,:)=ptr_FMSs%Formula
            
!            ptr_pcp => InputProblem01%MixtureInput%Solutes
!            nullify(ListPerformance)            
!            do while(associated(ptr_pcp))
!                ms(1,:,:) = ptr_pcp%Formula                
!                call Evaluate_Mixture(MolecularDesign,.True.,ms,NCOM,InputProblem01%T,ptr_pcp%BoilingPoint,&
!                    ptr_FMSs%tsolv,nousari,nousari,nousarl,ListPerformance,nousari)
!                    ptr_pcp => ptr_pcp%next
!            enddo
!            call Incorporate_Structure (ptr_FMSs%Formula,ptr_FMSs%tsolv,FMSs,ListProperties,ListPerformance)
        enddo     
    endif !(MolecularDesign)    
    
    close (unit=1)    
    return 
endsubroutine rec_inf
      
subroutine seleccion_grupos (IMPRIM,KGRUP,IFIN2,NGR1,NGR2,NGR3,NGR4,GruposSel,family)
      
    USE GRPS
    use input
    use CONSTANTES
    use SubGrupos
 !   use prop_pdb, only:check_inter_param
    PARAMETER (NG=70,INS=20,INS1=10,NAA=4)
    
    IMPLICIT real*8 (A-H,O-Z)
    
!---Commons      
    COMMON/NGX/NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,&
    		    NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1
    integer::NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,&
             NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1	

    common/PUNSUB/NPUNT(NMG),NGRUP(NMG),NUM
    integer::NPUNT,NGRUP,NUM
    common/PUNGRU/NPINT(NG),NINTT(NG),NUMINT
    integer::NPINT,NINTT,NUMINT
!---Variables
    CHARACTER*8 FS(NMG)
    CHARACTER*21 family(8)
    integer::GruposSel(NSEL)
    DIMENSION MAINGR(NMG)      
    DIMENSION isv(NMG),isd(NMG),idv(NMG),i3v(NMG),&
    i4v(NMG),icyc(NMG),iar(NMG),igrp(NMG),ngr1(ins1),ngr2(ins1),&
    ngr3(ins1),ngr4(ins1),ngr5(ins1),MGR(NSEL),NGR(NSEL),ngru1(nsel),&
    NGRU2(NSEL),KGRUP(NSEL)
    integer::mdum(nsel),GRUPOS(NSEL),lgrup(nsel),mgrup(nsel)
    LOGICAL PROB,ESC,entro,inter

    entro=.False.
    idev = 6
	idevr = 5
    call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,imasv,imasd,imadv,ima3v,ima4v,imacy,imaar,imagr)
!   isv:    subgrupos terminales de valencia simple.
!   isd:    subgrupos terminales de valencia simple sin subgrupos de valencia  dual  en  el mismo grupo principal.
!   idv:    subgrupos de valencia dual para compuestos alifaticos.
!   i3v:    subgrupos tri valentes para compuestos alifaticos.
!   i4v:    subgrupos tetra valentes para compuestos alifaticos.
!   icyc:   subgrupos de valencia dual para compuestos ciclicos.
!   iar:    subgrupos para compuestos aromaticos.
!   igrp:   subgrupos que son compuestos moleculares.
!   ima*:   N�mero de grupos en el vector *
                                                            
1320 GruposSel(:)=0

!==================================
!   Intermediate groups
!==================================     
    if (ifam.eq.1) then
!----------------------------------
!   Aromatic groups
!----------------------------------
        call sel_gru_fam (ipareq,family(ifam),ifam,iar,imaar,&
                          imprim,kgrup,ifin2,mdum,ndum,entro,mgr,MDV,&
                          ngru1,mgru1,prob) 
	    do i=1,MDV
	      call insert_group (mgr(i),GruposSel,NSEL)
	      !call CR_PUNTF(mgr(i),NPUNT,NGRUP,NPINT,NINTT)
	    enddo
	    MDV1=MDV
	  
		NGJ2=0
		NGK2=0
		NGK2J2=0
		NGK1J3=0
		call sel_gru_fam (ipareq,family(ifam),ifam,idv,imadv,&
                          imprim,kgrup,ifin2,mdum,ndum,entro,ngr,MDV2,&
                          ngru1,mgru1,prob)
!c
		do i=1,MDV2
			call car_combprop (ngr(i),ipareq,mm,mj,mk,mi,mh)
			if (mj.eq.0) then
				NGK2=NGK2+1
				ngr1(NGK2)=ngr(i)
			else if (mk.eq.2) then
				NGK2J2=NGK2J2+1
				ngr2(NGK2J2)=ngr(i)
			else if (mk.eq.1) then
				NGK1J3=NGK1J3+1
				ngr3(NGK1J3)=ngr(i)
			else
				NGJ2=NGJ2+1
				ngr4(NGJ2)=ngr(i)
	 		end if
	 		!call CR_PUNTF(ngr(i),NPUNT,NGRUP,NPINT,NINTT)
	 		call insert_group (ngr(i),GruposSel,NSEL)
		end do
		if (NGK2.ne.0) then
			do i=1,NGK2
				mgr(MDV+i)=ngr1(i)
			end do
		end if
		MDV=MDV+NGK2
		if (NGK2J2.ne.0) then
			do i=1,NGK2J2
				mgr(MDV+i)=ngr2(i)
			end do
		end if
		MDV=MDV+NGK2J2
		if (NGK1J3.ne.0) then
			do i=1,NGK1J3
				mgr(MDV+i)=ngr3(i)
			end do
		end if
		MDV=MDV+NGK1J3
		if (NGJ2.ne.0) then
			do i=1,NGJ2
				mgr(MDV+i)=ngr4(i)
			end do
		end if
		MDV=MDV+NGJ2
    else if (ifam.eq.2) then
!----------------------------------
!   Molecular groups
!----------------------------------
		call sel_gru_fam (ipareq,family(ifam),ifam,igrp,&
                          imagr,imprim,kgrup,ifin2,mdum,ndum,entro,mgr,&
                          MDV,ngru1,mgru1,prob)
        do i=1,MDV
           ! call CR_PUNTF(mgr(i),NPUNT,NGRUP,NPINT,NINTT)
            call insert_group (mgr(i),GruposSel,NSEL)
        enddo
    else if ((ifam.ge.3).and.(ifam.le.5)) then
!----------------------------------
!   Alifatic groups
!----------------------------------

!------Grupos tetravalentes
		NGJ4=0
		NGK4=0
		call sel_gru_fam (ipareq,family(8),ifam,i4v,ima4v,&
                          imprim,kgrup,ifin2,mdum,ndum,entro,ngr,igr4,&
                          ngru1,mgru1,prob)

		do i=1,igr4
			call car_combprop (ngr(i),ipareq,mm,mj,mk,mi,mh)
			if (mj.eq.4) then
				NGJ4=NGJ4+1
				ngr2(NGJ4)=ngr(i)
			else
				NGK4=NGK4+1
				ngr1(NGK4)=ngr(i)
			end if
			!call CR_PUNTF(ngr(i),NPUNT,NGRUP,NPINT,NINTT)
			call insert_group (ngr(i),GruposSel,NSEL)
		end do
		if (NGK4.ne.0) then
			do i=1,NGK4
				mgr(i)=ngr1(i)
			end do
		end if
		if (NGJ4.ne.0) then
			do i=1,NGJ4
				mgr(NGK4+i)=ngr2(i)
			end do
		end if
		MDV=NGK4+NGJ4
		
!------Grupos trivalentes
		NGJ3=0
		NGK3=0
		NGK3J2=0
		NGK2J3=0
		NGK1J4=0
		call sel_gru_fam (ipareq,family(7),ifam,i3v,ima3v,&
                          imprim,kgrup,ifin2,mdum,ndum,entro,ngr,igr3,&
                          ngru1,mgru1,prob)
!c
		do i=1,igr3
			call car_combprop (ngr(i),ipareq,mm,mj,mk,mi,mh)
			if (mj.eq.0) then
				NGK3=NGK3+1
				ngr1(NGK3)=ngr(i)
			else if (mk.eq.3) then !chequear si puede realmente pasar por ac�
				NGK3J2=NGK3J2+1
				ngr2(NGK3J2)=ngr(i)
			else if (mk.eq.2) then
				NGK2J3=NGK2J3+1
				ngr3(NGK2J3)=ngr(i)
			else if (mk.eq.1) then
				NGK1J4=NGK1J4+1
				ngr4(NGK1J4)=ngr(i)
			else
				NGJ3=NGJ3+1
				ngr5(NGJ3)=ngr(i)
			end if
			!call CR_PUNTF(ngr(i),NPUNT,NGRUP,NPINT,NINTT)
			call insert_group (ngr(i),GruposSel,NSEL)
		end do
		if (NGK3.ne.0) then
			do i=1,NGK3
				mgr(MDV+i)=ngr1(i)
			end do
		end if
		MDV=MDV+NGK3
		if (NGK3J2.ne.0) then
			do i=1,NGK3J2
				mgr(MDV+i)=ngr2(i)
			end do
		end if
		MDV=MDV+NGK3J2
		if (NGK2J3.ne.0) then
			do i=1,NGK2J3
				mgr(MDV+i)=ngr3(i)
			end do
		end if
		MDV=MDV+NGK2J3
		if (NGK1J4.ne.0) then
			do i=1,NGK1J4
				mgr(MDV+i)=ngr4(i)
			end do
		end if
		MDV=MDV+NGK1J4
		if (NGJ3.ne.0) then
			do i=1,NGJ3
				mgr(MDV+i)=ngr5(i)
			end do
		end if
		MDV=MDV+NGJ3
!c--------Grupos divalentes
		NGJ2=0
		NGK2=0
		NGK2J2=0
		NGK1J3=0
		call sel_gru_fam (ipareq,family(ifam),ifam,idv,imadv,&
                          imprim,kgrup,ifin2,mdum,ndum,entro,ngr,MDV2,&
                          ngru1,mgru1,prob)

		do i=1,MDV2
			call car_combprop (ngr(i),ipareq,mm,mj,mk,mi,mh)
			if (mj.eq.0) then
				NGK2=NGK2+1
				ngr1(NGK2)=ngr(i)
			else if (mk.eq.2) then
				NGK2J2=NGK2J2+1
				ngr2(NGK2J2)=ngr(i)
			else if (mk.eq.1) then
				NGK1J3=NGK1J3+1
				ngr3(NGK1J3)=ngr(i)
			else
				NGJ2=NGJ2+1
				ngr4(NGJ2)=ngr(i)
	 		end if
!	 		call CR_PUNTF(ngr(i),NPUNT,NGRUP,NPINT,NINTT)
	 		call insert_group (ngr(i),GruposSel,NSEL)
		end do
		if (NGK2.ne.0) then
			do i=1,NGK2
				mgr(MDV+i)=ngr1(i)
			end do
		end if
		MDV=MDV+NGK2
		if (NGK2J2.ne.0) then
			do i=1,NGK2J2
				mgr(MDV+i)=ngr2(i)
			end do
		end if
		MDV=MDV+NGK2J2
		if (NGK1J3.ne.0) then
			do i=1,NGK1J3
				mgr(MDV+i)=ngr3(i)
			end do
		end if
		MDV=MDV+NGK1J3
		if (NGJ2.ne.0) then
			do i=1,NGJ2
				mgr(MDV+i)=ngr4(i)
			end do
		end if
		MDV=MDV+NGJ2
!c
    else  ! Ahora no puede pasar mas por ac�
!----------------------------------
!   Cyclic groups
!---------------------------------- 
		call sel_gru_fam (ipareq,family(ifam),ifam,icyc,&
                         imacy,imprim,kgrup,ifin2,mdum,ndum,entro,mgr,&
                         MDV,ngru1,mgru1,prob)
        do i=1,MDV
           ! call CR_PUNTF(mgr(i),NPUNT,NGRUP,NPINT,NINTT)
            call insert_group (mgr(i),GruposSel,NSEL)
        enddo
    end if

 
    if (prob) then
		!go to 1200
    end if

	nalar = 0
	do 1340 i=1,naa
	    ialar(i) = 0
1340 continue

!---Ingreso de los grupos seleccionados al VGP
    if (MDV.gt.ins) then
!	    call limp (idev)         
		write (idev,460)
		CALL PAUSA 
		MDV = ins
    end if

    ifin01 = ifin2

    call asignar_grupos_principales (mgr,MDV,grupos)
    call armar_grupos_distintos (grupos,MDV,kgrup,ifin01,ifin3)

!---Confirmacion de los grupos seleccionados
3060 write (idev,1310) (Obtain_SubGroup_Name(mgr(i)),i=1,MDV)
    write (idev,970)
    read (idevr,980,err=3060) icomp
    if (icomp.eq.1) then
        go to 1320
    end if
    i=1
    do while(GruposSel(i)/=0)
        call CR_PUNTF(GruposSel(i),NPUNT,NGRUP,NPINT,NINTT)
        i=i+1
    enddo
    
!==============================================
!   Single valence groups
!==============================================
    if (ifam.ne.2) then
	    write (idev,1610) family(6)
		if (imprim.eq.2) then
			write (2,*)
			write (2,1610) family(6)
		end if
1330	GruposSel(:)=0
        NGK1=0
		NGK1J2=0
		NGM1=0
		MSV=0
	    call sel_gru_fam (ipareq,family(6),ifam,isv,imasv,&
                          imprim,kgrup,ifin2,mgr,MDV,entro,mval,nival,&
                          ngru2,mgru2,prob)
!
		if (nival.gt.ins1) then
!			call limp (idev)         
			write (idev,"(' ',/,' * Warning: Only 10 groups will be considered ',&
             'because memory space reasons. Please select only 10 groups')")
			CALL PAUSA 
			go to 1330
!			MSV = ins1
		end if
!
		do i=1,nival
			call car_combprop (mval(i),ipareq,mm,mj,mk,mi,mh)
			if (mk.eq.0) then
				NGM1=1
				ngr3(1)=mval(i)
			else
				if (mj.eq.0) then
					NGK1=NGK1+1
					ngr1(NGK1)=mval(i)
				else
					NGK1J2=NGK1J2+1
					ngr2(NGK1J2)=mval(i)
				end if
			end if
			call insert_group (mval(i),GruposSel,NSEL)
		end do
		
		if (NGK1.ne.0) then
			do i=1,NGK1
				mval(i)=ngr1(i)
			end do
		end if
		MSV=NGK1
		if (NGK1J2.ne.0) then
			do i=1,NGK1J2
				mval(MSV+i)=ngr2(i)
			end do
		end if
		MSV=MSV+NGK1J2
		if (NGM1.ne.0) then
			mval(MSV+1)=ngr3(1)
		end if
		MSV=MSV+NGM1
!C
		if (prob) then
	!		go to 1200
		end if

!c
!c------Ingreso de los subgrupos de valencia uno al VGP
		ifin0 = ifin3
		call asignar_grupos_principales (mval,MSV,grupos)
		call armar_grupos_distintos (grupos,MSV,kgrup,ifin0,ifin4)
!c
!c------Confirmacion de los grupos seleccionados
        
3070    write (idev,1310) (Obtain_SubGroup_Name(mval(i)),i=1,MSV)
		write (idev,970)
		read (idevr,980,err=3070) icomp
		if (icomp.eq.1) then
			go to 1330
		end if
        i=1
        do while(GruposSel(i)/=0)
            call CR_PUNTF(GruposSel(i),NPUNT,NGRUP,NPINT,NINTT)
            i=i+1
        enddo	
!c
!c------Control de los parametros de interaccion de los conjuntos de 
!c------grupos seleccionados para asegurar la generacion de solventes
        
		!if ((ifam.eq.3).or.(ifam.eq.4)) then
		if ((ifam.eq.3)) then
			esc = .false.
			call asignar_grupos_principales (mgr,MDV,grupos)
			call armar_grupos_distintos (grupos,MDV,lgrup,0,ifinl)
			call asignar_grupos_principales (mval,MSV,grupos)
			call armar_grupos_distintos (grupos,MSV,mgrup,0,ifinm)
			call inter2 (ipareq,0,lgrup,ifinl,mgrup,ifinm,inter)
			if (.not.inter) then
				write (idev,1540)
				esc = .true.
			end if
!			if (ifam.eq.4) then
!				call asignar_grupos_principales (maingr,iaux,laux,grupos)
!				call armar_grupos_distintos (grupos,laux,maux,0,ifaux)
!				call inter1 (ipareq,0,maux,ifaux,mgrup,
!     *                         ifinm,igrupi,ifarsv)
!				if (ifarsv.eq.0) then
!					write (idev,1510)
!					esc = .true.
!				end if
!				call inter1 (ipareq,0,maux,ifaux,lgrup,
!     *                         ifinl,igrupi,ifardv)
!				if (ifardv.eq.0) then
!					write (idev,1500)
!					esc = .true.
!				end if
!				if (ifarsv.eq.0) then
!					write (idev,1550)
!				end if
!			else if (ifam.eq.3) then
			if (.not.inter) then
		    	write (idev,1550)
			end if
!			end if
			if (esc) then
 1840				write (idev,1520)
				write (idev,1530)
				read (idevr,980,err=1840) icomp
				if (icomp.eq.1) then
			!		go to 1200
				else if (icomp.eq.2) then
					entro=.true.
					go to 1320
				else if (icomp.eq.3) then
					go to 1830
				else if (icomp.eq.4) then
					go to 10000
				else 
					go to 1840
				end if
1830				continue
			end if
		end if
!c
!c------Control de los parametros de interaccion entre los dos conjuntos de
!c------grupos ingresados (todos con todos)
!c        call limp (idev)         
1850		write (idev,1620)
		if (imprim.eq.2) then
			write (2,*)
			write (2,5010)
		end if
		if (imprim.eq.0) then
			imprim = 1
		end if
	  if(model /= 3)then	
	    if (check_inter_param (ipareq,imprim,kgrup,ifin01+1)) then
			CALL PAUSA 
1820    	write (idev,1810)
			read (idevr,980,err=1820) icomp
			if (icomp.eq.1) then
				go to 1850
			else if (icomp.eq.2) then
			!	go to 1200
			else if (icomp.eq.3) then
				entro=.true.
				go to 1320
			else if (icomp.eq.4) then
				go to 1860
			else if (icomp.eq.5) then
				go to 10000
			else 
				go to 1820
			end if
		end if
    endif

!c------Asignacion de los grupos con subgrupos sin valencia dual
 1860		msv1 = 0
		do 640 i=1,MSV
			do 650 j=1,imasd
				if (mval(i).eq.isd(j)) then
					msv1 = msv1 + 1
					msd(msv1) = isd(j)
				end if
 650			continue
 640		continue
      else
		MSV = 0
		msv1 = 0
      end if
      NGSV(1:MSV)=MVAL(1:MSV)


      NGDV(1:MDV)=MGR(1:MDV)

!-----FORMATOS
 461  format (' ',/,' * Warning: Only 10 groups will be considered ',&
             'because memory space reasons. Please select only 10 groups')
 460  FORMAT(' ',/,' * Warning: Only 20 groups will be considered because memory space reasons.')
 620  format (20i3)
 626  format (3i3)
 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 980  format (i2)
1310  format (' ',/,' ','You chose this groups: ',6a8,/,' ',26x,6a8,/,' ',26x,6a8)
1520  format (' ',/,'   What do you want to do?')
1530  format ('                 Change the family of solvents : 1',/,&
              '                 Change the chosen groups      : 2',/,&
              '                 Continue                      : 3',/,&
              '                 Exit program                  : 4',/,&
              63x,'> ',$)
1540  format (' ',/,' * No interaction parameters between the dual val',&
             'ence and the',/,'   single valence groups already entered.')
1550  format (' ',/,' * Imposibility of forming the family of solvents chosen.')
1610  format (' ',/,' *** ',a21,' preliminary selection ***')
1620  format (' ',/,' * Looking for chosen groups without',' interaction parameters *')
1810  format (' ',/,' The former pair of groups don''t have interactio',&
             'n parameters.',&
                 /,' If you decide to continue, the solvents ',&
             'generated with these pairs',&
                 /,' of groups won''t be taken into account.',//,&
                   ' What do you want to do?',/,&
             '                 Look the pairs of groups again: 1',/,&
             '                 Change the family of solvents : 2',/,&
             '                 Change the chosen groups      : 3',/,&
             '                 Continue                      : 4',/,&
             '                 Exit program                  : 5',/,&
             63x,'> ',$)
5010  format (' * Selected pairs of groups without interaction',' parameters.')
10000 RETURN
      endsubroutine SELECCION_GRUPOS
!==========================================================================
SUBROUTINE CHECK_PROP()
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------   
  use CONSTANTES
  use Subgrupos
!  use prop_pdb, only:checkdatos
  IMPLICIT real*8 (A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
  INTEGER GRUPO1(NMG),comp(10,2),idr(nsel),nyr(nsel),ident(8),model
  CHARACTER*37 nomco
  CHARACTER*35 nombre(8),formula(8)
  CHARACTER*8 fs(NMG)
  CHARACTER car*3
  COMMON/NOM/FS
!Sentencias----------------------------------------------------------------         
    call ab_ban1(model)					 

400 IPAREQ=4
	nomco = 'COMPOUND:                    '
    ngrmax = Size_LSubGroups()
    DO J=1,ngrmax
	  call carac (j,ipareq,fs(j),car)
	  grupo1(j) = j
	END DO
	write(6,541)
1960 call leer_comp (ipareq,nomco,comp)
!c-----formaci�n de comp (para el compuesto a checkear) a partir de idr y nyr 
	k=1
	do 5 i=1,19,2
	  j =	(i+1)/2
        comp(k,1)=idr(j)
	  comp(k,2)=nyr(j)
        k=k+1
   5	continue
      call checkdatos (comp,nisom,ident,nombre,formula)
	if (nisom.ne.0) then
		write (6,9)
	      do 10 i=1,nisom
			write (6,7) ident(i),nombre(i),formula(i)
  10		  continue
 		  write (6,8)
		read (5,*) nr
		if (nr.ne.0) call VER_COM(NR)
      else
		write(6,6)
		call pausa
	end if
	CALL ci_ban1
	close(unit=7)
      
!-----FORMATOS
   6    format (/,' There is no isomers corresponding to the UNIFAC',&
     	/,' structure you have entered in the data bank "PROP.PDB".')
   7	format (x,i4,4x,a35,a35) 
   8	format (/,&
     	'  If you wish, enter a Reference Number to view ',/,&
     	'  its properties from the data bank ',/,&
     	'  (or choose 0 to continue)',/,&
     	60x,'>')     
   9	format (/,&
     	'  The data bank "PROP.PDB" contains information about the',/,&
     	'  following isomers for the UNIFAC structure you have entered:'&
     	,//,&
     	'  R.N.	Name				   Formula')
 541  format (//,' Give the group composition of the compound: ',&
     		' id1,ny1,id2,ny2,etc.',/,&
             10x,'id1,id2: subgroup identification number',/,&
             10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
 999  return
endsubroutine check_prop
!==========================================================================
SUBROUTINE CHECK_DATABASE()
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------      
 ! use prop_pdb,only:check_intrcn,CHECK_GRUPOS
!Variables INTERNAS
  INTEGER OP
!Sentencias----------------------------------------------------------------      
      WRITE (6,100)
 10   WRITE (6,110)
      READ (5,200,ERR=10) OP
      IF (OP<1.OR.OP>3) GOTO 10
      
      IF (OP.EQ.1) THEN
        CALL CHECK_PROP
      ELSE IF (OP.EQ.2) THEN
        CALL CHECK_INTRCN
      ELSE IF (OP.EQ.3) THEN
        CALL CHECK_GRUPOS
      ENDIF
!-----FORMATOS  
 100  FORMAT(///,' Select the database:',///,&
           6X,'1 - Prop.pdb',//,&
           6X,'2 - Intrcn.mds',//,&
           6X,'3 - Gruposram.mds',//) 
 110  FORMAT(/,60X,'> ',$)  
 200  FORMAT(I1) 
      return
      endsubroutine check_database
      
! ====================================================================
subroutine enter_mixture()
!------------------------------------------------------------------------------------
!   Descripci�n
!   - Variables de entrada
!       none
!   - Variables de salida
!       MOP:    separation operation (LL extraction(1); extractive destilation (2))
!       tem:    operation temperature
!------------------------------------------------------------------------------------
   ! use PROP_PDB, only:check_inter_param,load_pure_comp_prop
    use CONSTANTES
    use SubGrupos
    use Input

    implicit none
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
!Variables INTERNAS	
    type(Compound),pointer::recorreSolutes
    type(Compound),pointer::ptr_pcp
    integer::return
    integer,dimension(NMG)::MGSolutes,MGMixture
    integer::i, idev, idevr
    integer::icant, ngrmax, icomp, isolu
    character*37::nomco
!FUNCIONES
    integer::count_groups    
!------------------------------------------------------------------------------------
    
!SENTENCIAS  

!---Inicializaci�n de variables  
    idev = 6
	idevr = 5
    icant = 0
    allocate(InputProblem01%MixtureInput)
    nullify(InputProblem01%MixtureInput%PCR)    
!---Selecci�n de modelo termodin�mico y apertura del banco de datos    
    call ab_ban1(model) 

!---Selecci�n de operaci�n de separaci�n
    call operation_separation(InputProblem01%mop)

!---Selecci�n de tabla de par�metros UNIFAC
    call tabla_parametros(ipareq) 

!---temperatura de operaci�n
    if(InputProblem01%mop.EQ.1)then !Liquid-Liquid extraction
		write(idev,"(1X,/,' Give operation temperature (K)',19X,' > ',$)")
		read(idevr,*) InputProblem01%T
		if(model == 3)then
		    write(idev,"(1X,/,' Give operation pressure (atm)',20X,' > ',$)")
		    read(idevr,*) InputProblem01%P
        endif
    elseif(InputProblem01%mop.eq.2)then !Extractive destilation
		InputProblem01%T = 1.0
    endif

960 InputProblem01%MixtureInput%SolutesNumber = 1
    if(InputProblem01%inv)then 
        nomco = 'Solvent:                             '
    else      
        if (InputProblem01%mop.eq.1) then !Separation operation
         nomco = 'Component to be recovered '
         InputProblem01%MixtureInput%SolutesNumber = 0
         do while(InputProblem01%MixtureInput%SolutesNumber == 0)
            write(idev,"(//,' How many solutes there are in the mixture?:',6X,' > ',$)")
            read(5,"(I2)") InputProblem01%MixtureInput%SolutesNumber
         enddo
        else
         nomco = 'Less  volatile  component:           '
        end if
    endif
    
    nullify(InputProblem01%MixtureInput%Solutes)
    InputProblem01%MainGroups(:) = 0
    InputProblem01%SubGroups(:) = 0
    MGSolutes(:) = 0
        
!---carga de solutos (o solvente si inv == .True.)
    ReadComp:&
    do isolu=1,InputProblem01%MixtureInput%SolutesNumber 
        nullify(ptr_pcp)
        allocate(ptr_pcp)
    !---Carga de componente    
        write (idev,"(//,' Give the group composition of the following ',&
                'components: id1,ny1,id2,ny2,etc.',/,' Where: ',/,&
                10x,'id1,id2: subgroup identification number',/,&
                10x,'ny1,ny2: number of subgroups id1,id2,etc',/)") 

        if(.not.InputProblem01%inv .and. InputProblem01%mop==1)nomco(36:)= '#'//char(48+isolu)//':'
        call leer_comp (ipareq,nomco,ptr_pcp%Formula)
    !---Control de parametros de interaccion         
        i=1 
        do while(ptr_pcp%Formula(i,1)/=0)
            call Load_Main_Groups(ptr_pcp%Formula(i,1),MGSolutes) 
            i=i+1
        enddo
        if (count_groups(ptr_pcp%Formula(:,1)) > 1) then 
            write (idev,"(' ','*** Checking interaction parameters ***',/)")
            if(model /= 3)then
            if (check_inter_param(ipareq,1,MGSolutes,1)) then
                write (idev,"(' ',/,' ','* No interaction parameters for some of the ','subgroups',/,'   already entered.')")
                write (idev,1010)
                read (idevr,510) return
                go to 960 ! Vuelve arriba a preguntar cantidad de solutos 
            end if
            endif
        end if

    !---Confirmacion del componente ingresado
 3000   write (idev,1210) nomco,(Obtain_SubGroup_Name(ptr_pcp%Formula(i,1)),&
                          ptr_pcp%Formula(i,2),i=1,count_groups(ptr_pcp%Formula(:,1))) 
        write (idev,970)
        read (idevr,"(i2)",err=3000) icomp
        
        if (icomp.eq.1) goto 960 ! Si se responde no OK, vuelve arriba a preguntar cantidad de solutos 

    !---B�squeda de propiedades de componente puro en prop.pdb
        call load_pure_comp_prop (ptr_pcp)
        
    !--Incorporporaci�n del componente
		call incorporate_compound(ptr_pcp,InputProblem01%MixtureInput%Solutes)  
    enddo & 
    ReadComp
    
!---Lectura del Componente Principal del Refinado (CPR).    
!   icant == 1 cuando el CPR ya fue ingresado, pero se volvi� 
!   a cambiar el o los CAR porque no estaban todos los par�metros de interacci�n
    nullify(ptr_pcp)
    allocate(ptr_pcp)
    MGMixture(:) = 0    
    if(icant == 0) then 
1960    nomco = 'Principal component in the raffinate:'    
        if (InputProblem01%mop.eq.2) nomco = 'More  volatile  component:           '
        call leer_comp (ipareq,nomco,ptr_pcp%Formula)
    endif

!---Chequea par�metros de interacci�n entre grupos de solutos y PCR
    i=1 
    MGMixture(:) = MGSolutes(:)
    do while(ptr_pcp%Formula(i,1)/=0)
        call Load_Main_Groups(ptr_pcp%Formula(i,1),MGMixture) 
        i=i+1
    enddo        
    if(model /= 3)then
      if (check_inter_param(ipareq,1,MGMixture,1) ) then
3020    write (idev,"(' ',/,' ','* No interaction parameters for some of the ','subgroups',/,'   already entered.')")
	    write (idev,1020)
	    read (idevr,980,err=3020) icomp
	    if (icomp.eq.1) then
		    ICANT=1
		    go to 960
	    else
		    ICANT=0
		    go to 1960
	    end if
      end if
    endif

!---Confirmacion del CPR
3030 write (idev,1210) nomco,(Obtain_SubGroup_Name(ptr_pcp%Formula(i,1)),&
                      ptr_pcp%Formula(i,2),i=1,count_groups(ptr_pcp%Formula(:,1)))
    write (idev,970)
    read (idevr,980,err=3030) icomp
    if (icomp.eq.1) then
        go to 1960
    end if
    call load_pure_comp_prop (ptr_pcp)      

!---Lectura de las composiciones azeotropicas 
  	if (InputProblem01%mop.eq.2) call azeotrope_input()

!---Incorporar PCR a la mezcla
    call incorporate_compound(ptr_pcp,InputProblem01%MixtureInput%PCR) 

!---incorporar los grupos de la mezcla
    InputProblem01%MixtureInput%SubGroups(:) = 0
    InputProblem01%MixtureInput%MainGroups(:) = 0

    !Solutos
    ptr_pcp => InputProblem01%Mixtureinput%Solutes
    do while(associated(ptr_pcp))
        i = 1
        do while(ptr_pcp%Formula(i,1) /= 0)
        !---Incorpora subgrupo 
            call insert_group (ptr_pcp%Formula(i,1), InputProblem01%MixtureInput%SubGroups, &
                                size(InputProblem01%MixtureInput%SubGroups))
        !---Incorpora grupo principal    
            call insert_group (Obtain_MainGroup_Number(ptr_pcp%Formula(i,1)), &
                               InputProblem01%MixtureInput%MainGroups, size(InputProblem01%MixtureInput%MainGroups))
            i = i + 1
        enddo
        ptr_pcp => ptr_pcp%next
    enddo

    !PCR
    ptr_pcp => InputProblem01%Mixtureinput%PCR
    i = 1
    do while(ptr_pcp%Formula(i,1) /= 0)
    !---Incorpora subgrupo 
        call insert_group (ptr_pcp%Formula(i,1), InputProblem01%MixtureInput%SubGroups, size(InputProblem01%MixtureInput%SubGroups))
    !---Incorpora grupo principal    
        call insert_group (Obtain_MainGroup_Number(ptr_pcp%Formula(i,1)), &
                           InputProblem01%MixtureInput%MainGroups, size(InputProblem01%MixtureInput%MainGroups))
        i = i + 1
    enddo    
          
!	formatos
2900  format (' ',//,70x,'<ret>',$)
 510  format (a)
1010  format ('   Change the info. about the component:  <ret>. ',$)
1020  format ('   Do you want to change the components?',/,&
             '               Component to be recovered: 1',/,&
             '               Principal component in the raffinate: 2',/,60x,'> ',$)                                                  
 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 980  format (i2)
 981  FORMAT (I1)
1210  format (' ',a37,10(a8,i2))
 535  format (1x,/,' Give the density (gr/ml) of the ',a37,/,64x,'> ',$)
      return 

endsubroutine enter_mixture 


subroutine operation_separation(mop)
!--------------------------------------------------
!   Selecci�n de la operaci�n de separaci�n
!   mop = 
!           1- Liquid-liquid extraction
!           2- Etractive distillation
!--------------------------------------------------
    implicit none
!Variables de SALIDA
    integer,intent(out)::mop    

    mop=0
60  do while(mop /= 1 .and. mop /= 2)
        write (6,"(' Choose the separation operation:',//,&
               29x,'Liquid-Liquid Extraction: 1',/,&
               29x,'Extractive  Distillation: 2')")
        write (6,"(1x,/,60x,'> ',$)")
        read (5,*,err=60)mop
    enddo

endsubroutine operation_separation


subroutine azeotrope_input()
!------------------------------------------------------
!   Pregunta al usuario si el sistema forma alg�n aze�-
!   tropo y pide los valores de temperatura y composi-
!   ci�n del mismo.
!------------------------------------------------------
    use Input
    implicit none 
!Variables INTERNAS
    integer::i1, icomp, idev, idevr
    real*8::Tazeo, Xazeo
    character*1::opt,iopt

!SENTENCIAS    
    idev = 6
    idevr = 5

80	write (idev,"(1x,/,' Does the system form an azeotrope?',//,51x,'Yes: 1',/,51x,'No : 2',//,60x,'> ',$)")
	read (idevr,"(a)",err=80) opt
	i1 = index ('12',opt)
	if (i1.eq.0) goto 80
	
	if (opt.eq.'1') then !si forma aze�tropo
110     write (idev,"(1x,/,' Give the azeotropic temperature (K)              ',' > ',$)") !90  format 
		read (idevr,*,err=110)	tazeo
 
120 	write (idev,"(1x,/,' Give the azeotropic composition',/,' (molar fraction of the less volatile component)   > ',$)") ! 100  format 
		read (idevr,*,err=120) xazeo
 
160	    write (idev,"(' Azeotropic temperature:    ',f8.3,' K')") tazeo
		write (idev,"(' Azeotropic composition:    ',f8.3)") xazeo
		write (idev,970)
		read (idevr,"(i2)",err=160) icomp
		
		if (icomp.eq.1) then
170	        write (idev,"(' ',/,' Which variable do you want to change?',//,21x,'tazeo: 1',/,21x,'xazeo: 2',//,42x,'> ',$)") 
			read (idevr,"(a)",err=160) iopt
			i1 = index ('12',iopt)
			if (i1.eq.0) goto 170

			if (iopt.eq.'1') then
180			    write (idev,"(1x,/,' Give the azeotropic temperature (K)              ',' > ',$)") 
				read (idevr,*,err=180) tazeo
				go to 160
			else if (iopt.eq.'2') then
190			    write (idev,"(1x,/,' Give the azeotropic composition',/,' (molar fraction of the less volatile component)   > ',$)") 
				read (idevr,*,err=190) xazeo
				go to 160
			end if
		end if
	else
		tazeo = 1.0
		xazeo = 0.0
	end if 
 
970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 
 
 endsubroutine azeotrope_input
 
endmodule

!=======================================================
subroutine Load_Main_Groups(subNumber,MainGroups)
!-------------------------------------------------------
!   Esta subrutina carga el subgrupo subNumber en el 
!   vector MainGroups si �ste a�n no est� presente en 
!   el mismo
!-------------------------------------------------------

use CONSTANTES
use SubGrupos
implicit none
!Variables de ENTRADA
    integer,intent(in)::subNumber
!Variables de ENTRADA/SALIDA
    integer,intent(inout)::MainGroups(NMG)
!Variables internas
    integer::i,j,mainNumber
    logical::presente

!SENTENCIAS

    mainNumber = Obtain_MainGroup_Number(subNumber)
    call buscaras(mainNumber,MainGroups,NMG,presente)
    if(.not.presente)then
        j=1
        do while(MainGroups(j)/=0)
            j=j+1
        enddo
        MainGroups(j)=mainNumber
    endif
    i=i+1
endsubroutine Load_Main_Groups

subroutine insert_group (grupo,vector,n)
!-----------------------------------------------------------------
!   Esta subrutina agrega el grupo "grupo" al vector "vector" de
!   dimensi�n "n"
!-----------------------------------------------------------------
implicit none
!Variables de ENTRADA
integer,intent(in)::grupo,n
!Variables de ENTRADA/SALIDA
integer,dimension(n),intent(inout)::vector
!Variables INTERNAS
integer::i
logical::logic

!SENTENCIAS
    call buscaras(grupo,vector,n,logic)
    if(.not.logic)then
        i=1
        do while(vector(i)/=0)!busca primera posici�n vac�a
            i=i+1
        enddo
        vector(i)=grupo
    endif
    
endsubroutine

subroutine write_input_file (MolecularDesign)
!-----------------------------------------------------------------
!   Esta subrutina 
!   
!-----------------------------------------------------------------    
    use Input
    use StructuresDesign
    use GRPS
    use PropertiesData
    use blockdatas  
    implicit none
!Variables de ENTRADA
    logical,intent(in)::MolecularDesign
!Variables INTERNAS
    integer::i,j,qLimits,idev
    integer,dimension(NSEL):: mgr,mst
    real*8::temp,tazeo,sclli,dtc2s
    type(Compound),pointer::recorreSolutes
    type(FinalStructure),pointer::recorreSolvents
    common/NGX/NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1
    integer::NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2,NGK1,NGK1J2,NGM1

    idev=1
!SENTENCIAS
    write (idev,"(a70,' !Problem Title')") InputProblem01%ProblemTitle !l�nea 1
    write (idev,"(3i3,' !ipareq, mop, npepis')") ipareq,InputProblem01%mop,npepis !l�nea 2

!---Datos de la mezcla    
   !Solutos
    write (idev,"(I2,' !Solutes number')") InputProblem01%MixtureInput%SolutesNumber !l�nea 3
    recorreSolutes => InputProblem01%MixtureInput%Solutes
    do while(associated(recorreSolutes))
        write(idev,"(20i3,' !Solute formula')") (recorreSolutes%Formula(j,1),&
                    recorreSolutes%Formula(j,2),j=1,recorreSolutes%GroupsNumber) !linea 4
        write(idev,"(F9.2,' !Boilling point')")recorreSolutes%BoilingPoint !l�nea 5
        recorreSolutes => recorreSolutes%next
    enddo
    
   !Ppal comp en el refinado    
    write(idev,"(20i3,' !PCR formula')") (InputProblem01%MixtureInput%PCR%Formula(j,1),&
                                InputProblem01%MixtureInput%PCR%Formula(j,2)&
                    ,j=1,InputProblem01%MixtureInput%PCR%GroupsNumber) !l�nea 6
   !Operation temperature
    write(idev,*) InputProblem01%T 
!---Datos para el dise�o molecular    
    if(MolecularDesign)then
        write (idev,"(6i3,' !MDV,MSV,msv1,mjk,idato(3),MDV2')") MDV,MSV,msv1,mjk,idato(3),MDV2 !linea 7
        InputProblem01%kilout=IDATO(3)
   !    Carca de grupos
   !    Grupos intermedios
        if (MDV.ne.0) then
	        write (idev,"(11i3,' !NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2')")&
            	        NGK4,NGJ4,NGK3,NGK3J2,NGK2J3,NGK1J4,NGJ3,NGK2,NGK2J2,NGK1J3,NGJ2 !linea 8
	    write (idev,"(20i3,' !(mgr(i),i=1,MDV)')") (NGDV(i),i=1,MDV) !linea 9
        end if
   !    Grupos terminales    
        if(MSV.ne.0)then
	        write (idev,"(3i3,' !NGK1,NGK1J2,NGM1')") NGK1,NGK1J2,NGM1 !l�nea 10
	    write (idev,"(20i3,' !(mval(i),i=1,MSV)')") (NGSV(i),i=1,MSV) !l�nea 11
        end if
        !Grupos terminales que no comparten un mismo grupo ppal con subgrupos de valencia dual.    
        if (msv1.ne.0) then
	        write (idev,"(20i3,' !(mval(i),i=1,MSV)')") (msd(i),i=1,msv1) !l�nea 12
	        NGSV1(1:MSV1)=MSD(1:MSV1) 
        end if
        
!-------Boundaries    
        if (InputProblem01%mop.eq.1) then
            i=1
            do while (lugar(i,1)/=0)
                i=i+1         
            enddo
            qLimits=i-1
          
            write(idev,"(i3,' !qLimits')")qLimits
      
            do i=1,qLimits
                if(lugar(i,2)==0) then
                    write(idev,"(2i3,' !lugar(i,1), lugar(i,2)')")lugar(i,1), lugar(i,2)
                    write(idev,"(F9.2,' !RunValue')")limits(lugar(i,1))%RunValue
                elseif(lugar(i,2)==1)then
                    write(idev,"(2i3,' !lugar(i,1), lugar(i,2)')")lugar(i,1), lugar(i,2)
                    write(idev,"(F9.2,' !LowerBound')")limits(lugar(i,1))%LowerBound    
                elseif(lugar(i,2)==2)then
                    write(idev,"(2i3,' !lugar(i,1), lugar(i,2)')")lugar(i,1), lugar(i,2)
                    write(idev,"(F9.2,' !UpperBound')")limits(lugar(i,1))%UpperBound
                endif   
            enddo
	        write (idev,*) DENS2,ifam,&
	                    idato(1),idato(2),nalar,(ialar(i),i=1,nalar) !l�nea 13
            NSOL=IDATO(1)
            IS=IDATO(2)
        else !if (InputProblem01%mop.eq.2)
	        write (1,*) dato2(1),dato2(2),dato2(3),dato2(4),dato2(5),&
        			bp1,DENS2,tazeo,xazeo,a1,a2,a3,b1,b2,b3,ifam,idato(1),&
        			idato(2),nalar,(ialar(i),i=1,nalar) !l�nea 13
            SCLLI=DATO2(1) 
            DTC2S=DATO2(5)
            NSOL=IDATO(1)
            IS=IDATO(2)
        endif
        
        if (mjk.eq.1) then
            write (idev,"(20i3),' !ngm,(ngm1V(i),i=1,ngm)'") ngm,(ngm1V(i),i=1,ngm) !l�nea 14
	        do i=1,ngm
	    	    write (idev,"(20i3,' !mst(i),JST1(i),KST1(i),IST1(i),hst(i)')") mst(i),JST1(i),KST1(i),IST1(i),hst(i) !l�nea 15
            enddo
        end if
    else ! not MolecularDesign
        write(idev,"(I2,' !Solvents number')") Size_LFMSs(FMSs)
        recorreSolvents => FMSs
        do while(associated(recorreSolvents))
            write(idev,"(20i3,' !Solute formula')") (recorreSolvents%Formula(j,1),&
            recorreSolvents%Formula(j,2),j=1,recorreSolvents%GroupsNumber)
            write(idev,"(i3,' !tIPO DE SOLVENTE')") recorreSolvents%tsolv
            recorreSolvents => recorreSolvents%next
        enddo        
        
    endif !(MolecularDesign)    
    endsubroutine write_input_file
    
    
    
    
    
    
    
    
    
    

    
    
    