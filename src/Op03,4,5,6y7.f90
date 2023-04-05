subroutine distcoef
	!USE PROPI
	use PureProp
	use SubGrupos
	use constantes
	use Input
!	use PROP_PDB, only:checkdatos,elegir_CXR,pvomegaybp
	use input_data,only:Leer_Comp,ingresar_componentes,checkdatos,elegir_CXR,pvomegaybp
	
   PARAMETER(NMODEL=3,NGPM=30,NGA=70,NSCM=10,ins=20)

      implicit real*8(A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
    use input_data, only:Model_S
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
    type(Compound),pointer::ptr_pcp  
  !  external::Size_LSubGroups
      DIMENSION pmg(ins),X1(10),X2(10)
	common/solub/X1,X2
      COMMON/US/MS(NCOM,DiffStructGroups,2),nms(NCOM)
	COMMON/NOM/FS
      dimension idr(nsel),nyr(nsel),nr(2),pm(2),dens(2)
      integer grupo1(NMG),maingr(NMG),h1
	integer comp(20),ident(8)
	character*35 nombre(8),formula(8) 
      character*8 fs(NMG)
      character*37 nomco
      character car*3
      character*17 tabla(nmodel)
      LOGICAL NOACE,PROB
      DIMENSION X(NCOM),ACT(NCOM)

!SENTENCES
      call ab_ban1(model) 							 
2900  CALL TABLA_PARAMETROS(IPAREQ)
     
     
     ngrmax = Size_LSubGroups()
      DO J=1,ngrmax
	   call carac (j,ipareq,fs(j),car)
			grupo1(j) = j
			maingr(j) = mainsg (j,ipareq)
	END DO
	NC=-2
	write(6,*)' First, enter the system',' (principal component of each phase)'
	call ingresar_componentes (NC,prob)
	nrcpr=0
      if (ipareq.eq.2) then
!c	formación de comp (para el CPR) a partir de idf y nyf 
	do k=1,2
        call checkdatos (ms(k,:,:),nisom,ident,nombre,formula)
	  if (nisom.ne.0) then
            call elegir_CXR (nisom,ident,nombre,formula,nrcpr)
	      nr(k)=nrcpr
	      nrcpr=0
        end if
	end do
	end if
!c----------------------------------------
	nomco = 'SOLUTE:                      '
	write(6,541)
	call leer_comp (ipareq,nomco,comp)
		nms(3)=msol
		do 12 j=1,msol
			ms(3,j,1) = idr(j)
			ms(3,j,2) = nyr(j)
  12		continue
      CALL CR_PUNT2 (3,IPAREQ)
	CALL Store_Pr (IPAREQ)
      CALL Store_In ()
	WRITE(6,5556)
	READ(5,*) Top
      tabla(1)='    liquid-liquid'
	tabla(2)='     liquid-vapor'
	tabla(3)=' ifinite dilution'
	call PARIN (NC,NG,6,NOACE,MS)
      if (NOACE) then
       write (6,1490) tabla(ipareq)
          CALL PAUSA 
	 GO TO 99
	end if
	call PARAM (NC,NG,Top)
	call solbin (NG,Top)
	X(1)=X1(1)
	X(2)=X1(2)
	X(3)=0.
	call Gamma(IPAREQ,3,Top,X,ACT)
	ACT1=ACT(3)
	X(1)=X2(1)
	X(2)=X2(2)
	X(3)=0.
	call Gamma(IPAREQ,3,Top,X,ACT)
	ACT2=ACT(3)
!c-----Peso Molecular y Densidad de los 2 principales comp. de las fases--
	do k=1,2
	pm(k)=0
	if (nr(k).ne.0) then
		do 961 i=1,nms(k)
          CALL Leer_Pr (ms(k,I,1),IPAREQ,fs(ms(k,I,1)),M1,J1,K1,I1,H1,&
     					PMG(i),RPAR,QPAR,DTC,DPC,DelV,dvi,dn,DTCA)
		PM(k) = PM(k) + PMG(I)*ms(k,I,2)
 961		continue
		call pvomegaybp (nr(k),ptr_pcp)
		Dens(k)=0.001*PM(k)/ptr_pcp%VLIQ
	else
		write (6,693) k
		write (6,694) 
		read (5,*) icomp
        call Calc_Prop_Pure(MS(K,:,:),icomp,Top,MWT=PM(k),Denst=Dens(k))
	end if
	end do
!c----------------------------------------
	pmA=X1(1)*PM(1)+X1(2)*PM(2)
	pmB=X2(1)*PM(1)+X2(2)*PM(2)
	densA=X1(1)*dens(1)+X1(2)*dens(2)
	densB=X2(1)*dens(1)+X2(2)*dens(2)
      write (1,*) ' Estimation of the partition',' coefficient for: '
	write (1,*) ' Principal component in phase A: '
      write (1,420) (ms(1,i,1),ms(1,i,2),i=1,nms(1))
	write (1,*) ' Principal component in phase B: '
      write (1,420) (ms(2,i,1),ms(2,i,2),i=1,nms(1))
	write (1,*) ' Solute: '
      write (1,420) (ms(3,i,1),ms(3,i,2),i=1,nms(1))
	write (1,*) ' Temperature: ', Top
	write (1,*) ' Unifac table used: ', ipareq
	WRITE(6,*)' Densidad de las fases: '
	write(6,*)' DensA= ',DensA,'  DensB= ',DensB
	WRITE(1,*)' Densidad de las fases: '
	write(1,*)' DensA= ',DensA,'  DensB= ',DensB
	P=(ACT2*densA*pmB)/(ACT1*densB*pmA)
	WRITE(6,*) 'P= ',P
	WRITE(6,*) 'logP= ',dlog10(P)
	WRITE(1,*) 'P= ',P
	WRITE(1,*) 'logP= ',dlog10(P)
	call pausa
  99	close(unit=7)
 100	CALL ci_ban1
 541  format (' Now give the group composition of the solute: ',&
     		' id1,ny1,id2,ny2,etc.',/,&
             10x,'id1,id2: subgroup identification number',/,&
             10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
 420  format (20i3)
1490  format (' ',/,' * No interaction parameters available for ',&
             'this system',/,'   in the ',a17,&
             ' UNIFAC table')
5556  FORMAT (1X,/,' Give the system temperature (K)      ',' > ',$)
 693	format (' ','Give the structure type of the component ',i,&
     		':',//,' ','Aliphatic:',20x,'0',&
     		 /,' ','Aromatic:',21x,'1',&
     		 /,' ','Single group solvent:',9x,'2',&
     		 /,' ','Cyclic:',23x,'3')
 694  FORMAT (' ',//,70X,'>',$)
 	return
	end

subroutine log_POW()
      use input_data,only:Leer_Comp  
      use SubGrupos
      use CONSTANTES
     
!-------------------------------------------------
!   Descripción
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
      implicit real*8(A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
!...Integers
      integer maingr(NMG),grupo1(NMG)
      integer::idr(nsel),nyr(nsel),MS(NCOM,DiffStructGroups,2),comp(DiffStructGroups,2)
!...Characters
      character*3 car
      character*6 nomb
      character*8 fs(NMG)
      character*37 nomco
!...Commons
      common/PAREQ/IPAREQ	
!...Sentencias 
      call ab_ban1(model)
      call TABLA_PARAMETROS(IPAREQ)							 
	nomco = 'COMPOUND:                    '
      ngrmax = Size_LSubGroups()
      do J=1,ngrmax
	    call carac (j,ipareq,fs(j),car)
	    grupo1(j) = j
		maingr(j) = mainsg (j,ipareq)
	  enddo
	  write(6,541)
 541  format (' Give the group composition of the compound: ',&
       		' id1,ny1,id2,ny2,etc.',/,&
              10x,'id1,id2: subgroup identification number',/,&
              10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
 	  call leer_comp (ipareq,nomco,comp)
	  do 12 j=1,msol
	    ms(3,j,1) = idr(j)
		ms(3,j,2) = nyr(j)
  12  continue
      call OW_Partition_Coefficient(POW,MS(3,:,:))
3100  write (6,"(/,' Give output file name (maximum ','6 characters).')")
      read (5,'(a)',err=3100) nomb
      open (unit=1,file=nomb//'.pow',form='formatted')
	  write (1,*) ' ECOFAC 1.0 '
	  write (1,*) ' Case: '
      write (1,"(6x,a6)") nomb
      write (1,*) ' Estimation of n-octanol/water partition',' coefficient for: '
      write (1,"(20i3)") (idr(i),nyr(i),i=1,msol)
	  write (1,*) ' Temperature: 298.15 K '
	  write (1,*) ' Unifac table used: ', ipareq
	  write(1,*) 'Pow= ',Pow
	  write(1,*) 'logPow= ',dlog10(Pow)
	  write(6,*) 'Pow= ',Pow
	  write(6,*) 'logPow= ',dlog10(Pow)
	  write(6,*) ' Find results in the output file .pow'
	  close (unit=1)
	  call pausa      
!      
      endsubroutine     
      
      

      subroutine solubilities
!-------------------------------------------------
!   Descripción
!-------------------------------------------------
      use input_data,only:ingresar_componentes
      use SubGrupos
      use constantes
      parameter(nmodel=3)
      implicit real*8(A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface  
!...Integers
      integer MS(NCOM,DiffStructGroups,2),nms(NCOM)
!...Characters
      character*3 car
      character*6 nomb
      character*8 fs(NMG)
      character*17 tabla(nmodel)
!...Logical
      logical::prob,NOACE
!...Commons
      common/PAREQ/IPAREQ	
      common/NOM/FS 
      
3201  write (6,"(/,' Give output file name (maximum ','6 characters).')")
      read (5,"(a)",err=3201) nomb
      open (unit=1,file=nomb//'.sol',form='formatted')
	  write (1,*) ' ECOFAC 1.0 '
	  write (1,*) ' Case: '
      write (1,"(6x,a6)") nomb
      call ab_ban1(model)
      call TABLA_PARAMETROS(IPAREQ)							 
      ngrmax = Size_LSubGroups()
      do J=1,ngrmax
	    call carac (j,ipareq,fs(j),car)
	  enddo
	  NC=-2
	  call ingresar_componentes (NC,prob)
      call CR_PUNT2 (NC,IPAREQ)
	  call Store_Pr(IPAREQ)
      call Store_In()
	  write(6,"(1X,/,' Give the system temperature (K)      ',' > ',$)")
	  read(5,*) T
      tabla(1)='    liquid-liquid'
	  tabla(2)='     liquid-vapor'
	  tabla(3)=' ifinite dilution'
	  call PARIN (NC,NG,IDEV,NOACE,MS)
      if (NOACE) then
        write (6,"(' ',/,' * No interaction parameters available for ',&
                'this system',/,'   in the ',a17,&
                ' UNIFAC table')") tabla(ipareq)
        CALL PAUSA 
	    GO TO 1000
	  end if
	  call PARAM (NC,NG,T)
      write (1,*) ' Estimation of mutual solubilities for: '
	  write (1,*) ' Principal component in phase A: '
      write (1,"(20i3)") (ms(1,i,1),ms(1,i,2),i=1,nms(1))
	  write (1,*) ' Principal component in phase B: '
      write (1,"(20i3)") (ms(2,i,1),ms(2,i,2),i=1,nms(1))
	  write (1,*) ' Temperature: ', T
	  write (1,*) ' Unifac table used: ', ipareq
	  call solbin (NG,T)
	  write(6,*) ' Find results in the output file .sol'
	  close (unit=1)
	  call ci_ban1
1000  continue	
      end subroutine
      
!
!
      subroutine activity_coefficients
!-------------------------------------------------
!   Descripción
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
      use input_data,only:ingresar_componentes
      use SubGrupos
      use constantes

      implicit real*8(A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface  
!...Integers

!...real*8
      dimension X(NCOM), ACT(NCOM)
!...Characters
      character*3 car
      character*6 nomb
      character*8 fs(NMG)
!...Logical
      logical prob
!...Commons
      common/PAREQ/IPAREQ
      common/US/MS(NCOM,DiffStructGroups,2),nms(NCOM)
      common/NOM/FS
!
!...Sentencias
3202  write (6,"(/,' Give output file name (maximum ','6 characters).')")
      read (5,"(a)",err=3202) nomb
      open (unit=1,file=nomb//'.gam',form='formatted')
	  write (1,*) ' ECOFAC 1.0 '
	  write (1,*) ' Case: '
      write (1,"(6x,a6)") nomb
      write (1,*) ' Estimation of activity coefficients: '
      call ab_ban1(model)						 
      CALL TABLA_PARAMETROS(IPAREQ)
      ngrmax = Size_LSubGroups()
      do J=1,ngrmax
	    call carac (j,ipareq,fs(j),car)
	  enddo
	  call ingresar_componentes (NC,prob)
	  if (prob) goto 888
      CALL CR_PUNT2 (NC,IPAREQ)
	  CALL Store_Pr (IPAREQ)
      CALL Store_In ()
  27  write(6,"(1x,/,' Give the solution composition',/,' (molar fraction of each component)   > ')")
	  do i=1,nc
	    write(6,*)' X',i,'= '
		read(5,*) X(i)
	  end do
	  sumx=0
	  do i=1,nc
	    sumx=sumx+X(i)
	  enddo
	  if (sumx.ne.1) then
	    write(6,*) ' The sum of the molar fractions is not equal to 1'
		goto 27
	  endif
      write(6,"(1X,/,' Give the solution temperature (K)      ',' > ',$)")
	  read(5,*) T
	  call Gamma(IPAREQ,NC,T,X,ACT)
	  write(1,*)'  '
	  do i=1,nc
	    write (1,*) ' Component number: ',i
        write (1,"(20i3)") (ms(i,j,1),ms(i,j,2),j=1,nms(1))
        write(1,*)' X ',i,'= ',X(i)
		write(6,*)' Gamma ',i,'= ',ACT(i)
		write(1,*)' Gamma ',i,'= ',ACT(i)
		write(1,*)'  '
	  end do
	  write (1,*) ' Temperature: ', T
	  write (1,*) ' Unifac table used: ', ipareq
	  write(6,*) ' Find results in the output file .gam'
	  close (unit=1)
      CALL PAUSA 
 888  CALL ci_ban1
      endsubroutine      
      