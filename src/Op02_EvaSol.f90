subroutine Evasol()
    use CONSTANTES
    use input_data, only:enter_mixture
    implicit none


!VARIABLES INTERNAS
    integer::i,j,s
    integer::idev,idevr,nopt,nsolut,mop,msol,mraf
    integer::idr(NMG),nyr(NMG),idf(NMG),nyf(NMG),idrs(NMG,NSEL),nyrs(NMG,NSEL)
    real*8::tem,bp1,a1,a2,a3,dens2,b1,b2,b3,tazeo,xazeo
    character*35::ifile

!...SENTENCIAS
    idev = 6
    idevr= 5
!...LECTURA DE SOLUTOS
    nopt = 0
    do while (nopt<1.or.nopt>2)
        write (idev,"(////' ','Component read from:',&
                              //,11x,' keyboard:           :  1',&
                               /,11x,' input file          :  2',&
     	                      //,60x,'> ',$)")
        read (idevr,*) nopt
    enddo
!    if (nopt==1) then !por teclado
!        write (idev,"(///' ','How many component to be recovered do you want evaluate?',//,60x,'> ',$)") 
!        read(idevr,*) NSOLUT
!    else !por archivo
! 112    write (idev,*) "Give the file name (maximun 35 characters):"
!        read (idevr,"(a)",err=112) IFILE
!        open(unit=20,file=IFILE,status='old',access='direct', form='formatted',RECL=62)
!        Do j=1, 100
!            read  (20,"(20i3)",rec=j,err=113) (idrs(j,i),nyrs(j,i),i=1,10)
!            write (6,"(20i3)") (idrs(j,i),nyrs(j,i),i=1,10)
!        enddo
! 113    NSOLUT = j-1   
!    endif
  !  MS(:,:,:) = 0
    call enter_mixture()
    


endsubroutine Evasol