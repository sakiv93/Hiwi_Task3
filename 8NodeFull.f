cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  UEL FOR PHASE FIELD MODELLING OF FERROELECTRIC MATERIALS OF 1st     c
c                            ORDER                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)


      INCLUDE 'ABA_PARAM.INC'


      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),SVARS(*),
     1 ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),
     2 A(NDOFEL),TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


        double precision, dimension(24,24)::Kuu
        double precision, dimension(24,24)::Kup
        double precision, dimension(24,24)::Kpu
        double precision, dimension(8,8)::Kphi
        double precision, dimension(24,24)::Kppv
        double precision, dimension(24,24)::Kele
        double precision, dimension(8,24)::Kelea
        double precision, dimension(24,24)::Kpups
        double precision, dimension(24,24)::Kpps
        double precision, dimension(24,24)::Kpp
        double precision, dimension(8,24)::Kphip
        double precision, dimension(24,8)::Kpphi
        double precision, dimension(56,56)::Ke
        double precision, dimension(56,56)::Ke2
        double precision, dimension(24,24)::Kup_r
        double precision, dimension(24,24)::Kpps_r
        double precision, dimension(56,56)::Ke_r
        double precision, dimension(56,56)::Ke_r2


        double precision, dimension(6,3)::qmat
        double precision, dimension(6,3)::qmat2
        double precision, dimension(3,3)::bata
        double precision, dimension(3,3)::batam
        double precision, dimension(3,3)::alfamat
        double precision, dimension(3,3)::alfamatr
        double precision, dimension(6,6)::cmat
        double precision, dimension(3,3)::kapmat
        double precision, dimension(9,9)::grdmat
        double precision, dimension(3,3)::xjaci
        double precision, dimension(3,3)::xjacm
        double precision, dimension(4,8)::shp
        double precision               ::djacb
        double precision              ::interval

        double precision, dimension(17)::B
c        double precision, dimension(56)::ue
c        double precision, dimension(56)::ve

        double precision, dimension(6,24)::Bu
        double precision, dimension(3,24)::Nu
        double precision, dimension(3,8)::Bphi
        double precision, dimension(9,24)::Bp
        double precision, dimension(3,8)::wg
        double precision, dimension(3,3)::Tra
        double precision, dimension(6,6)::Trb
        double precision, dimension(9,9)::Trc
        double precision, dimension(3,3)::Tran
        double precision, dimension(3,8)::pve
        double precision, dimension(3,8)::pe
        double precision, dimension(3)::xyz0
        double precision, dimension(3,8)::uelem
        double precision, dimension(3,3)::qmatpS
        double precision, dimension(3,8)::Evor
        double precision                ::vor

        double precision, dimension(6,24)::temp
        double precision, dimension(3,8)::temp_e
        double precision, dimension(9,24)::temp_p
        double precision, dimension(6,24)::temp_pu
        double precision, dimension(3,24)::temp_up
        double precision, dimension(3,24)::temp_pps
        double precision, dimension(3,24)::temp_ppv
        double precision, dimension(3,24)::temp_ele
        double precision, dimension(3,3)::temp_batam
        double precision, dimension(3,24)::temp_pups
        double precision, dimension(6,24)::tempup_r
        double precision, dimension(3,24)::temppps_r


        double precision, dimension(6)::strain
        double precision, dimension(6)::stress
        double precision, dimension(3)::elefld
        double precision, dimension(3)::eledispl
        double precision, dimension(9)::polargrd
        double precision, dimension(24)::mdispl
        double precision, dimension(8)::edispl
        double precision, dimension(24)::pdispl
        double precision, dimension(24)::mforc
        double precision, dimension(8)::eforc
        double precision, dimension(24)::pforc
        double precision, dimension(3)::mut
        double precision, dimension(9)::pol

        double precision :: s
        double precision :: t
        double precision :: y

        integer :: i,j,m,n,ig,ii,jj,ki,pj

        dt = 0.1
        ba = 1
        vor = 0.d0

cc....... create applied vortex field

c      do i=1,8
c        Evor(1,i)=-vor*COORDS(2,i)
c        Evor(2,i)= vor*COORDS(1,i)
c        Evor(3,i)= 0.d0
c      enddo

c       CREATING GAUSS POINTS

        do i=1,4
        wg(3,i  )=-0.57735
        wg(3,i+4)= 0.57735
      enddo

      do i=1,2
        wg(2,i  )=-0.57735
        wg(2,i+2)= 0.57735
        wg(2,i+4)=-0.57735
        wg(2,i+6)= 0.57735
      enddo

        wg(1,1)=-0.57735
        wg(1,2)= 0.57735
        wg(1,3)= 0.57735
        wg(1,4)=-0.57735
        wg(1,5)=-0.57735
        wg(1,6)= 0.57735
        wg(1,7)= 0.57735
        wg(1,8)=-0.57735


c       SETTING THE INITIAL VALUE AS ZERO FOR K MATRIX

        do i=1,24
            do j=1,24
               Kuu(i,j)=0.d0
               Kpu(i,j)=0.d0
               Kup(i,j)=0.d0
               Kup_r(i,j)=0.d0
               Kpp(i,j)=0.d0
               Kpps(i,j)=0.d0
               Kpps_r(i,j)=0.d0
               Kppv(i,j)=0.d0
               Kele(i,j)=0.d0
               Kpups(i,j)=0.d0
            enddo
         enddo

         do i = 1,8
           do j = 1,8
               Kphi(i,j) = 0.d0
            enddo
         enddo

         do i = 1,24
           do j = 1,8
               Kpphi(i,j) = 0.d0
            enddo
         enddo


        do i = 1,8
           do j = 1,24
               Kphip(i,j) = 0.d0
               Kelea(i,j) = 0.d0
            enddo
         enddo


        do i = 1,3
           do j = 1,3
               bata(i,j) = 0.d0
               Tran(i,j) = 0.d0
            enddo
         enddo


         write(6,*) '************** Kuu***************'
         do i=1,24
            do j=1,24
                write(6,*) Kuu(i,j)
            enddo
        enddo

        write(6,*) '***************** Kpu**************'
         do i=1,24
            do j=1,24
                write(6,*) Kpu(i,j)
            enddo
        enddo

        write(6,*) '***************** Kup***************'
         do i=1,24
            do j=1,24
                write(6,*) Kup(i,j)
            enddo
        enddo

                write(6,*) '***************** Kup_r***************'
         do i=1,24
            do j=1,24
                write(6,*) Kup_r(i,j)
            enddo
        enddo

        write(6,*) '****************** Kpp****************'
         do i=1,24
            do j=1,24
                write(6,*) Kpp(i,j)
            enddo
        enddo

        write(6,*) '***************** Kpps****************'
         do i=1,24
            do j=1,24
                write(6,*) Kpps(i,j)
            enddo
        enddo

                write(6,*) '***************** Kpps_r****************'
         do i=1,24
            do j=1,24
                write(6,*) Kpps_r(i,j)
            enddo
        enddo

        write(6,*) '******************* kppv**************'
         do i=1,24
            do j=1,24
                write(6,*) Kppv(i,j)
            enddo
        enddo

        write(6,*) '****************** Kele***************'
         do i=1,24
            do j=1,24
                write(6,*) Kele(i,j)
            enddo
        enddo

        write(6,*) '******************* Kpups*****************'
         do i=1,24
            do j=1,24
                write(6,*) Kpups(i,j)
            enddo
        enddo

                 write(6,*) '******************** Kphi*****************'
         do i=1,8
            do j=1,8
                write(6,*) Kphi(i,j)
            enddo
        enddo

                write(6,*) '******************** Kpphi*****************'
         do i=1,24
            do j=1,8
                write(6,*) Kpphi(i,j)
            enddo
        enddo

        write(6,*) '******************** Kphip*********************'
         do i=1,8
            do j=1,24
                write(6,*) Kphip(i,j)
            enddo
        enddo

        write(6,*) '********************* Kelea********************'
         do i=1,8
            do j=1,24
                write(6,*) Kelea(i,j)
            enddo
        enddo

                 write(6,*) '*********************** Tran**************'
         do i=1,3
            do j=1,3
                write(6,*) Tran(i,j)
            enddo
        enddo

        write(6,*) '*************************** bata*******************'
         do i=1,3
            do j=1,3
                write(6,*) bata(i,j)
            enddo
        enddo



c      INITIALIZING THE INITIAL CONDITION

         do j = 1,3
           do k = 1,8
               pe(j,k) = u(7*(k-1)+4+j)
               pve(j,k) = v(7*(k-1)+4+j)
               uelem(j,k) = u(7*(k-1)+j)
            enddo
         enddo




        do j=1,3
           xyz0(j)=0
            do k=1,8
                xyz0(j)=xyz0(j) + COORDS(j,k)/8.0
            enddo
        enddo

        write(6,*) '******************* xyz********************'
         do i=1,3
                write(6,*) xyz0(i)
            enddo


         write(6,*) '********************** pe************************'
         do i=1,3
            do j=1,8
                write(6,*) pe(i,j)
            enddo
        enddo

        write(6,*) '*********************** pve************************'
         do i=1,3
            do j=1,8
                write(6,*) pve(i,j)
            enddo
        enddo

        write(6,*) '************************ uelem*********************'
         do i=1,3
            do j=1,8
                write(6,*) uelem(i,j)
            enddo
        enddo


            call trans(Tra,Trb,Trc,xyz0)
            call material(cmat,kapmat,B,grdmat,xyz0)


c       PRINTING THE PROPERTIES FROM INPUT FILE

        write(6,*) '********PROPS************'
        write(6,*)  'PROPS 1', PROPS(1)
        write(6,*)  'PROPS 2', PROPS(2)
        write(6,*)  'PROPS 3', PROPS(3)
        write(6,*)  'PROPS 4', PROPS(4)
        write(6,*)  'PROPS 5', PROPS(5)
        write(6,*)  'PROPS 6', PROPS(6)
        write(6,*)  'PROPS 7', PROPS(7)
        write(6,*)  'PROPS 8', PROPS(8)
        write(6,*)  'PROPS 9', PROPS(9)
        write(6,*)  'PROPS 10', PROPS(10)
        write(6,*)  'PROPS 11', PROPS(11)
        write(6,*)  'PROPS 12', PROPS(12)
        write(6,*)  'PROPS 13', PROPS(13)
        write(6,*)  'PROPS 14', PROPS(14)
        write(6,*)  'PROPS 15', PROPS(15)
        write(6,*)  'PROPS 16', PROPS(16)
        write(6,*)  'PROPS 17', PROPS(17)


c       ASSIGNING PROPERTIES FROM INPUT FILE

        do i = 1,17
           B(i) = PROPS(i)
        enddo

        write(6,*) '******************PROPS***************************'
        do i=1,17
            write(6,*) B(i)
        enddo



c       THE SUMMATION OF GAUSSIAN POINTS

        do ig = 1,8

         s = wg(1,ig)
         t = wg(2,ig)
         y = wg(3,ig)

         call sfr(shp,s,t,y,nnode)
         call jacob(shp,djacb,coords,nnode,mcrd,xjacm,xjaci)
         call non_material(pe,qmat,alfamat,B,alfamatr,shp
     &          ,xyz0,uelem,qmatpS)



c        B MATRIX INITIALIZATION

        do i = 1,6
           do j = 1,24
               Bu(i,j)=0
           enddo
         enddo


c       FORMULATING Bu MATRIX

        do j = 1,8
            Bu(1,j*3-2) = shp(1,j)
            Bu(2,j*3-1) = shp(2,j)
            Bu(3,j*3  ) = shp(3,j)
            Bu(4,j*3-1) = shp(3,j)
            Bu(4,j*3  ) = shp(2,j)
            Bu(5,j*3-2) = shp(3,j)
            Bu(5,j*3  ) = shp(1,j)
            Bu(6,j*3-2) = shp(2,j)
            Bu(6,j*3-1) = shp(1,j)
         enddo

        write(6,*) '********************** Bu*************************'
        do i=1,6
            do j=1,24
                write(6,*) Bu(i,j)
            enddo
        enddo

c        FORMULATING Kuu MATRIX

            temp = matmul(cmat,Bu)
            Kuu = Kuu + (djacb) * (matmul(transpose(Bu),temp))
            temp = 0.0

            write(6,*) '*******************Kuu*************************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kuu(i,j)
                enddo
            enddo

c       N MATRIX INITIALIZATION

        do i = 1,3
           do j = 1,24
               Nu(i,j) = 0
            enddo
         enddo

c       FORMULATING N MATRIX

        do j = 1,8
            Nu(1,j*3-2) = shp(4,j)
            Nu(2,j*3-1) = shp(4,j)
            Nu(3,j*3  ) = shp(4,j)
         enddo


          write(6,*) '********************* Nu************************'
        do i=1,3
            do j=1,24
                write(6,*) Nu(i,j)
            enddo
        enddo

c       Kpu MATRIX

            temp_pu = matmul(qmat,Nu)
            Kpu = Kpu + (djacb) * (matmul(transpose(Bu),temp_pu))
            temp_pu = 0.0

            write(6,*) '************************Kpu********************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kpu(i,j)
                enddo
            enddo

c       Kup MATRIX

            temp_up = matmul(transpose(qmat),Bu)
            Kup = Kup + (djacb) * (matmul(transpose(Nu),temp_up))
            temp_up = 0.0

            write(6,*) '******Kup************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kup(i,j)
                enddo
            enddo

c       qmat2 INITIALIZATION

        do i = 1,6
          do j = 1,3
            qmat2(i,j) = 0.5*qmat(i,j)
          enddo
        enddo

c       Kup_r MATRIX (NEEDED FOR RESIDUAL FORMATION)

            tempup_r = matmul(qmat2,Nu)
            Kup_r = Kup_r + (djacb) * (matmul(transpose(Bu),tempup_r))
            tempup_r = 0.0

            write(6,*) '******************Kup_r************************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kup_r(i,j)
                enddo
            enddo

c       Bphi MATRIX INITIALIZATION

        do i = 1,3
          do j = 1,8
            Bphi(i,j) = 0
          enddo
        enddo



c       FORMULATING Bphi MATRIX

        do j = 1,8
          Bphi(1,j) = shp(1,j)
          Bphi(2,j) = shp(2,j)
          Bphi(3,j) = shp(3,j)
        enddo


        write(6,*) '*******************Bphi****************************'
            do i = 1,3
                do j = 1,8
                write(6,*) Bphi(i,j)
                enddo
            enddo

c       Kphi MATRIX

            temp_e = matmul(kapmat,Bphi)
            Kphi = Kphi + (djacb) * (matmul(transpose(Bphi),temp_e))
            temp_e = 0.0

            write(6,*) '********************Kphi***********************'
            do i = 1,8
                do j = 1,8
                write(6,*) Kphi(i,j)
                enddo
            enddo

c       Kphip,Kpphi, Kelea MATRIX

        Kphip = Kphip + (djacb)* (matmul(transpose(Bphi),Nu))
        Kpphi = Kpphi + (djacb) * (matmul(transpose(Nu),Bphi))
        Kelea = Kelea + (djacb) * kapmat(1,1) *
     &          (matmul(transpose(Bphi),Nu))

        write(6,*) '***********************Kphip**********************'
            do i = 1,8
                do j = 1,24
                write(6,*) Kphip(i,j)
                enddo
            enddo

            write(6,*) '********************Kpphi**********************'
            do i = 1,24
                do j = 1,8
                write(6,*) Kpphi(i,j)
                enddo
            enddo

            write(6,*) '***********************kelea*******************'
            do i = 1,8
                do j = 1,24
                write(6,*) kelea(i,j)
                enddo
            enddo

c       Kpps MATRIX

        temp_pps = matmul(alfamat,Nu)
        Kpps = Kpps + (djacb) * (matmul(transpose(Nu),temp_pps))
        temp_pps = 0.0

        write(6,*) '*****************Kpps******************************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kpps(i,j)
                enddo
            enddo


c       Kpps_r RESIDUAL MATRIX

        temppps_r = matmul(alfamatr,Nu)
        Kpps_r = Kpps_r + (djacb) * (matmul(transpose(Nu),temppps_r))
        temppps_r = 0.0

        write(6,*) '*****************Kpps_r***************************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kpps_r(i,j)
                enddo
            enddo

c       bata MATRIX

         bata(1,1) = 1.0/(ba*dt)
         bata(2,2) = 1.0/(ba*dt)
         bata(3,3) = 1.0/(ba*dt)

         write(6,*) '***********************bata**********************'
            do i = 1,3
                do j = 1,3
                write(6,*) bata(i,j)
                enddo
            enddo

        do i = 1,3
          do j = 1,3
            batam(i,j) = 0
          enddo
        enddo

        temp_batam = matmul(bata,Tra)
        batam = batam + (matmul(transpose(Tra),temp_batam))
        temp_batam = 0.0

        write(6,*) '*****************batam*****************************'
            do i = 1,3
                do j = 1,3
                write(6,*) batam(i,j)
                enddo
            enddo

        do i = 1,3
          do j = 1,3
            bata(i,j) = batam(i,j)
          enddo
        enddo

c       Kppv MATRIX

        temp_ppv = matmul(bata,Nu)
        Kppv = Kppv + (djacb) * (matmul(transpose(Nu),temp_ppv))
        temp_ppv = 0.0

        write(6,*) '******Kppv************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kppv(i,j)
                enddo
            enddo

c       Kele MATRIX

         Tran(1,1) = 1.0
         Tran(2,2) = 1.0
         Tran(3,3) = 1.0


        temp_ele = matmul(Tran,Nu)
        Kele = Kele + (djacb) * (matmul(transpose(Nu),temp_ele))
        temp_ele = 0.0

                write(6,*) '******************Kele*********************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kele(i,j)
                enddo
            enddo

c     Bp MATRIX INITIALIZATION

        do i = 1,9
          do j = 1,24
            Bp(i,j) = 0
          enddo
        enddo

c       FORMULATING Bp MATRIX

        do j = 1,8
          Bp(1,j*3-2) = shp(1,j)
          Bp(2,j*3-1) = shp(2,j)
          Bp(3,j*3  ) = shp(3,j)
          Bp(4,j*3-2) = shp(2,j)
          Bp(5,j*3-1) = shp(1,j)
          Bp(6,j*3-2) = shp(3,j)
          Bp(7,j*3  ) = shp(1,j)
          Bp(8,j*3-1) = shp(3,j)
          Bp(9,j*3  ) =s hp(2,j)
        enddo

         write(6,*) '************************Bp************************'
            do i = 1,9
                do j = 1,24
                write(6,*) Bp(i,j)
                enddo
            enddo

c       Kpp MATRIX

        temp_p = matmul(grdmat,Bp)
        Kpp = Kpp + (djacb) * (matmul(transpose(Bp),temp_p))
        temp_p = 0.0

        write(6,*) '***********************Kpp*************************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kpp(i,j)
                enddo
            enddo

c       Kpups MATRIX

       temp_pups = matmul(qmatpS,Nu)
       Kpups = Kpups + (djacb) * (matmul(transpose(Nu),temp_pups))
       temp_pups = 0.0

       write(6,*) '**************************Kpups*********************'
            do i = 1,24
                do j = 1,24
                write(6,*) Kpups(i,j)
                enddo
            enddo




            enddo


            do m = 1,24
               do n = 1,6
                  mforc(m) = mforc(m) - (Bu(n,m) * stress(n) * djacb)
               enddo
            enddo

             write(6,*) '*********************mforc*******************'
            do i = 1,24
                write(6,*) mforc(i)
            enddo

            do m =1,8
               do n = 1,3
                eforc(m) = eforc(m) - (Bphi(n,m) * eledispl(n) * djacb)
               enddo
            enddo

             write(6,*) '******************eforc**********************'
            do i = 1,8
                write(6,*) eforc(i)
            enddo

            do m = 1,24
               do n = 1,9
                  pforc(m) = pforc(m) - (Bp(n,m) * pol(n) * djacb)
               enddo
            enddo

            write(6,*) '************************pforc*****************'
            do i = 1,24
                write(6,*) pforc(i)
            enddo

c       ASSIGNING DISPLACEMENTS

        mdispl(1) = u(1)
        mdispl(2) = u(2)
        mdispl(3) = u(3)
        mdispl(4) = u(8)
        mdispl(5) = u(9)
        mdispl(6) = u(10)
        mdispl(7) = u(15)
        mdispl(8) = u(16)
        mdispl(9) = u(17)
        mdispl(10) = u(22)
        mdispl(11) = u(23)
        mdispl(12) = u(24)
        mdispl(13) = u(29)
        mdispl(14) = u(30)
        mdispl(15) = u(31)
        mdispl(16) = u(36)
        mdispl(17) = u(37)
        mdispl(18) = u(38)
        mdispl(19) = u(43)
        mdispl(20) = u(44)
        mdispl(21) = u(45)
        mdispl(22) = u(50)
        mdispl(23) = u(51)
        mdispl(24) = u(52)


        edispl(1) = u(4)
        edispl(2) = u(11)
        edispl(3) = u(18)
        edispl(4) = u(25)
        edispl(5) = u(32)
        edispl(6) = u(39)
        edispl(7) = u(46)
        edispl(8) = u(53)

        pdispl(1) = u(5)
        pdispl(2) = u(6)
        pdispl(3) = u(7)
        pdispl(4) = u(12)
        pdispl(5) = u(13)
        pdispl(6) = u(14)
        pdispl(7) = u(19)
        pdispl(8) = u(20)
        pdispl(9) = u(21)
        pdispl(10) = u(26)
        pdispl(11) = u(27)
        pdispl(12) = u(28)
        pdispl(13) = u(33)
        pdispl(14) = u(34)
        pdispl(15) = u(35)
        pdispl(16) = u(40)
        pdispl(17) = u(41)
        pdispl(18) = u(42)
        pdispl(19) = u(47)
        pdispl(20) = u(48)
        pdispl(21) = u(49)
        pdispl(22) = u(54)
        pdispl(23) = u(55)
        pdispl(24) = u(56)

c       COMPUTING SDV'S

        strain = matmul(Bu,mdispl)
        elefld = -matmul(Bphi,edispl)
        polargrd = matmul(Bp,pdispl)
        stress = matmul(cmat,strain) - matmul(qmat,elefld)
        mut = matmul(kapmat,elefld)
        eledispl = matmul(transpose(qmat),strain) + mut
        pol = matmul(grdmat,polargrd)





c       ASSIGNING SDV'S

            svars(1) = strain(1)
            svars(2) = strain(2)
            svars(3) = strain(3)
            svars(4) = strain(4)
            svars(5) = strain(5)
            svars(6) = strain(6)
            svars(7) = stress(1)
            svars(8) = stress(2)
            svars(9) = stress(3)
            svars(10) = stress(4)
            svars(11) = stress(5)
            svars(12) = stress(6)
            svars(13) = elefld(1)
            svars(14) = elefld(2)
            svars(15) = elefld(3)
            svars(16) = eledispl(1)
            svars(17) = eledispl(2)
            svars(18) = eledispl(3)
            svars(19) = polargrd(1)
            svars(20) = polargrd(2)
            svars(21) = polargrd(3)
            svars(22) = polargrd(4)
            svars(23) = polargrd(5)
            svars(24) = polargrd(6)
            svars(25) = polargrd(7)
            svars(26) = polargrd(8)
            svars(27) = polargrd(9)
            svars(28) = pol(1)
            svars(29) = pol(2)
            svars(30) = pol(3)
            svars(31) = pol(4)
            svars(32) = pol(5)
            svars(33) = pol(6)
            svars(34) = pol(7)
            svars(35) = pol(8)
            svars(36) = pol(9)

            write(6,*) '**********************SVARS*******************'
            do m=1,36
                write(6,*) svars(m)
            enddo


c       Ke2 MATRIX

            do i = 1,56
               do j = 1,56
                  Ke2(i,j) = 0
               enddo
            enddo

           write(6,*) '********************Ke2*************************'
            do m = 1,56
               do n = 1,56
                write(6,*) Ke2(m,n)
               enddo
            enddo

            do i = 1,24
                do j = 1,24
                   Ke2(i,j) = Kuu(i,j)
                enddo
            enddo

            do i = 1,24
               do j = 33,56
                  Ke2(i,j) = -Kup(i,j-32)
               enddo
            enddo

            do i = 25,32
               do j = 25,32
                  Ke2(i,j) = -Kphi(i-24,j-24)
               enddo
            enddo

            do i = 25,32
               do j = 33,56
                  Ke2(i,j) = Kphip(i-24,j-32)
               enddo
            enddo

            do i =25,32
               do j = 33,56
                  Ke2(j,i) = Ke2(i,j)
               enddo
            enddo

            do i = 33,56
               do j = 1,24
                  Ke2(i,j) = -Kpu(i-32,j)
               enddo
            enddo

            do i = 33,56
               do j = 33,56
                  ii=i-32
                  jj=j-32
                Ke2(i,j)=Kpp(ii,jj)+Kpps(ii,jj)+Kppv(ii,jj)-Kpups(ii,jj)
               enddo
            enddo

c       REORDER THE MATRIX Ke2, NEW STIFFNESS MATRIX IS AMATRX

             do i = 1,56
               do j = 1,56
                  amatrx(i,j) = 0.d0
               enddo
            enddo



            do i = 1,8
               do ki = 1,3
                  m = (i-1) * 7 + ki
                  ii=3*(i-1) + ki

            do j = 1,8
               do pj = 1,3
                  n = (j-1) * 7 + pj
                  jj = 3 * (j-1) + pj
                  amatrx(m,n) = -Ke2(ii,jj)
               enddo

            pj = 4
                  n = (j-1) * 7 + pj
                  jj = 24 + j
                  amatrx(m,n) = -Ke2(ii,jj)

            do pj = 5,7
                  n = (j-1) * 7 + pj
                  jj= 32 + 3 * (j-1) + (pj-4)
                  amatrx(m,n) = -Ke2(ii,jj)
                 enddo
               enddo
            enddo

            ki = 4
                  m = (i-1) * 7 + ki
                  ii = 24 + i

            do j = 1,8
               do pj = 1,3
                  n = (j-1) * 7 + pj
                  jj = 3 * (j-1) + pj
                  amatrx(m,n) = -Ke2(ii,jj)
               enddo

            pj = 4
                  n = (j-1) * 7 + pj
                  jj = 24 + j
                  amatrx(m,n) = -Ke2(ii,jj)

            do pj = 5,7
                  n = (j-1) * 7 + pj
                  jj = 32 + 3 * (j-1) + (pj-4)
                  amatrx(m,n) = -Ke2(ii,jj)
                 enddo
             enddo

            do ki = 5,7
                  m = (i-1) * 7 + ki
                  ii = 32 + 3 * (i-1) + (ki-4)

            do j = 1,8
               do pj = 1,3
                  n = (j-1) * 7 + pj
                  jj = 3 * (j-1) + pj
                  amatrx(m,n) = -Ke2(ii,jj)
               enddo

            pj = 4
                  n = (j-1) * 7 + pj
                  jj = 24 + j
                  amatrx(m,n) = -Ke2(ii,jj)

            do pj = 5,7
                  n = (j-1) * 7 + pj
                  jj = 32 + 3 * (j-1) + (pj-4)
                  amatrx(m,n) = -Ke2(ii,jj)
                enddo
              enddo
            enddo
          enddo



           write(6,*) '********************amatrx*********************'
            do m = 1,56
               do n = 1,56
                write(6,*) amatrx(m,n)
               enddo
            enddo



c       TO FORM RESIDUAL MATRIX

            do i = 1,56
               do j = 1,56
                  Ke_r(i,j) = Ke2(i,j)
               enddo
            enddo

        write(6,*) '********************Ke_r*************************'
            do m = 1,56
               do n = 1,56
                write(6,*) Ke_r(m,n)
               enddo
            enddo

            do i = 1,24
               do j = 33,56
                  Ke_r(i,j )= -Kup_r(i,j-32)
               enddo
            enddo

            do i = 33,56
               do j = 33,56
                   Ke_r(i,j)= Kpp(i-32,j-32) + Kpps_r(i-32,j-32)
               enddo
            enddo

c       REORDER THE MATRIX Ke_r, THE NEW STIFFNESS MATRIX IS Ke_r2

            do i = 1,56
               do j = 1,56
                  Ke_r2(i,j) = 0
               enddo
            enddo

            do i = 1,8
               do ki = 1,3
                  m =( i-1) * 7 + ki
                  ii = 3 * (i-1) + ki

            do j = 1,8
               do pj = 1,3
                  n = (j-1) * 7 + pj
                  jj = 3 * (j-1) + pj
                  Ke_r2(m,n) = Ke_r(ii,jj)
               enddo

            pj = 4
                  n = (j-1) * 7 + pj
                  jj = 24 + j
                  Ke_r2(m,n) = Ke_r(ii,jj)

             do pj = 5,7
                  n = (j-1) * 7 + pj
                  jj = 32 + 3 * (j-1) + (pj-4)
                  Ke_r2(m,n) = Ke_r(ii,jj)
                enddo
              enddo
            enddo

            ki = 4
                  m = (i-1) * 7 + ki
                  ii = 24 + i

            do j = 1,8
               do pj = 1,3
                  n = (j-1) * 7 + pj
                  jj = 3 * (j-1) + pj
                  Ke_r2(m,n) = Ke_r(ii,jj)
               enddo

            pj = 4
                  n = (j-1)*7+pj
                  jj = 24 + j
                  Ke_r2(m,n) = Ke_r(ii,jj)

            do pj = 5,7
                  n = (j-1) * 7 + pj
                  jj = 32 + 3 * (j-1) + (pj-4)
                  Ke_r2(m,n) = Ke_r(ii,jj)
               enddo
            enddo

            do ki = 5,7
                  m = (i-1) * 7 + ki
                  ii = 32 + 3 * (i-1) + (ki-4)

            do j = 1,8
               do pj = 1,3
                  n = (j-1) * 7 + pj
                  jj = 3 * (j-1) + pj
                  Ke_r2(m,n) = Ke_r(ii,jj)
               enddo

            pj = 4
                  n = (j-1 )* 7 + pj
                  jj = 24 + j
                  Ke_r2(m,n) = Ke_r(ii,jj)

            do pj = 5,7
                  n = (j-1) * 7 + pj
                  jj = 32 + 3 * (j-1) + (pj-4)
                  Ke_r2(m,n) = Ke_r(ii,jj)
                enddo
              enddo
            enddo
          enddo

            write(6,*) '**********************Ke_r2*******************'
            do m = 1,56
               do n = 1,56
                write(6,*) Ke_r2(m,n)
               enddo
            enddo




cc........ calculate Rhs

c      do i=1,56
c            rhs(i,1)=0.0d0
c        enddo

c      do i=1,56
c        do j=1,56
c          rhs(i,1)=rhs(i,1)+Ke_r2(i,j)*u(j)
c        enddo
c      enddo

c      write(6,*)  '*****************************rhs*******************'
c        do i = 1,56
c          write(6,*) rhs(i,1)
c        enddo

cc......, including the rate effect

c      do j=1,8
c        do i=1,3
c          interval=0
c          do ii=1,8
c            do k=1,3
c              interval=interval+Kppv((j-1)*3+i,(ii-1)*3+k)*pve(k,ii)
c            enddo
c          enddo

c          rhs((j-1)*7+4+i,1)=rhs((j-1)*7+4+i,1)+interval*dt

c        enddo
c      enddo

c      write(6,*)  '************************rhs************************'
c        do i = 1,56
c          write(6,*) rhs(i,1)
c        enddo


cc...... including the applied vortex field

c      do j=1,8
c        do i=1,3
c          interval=0
c          do ii=1,8
c            do k=1,3
c            interval=interval+Kele((j-1)*3+i,(ii-1)*3+k)*Evor(k,ii)
c            enddo
c          enddo

c          rhs((j-1)*7+4+i,1)=rhs((j-1)*7+4+i,1)-interval

c        enddo
c      enddo

c      write(6,*)  '*******************************rhs*****************'
c        do i = 1,56
c          write(6,*) rhs(i,1)
c        enddo



cc   coupling between Ea and Ein

c      do j=1,8

c        interval=0
c          do ii=1,8
c            do k=1,3
c            interval=interval+Kelea(j,(ii-1)*3+k)*Evor(k,ii)
c            enddo
c          enddo

c          rhs((j-1)*7+4,1)=rhs((j-1)*7+4,1)-interval

c      enddo

c      write(6,*)  '*****************************rhs*******************'
c        do i = 1,56
c          write(6,*) rhs(i,1)
c        enddo



        write(6,*)'***********NODAL DISPLACEMENT U*********'
        write(6,*) 'MECHDISPL 1      =', u(1)
        write(6,*) 'MECHDISPL 2      =', u(2)
        write(6,*) 'MECHDISPL 3      =', u(3)
        write(6,*) 'ELECPOTENTIAL 4  =', u(4)
        write(6,*) 'POLARIZATION 5   =', u(5)
        write(6,*) 'POLARIZATION 6   =', u(6)
        write(6,*) 'POLARIZATION 7   =', u(7)
        write(6,*) 'MECHDISPL 8      =', u(8)
        write(6,*) 'MECHDISPL 9      =', u(9)
        write(6,*) 'MECHDISPL 10     =', u(10)
        write(6,*) 'ELECPOTENTIAL 11 =', u(11)
        write(6,*) 'POLARIZATION 12  =', u(12)
        write(6,*) 'POLARIZATION 13  =', u(13)
        write(6,*) 'POLARIZATION 14  =', u(14)
        write(6,*) 'MECHDISPL 15     =', u(15)
        write(6,*) 'MECHDISPL 16     =', u(16)
        write(6,*) 'MECHDISPL 17     =', u(17)
        write(6,*) 'ELECPOTENTIAL 18 =', u(18)
        write(6,*) 'POLARIZATION 19  =', u(19)
        write(6,*) 'POLARIZATION 20  =', u(20)
        write(6,*) 'POLARIZATION 21  =', u(21)
        write(6,*) 'MECHDISPL 22     =', u(22)
        write(6,*) 'MECHDISPL 23     =', u(23)
        write(6,*) 'MECHDISPL 24     =', u(24)
        write(6,*) 'ELECPOTENTIAL 25 =', u(25)
        write(6,*) 'POLARIZATION 26  =', u(26)
        write(6,*) 'POLARIZATION 27  =', u(27)
        write(6,*) 'POLARIZATION 28  =', u(28)
        write(6,*) 'MECHDISPL 29     =', u(29)
        write(6,*) 'MECHDISPL 30     =', u(30)
        write(6,*) 'MECHDISPL 31     =', u(31)
        write(6,*) 'ELECPOTENTIAL 32 =', u(32)
        write(6,*) 'POLARIZATION 33  =', u(33)
        write(6,*) 'POLARIZATION 34  =', u(34)
        write(6,*) 'POLARIZATION 35  =', u(35)
        write(6,*) 'MECHDISPL 36     =', u(36)
        write(6,*) 'MECHDISPL 37     =', u(37)
        write(6,*) 'MECHDISPL 38     =', u(38)
        write(6,*) 'ELECPOTENTIAL 39 =', u(39)
        write(6,*) 'POLARIZATION 40  =', u(40)
        write(6,*) 'POLARIZATION 41  =', u(41)
        write(6,*) 'POLARIZATION 42  =', u(42)
        write(6,*) 'MECHDISPL 43     =', u(43)
        write(6,*) 'MECHDISPL 44     =', u(44)
        write(6,*) 'MECHDISPL 45     =', u(45)
        write(6,*) 'ELECPOTENTIAL 46 =', u(46)
        write(6,*) 'POLARIZATION 47  =', u(47)
        write(6,*) 'POLARIZATION 48  =', u(48)
        write(6,*) 'POLARIZATION 49  =', u(49)
        write(6,*) 'MECHDISPL 50     =', u(50)
        write(6,*) 'MECHDISPL 51     =', u(51)
        write(6,*) 'MECHDISPL 52     =', u(52)
        write(6,*) 'ELECPOTENTIAL 53 =', u(53)
        write(6,*) 'POLARIZATION 54  =', u(54)
        write(6,*) 'POLARIZATION 55  =', u(55)
        write(6,*) 'POLARIZATION 56  =', u(56)

        write(7,*)'***********NODAL DISPLACEMENT U*********'
        write(7,*) 'MECHDISPL 1      =', u(1)
        write(7,*) 'MECHDISPL 2      =', u(2)
        write(7,*) 'MECHDISPL 3      =', u(3)
        write(7,*) 'ELECPOTENTIAL 4  =', u(4)
        write(7,*) 'POLARIZATION 5   =', u(5)
        write(7,*) 'POLARIZATION 6   =', u(6)
        write(7,*) 'POLARIZATION 7   =', u(7)
        write(7,*) 'MECHDISPL 8      =', u(8)
        write(7,*) 'MECHDISPL 9      =', u(9)
        write(7,*) 'MECHDISPL 10     =', u(10)
        write(7,*) 'ELECPOTENTIAL 11 =', u(11)
        write(7,*) 'POLARIZATION 12  =', u(12)
        write(7,*) 'POLARIZATION 13  =', u(13)
        write(7,*) 'POLARIZATION 14  =', u(14)
        write(7,*) 'MECHDISPL 15     =', u(15)
        write(7,*) 'MECHDISPL 16     =', u(16)
        write(7,*) 'MECHDISPL 17     =', u(17)
        write(7,*) 'ELECPOTENTIAL 18 =', u(18)
        write(7,*) 'POLARIZATION 19  =', u(19)
        write(7,*) 'POLARIZATION 20  =', u(20)
        write(7,*) 'POLARIZATION 21  =', u(21)
        write(7,*) 'MECHDISPL 22     =', u(22)
        write(7,*) 'MECHDISPL 23     =', u(23)
        write(7,*) 'MECHDISPL 24     =', u(24)
        write(7,*) 'ELECPOTENTIAL 25 =', u(25)
        write(7,*) 'POLARIZATION 26  =', u(26)
        write(7,*) 'POLARIZATION 27  =', u(27)
        write(7,*) 'POLARIZATION 28  =', u(28)
        write(7,*) 'MECHDISPL 29     =', u(29)
        write(7,*) 'MECHDISPL 30     =', u(30)
        write(7,*) 'MECHDISPL 31     =', u(31)
        write(7,*) 'ELECPOTENTIAL 32 =', u(32)
        write(7,*) 'POLARIZATION 33  =', u(33)
        write(7,*) 'POLARIZATION 34  =', u(34)
        write(7,*) 'POLARIZATION 35  =', u(35)
        write(7,*) 'MECHDISPL 36     =', u(36)
        write(7,*) 'MECHDISPL 37     =', u(37)
        write(7,*) 'MECHDISPL 38     =', u(38)
        write(7,*) 'ELECPOTENTIAL 39 =', u(39)
        write(7,*) 'POLARIZATION 40  =', u(40)
        write(7,*) 'POLARIZATION 41  =', u(41)
        write(7,*) 'POLARIZATION 42  =', u(42)
        write(7,*) 'MECHDISPL 43     =', u(43)
        write(7,*) 'MECHDISPL 44     =', u(44)
        write(7,*) 'MECHDISPL 45     =', u(45)
        write(7,*) 'ELECPOTENTIAL 46 =', u(46)
        write(7,*) 'POLARIZATION 47  =', u(47)
        write(7,*) 'POLARIZATION 48  =', u(48)
        write(7,*) 'POLARIZATION 49  =', u(49)
        write(7,*) 'MECHDISPL 50     =', u(50)
        write(7,*) 'MECHDISPL 51     =', u(51)
        write(7,*) 'MECHDISPL 52     =', u(52)
        write(7,*) 'ELECPOTENTIAL 53 =', u(53)
        write(7,*) 'POLARIZATION 54  =', u(54)
        write(7,*) 'POLARIZATION 55  =', u(55)
        write(7,*) 'POLARIZATION 56  =', u(56)



        return
        end subroutine uel

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc        NON-LINEAR MATERIAL PROPERTIES (NON-CONSTANT MATRIX)          c
cc              IN LOCAL OR MATERIAL CO-ORDINATES                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        SUBROUTINE non_material(pe,qmat,alfamat,B,alfamatr,shp,
     &                  xyz0,uelem,qmatpS)

      INCLUDE 'ABA_PARAM.INC'

      double precision, dimension(6,3)::qmat
      double precision, dimension(3,3)::alfamat
      double precision, dimension(3,8)::pe
      double precision, dimension(3)::p
      double precision, dimension(4,8)::shp
      double precision, dimension(3)::pa
      double precision, dimension(3)::xyz0
      double precision, dimension(3,3)::alfamatr
      double precision, dimension(3,3)::alfamat1
      double precision, dimension(3,3)::alfamatr1
      double precision, dimension(6,3)::qmat1
      double precision, dimension(3,8)::uelem
      double precision, dimension(3,3)::uxdx
      double precision, dimension(3,3)::uada
      double precision, dimension(3,3)::eps
      double precision, dimension(3,3)::qmatpS
      double precision, dimension(3,3)::qmatpS1
      double precision, dimension(3,3)::Tra
      double precision, dimension(6,6)::Trb
      double precision, dimension(9,9)::Trc
      double precision, dimension(17)::B
      double precision, dimension(6,3)::temp_qmat
      double precision, dimension(3,3)::temp_alfa
      double precision, dimension(3,3)::temp_alfar
      double precision              ::interval

      integer :: i,j,k,m,n

      call trans(Tra,Trb,Trc,xyz0)

c       POLARIZATION COMPONENTS IN MATERIAL CO-ORDINATE SYSTEM

        p(1) = 0.0
        p(2) = 1.0
        p(3) = 0.0

            do j = 1,3
               do i = 1,8
                  p(j) = p(j) + shp(4,i) * pe(j,i)
               enddo
            enddo


            write(6,*) '********************p************************'
            do i=1,3
                write(6,*) p(i)
            enddo

c       DEFINING PARAMETER pa IN MATERIAL CO-ORDINATE SYSTEM

            do i = 1,3
               pa(i) = 0
            enddo

            do i = 1,3
               do j = 1,3
                  pa(i) = pa(i) + Tra(i,j) * p(j)
               enddo
            enddo

            write(6,*) '**********************pa**********************'
            do i=1,3
                write(6,*) pa(i)
            enddo

            do j = 1,3
               p(j) = pa(j)
            enddo

c       STRAIN COMPONENTS IN MATERIAL CO-ORDINATE SYSTEM

            do i = 1,3
               do j = 1,3
                  uxdx(i,j)=0.0
                  uada(i,j)=0.0
               enddo
            enddo

            do i = 1,3
               do j=1,3
                  do k=1,8
                     uxdx(i,j)=uxdx(i,j) +shp(j,k)*uelem(i,k)
                  enddo
               enddo
            enddo

            write(6,*) '**********************uxdx********************'
            do i = 1,3
               do j = 1,8
                  write(6,*) uxdx(i,j)
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  do m=1,3
                     do n=1,3
                        uada(i,j)=uada(i,j) +Tra(i,m)*Tra(j,n)*uxdx(m,n)
                     enddo
                   enddo
               enddo
            enddo

             write(6,*) '*******************uada**********************'
            do i = 1,3
               do j = 1,8
                  write(6,*) uada(i,j)
               enddo
            enddo

            eps(1,1)=uada(1,1)
            eps(2,2)=uada(2,2)
            eps(3,3)=uada(3,3)
            eps(1,2)=(uada(1,2)+uada(2,1))/2.0
            eps(1,3)=(uada(1,3)+uada(3,1))/2.0
            eps(2,3)=(uada(2,3)+uada(3,2))/2.0
            eps(2,1)=eps(1,2)
            eps(3,1)=eps(1,3)
            eps(3,2)=eps(2,3)

             write(6,*) '********************eps**********************'
            do i = 1,3
               do j = 1,3
                  write(6,*) eps(i,j)
               enddo
            enddo

c       qmat STIFFNESS MATRIX

            do i = 1,6
               do j = 1,3
                  qmat(i,j) = 0
               enddo
            enddo

             do i = 1,3
               do j = 1,3
                  if (i.ne.j) then
                     qmat(i,j) = 2*B(10)*p(j)
                  else
                     qmat(i,j) = 2*B(9)*p(j)
                  endif
                enddo
            enddo

            qmat(6,1) = B(11)*p(2)
            qmat(5,1) = B(11)*p(3)
            qmat(6,2) = B(11)*p(1)
            qmat(4,2) = B(11)*p(3)
            qmat(5,3) = B(11)*p(1)
            qmat(4,3) = B(11)*p(2)

            write(6,*) '***********************qmat*******************'
            do i = 1,6
               do j = 1,3
                  write(6,*) qmat(i,j)
               enddo
            enddo

            do i = 1,6
               do j = 1,3
                  qmat1(i,j) = 0
               enddo
            enddo

            do m=1,6
          do n=1,3
            do i=1,6
              interval=0
              do j=1,3
                interval=interval+qmat(i,j)*Tra(j,n)
              enddo
              qmat1(m,n)=qmat1(m,n)+Trb(i,m)*interval
            enddo
          enddo
        enddo

            write(6,*) '*******************qmat1**********************'
            do i = 1,6
               do j = 1,3
                  write(6,*) qmat1(i,j)
               enddo
            enddo

            do i = 1,6
               do j = 1,3
                  qmat(i,j) = qmat1(i,j)
               enddo
            enddo



c       qmatpS STIFFNESS MATRIX

            do i = 1,3
               do j=1,3
                  if (i.ne.j) then
                     qmatpS(i,j) = 2*B(11)*eps(i,j)
                  endif
               enddo
            enddo

        qmatpS(1,1) = 2*B(9)*eps(1,1)+2*B(10)*eps(2,2)+2*B(10)*eps(3,3)
        qmatpS(2,2) = 2*B(10)*eps(1,1)+2*B(9)*eps(2,2)+2*B(10)*eps(3,3)
        qmatpS(3,3) = 2*B(10)*eps(1,1)+2*B(10)*eps(2,2)+2*B(9)*eps(3,3)

         write(6,*) '*******************qmatpS************************'
            do i = 1,3
               do j = 1,3
                  write(6,*) qmatpS(i,j)
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  qmatpS1(i,j)=0
               enddo
            enddo

            do m = 1,3
            do n = 1,3
            do i = 1,3
            do j = 1,3
               qmatpS1(m,n) = qmatpS1(m,n)+qmatpS(i,j)*Tra(i,m)*Tra(j,n)
            enddo
            enddo
            enddo
            enddo

            write(6,*) '******qmatps1******'
            do i = 1,3
               do j = 1,3
                  write(6,*) qmatpS1(i,j)
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  qmatpS(i,j) = qmatpS1(i,j)
               enddo
            enddo

c            write(6,*) '******qmatpS******'
c            do i = 1,3
c               do j = 1,3
c                  write(6,*) qmatpS(i,j)
c               enddo
c            enddo

c          alfamat STIFFNESS MATRIX

c            do i = 1,3
c               do j = 1,3
c                  alfamat(i,j) = 0
c               enddo
c            enddo


        alfamat(1,1)=2*B(12)+12*B(13)*p(1)**2+2*B(14)*(p(2)**2+p(3)**2)
     &    +30*B(15)*p(1)**4+2*B(17)*p(2)**2*p(3)**2
     &    +B(16)*(12*p(1)**2*(p(2)**2+p(3)**2)+2*(p(2)**4+p(3)**4))


        alfamat(2,2)=2*B(12)+12*B(13)*p(2)**2+2*B(14)*(p(1)**2+p(3)**2)
     &    +30*B(15)*p(2)**4+2*B(17)*p(1)**2*p(3)**2
     &    +B(16)*(12*p(2)**2*(p(1)**2+p(3)**2)+2*(p(1)**4+p(3)**4))


        alfamat(3,3)=2*B(12)+12*B(13)*p(3)**2+2*B(14)*(p(2)**2+p(1)**2)
     &    +30*B(15)*p(3)**4+2*B(17)*p(2)**2*p(1)**2
     &    +B(16)*(12*p(3)**2*(p(2)**2+p(1)**2)+2*(p(2)**4+p(1)**4))


        alfamat(1,2)=4*B(14)*p(1)*p(2)+8*B(16)*(p(1)**3*p(2)+
     &     p(1)*p(2)**3)+4*B(17)*p(2)*p(1)*p(3)**2

        alfamat(1,3)=4*B(14)*p(1)*p(3)+8*B(16)*(p(1)**3*p(3)+
     &    p(1)*p(3)**3)+4*B(17)*p(3)*p(1)*p(2)**2

        alfamat(2,3)=4*B(14)*p(2)*p(3)+8*B(16)*(p(2)**3*p(3)+
     &    p(2)*p(3)**3)+4*B(17)*p(3)*p(2)*p(1)**2

            alfamat(2,1) = alfamat(1,2)
            alfamat(3,1) = alfamat(1,3)
            alfamat(3,2) = alfamat(2,3)

            write(6,*) '**********************alfamat******************'
            do i = 1,3
               do j = 1,3
                  write(6,*) alfamat(i,j)
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  alfamat1(i,j) = 0
               enddo
            enddo

            do m=1,3
        do n=1,3
          do i=1,3
            interval=0
            do j=1,3
              interval=interval+alfamat(i,j)*Tra(j,n)
            enddo
            alfamat1(m,n)=alfamat1(m,n)+Tra(i,m)*interval
          enddo
        enddo
      enddo

            write(6,*) '****************alfamat1**********************'
            do i = 1,3
               do j = 1,3
                  write(6,*) alfamat1(i,j)
               enddo
            enddo

            do i = 1,3
               do j=1,3
                  alfamat(i,j) = alfamat1(i,j)
               enddo
            enddo

c             write(6,*) '******alfamat******'
c            do i = 1,3
c               do j = 1,3
c                  write(6,*) alfamat(i,j)
c               enddo
c            enddo

c      alfamatr STIFFNESS MATRIX

c            do i = 1,3
c               do j = 1,3
c                  alfamatr(i,j)=0
c               enddo
c            enddo

        alfamatr(1,1)=2*B(12)+4*B(13)*p(1)**2+2*B(14)*(p(2)**2+p(3)**2)
     &    +6*B(15)*p(1)**4+2*B(17)*p(2)**2*p(3)**2
     &    +B(16)*(4*p(1)**2*(p(2)**2+p(3)**2)+2*(p(2)**4+p(3)**4))

        alfamatr(2,2)=2*B(12)+4*B(13)*p(2)**2+2*B(14)*(p(1)**2+p(3)**2)
     &    +6*B(15)*p(2)**4+2*B(17)*p(1)**2*p(3)**2
     &    +B(16)*(4*p(2)**2*(p(1)**2+p(3)**2)+2*(p(1)**4+p(3)**4))

        alfamatr(3,3)=2*B(12)+4*B(13)*p(3)**2+2*B(14)*(p(2)**2+p(1)**2)
     &    +6*B(15)*p(3)**4+2*B(17)*p(2)**2*p(1)**2
     &    +B(16)*(4*p(3)**2*(p(2)**2+p(1)**2)+2*(p(2)**4+p(1)**4))

            alfamatr(1,2) = 0
            alfamatr(1,3) = 0
            alfamatr(2,3) = 0
            alfamatr(2,1) = 0
            alfamatr(3,1) = 0
            alfamatr(3,2) = 0

            write(6,*) '******************alfamatr*********************'
            do i = 1,3
               do j = 1,3
                  write(6,*) alfamatr(i,j)
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  alfamatr1(i,j)=0
               enddo
            enddo

            do m=1,3
        do n=1,3
          do i=1,3
            interval=0
            do j=1,3
              interval=interval+alfamatr(i,j)*Tra(j,n)
            enddo
            alfamatr1(m,n)=alfamatr1(m,n)+Tra(i,m)*interval
          enddo
        enddo
      enddo


            write(6,*) '*****************alfamatr1*********************'
            do i = 1,3
               do j = 1,3
                  write(6,*) alfamatr1(i,j)
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  alfamatr(i,j) = alfamatr1(i,j)
               enddo
            enddo

c             write(6,*) '******alfamatr******'
c            do i = 1,3
c               do j = 1,3
c                  write(6,*) alfamatr(i,j)
c               enddo
c            enddo

            return
            end subroutine non_material

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc          MATERIAL PROPERTIES (CONSTANT MATRIX)                       c
cc               IN GLOBAL CO-ORDINATES                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        SUBROUTINE material(cmat,kapmat,B,grdmat,xyz0)

        INCLUDE 'ABA_PARAM.INC'

        double precision, dimension(6,6)::cmat
        double precision, dimension(3,3)::kapmat
        double precision, dimension(9,9)::grdmat
        double precision, dimension(3,3)::Tra
        double precision, dimension(6,6)::Trb
        double precision, dimension(9,9)::Trc
        double precision, dimension(3)::xyz0
        double precision, dimension(6,6)::cmat1
        double precision, dimension(3,3)::kapmat1
        double precision, dimension(9,9)::grdmat1
        double precision, dimension(17)::B
        double precision               ::interval
c        double precision, dimension(6,6)::temp
c        double precision, dimension(3,3)::temp_kap
c        double precision, dimension(9,9)::temp_grd

        integer :: i,j,m,n

        call trans(Tra,Trb,Trc,xyz0)

c       cmat STIFFNESS MATRIX

            do i = 1,6
               do j = 1,6
                  cmat(i,j) = 0
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  if (i.ne.j) then
                     cmat(i,j) = B(2)
                  else
                     cmat(i,j) = B(1)
                  endif
               enddo
                     cmat(3+i,3+i) = B(3)
            enddo

             write(6,*) '******************* cmat**********************'
            do i=1,6
                do j=1,6
                write(6,*) cmat(i,j)
            enddo
            enddo

            do i = 1,6
               do j = 1,6
                  cmat1(i,j) = 0
               enddo
            enddo

             do m=1,6
        do n=1,6
          do i=1,6
            interval=0
            do j=1,6
              interval=interval+cmat(i,j)*Trb(j,n)
            enddo
            cmat1(m,n)=cmat1(m,n)+Trb(i,m)*interval
          enddo
        enddo
      enddo

            write(6,*) '*******************cmat1***********************'
            do i = 1,6
               do j = 1,6
                  write(6,*) cmat1(i,j)
               enddo
            enddo

            do i = 1,6
               do j = 1,6
                  cmat(i,j) = cmat1(i,j)
               enddo
            enddo

c             write(6,*) '*************************cmat****************'
c            do i = 1,6
c               do j = 1,6
c                  write(6,*) cmat(i,j)
c               enddo
c            enddo

c       kapmat STIFFNESS MATRIX

            do i = 1,3
               do j = 1,3
                  kapmat(i,j) = 0
               enddo
            enddo

            do i = 1,3
               kapmat(i,i) = B(4)
            enddo

            write(6,*) '*******************kapmat**********************'
            do i = 1,3
               do j = 1,3
                  write(6,*) kapmat(i,j)
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  kapmat1(i,j) = 0
               enddo
            enddo

            do m=1,3
        do n=1,3
          do i=1,3
            interval=0
            do j=1,3
              interval=interval+kapmat(i,j)*Tra(j,n)
            enddo
            kapmat1(m,n)=kapmat1(m,n)+Tra(i,m)*interval
          enddo
        enddo
      enddo

            write(6,*) '********************kapmat1********************'
            do i = 1,3
               do j = 1,3
                  write(6,*) kapmat1(i,j)
               enddo
            enddo


            do i = 1,3
               do j = 1,3
                  kapmat(i,j) = kapmat1(i,j)
               enddo
            enddo


c       grdmat STIFFNESS MATRIX

            do i = 1,9
               do j = 1,9
                  grdmat(i,j) = 0
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  if (i.ne.j) then
                     grdmat(i,j) = B(6)
                  else
                     grdmat(i,j) = B(5)
                  endif
               enddo

            grdmat(3+2*i-1,3+2*i-1) = B(7) + B(8)
            grdmat(3+2*i  ,3+2*i  ) = B(7) + B(8)
            grdmat(3+2*i-1,3+2*i  ) = B(7) - B(8)
            grdmat(3+2*i  ,3+2*i-1) = B(7) - B(8)
            enddo

            write(6,*) '******grdmat******'
            do i = 1,9
               do j = 1,9
                  write(6,*) grdmat(i,j)
               enddo
            enddo


            do i = 1,9
               do j = 1,9
                  grdmat1(i,j) = 0
               enddo
            enddo

            do m=1,9
        do n=1,9
          do i=1,9
            interval=0
            do j=1,9
              interval=interval+grdmat(i,j)*Trc(j,n)
            enddo
            grdmat1(m,n)=grdmat1(m,n)+Trc(i,m)*interval
          enddo
        enddo
      enddo

            write(6,*) '*******************grdmat1*********************'
            do i = 1,9
               do j = 1,9
                  write(6,*) grdmat1(i,j)
               enddo
            enddo

            do i = 1,9
               do j = 1,9
                  grdmat(i,j) = grdmat1(i,j)
               enddo
            enddo

c             write(6,*) '******grdmat1******'
c            do i = 1,9
c               do j = 1,9
c                  write(6,*) grdmat(i,j)
c               enddo
c            enddo

            return
            end subroutine material

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             JACOBI MATRIX AND CARTESIAN DERIVATIVES                  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        SUBROUTINE jacob(shp,djacb,coords,nnode,mcrd,xjacm,xjaci)

        INCLUDE 'ABA_PARAM.INC'

        dimension::coords(mcrd,nnode)
        dimension::shp(4,8)
        dimension::xjaci(3,3)
        dimension::xjacm(3,3)
        dimension::cartd(3,8)
        dimension::a(3,3)

        integer :: i,j,k

c       xjacm JACOBI MATRIX

             do i = 1,3
               do j = 1,3
                  xjacm(i,j) = 0.0d0
                  xjaci(i,j) = 0.0d0
                do k = 1,8
                  xjacm(i,j)=xjacm(i,j)+shp(i,k)*coords(j,k)
             enddo
           enddo
         enddo

         write(6,*) '**************** xjacm****************************'
        do i=1,3
            do j=1,3
                write(6,*) xjacm(i,j)
            enddo
        enddo


c       DETERMINANT

        djacb1 = xjacm(1,1)*(xjacm(2,2)*xjacm(3,3)-xjacm(2,3)*xjacm(3,2))
        djacb2 =-xjacm(1,2)*(xjacm(2,1)*xjacm(3,3)-xjacm(2,3)*xjacm(3,1))
        djacb3 = xjacm(1,3)*(xjacm(2,1)*xjacm(3,2)-xjacm(2,2)*xjacm(3,1))
        djacb= djacb1 + djacb2 + djacb3


c             write(6,*) '**************** djacb***********************'
c        do i=1,3
c            do j=1,3
c                write(6,*) djacb(i,j)
c            enddo
c        enddo


c       INVERSE OF JACOBIAN

            do i = 1,3
               do j = 1,3
                  a(i,j) = 0.0d0
               enddo
            enddo

            a(1,1) =  xjacm(2,2)*xjacm(3,3)-xjacm(2,3)*xjacm(3,2)
            a(1,2) =-(xjacm(2,1)*xjacm(3,3)-xjacm(2,3)*xjacm(3,1))
            a(1,3) =  xjacm(2,1)*xjacm(3,2)-xjacm(2,2)*xjacm(3,1)
            a(2,1) =-(xjacm(1,2)*xjacm(3,3)-xjacm(1,3)*xjacm(3,2))
            a(2,2) =  xjacm(1,1)*xjacm(3,3)-xjacm(1,3)*xjacm(3,1)
            a(2,3) =-(xjacm(1,1)*xjacm(3,2)-xjacm(1,2)*xjacm(3,1))
            a(3,1) =  xjacm(1,2)*xjacm(2,3)-xjacm(1,3)*xjacm(2,2)
            a(3,2) =-(xjacm(1,1)*xjacm(2,3)-xjacm(1,3)*xjacm(2,1))
            a(3,3) =  xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)

             write(6,*) '*********************** a*********************'
        do i=1,3
            do j=1,3
                write(6,*) a(i,j)
            enddo
        enddo


            xjaci(1,1)= a(1,1)/djacb
            xjaci(1,2)= a(2,1)/djacb
            xjaci(1,3)= a(3,1)/djacb
            xjaci(2,1)= a(1,2)/djacb
            xjaci(2,2)= a(2,2)/djacb
            xjaci(2,3)= a(3,2)/djacb
            xjaci(3,1)= a(1,3)/djacb
            xjaci(3,2)= a(2,3)/djacb
            xjaci(3,3)= a(3,3)/djacb

             write(6,*) '******************** xjaci********************'
        do i=1,3
            do j=1,3
                write(6,*) xjaci(i,j)
            enddo
        enddo

cc         CARTESIAN DERIVATIVES

            do  i = 1,3
               do  j = 1,8
                  cartd(i,j) = 0.0d0
               do  k = 1,3
                  cartd(i,j)=cartd(i,j)+xjaci(i,k)*shp(k,j)
                 enddo
                enddo
              enddo

              write(6,*) '****************** cartd*********************'
        do i=1,3
            do j=1,8
                write(6,*) cartd(i,j)
            enddo
        enddo

            do  i=1,3
               do  j =1,8
                  shp(i,j) = cartd(i,j)
               enddo
            enddo

c            write(6,*) '*** shp***'
c        do i=1,3
c            do j=1,8
c                write(6,*) shp(i,j)
c            enddo
c        enddo

            return
            end subroutine jacob

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc            SHAPE FUNCTION AND THEIR DERIVATIVES                      c
cc              GLOBAL CO-ORDINATE SYSTEN x,y,z                         c
cc              LOCAL CO-ORDINATE SYSTEM s,t,y                          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        SUBROUTINE sfr(shp,s,t,y,nnode)

        INCLUDE 'ABA_PARAM.INC'

        dimension::shp(4,nnode)
        real :: od,ed

            od = 1.0d0
            ed = 8.0d0

c          SHAPE FUNCTION

            shp(4,1) = (od-s)*(od-t)*(od-y)/ed
            shp(4,2) = (od+s)*(od-t)*(od-y)/ed
            shp(4,3) = (od+s)*(od+t)*(od-y)/ed
            shp(4,4) = (od-s)*(od+t)*(od-y)/ed
            shp(4,5) = (od-s)*(od-t)*(od+y)/ed
            shp(4,6) = (od+s)*(od-t)*(od+y)/ed
            shp(4,7) = (od+s)*(od+t)*(od+y)/ed
            shp(4,8) = (od-s)*(od+t)*(od+y)/ed

c           DERIVATIVES

            shp(1,1) = -(od-t)*(od-y)/ed
            shp(1,2) =  (od-t)*(od-y)/ed
            shp(1,3) =  (od+t)*(od-y)/ed
            shp(1,4) = -(od+t)*(od-y)/ed
            shp(1,5) = -(od-t)*(od+y)/ed
            shp(1,6) =  (od-t)*(od+y)/ed
            shp(1,7) =  (od+t)*(od+y)/ed
            shp(1,8) = -(od+t)*(od+y)/ed

            shp(2,1) = (od-s)*(-od)*(od-y)/ed
            shp(2,2) = (od+s)*(-od)*(od-y)/ed
            shp(2,3) = (od+s)*( od)*(od-y)/ed
            shp(2,4) = (od-s)*( od)*(od-y)/ed
            shp(2,5) = (od-s)*(-od)*(od+y)/ed
            shp(2,6) = (od+s)*(-od)*(od+y)/ed
            shp(2,7) = (od+s)*( od)*(od+y)/ed
            shp(2,8) = (od-s)*( od)*(od+y)/ed

            shp(3,1) = (od-s)*(od-t)*(-od)/ed
            shp(3,2) = (od+s)*(od-t)*(-od)/ed
            shp(3,3) = (od+s)*(od+t)*(-od)/ed
            shp(3,4) = (od-s)*(od+t)*(-od)/ed
            shp(3,5) = (od-s)*(od-t)*(+od)/ed
            shp(3,6) = (od+s)*(od-t)*(+od)/ed
            shp(3,7) = (od+s)*(od+t)*(+od)/ed
            shp(3,8) = (od-s)*(od+t)*(+od)/ed


            write(6,*) '********************* shp**********************'
        do i=1,4
            do j=1,8
                write(6,*) shp(i,j)
            enddo
        enddo

            return
            end subroutine sfr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc              CO-ORDINATE TRANSFORM MATRIX                            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        SUBROUTINE trans(Tra,Trb,Trc,xyz0)

      INCLUDE 'ABA_PARAM.INC'

        double precision, dimension(3,3)::Tra
        double precision, dimension(6,6)::Trb
        double precision, dimension(9,9)::Trc
        double precision::arfa
        double precision, dimension(3)::xyz0
        double precision::pi
        double precision::deg
        double precision::beta
        double precision::theta

        integer :: i,j,m,n,k,l

            deg = 0.3
            pi  = 3.141592653589793
            arfa = deg*pi/180.0
            beta = 0
            theta = 0

c        TRANSFORM MATRIX FOR DISPLACEMENT VECTORS

          Tra(1,1)= cos(arfa)*cos(theta)-sin(arfa)*cos(beta)*sin(theta)
          Tra(1,2)= sin(arfa)*cos(theta)+cos(arfa)*cos(beta)*sin(theta)
          Tra(1,3)= sin(beta)*sin(theta)
          Tra(2,1)=-cos(arfa)*sin(theta)-sin(arfa)*cos(beta)*cos(theta)
          Tra(2,2)=-sin(arfa)*sin(theta)+cos(arfa)*cos(beta)*cos(theta)
          Tra(2,3)= sin(beta)*cos(theta)
          Tra(3,1)= sin(arfa)*sin(beta)
          Tra(3,2)=-cos(arfa)*sin(beta)
          Tra(3,3)= cos(beta)

c        TRANSFORM MATRIX FOR PLOARIZATION GRADIENT VECTOR

            do i = 1,9
               do j = 1,9
               Trc(i,j) = 0
               enddo
            enddo

            do j = 1,9
               if(j.eq.1.or.j.eq.2.or.j.eq.3)then
               n = j
               l = j
            else if(j.eq.4)then
               n = 1
               l = 2
            else if(j.eq.5)then
               n = 2
               l = 1
            else if(j.eq.6)then
               n = 1
               l = 3
            else if(j.eq.7)then
               n = 3
               l = 1
            else if(j.eq.8)then
               n = 2
               l = 3
            else if(j.eq.9)then
               n = 3
               l = 2
            end if

            do i = 1,9
               if(I.eq.1.or.I.eq.2.or.I.eq.3)then
               m = i
               k = i
            else if(I.eq.4)then
               m = 1
               k = 2
            else if(I.eq.5)then
               m = 2
               k = 1
            else if(I.eq.6)then
               m = 1
               k = 3
            else if(I.eq.7)then
               m = 3
               k = 1
            else if(I.eq.8)then
               m = 2
               k = 3
            else if(I.eq.9)then
               m = 3
               k = 2
            endif

               Trc(i,j)=Tra(m,n)*Tra(k,l)
               enddo
            enddo

c            TRANSFORM MATRIX FOR STRESS/STRAIN VECTOR

             do i = 1,6
               do j = 1,6
                  Trb(i,j) = 0
               enddo
            enddo

            do i = 1,3
               do j = 1,3
                  Trb(i,j)=Tra(i,j)*Tra(i,j)
               enddo
            enddo

            do i = 1,3
               Trb(i,4)=Tra(i,2)*Tra(i,3)
               Trb(i,5)=Tra(i,1)*Tra(i,3)
               Trb(i,6)=Tra(i,1)*Tra(i,2)
            enddo

            do j = 1,3
               Trb(4,j)=2*Tra(2,j)*Tra(3,j)
               Trb(5,j)=2*Tra(1,j)*Tra(3,j)
               Trb(6,j)=2*Tra(1,j)*Tra(2,j)
            enddo

            Trb(4,4)=Tra(2,2)*Tra(3,3)+Tra(2,3)*Tra(3,2)
            Trb(4,5)=Tra(2,1)*Tra(3,3)+Tra(2,3)*Tra(3,1)
            Trb(4,6)=Tra(2,1)*Tra(3,2)+Tra(2,2)*Tra(3,1)
            Trb(5,4)=Tra(1,2)*Tra(3,3)+Tra(1,3)*Tra(3,2)
            Trb(5,5)=Tra(1,1)*Tra(3,3)+Tra(1,3)*Tra(3,1)
            Trb(5,6)=Tra(1,1)*Tra(3,2)+Tra(1,2)*Tra(3,1)
            Trb(6,4)=Tra(1,2)*Tra(2,3)+Tra(1,3)*Tra(2,2)
            Trb(6,5)=Tra(1,1)*Tra(2,3)+Tra(1,3)*Tra(2,1)
            Trb(6,6)=Tra(1,1)*Tra(2,2)+Tra(1,2)*Tra(2,1)

            write(6,*) '******************** Tra***********************'
         do i=1,3
            do j=1,3
                write(6,*) Tra(i,j)
            enddo
        enddo

        write(6,*) '*********************** Trb************************'
        do i=1,6
            do j=1,6
                write(6,*) Trb(i,j)
            enddo
        enddo

         write(6,*) '*********************** Trc***********************'
        do i=1,9
            do j=1,9
                write(6,*) Trc(i,j)
            enddo
        enddo




            return
            end subroutine trans
