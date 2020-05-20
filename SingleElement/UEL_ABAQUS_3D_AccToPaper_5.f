ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c....... The subroutine is to form the element matrix                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Literature reference:
c   Three-dimensional finite element modeling       c 
c   of polarization switching in a ferroelectric    c
c   single domain with an impermeable notch         c

        SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)


      INCLUDE 'ABA_PARAM.INC'


      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),SVARS(*),
     1 ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),
     2 A(NDOFEL),TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
        double precision, dimension(3,3)::Tra                           !Transform matrix for coordinate or displacment vectors 
        double precision, dimension(6,6)::Trb                           !Transform matrix for stress or strain
        double precision, dimension(9,9)::Trc                           !Transform matrix for polarization gradient
        double precision, dimension(3,8)::wg                            !Gauss Points with x and y co-ordinates
        double precision                ::xy0(3)                        !Input to transposition matrix
        
        
        double precision, dimension(24,24)::Kuu                         !K_uu as in literature
        double precision, dimension(24,24)::Kup                         !K_up as in literature
        double precision, dimension(24,24)::Kpu                         !K_pu as in literature
        double precision, dimension(24,24)::Kpp                         !K_pp as in literature
        double precision, dimension(24,24)::Kppd                        !K_pp' as in literature
        double precision, dimension(8,24) ::Kphip                       !K_phip as in literature
        double precision, dimension(8,8)  ::Kphi                        !K_phiphi as in literature
        double precision, dimension(56,56)::Ke                          
        double precision, dimension(56,56)::Ke2                         !Ke2 for d_j=[u1 v1 ... u4 v4 phi1 ... phi4 Px1 Py1 ... Px4 Py4]'
        
        
        double precision, dimension(3,24)::Nu                           !Same for Np
        double precision, dimension(6,24)::Bu                           !B_u as in literature
        double precision, dimension(3,8) ::Bphi                         !B_phi as in literature
        double precision, dimension(9,24)::Bp                           !B_p as in literature
        
        double precision, dimension(6,6)::cmat                          ! c : Elastic constants
        double precision, dimension(3,3)::kapmat                        ! kappa
        double precision, dimension(9,9)::grdmat                        ! G
        double precision, dimension(6,3)::qmat                          ! Q(P_i) Electrostrictive coefficients as function of polarization
        double precision, dimension(3,3)::alfamat                       !alpha(P_i) as in literature
        
        
        double precision, dimension(3,8)::uelem                         !Corresponds to U(NDOFEL) Displacement DOF
        double precision, dimension(3,8)::pe                            !Corresponds to U(NDOFEL) Polarization DOF           
        
        ! xe=COORD, ue=U, ve=V, re=RHS
        
        double precision,  dimension(56) ::Re
        double precision,  dimension(24) ::mdispl                       !Mechanical displacements 
        double precision,  dimension(8)  ::edispl                       !Electrical potentials
        double precision,  dimension(24) ::pdispl                       !Polarizations
        
        double precision,  dimension(6) ::strain                        !strains
        double precision,  dimension(3) ::elefld                        !Electrical fields
        double precision,  dimension(9) ::polargrd                      !Polarization gradient
        double precision,  dimension(6) ::stress                        !Stresses
        double precision,  dimension(3) ::mut
        double precision,  dimension(3) ::eledispl                      !Electrical displacements
        double precision,  dimension(9) ::pol                           !???  ???!

                        
        double precision, dimension(3,3)::xjaci                         !Input to Jacobian Matrix
        double precision, dimension(3,3)::xjacm                         !Input to Jacobian Matrix
        double precision, dimension(4,8)::shp                           !Shape function shp(3,:) and shapefunction derivative matrix shp(1:2,:)
        double precision                ::djacb                         !Determinant of Jacobian Matrix
        
        double precision s,t,y,interval
        integer i,j,k,nel,ndm,ig,ii,jj,m,n,ki,pj
        
        double precision,  dimension(17) ::B

        nel=NNODE
        ndm=MCRD

c.......Assigning properties to a variable
        DO i=1,17
            B(i) = PROPS(i)
        END DO

c....... create Gaussian points     

        wg(1,1)=-0.57735
        wg(1,2)= 0.57735
        wg(1,3)= 0.57735
        wg(1,4)=-0.57735
        wg(1,5)=-0.57735
        wg(1,6)= 0.57735
        wg(1,7)= 0.57735
        wg(1,8)=-0.57735
        
        do i=1,2
            wg(2,i  )=-0.57735
            wg(2,i+2)= 0.57735
            wg(2,i+4)=-0.57735
            wg(2,i+6)= 0.57735
        enddo
      
        do i=1,4
            wg(3,i  )=-0.57735
            wg(3,i+4)= 0.57735
        enddo

      write(6,*) "Gauss Points"
      do i = 1, size(wg,1)
            write(6,'(20G12.4)') wg(i,:)
      end do
        
c     Gaussian integration weight equals one for two points in one direction
c     so it does not need weight vector
c     set the initial value as zero for matrix

      do i=1,24
      do j=1,24
          Kuu(i,j) =0.0d0
          Kpu(i,j) =0.0d0
          Kup(i,j) =0.0d0
          Kpp(i,j) =0.0d0
          Kppd(i,j)=0.0d0
      enddo
      enddo
        
      do i=1,8
      do j=1,8
          Kphi(i,j)=0.0d0         
      enddo
      enddo

      do i=1,8
      do j=1,24
          Kphip(i,j)=0.0d0
      enddo
      enddo
c........ set previous values of U matrix 

      do j=1,3 
        do k=1,8
          pe(j,k) =U(7*(k-1)+4+j) ![5,12,19,26,33,40,47,54,6,13,20,27,34,
                                  !41,48,55,7,14,21,28,35,42,49,56] Polarization DOF
          uelem(j,k)=U(7*(k-1)+j) ![1,8,15,22,29,36,43,50,2,9,16,23,30
                                  !37,44,51,3,10,17,24,31,38,45,52] Displacement DOF
        enddo
      enddo
      
      write(6,*) "pe Matrix"
      do i = 1, size(pe,1)
            write(6,'(20G12.4)') pe(i,:)
      end do
      
      write(6,*) "uelem Matrix"
      do i = 1, size(uelem,1)
            write(6,'(20G12.4)') uelem(i,:)
      end do

      
c....................(End Extended for 3D).....................c
      
      xy0(1)=0.0  !!Not used 
      xy0(2)=0.0  !!Not used 
      xy0(3)=0.0  !!Not used
      do j=1,3        
        do k=1,8
          xy0(j)=xy0(j) +COORDS(j,k)/4.0
        enddo
      enddo
    
        call trans_matrix01(Tra,Trb,Trc,xy0)
        call material_matrix01(cmat,kapmat,grdmat,xy0,B)
                                                                        
c....................(Start Extended for 3D).....................c
c....... The summation of Gaussian points
      do ig=1,8    !! Running loop to summate over integration points                    
        s=wg(1,ig) !! x direction
        t=wg(2,ig) !! y direction
        y=wg(3,ig) !! z direction

        call sfr23dc01(shp,s,t,y)
        call jaco3dc01(shp,djacb,COORDS,nel,ndm,xjacm,xjaci)
        call non_material_matrix01(pe,qmat,alfamat,shp,
     $       xy0,uelem,B)
     
c....... stiffness matrix K_uu
        !Initiate Bu Matrix
        do i=1,6
          do j=1,24
            Bu(i,j)=0.0d0
          enddo
        enddo
        !Stress matrix = [S_11,S_22,S_33,S_23,S_31,S_12]

        do j=1,8
            Bu(1,j*3-2)=shp(1,j) !Derivative of shape function wrt x
            Bu(2,j*3-1)=shp(2,j) !Derivative of shape function wrt y
            Bu(3,j*3-0)=shp(3,j) !Derivative of shape function wrt z
            Bu(4,j*3-1)=shp(3,j) !Derivative of shape function wrt z
            Bu(4,j*3-0)=shp(2,j) !Derivative of shape function wrt y
            Bu(5,j*3-2)=shp(3,j) !Derivative of shape function wrt z
            Bu(5,j*3-0)=shp(1,j) !Derivative of shape function wrt x
            Bu(6,j*3-2)=shp(2,j) !Derivative of shape function wrt y
            Bu(6,j*3-1)=shp(1,j) !Derivative of shape function wrt x   
        enddo
        
      write(6,*) "djacb"
            write(6,*) djacb

      write(6,*) "Bu Matrix"
      do i = 1, size(Bu,1)
            write(6,'(20G12.4)') Bu(i,:)
      end do
        
       !Kuu= S_V(B_u^T*C*B_u)dV !!!S_V indicates integral over volume
        do m=1,24
        do n=1,24
          do i=1,6 
            interval=0.0d0
            do j=1,6
              interval=interval+cmat(i,j)*Bu(j,n)
            enddo
            Kuu(m,n)=Kuu(m,n)+Bu(i,m)*interval*djacb  ! djacb:Determinent of Jacobian   
          enddo
        enddo
        enddo
        
      write(6,*) "Kuu Matrix"
      do i = 1, size(Kuu,1)
            write(6,'(20G12.4)') Kuu(i,:)
      end do

c....... stiffness matrix K_up and K_pu     
        do i=1,3
          do j=1,24
            Nu(i,j)=0.0d0               !Same for Np : Term mentioned in literature
          enddo
        enddo

        do j=1,8
            Nu(1,j*3-2)=shp(4,j)    !Shape functions are stored in shp(4,:)
            Nu(2,j*3-1)=shp(4,j)
            Nu(3,j*3-0)=shp(4,j)
        enddo
        
      write(6,*) "Nu Matrix"
      do i = 1, size(Nu,1)
            write(6,'(20G12.4)') Nu(i,:)
      end do

      ! Kpu= S_V(N_p^T*Q'*B_u)dV
        do n=1,24
        do m=1,24
          do i=1,6 
                interval=0.0d0
            do j=1,3
                interval=interval+Nu(j,n)*qmat(i,j)  !??? Here how Q' is qmat
            enddo
            Kpu(n,m)=Kpu(n,m)+Bu(i,m)*interval*djacb
          enddo
        enddo
        enddo
        
      write(6,*) "Kpu Matrix"
      do i = 1, size(Kpu,1)
            write(6,'(20G12.4)') Kpu(i,:)
      end do

      ! Here, Kup should be Kup^S in manus and is the symetric matrix of Kpu.
        do m=1,24
        do n=1,24
            Kup(m,n)=Kpu(n,m) !symmetric components  
        enddo
        enddo
        
      write(6,*) "Kup Matrix"
      do i = 1, size(Kup,1)
            write(6,'(20G12.4)') Kup(i,:)
      end do

c....... stiffness matrix K_phi and Kphip     
        do i=1,3
          do j=1,8
            Bphi(i,j)=0.0d0
          enddo
        enddo

        do j=1,8
          Bphi(1,j)=shp(1,j)
          Bphi(2,j)=shp(2,j)
          Bphi(3,j)=shp(3,j)
        enddo
        
      write(6,*) "Bphi Matrix"
      do i = 1, size(Bphi,1)
            write(6,'(20G12.4)') Bphi(i,:)
      end do

        ! Kphi= S_V(B_phi^T*kapmat*B_phi)dV
        do m=1,8
        do n=1,8
           do i=1,3 
             interval=0.0d0
             do j=1,3
             interval=interval+kapmat(i,j)*Bphi(j,n)
             enddo
             Kphi(m,n)=Kphi(m,n)+Bphi(i,m)*interval*djacb     
           enddo
        enddo
        enddo
        
      write(6,*) "Kphi Matrix"
      do i = 1, size(Kphi,1)
            write(6,'(20G12.4)') Kphi(i,:)
      end do

        ! Kphip= S_V(B_phi^T*N_p)dV
        do m=1,8
        do n=1,24                 
           do i=1,3 
             Kphip(m,n)=Kphip(m,n)+Bphi(i,m)*Nu(i,n)*djacb
           enddo
        enddo
        enddo
        
      write(6,*) "Kphip Matrix"
      do i = 1, size(Kphip,1)
            write(6,'(20G12.4)') Kphip(i,:)
      end do

c....... stiffness matrix K_pp           
        do i=1,9
          do j=1,24
            Bp(i,j)=0.0d0
          enddo
        enddo
        
        !Polarization gradient matrix
        ![Bp]*[P']=[P(x,x),P(y,y),P(z,z),P(x,y),P(y,x),P(y,z),P(z,y),P(x,z),P(z,x)]^T
        !Formulated according with Plarization gradient matrix as above 
        do j=1,8
            Bp(1,j*3-2)=shp(1,j) !Derivative of shape function wrt x
            Bp(2,j*3-1)=shp(2,j) !Derivative of shape function wrt y
            Bp(3,j*3-0)=shp(3,j) !Derivative of shape function wrt z
            Bp(4,j*3-2)=shp(2,j) !Derivative of shape function wrt y
            Bp(5,j*3-1)=shp(1,j) !Derivative of shape function wrt x
            Bp(6,j*3-1)=shp(3,j) !Derivative of shape function wrt z
            Bp(7,j*3-0)=shp(2,j) !Derivative of shape function wrt y
            Bp(8,j*3-2)=shp(3,j) !Derivative of shape function wrt z
            Bp(9,j*3-0)=shp(1,j) !Derivative of shape function wrt x

        enddo
        
      write(6,*) "Bp Matrix"
      do i = 1, size(Bp,1)
            write(6,'(20G12.4)') Bp(i,:)
      end do

      ! Kpp=S_V(B_p^T*G*B_p)dV 
        do m=1,24
        do n=1,24
           do i=1,9 
             interval=0.0d0 
             do j=1,9
               interval=interval+grdmat(i,j)*Bp(j,n)
             enddo
             Kpp(m,n)=Kpp(m,n)+Bp(i,m)*interval*djacb
           enddo
        enddo
        enddo
        
      write(6,*) "Kpp Matrix"
      do i = 1, size(Kpp,1)
            write(6,'(20G12.4)') Kpp(i,:)
      end do
        
        ! Kpp'=S_V(N_p*alpha^S*N_p)dV !??? Why this matrix ???!
        do m=1,24
        do n=1,24
           do i=1,3 
             interval=0.0d0  
             do j=1,3
             interval=interval+alfamat(i,j)*Nu(j,n)
             enddo
             Kppd(m,n)=Kppd(m,n)+Nu(i,m)*interval*djacb 
           enddo         
        enddo
        enddo
      enddo ! ****************end of ig*****************  
      
c   Ke2 for d_j=
c   [u11 u12 u13 ... u81 u82 u83 phi1 ... phi8 P11 Py12 p13 ... P81 P82 p83]'
c   with u12 is dispalment of 1st node in y-direction
      do i=1,56
        do j=1,56
           Ke2(i,j)=0.0d0 ! element stiffness matrix S_ij
        enddo
      enddo

      do i=1,24
        do j=1,24
           Ke2(i,j)=Kuu(i,j)
        enddo
      enddo

      do i=1,24
        do j=33,56
           Ke2(i,j)=-Kup(i,j-32)
        enddo
      enddo

      do i=25,32                  ! Matrix 8X8
        do j=25,32
           Ke2(i,j)=-Kphi(i-24,j-24)
        enddo
      enddo

      do i=25,32                  ! Matrix 8X24  !????? Minus sign not included why ?????!
        do j=33,56
           Ke2(i,j)=-Kphip(i-24,j-32)  !********** Gave negative sign **********  
        enddo
      enddo

      do i=33,56                  ! Matrix 24X8  !????? Minus sign not included why ?????!
        do j=25,32
           Ke2(i,j)=Ke2(j,i)   
        enddo
      enddo

      do i=33,56                         ! Matrix 8X8
        do j=1,24
           Ke2(i,j)=-Kpu(i-32,j)   
        enddo
      enddo

      do i=33,56                         ! Matrix 8X8
      do j=33,56
        ii=i-32
        jj=j-32
        Ke2(i,j)=Kpp(ii,jj)+Kppd(ii,jj)         
      enddo
      enddo  
      ! end of ke2
      
      write(6,*) "Ke2 Matrix"
      do i = 1, size(Ke2,1)
            write(6,'(20G12.4)') Ke2(i,:)
      end do
      
c........  reorder the matrix Ke2, new stiff matrix is Ke
c   Ke for d_j=
c   [u11 u12 u13 phi1 P11 P12 p13 ... u81 u82 u83 phi8 P81 P82 p83]'

        do i = 1,56
           do j = 1,56
              Ke(i,j) = 0.0d0
           enddo
        enddo



        do i = 1,8
            do ki = 1,3
                m = (i-1)*7 + ki    ![1 2 3 8 9 10 ...... 50 51 52]
                ii= (i-1)*3 + ki    ![1 2 3 4 5 6  ...... 22 23 24]
                do j = 1,8
                    do pj = 1,3
                      n = (j-1)*7 + pj  ![1 2 3 8 9 10 ...... 50 51 52]
                      jj= (j-1)*3 + pj  ![1 2 3 4 5 6  ...... 22 23 24]
                      Ke(m,n) = Ke2(ii,jj)
                    enddo

                    pj = 4
                    n = (j-1) * 7 + pj      ![4 11 18 ..... 53]
                    jj = 24 + j             ![25 26 27..... 32]
                    Ke(m,n) = Ke2(ii,jj)

                    do pj = 5,7                 
                    n = (j-1) * 7 + pj          ![5 6 7 12 13 14 ...... 54 55 56]
                    jj= 32 + 3 * (j-1) + (pj-4) ![33 34 35 36 37 38 ...... 54 55 56]
                    Ke(m,n) = Ke2(ii,jj)
                    enddo
                enddo
            enddo

            ki = 4
            m = (i-1) * 7 + ki  ![ 4 11 18  ..... 53]
            ii = 24 + i         ![25 26 27  ..... 32]

            do j = 1,8
                do pj = 1,3
                    n = (j-1) * 7 + pj      ![1 2 3 8 9 10 ...... 50 51 52]
                    jj = 3 * (j-1) + pj     ![1 2 3 4 5 6  ...... 22 23 24]
                    Ke(m,n) = Ke2(ii,jj)
                enddo
                pj = 4
                n = (j-1) * 7 + pj      ![4 11 18 ..... 53]
                jj = 24 + j             ![25 26 27  ..... 32]
                Ke(m,n) = Ke2(ii,jj)

                do pj = 5,7
                  n = (j-1) * 7 + pj            ![4 11 18 ..... 53]
                  jj = 32 + 3 * (j-1) + (pj-4)  ![33 34 35 36 37 38 ...... 54 55 56]
                  Ke(m,n) = Ke2(ii,jj)
                enddo
            enddo

            do ki = 5,7
                m = (i-1) * 7 + ki              ![4 11 18 ..... 53]
                ii = 32 + 3 * (i-1) + (ki-4)    ![33 34 35 36 37 38 ...... 54 55 56]
            do j = 1,8
            do pj = 1,3
                n = (j-1) * 7 + pj              ![4 11 18 ..... 53]
                jj = 3 * (j-1) + pj             ![1 2 3 4 5 6  ...... 22 23 24]
                Ke(m,n) = Ke2(ii,jj)
            enddo
            pj = 4
            n = (j-1) * 7 + pj                  ![ 4 11 18  ..... 53]
            jj = 24 + j                         ![25 26 27  ..... 32]
            Ke(m,n) = Ke2(ii,jj)

            do pj = 5,7
            n = (j-1) * 7 + pj              ![4 11 18 ..... 53]
            jj = 32 + 3 * (j-1) + (pj-4)    ![33 34 35 36 37 38 ...... 54 55 56]
            Ke(m,n) = Ke2(ii,jj)
            enddo
          enddo
        enddo
      enddo
c........  end of Ke

      write(6,*) "Ke Matrix"
      do i = 1, size(Ke,1)
            write(6,'(20G12.4)') Ke(i,:)
      end do

c........ calculate Re
      do i=1,56 
        Re(i)=0.0d0
      enddo

      do i=1,56 
        do j=1,56
          Re(i)=Re(i)+Ke(i,j)*U(j) 
        enddo
      enddo
      
      write(6,*) "Re"
      do i = 1, 56
        write(6,*) Re(i)
      end do
      
      
c       ASSIGNING SDV'S !! Have to study and understand about SDV s and on how to calculate them
      
c       ASSIGNING DISPLACEMENTS

        mdispl(1)   = u(1)
        mdispl(2)   = u(2)
        mdispl(3)   = u(3)
        mdispl(4)   = u(8)
        mdispl(5)   = u(9)
        mdispl(6)   = u(10)
        mdispl(7)   = u(15)
        mdispl(8)   = u(16)
        mdispl(9)   = u(17)
        mdispl(10)  = u(22)
        mdispl(11)  = u(23)
        mdispl(12)  = u(24)
        mdispl(13)  = u(29)
        mdispl(14)  = u(30)
        mdispl(15)  = u(31)
        mdispl(16)  = u(36)
        mdispl(17)  = u(37)
        mdispl(18)  = u(38)
        mdispl(19)  = u(43)
        mdispl(20)  = u(44)
        mdispl(21)  = u(45)
        mdispl(22)  = u(50)
        mdispl(23)  = u(51)
        mdispl(24)  = u(52)

        edispl(1)   = u(4)
        edispl(2)   = u(11)
        edispl(3)   = u(18)
        edispl(4)   = u(25)
        edispl(5)   = u(32)
        edispl(6)   = u(39)
        edispl(7)   = u(46)
        edispl(8)   = u(53)

        pdispl(1)   = u(5)
        pdispl(2)   = u(6)
        pdispl(3)   = u(7)
        pdispl(4)   = u(12)
        pdispl(5)   = u(13)
        pdispl(6)   = u(14)
        pdispl(7)   = u(19)
        pdispl(8)   = u(20)
        pdispl(9)   = u(21)
        pdispl(10)  = u(26)
        pdispl(11)  = u(27)
        pdispl(12)  = u(28)
        pdispl(13)  = u(33)
        pdispl(14)  = u(34)
        pdispl(15)  = u(35)
        pdispl(16)  = u(40)
        pdispl(17)  = u(41)
        pdispl(18)  = u(42)
        pdispl(19)  = u(47)
        pdispl(20)  = u(48)
        pdispl(21)  = u(49)
        pdispl(22)  = u(54)
        pdispl(23)  = u(55)
        pdispl(24)  = u(56)


c       COMPUTING SDV'S

        strain = matmul(Bu,mdispl)
        elefld = -matmul(Bphi,edispl)
        polargrd = matmul(Bp,pdispl)
        stress = matmul(cmat,strain) - matmul(qmat,elefld)
        mut = matmul(kapmat,elefld)
        eledispl = matmul(transpose(qmat),strain) + mut
        pol = matmul(grdmat,polargrd) !??? What is this pol ???! 





c       ASSIGNING SDV'S !! Have to study and understand about SDV s and on how to calculate them

            svars(1)  = strain(1)
            svars(2)  = strain(2)
            svars(3)  = strain(3)
            svars(4)  = strain(4)
            svars(5)  = strain(5)
            svars(6)  = strain(6)
            svars(7)  = stress(1)
            svars(8)  = stress(2)
            svars(9)  = stress(3)
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
      do i=1,56 
        RHS(i,1)=-Re(i)
        do j=1,56
        AMATRX(i,j) = Ke(i,j)
        enddo
      enddo
!print the AMATRX into the .dat file.      
      write(6,*) "AMATRX"
      do i = 1, size(AMATRX,1)
          write(6,'(40G12.4)')  AMATRX(i,:)
      end do
      write(6,*) "RHS"
      do i = 1, 56
        write(6,*) RHS(i,1)
      end do
      return
      end subroutine UEL  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c       Jacobi-Matrix, Inverse Jacobian and Determinant of Jacobian  
c       Cartesian derivatives of shape functions     
c       xe: Co-ordinates of nodes in direction 1,2 and 3
c       In xe: Rows stores x,y,z info and coloumns store node info of each element
c       djacb : determinant of jacobian matrix
c       ndm: number of dimentions
c       xjacm : jacobian matrix
c       xjaci : inverse of jacobian matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jaco3dc01(shp,djacb,xe,nel,ndm,xjacm,xjaci) 
        !xe=COORDS  nel=NNODES  ndm=MCRD
        implicit double precision(a-h,o-z)
        
        double precision ndm,nel,djacb,djacb1,djacb2,djacb3
        double precision shp(4,8), xe(3,8), xjaci(3,3), xjacm(3,3)
        double precision cartd(3,8), a(3,3)
        integer idm,jdm,inode
        
        
        
c       Jacobi-matrix : xjacm
c       Jacobian matrix is initiated and values are stored
c       Jacobian inverse matrix is initiated with zeros
        do idm = 1,3                                              
         do jdm = 1,3
            xjacm(idm,jdm) = 0.0d0
            xjaci(idm,jdm) = 0.0d0
            do inode = 1,8
            xjacm(idm,jdm)=xjacm(idm,jdm)+shp(idm,inode)*xe(jdm,inode)
            enddo
         enddo
        enddo
        
         write(6,*) '**************** xjacm****************************'
        do idm=1,3
            do jdm=1,3
                write(6,*) xjacm(idm,jdm)
            enddo
        enddo

c       Determinant of Jacobian
        djacb1= xjacm(1,1)*(xjacm(2,2)*xjacm(3,3)-xjacm(2,3)*xjacm(3,2))
        djacb2=-xjacm(1,2)*(xjacm(2,1)*xjacm(3,3)-xjacm(2,3)*xjacm(3,1))
        djacb3= xjacm(1,3)*(xjacm(2,1)*xjacm(3,2)-xjacm(2,2)*xjacm(3,1))
        djacb= djacb1 + djacb2 + djacb3

c       Inverse of Jacobian : xjaci
        do idm = 1,3
         do jdm = 1,3
            a(idm,jdm) = 0.0d0
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

        xjaci(1,1)= a(1,1)/djacb
        xjaci(1,2)= a(2,1)/djacb
        xjaci(1,3)= a(3,1)/djacb
        xjaci(2,1)= a(1,2)/djacb
        xjaci(2,2)= a(2,2)/djacb
        xjaci(2,3)= a(3,2)/djacb
        xjaci(3,1)= a(1,3)/djacb
        xjaci(3,2)= a(2,3)/djacb
        xjaci(3,3)= a(3,3)/djacb

c       cartesian derivatives : cartd
        do idm = 1,3
         do i = 1,8
           cartd(idm,i) = 0.0d0
             do jdm = 1,3
                cartd(idm,i) = cartd(idm,i)+xjaci(idm,jdm)*shp(jdm,i)
             enddo
         enddo
        enddo
        
        !??? Why replacing shape function derivatives wrt master space 
        ! with shp derivatives wrt global coordinates ???!
        
        do idm=1,3                                        
         do i=1,8
            shp(idm,i) = cartd(idm,i)
         enddo
        enddo

        return
        end       
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c       compute shape functions and their derivatives for 2D elements     
c       global coordinate system x,y; local coordinate system s,t      
c       shp: Shape function matrix with shp(3,:) stores shape functions
c            and shp(1,:), shp(2,:) shape function derivatives for 2D
c       nel: Number of nodes per Elem                                  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        SUBROUTINE sfr23dc01(shp,s,t,y) !Shape_Function_And_SFDerivatives
        INCLUDE 'ABA_PARAM.INC'
        double precision shp(4,8), s, t, y
        REAL :: od,ed
        od=1.0d0
        ed=8.0d0

c       shape functions for 8-node linear element (Hexa element)
        shp(4,1) = (od-s)*(od-t)*(od-y)/ed
        shp(4,2) = (od+s)*(od-t)*(od-y)/ed
        shp(4,3) = (od+s)*(od+t)*(od-y)/ed
        shp(4,4) = (od-s)*(od+t)*(od-y)/ed
        shp(4,5) = (od-s)*(od-t)*(od+y)/ed
        shp(4,6) = (od+s)*(od-t)*(od+y)/ed
        shp(4,7) = (od+s)*(od+t)*(od+y)/ed
        shp(4,8) = (od-s)*(od+t)*(od+y)/ed

c       derivatives
c       shp(a,b) is equivalent to ShapeF b,a with b=1 to 4  and a=1 to 3

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
        
      write(6,*) "Shape Functions"
      do i = 1, size(shp,1)
          write(6,'(40G12.4)')  shp(i,:)
      end do

        RETURN
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     nonlinear material-properties (non-constant matrix)              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine non_material_matrix01(pe,qmat,alfamat,shp
     $   ,xy0,uelem,B)

      implicit none

      double precision qmat(6,3),alfamat(3,3),shp(4,8),p(3)
      double precision pa(3)
      double precision pe(3,8)
      double precision qmat1(6,3),alfamat1(3,3),alfamatr1(2,2)
      double precision a1,a11,a12,a111,a112,a123,q11,q12,q44
      double precision interval
      double precision Tra(3,3),Trb(6,6),Trc(9,9),xy0(3) !?
      double precision uelem(3,8),uxdx(3,3),uada(3,3),eps(3,3) !?
      double precision B(17)
      integer i,j,k,m,n
      
      call trans_matrix01(Tra,Trb,Trc,xy0)
      
c......dimensionless material parameters for PTO3
        a1   =-1.0d0
        a11  =-0.24d0  !1.4047  !stress controled
        a12  = 2.5d0    !-0.8893 !stress controled
        a111 = 0.49d0
        a112 = 1.2d0
        a123 = -7.0d0

        q11 = 63.3d0  !0.661513d+02
        q12 = 2.6d0   !0.26686d+01
        q44 = 41.6d0  !0.434728d+02
      
      p(1) = 0.0d0
      p(2) = 0.0d0
      p(3) = 0.0d0
      do j=1,3
        do i=1,8
          p(j)=p(j)+shp(4,i)*pe(j,i)
        enddo
      enddo
  
      pa(1)=0.0d0
      pa(2)=0.0d0
      pa(3)=0.0d0
      do i=1,3
        do j=1,3
          pa(i)=pa(i)+Tra(i,j)*p(j) ! in material coordinate system
        enddo 
      enddo
      
      do j=1,3
        p(j)=pa(j) ! in material coordinate system
      enddo
      
      
      do i=1,3
      do j=1,3
          uxdx(i,j)=0.0d0
          uada(i,j)=0.0d0
      enddo
      enddo

      ! uxdx stores cartesian derivatives of displacements (Cross check)
      ! This routine is called after jacobian routine in which shp derivatives are changed from master space to global space 
      ! shp(1,k),shp(2,k),shp(3,k) contain derivatives in x,y,z directions
      do i=1,3
      do j=1,3
        do k=1,8
          uxdx(i,j)=uxdx(i,j) +shp(j,k)*uelem(i,k)
        enddo
      enddo
      enddo
      ! uxdx is a tensor (derivative of a vector w.r.t vector)
      ! Q.Tensor.Q_T for coordinate transformation from global to material point
      do I=1,3
      do J=1,3
        do m=1,3
        do n=1,3
          uada(I,J)=uada(I,J) +Tra(I,m)*Tra(J,n)*uxdx(m,n)
        enddo
        enddo
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
        
        !Left here
      
c....... Stiffmatrix qmat  
      qmat(1,1) = 2*q11*p(1)
      qmat(1,2) = 2*q12*p(2)
      qmat(1,3) = 2*q12*p(3)
      qmat(2,1) = 2*q12*p(1)
      qmat(2,2) = 2*q11*p(2)
      qmat(2,3) = 2*q12*p(3)
      qmat(3,1) = 2*q12*p(1)
      qmat(3,2) = 2*q12*p(2)
      qmat(3,3) = 2*q11*p(3)
      qmat(4,1) = 0.0d0
      qmat(4,2) = q44*p(3)
      qmat(4,3) = q44*p(2)
      qmat(5,1) = q44*p(3)
      qmat(5,2) = 0.0d0
      qmat(5,3) = q44*p(2)
      qmat(6,1) = q44*p(2)
      qmat(6,2) = q44*p(1)
      qmat(6,3) = 0.0d0
      
      write(6,*) "qmat before transformation"
      do i = 1, size(qmat,1)
          write(6,'(40G12.4)')  qmat(i,:)
      end do

      do i=1,6
        do j=1,3
          qmat1(i,j)=0.0d0
        enddo
      enddo

      do m=1,6
        do n=1,3
          do i=1,6
            interval=0.0d0
            do j=1,3
              interval=interval+qmat(i,j)*Tra(j,n)
            enddo
            qmat1(m,n)=qmat1(m,n)+Trb(i,m)*interval
          enddo
        enddo
      enddo

      do i=1,6
        do j=1,3
          qmat(i,j)=qmat1(i,j)
        enddo
      enddo
      
      write(6,*) "qmat after transformation"
      do i = 1, size(qmat,1)
          write(6,'(40G12.4)')  qmat(i,:)
      end do
  
c....... Stiffmatrix alfamat
c       Equations for alpha from # Wang J, Zhang TY. Phys Rev B 2006;73:144107.
      alfamat(1,1)=2*(a1 + a11*p(1)**2 + a12*(p(2)**2+p(3)**2)
     $ + a111*p(1)**4 + a112*(p(1)**2*(p(2)**2+p(3)**2)+p(2)**4+p(3)**4)
     $ + a123*p(2)**2*p(3)**2)

      alfamat(2,2)=2*(a1 + a11*p(2)**2 + a12*(p(1)**2+p(3)**2)
     $ + a111*p(2)**4 + a112*(p(2)**2*(p(1)**2+p(3)**2)+p(1)**4+p(3)**4)
     $ + a123*p(1)**2*p(3)**2)
     
      alfamat(3,3)=2*(a1 + a11*p(3)**2 + a12*(p(1)**2+p(2)**2)
     $ + a111*p(3)**4 + a112*(p(3)**2*(p(1)**2+p(2)**2)+p(1)**4+p(2)**4)
     $ + a123*p(1)**2*p(2)**2)

      alfamat(1,2)= a12*p(1)*p(2) 
     $ + a112*(p(1)**3*p(2)+p(1)*p(2)**3)
     $ + a123*p(1)*p(2)*p(3)**2
     
      alfamat(2,1)=alfamat(1,2)
      
      alfamat(1,3)= a12*p(1)*p(3) 
     $ + a112*(p(1)**3*p(3)+p(1)*p(3)**3)
     $ + a123*p(1)*p(3)*p(2)**2
     
      alfamat(3,1)=alfamat(1,3)
      
      alfamat(2,3)= a12*p(2)*p(3) 
     $ + a112*(p(2)**3*p(3)+p(2)*p(3)**3)
     $ + a123*p(2)*p(3)*p(1)**2
     
      alfamat(3,2)=alfamat(2,3)
      
      write(6,*) "alfamat before transformation"
      do i = 1, size(alfamat,1)
          write(6,'(40G12.4)')  alfamat(i,:)
      end do

      do i=1,3
        do j=1,3
          alfamat1(i,j)=0.0d0
        enddo
      enddo

      do m=1,3
        do n=1,3
          do i=1,3
            interval=0.0d0
            do j=1,3
              interval=interval+alfamat(i,j)*Tra(j,n)
            enddo
            alfamat1(m,n)=alfamat1(m,n)+Tra(i,m)*interval
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
           alfamat(i,j)=alfamat1(i,j)
        enddo
      enddo 
      
      write(6,*) "alfamat after transformation"
      do i = 1, size(alfamat,1)
          write(6,'(40G12.4)')  alfamat(i,:)
      end do
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       material-properties (constant matrix)                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine material_matrix01(cmat,kapmat,grdmat,xy0,B)

      implicit none

      double precision cmat(6,6),kapmat(3,3),grdmat(9,9)
      double precision cmat1(6,6),interval,kapmat1(3,3),grdmat1(9,9)
      double precision c11,c12,c44,G11,G12,G44,G44p,Kapa
      double precision Tra(3,3),Trb(6,6),Trc(9,9),xy0(3) !?
      double precision B(17)

      integer i,j,m,n

      call trans_matrix01(Tra,Trb,Trc,xy0)

c     dimensionless material parameters for PTO3 ## Refer paper by wang.
c     Elastic constants : c11=PROPS(1),c12=PROPS(2),c44=PROPS(3)
      c11 = 1766.0d0
      c12 = 802.0d0
      c44 = 1124.0d0

c     Gradient coefficients : G11=PROPS(4),G12=PROPS(5),G44=PROPS(6),G44p=PROPS(7)
      G11 = 1.6d0 !6.4d0   !1.6d0
      G12 = 0.0d0
      G44 = 0.8d0 !3.2d0   !0.8d0
      G44p= 0.8d0 !3.2d0   !0.8d0
c     ????????????????????????????????????????????????????????????????
!     Kapa = 0.0015d0 !total polarization
      Kapa = 66*0.0015d0 !spontaneous polarization
c     ????????????????????????????????????????????????????????????????
     
c     Stiffmatrix : cmat  
        do i=1,6
            do j=1,6
            cmat(i,j)=0.0d0
            enddo
        enddo
        
        cmat(1,1) = c11
        cmat(2,2) = c11
        cmat(3,3) = c11
        cmat(1,2) = c12
        cmat(1,3) = c12
        cmat(2,1) = c12
        cmat(2,3) = c12
        cmat(3,1) = c12
        cmat(3,2) = c12
        cmat(4,4) = c44
        cmat(5,5) = c44
        cmat(6,6) = c44
        
      write(6,*) "cmat before transformation"
      do i = 1, size(cmat,1)
          write(6,'(40G12.4)')  cmat(i,:)
      end do

      do i=1,6
        do j=1,6
          cmat1(i,j)=0.0d0
        enddo
      enddo

      do m=1,6
        do n=1,6
          do i=1,6
            interval=0.0d0
            do j=1,6
              interval=interval+cmat(i,j)*Trb(j,n)
            enddo
            cmat1(m,n)=cmat1(m,n)+Trb(i,m)*interval
          enddo
        enddo
      enddo

      do i=1,6
        do j=1,6
          cmat(i,j)=cmat1(i,j)
        enddo
      enddo
      
      write(6,*) "cmat after transformation"
      do i = 1, size(cmat,1)
          write(6,'(40G12.4)')  cmat(i,:)
      end do
      
!??? How to form this matrix ???!
!??? Is this correrct way ???!

c....... Stiffmatrix kapmat
      do i=1,3
        do j=1,3
        kapmat(i,j)=0.0d0
        enddo
      enddo

        kapmat(1,1)=Kapa
        kapmat(2,2)=Kapa
        kapmat(3,3)=Kapa
        
      write(6,*) "kapmat before transformation"
      do i = 1, size(kapmat,1)
          write(6,'(40G12.4)')  kapmat(i,:)
      end do

      do i=1,3
        do j=1,3
          kapmat1(i,j)=0.0d0
        enddo
      enddo

      do m=1,3
        do n=1,3
          do i=1,3
            interval=0.0d0
            do j=1,3
              interval=interval+kapmat(i,j)*Tra(j,n)
            enddo
            kapmat1(m,n)=kapmat1(m,n)+Tra(i,m)*interval
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
          kapmat(i,j)=kapmat1(i,j)
        enddo
      enddo

      write(6,*) "kapmat after transformation"
      do i = 1, size(kapmat,1)
          write(6,'(40G12.4)')  kapmat(i,:)
      end do
      
c....... Stiffmatrix grdmat

!??? How to form this matrix ???!
      do i=1,9
        do j=1,9
        grdmat(i,j)=0.0d0
        enddo
      enddo
        !Nanostructures in Ferroelectric Films for Energy Applications Page.no. 96
        !Equation for Gradient energy
        grdmat(1,1) = G11
        grdmat(2,2) = G11
        grdmat(3,3) = G11
        grdmat(1,2) = G12
        grdmat(1,3) = G12
        grdmat(2,1) = G12
        grdmat(2,3) = G12
        grdmat(3,1) = G12
        grdmat(3,2) = G12
        grdmat(4,4) = G44+G44p
        grdmat(5,5) = G44+G44p 
        grdmat(6,6) = G44+G44p
        grdmat(7,7) = G44+G44p 
        grdmat(8,8) = G44+G44p 
        grdmat(9,9) = G44+G44p
        grdmat(4,5) = G44-G44p
        grdmat(5,4) = G44-G44p
        grdmat(6,7) = G44-G44p
        grdmat(7,6) = G44-G44p
        grdmat(8,9) = G44-G44p
        grdmat(9,8) = G44-G44p
        
      write(6,*) "grdmat before transformation"
      do i = 1, size(grdmat,1)
          write(6,'(40G12.4)')  grdmat(i,:)
      end do

      do i=1,9
        do j=1,9
          grdmat1(i,j)=0.0d0
        enddo
      enddo

      do m=1,9
        do n=1,9
          do i=1,9
            interval=0.0d0
            do j=1,9
              interval=interval+grdmat(i,j)*Trc(j,n)
            enddo
            grdmat1(m,n)=grdmat1(m,n)+Trc(i,m)*interval
          enddo
        enddo
      enddo

      do i=1,9
        do j=1,9
           grdmat(i,j)=grdmat1(i,j)
        enddo
      enddo
      
      write(6,*) "grdmat after transformation"
      do i = 1, size(grdmat,1)
          write(6,'(40G12.4)')  grdmat(i,:)
      end do

      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c....... coordinate transform matrix from x-axis to X-axis             c
c....... global coordinate system x,y                                  c 
c....... material coordinate system X,Y                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trans_matrix01(Tra,Trb,Trc,xy0)
c     Tra(3,3): Transform matrix for coordinate or displacment vectors 
c     Trb(6,6): Transform matrix for stress or strain vectors
c     Trc(9,9): Transform matrix for polarization gradient vector 
c     Tra =[cos(X,x) cos(X,y); cos(Y,x) cos(Y,y)]
c     [X Y Z]' = Tra*[x y z]';
c     [s_XX s_YY s_ZZ s_YZ s_XZ s_XY]' = Trb*[s_xx s_yy s_zz s_yz s_xz s_xy ]';
c     [P_X,X  P_Y,Y  P_Z,Z  P_X,Y P_Y,X P_Y,Z P_Z,Y P_X,Z P_Z,X]' = 
c           Trc*[P_x,x  P_y,y  P_z,z  P_x,y P_y,x P_y,z P_z,y P_x,z P_z,x ]';
      implicit none
      double precision Tra(3,3),Trb(6,6),Trc(9,9)
      double precision xi,theta,phi,xy0(3)
      double precision pi,deg
      integer i,j,m,n,k,l
      
      pi = 3.141592653589793
      deg = 0.0
      xi = 0.       ! A rotation through angle ψ (xi) about the z-axis.
      theta = 0.    ! A rotation of angle θ (theta) about the new y-axis.
      phi = 0.      ! a second rotation of ϕ (phi) about the now-tilted z-axis.


  
c~       ! Tra(I,j)
c~       Tra(1,1)= cos(arfa) !cos(X,x)
c~       Tra(1,2)= sin(arfa) !cos(X,y)
c~       Tra(2,1)=-sin(arfa) !cos(Y,x)
c~       Tra(2,2)= cos(arfa) !cos(Y,y)

c       TRANSFORM MATRIX FOR VECTORS According to ROE Convention
c       https://www.continuummechanics.org/transformmatrix.html

        Tra(1,1)= cos(xi)*cos(theta)*cos(phi)-sin(xi)*sin(phi)
        Tra(1,2)= sin(xi)*cos(theta)*cos(phi)+cos(xi)*sin(phi)
        Tra(1,3)=-sin(theta)*cos(phi)
        Tra(2,1)=-cos(xi)*cos(theta)*sin(phi)-sin(xi)*cos(phi)
        Tra(2,2)=-sin(xi)*cos(theta)*sin(phi)+cos(xi)*cos(phi)
        Tra(2,3)= sin(theta)*sin(phi)
        Tra(3,1)= cos(xi)*sin(theta) 
        Tra(3,2)= sin(xi)*sin(theta) 
        Tra(3,3)= cos(theta)
        
      write(6,*) "TRA"
      do i = 1, size(Tra,1)
          write(6,'(40G12.4)')  Tra(i,:)
      end do


      ! Multiply Q.Tensor.Q_T and reduce it to 6*6 matrix in voigt notation  
      ! Trb(I,j)
      Trb(1,1)= Tra(1,1)*Tra(1,1) 
      Trb(1,2)= Tra(1,2)*Tra(1,2) 
      Trb(1,3)= Tra(1,3)*Tra(1,3)
      Trb(1,4)= Tra(1,2)*Tra(1,3) 
      Trb(1,5)= Tra(1,1)*Tra(1,3) 
      Trb(1,6)= Tra(1,1)*Tra(1,2)
       
      Trb(2,1)= Tra(2,1)*Tra(2,1) 
      Trb(2,2)= Tra(2,2)*Tra(2,2) 
      Trb(2,3)= Tra(2,3)*Tra(2,3)
      Trb(2,4)= Tra(2,2)*Tra(2,3) 
      Trb(2,5)= Tra(2,1)*Tra(2,3) 
      Trb(2,6)= Tra(2,1)*Tra(2,2) 
       
      Trb(3,1)= Tra(3,1)*Tra(3,1) 
      Trb(3,2)= Tra(3,2)*Tra(3,2) 
      Trb(3,3)= Tra(3,3)*Tra(3,3)
      Trb(3,4)= Tra(3,2)*Tra(3,3) 
      Trb(3,5)= Tra(3,1)*Tra(3,3) 
      Trb(3,6)= Tra(3,1)*Tra(3,2)
      
      Trb(4,1)= 2*Tra(3,1)*Tra(2,1) 
      Trb(4,2)= 2*Tra(3,2)*Tra(2,2) 
      Trb(4,3)= 2*Tra(3,3)*Tra(2,3)
      Trb(4,4)= Tra(3,3)*Tra(2,2)+Tra(3,2)*Tra(2,3) 
      Trb(4,5)= Tra(3,1)*Tra(2,3)+Tra(3,3)*Tra(2,1) 
      Trb(4,6)= Tra(3,2)*Tra(2,1)+Tra(3,1)*Tra(2,2)
      
      Trb(5,1)= 2*Tra(3,1)*Tra(1,1) 
      Trb(5,2)= 2*Tra(3,2)*Tra(1,2) 
      Trb(5,3)= 2*Tra(3,3)*Tra(1,3)
      Trb(5,4)= Tra(3,3)*Tra(1,2)+Tra(3,2)*Tra(1,3) 
      Trb(5,5)= Tra(3,1)*Tra(1,3)+Tra(3,3)*Tra(1,1) 
      Trb(5,6)= Tra(3,2)*Tra(1,1)+Tra(3,1)*Tra(1,2) 
       
      Trb(6,1)= 2*Tra(2,1)*Tra(1,1) 
      Trb(6,2)= 2*Tra(2,2)*Tra(1,2) 
      Trb(6,3)= 2*Tra(2,3)*Tra(1,3)
      Trb(6,4)= Tra(2,3)*Tra(1,2)+Tra(2,2)*Tra(1,3) 
      Trb(6,5)= Tra(2,3)*Tra(1,1)+Tra(2,1)*Tra(1,3) 
      Trb(6,6)= Tra(2,2)*Tra(1,1)+Tra(2,1)*Tra(1,2)
      
      write(6,*) "TRB"
      do i = 1, size(Trb,1)
          write(6,'(40G12.4)')  Trb(i,:)
      end do
  

c      do i=1,3
c         write(20,102) Trb(i,1),Trb(i,2),Trb(i,3)
c      enddo

      ! Multiply Q.Tensor.Q_T and reduce it to 9*9 matrix for 9 tensor components
      ![(1,1),(2,2),(3,3),(1,2),(2,1),(2,3),(3,2),(1,3),(3,1)]^T 
      ! Trc(I,j)
      Trc(1,1)= Tra(1,1)*Tra(1,1) 
      Trc(1,2)= Tra(1,2)*Tra(1,2) 
      Trc(1,3)= Tra(1,3)*Tra(1,3) 
      Trc(1,4)= Tra(1,1)*Tra(1,2) 
      Trc(1,5)= Tra(1,2)*Tra(1,1) 
      Trc(1,6)= Tra(1,2)*Tra(1,3)
      Trc(1,7)= Tra(1,3)*Tra(1,2) 
      Trc(1,8)= Tra(1,1)*Tra(1,3) 
      Trc(1,9)= Tra(1,3)*Tra(1,1)
      
      Trc(2,1)= Tra(2,1)*Tra(2,1) 
      Trc(2,2)= Tra(2,2)*Tra(2,2) 
      Trc(2,3)= Tra(2,3)*Tra(2,3) 
      Trc(2,4)= Tra(2,1)*Tra(2,2) 
      Trc(2,5)= Tra(2,2)*Tra(2,1) 
      Trc(2,6)= Tra(2,2)*Tra(2,3)
      Trc(2,7)= Tra(2,3)*Tra(2,2) 
      Trc(2,8)= Tra(2,1)*Tra(2,3) 
      Trc(2,9)= Tra(2,3)*Tra(2,1)
      
      Trc(3,1)= Tra(3,1)*Tra(3,1) 
      Trc(3,2)= Tra(3,2)*Tra(3,2) 
      Trc(3,3)= Tra(3,3)*Tra(3,3) 
      Trc(3,4)= Tra(3,1)*Tra(3,2) 
      Trc(3,5)= Tra(3,2)*Tra(3,1) 
      Trc(3,6)= Tra(3,2)*Tra(3,3)
      Trc(3,7)= Tra(3,3)*Tra(3,2) 
      Trc(3,8)= Tra(3,1)*Tra(3,3) 
      Trc(3,9)= Tra(3,3)*Tra(3,1)
      
      Trc(4,1)= Tra(1,1)*Tra(2,1) 
      Trc(4,2)= Tra(1,2)*Tra(2,2) 
      Trc(4,3)= Tra(1,3)*Tra(2,3) 
      Trc(4,4)= Tra(1,1)*Tra(2,2) 
      Trc(4,5)= Tra(1,2)*Tra(2,1) 
      Trc(4,6)= Tra(1,2)*Tra(2,3)
      Trc(4,7)= Tra(1,3)*Tra(2,2) 
      Trc(4,8)= Tra(1,1)*Tra(2,3) 
      Trc(4,9)= Tra(1,3)*Tra(2,1)
      
      Trc(5,1)= Tra(2,1)*Tra(1,1) 
      Trc(5,2)= Tra(2,2)*Tra(1,2) 
      Trc(5,3)= Tra(2,3)*Tra(1,3) 
      Trc(5,4)= Tra(2,1)*Tra(1,2) 
      Trc(5,5)= Tra(2,2)*Tra(1,1) 
      Trc(5,6)= Tra(2,2)*Tra(1,3)
      Trc(5,7)= Tra(2,3)*Tra(1,2) 
      Trc(5,8)= Tra(2,1)*Tra(1,3) 
      Trc(5,9)= Tra(2,3)*Tra(1,1)
      
      Trc(6,1)= Tra(2,1)*Tra(3,1) 
      Trc(6,2)= Tra(2,2)*Tra(3,2) 
      Trc(6,3)= Tra(2,3)*Tra(3,3) 
      Trc(6,4)= Tra(2,1)*Tra(3,2) 
      Trc(6,5)= Tra(2,2)*Tra(3,1) 
      Trc(6,6)= Tra(2,2)*Tra(3,3)
      Trc(6,7)= Tra(2,3)*Tra(3,2) 
      Trc(6,8)= Tra(2,1)*Tra(3,3) 
      Trc(6,9)= Tra(2,3)*Tra(3,1)
      
      Trc(7,1)= Tra(3,1)*Tra(2,1) 
      Trc(7,2)= Tra(3,2)*Tra(2,2) 
      Trc(7,3)= Tra(3,3)*Tra(2,3) 
      Trc(7,4)= Tra(3,1)*Tra(2,2) 
      Trc(7,5)= Tra(3,2)*Tra(2,1) 
      Trc(7,6)= Tra(3,2)*Tra(2,3)
      Trc(7,7)= Tra(3,3)*Tra(2,2) 
      Trc(7,8)= Tra(3,1)*Tra(2,3) 
      Trc(7,9)= Tra(3,3)*Tra(2,1)
      
      Trc(8,1)= Tra(1,1)*Tra(3,1) 
      Trc(8,2)= Tra(1,2)*Tra(3,2) 
      Trc(8,3)= Tra(1,3)*Tra(3,3) 
      Trc(8,4)= Tra(1,1)*Tra(3,2) 
      Trc(8,5)= Tra(1,2)*Tra(3,1) 
      Trc(8,6)= Tra(1,2)*Tra(3,3)
      Trc(8,7)= Tra(1,3)*Tra(3,2) 
      Trc(8,8)= Tra(1,1)*Tra(3,3) 
      Trc(8,9)= Tra(1,3)*Tra(3,1)
      
      Trc(9,1)= Tra(3,1)*Tra(1,1) 
      Trc(9,2)= Tra(3,2)*Tra(1,2) 
      Trc(9,3)= Tra(3,3)*Tra(1,3) 
      Trc(9,4)= Tra(3,1)*Tra(1,2) 
      Trc(9,5)= Tra(3,2)*Tra(1,1) 
      Trc(9,6)= Tra(3,2)*Tra(1,3)
      Trc(9,7)= Tra(3,3)*Tra(1,2) 
      Trc(9,8)= Tra(3,1)*Tra(1,3) 
      Trc(9,9)= Tra(3,3)*Tra(1,1)
      
      write(6,*) "TRC"
      do i = 1, size(Trc,1)
          write(6,'(40G12.4)')  Trc(i,:)
      end do
      
c      do i=1,4
c         write(10,101) Trc(i,1),Trc(i,2),Trc(i,3),Trc(i,4)
c      enddo

101    format(4(2x,f6.3))
102    format(3(2x,f6.3))

!Have to reorder AMATRX
!Have to check how to formulate Bp matrix and accordingly grdmat
!Have to ask transformation matrix
!Have to modify material matrix

      return
      end
