!! Copyright (C) 2016 M. Moreira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!!

module run_cells_m

  use global_m
  use derivatives_m
  use sim_init_m
  use misc_m

  implicit none

  private

  public :: run_cells


  contains

    subroutine run_cells(sim_id, iseed, porosity, rinit)

      implicit none

      real, intent(in) :: porosity
      integer, intent(in) :: iseed
      character(3), intent(in) :: sim_id
      integer, intent(inout) :: rinit(3)
      dir_name = sim_id

      ! initializing parameters
      call  parameters_init(cell_radius, ntype, density, interface_width, tstep, dt, Lsize, dr, dir_name, iseed,&
           np_bndry, depletion_weight, adh11, adhs, chi(1), metcoef, output_period, periodic)

      ! number of points in the mesh
      np = 8*Lsize(1)*Lsize(2)*Lsize(3) ! number of points
      np_tt = 8*(Lsize(1)+2*np_bndry)*(Lsize(2)+2*np_bndry)*(Lsize(3)+2*np_bndry) ! number of points plus boundary points

      ALLOCATE(ncell(0:ntype))

      ncell(0:ntype) = 1 ! (np*density(1:ntype))/(4*cell_radius**2)
      tcell = 1!sum(ncell(1:ntype))

      ! partition mesh
      Lsize_part(1:3) = (/2*cell_radius, 2*cell_radius, 2*cell_radius/)
      np_part = 8*Lsize_part(1)*Lsize_part(2)*Lsize_part(3)
      np_part_bndry = 4*max(Lsize(1),Lsize(2),Lsize(3))
      np_part_tt = 8*(Lsize_part(1) + 2*np_part_bndry)*(Lsize_part(2) + 2*np_part_bndry)*(Lsize_part(3) + 2*np_part_bndry)
      ! allocating matrices and vectors

      ALLOCATE(lxyz(np_tt,1:3))
      ALLOCATE(lxyz_inv(-Lsize(1)-np_bndry:Lsize(1)+np_bndry, &
                        -Lsize(2)-np_bndry:Lsize(2)+np_bndry,&
                        -Lsize(3)-np_bndry:Lsize(3)+np_bndry))
      ALLOCATE(lxyz_part(np_part_tt,1:3))
      ALLOCATE(lxyz_inv_part(-Lsize_part(1)-np_part_bndry:Lsize_part(1)+np_part_bndry, &
                        -Lsize_part(2)-np_part_bndry:Lsize_part(2)+np_part_bndry, &
                        -Lsize_part(3)-np_part_bndry:Lsize_part(3)+np_part_bndry))


      ALLOCATE(aux(np,ntype))
      ALLOCATE(gg(1:np,1:3))
      ALLOCATE(s(0:np))
      ALLOCATE(seq(0:np))
      ALLOCATE(shfield(0:np))
      ALLOCATE(shfield_lapl(0:np))
      ALLOCATE(chem(0:np))
      ALLOCATE(mmps(0:np))
      ALLOCATE(mmps_lapl(1:np))
      ALLOCATE(sources(1:np))
      ALLOCATE(chem_lapl(1:np))
      ALLOCATE(gchem(1:np,1:3))
      ALLOCATE(r(1:1000))
      ALLOCATE(velocity(1:tstep))
      ALLOCATE(path(1:np))
      ALLOCATE(vinst(1:tstep/20))
      !call system('rm 001/*')
      mmps(:) = 0.d0
      mmps_lapl = 0.d0
!      ALLOCATE(gammaw(ntype))
!      gammaw(1:3) = (/ 1.0, 2.0 /)

      ! initializing space matrices
      call space_init(Lsize, lxyz, lxyz_inv, np_bndry, np, periodic)  ! auxiliar fields
      ! single cell fields Dirichlet boundary condition
      call space_init(Lsize_part, lxyz_part, lxyz_inv_part, np_part_bndry, np_part, .false.)

      call chemical_init(chem, sources, np, Lsize, lxyz, lxyz_inv)
      ! generating sphere points
      call gen_cell_points(cell_radius,sphere,np_sphere)


      ! adhesions
      adh1(1,1) = adh11
      adh1(2,2) = adh11
    !  adh1(1,2) = adh12
  !    adh1(2,1) = adh1(1,2)
      adh2(1) = 1.49*adh1(1,1)
      adh2(2) = 1.49*adh1(2,2)
      alpha_p = 0.10 ! production of mmps
      beta_d = 1.00 ! ecm degradation
      volume_target = (4.d0/3.d0)*M_PI*cell_radius**3
      vol_lagrangian = 1.00
      scoef = 1.0
      !metcoef = 1.d0
      chi(2) = 4.d0
      !chi(1) = 10.d0
      chemresponse(1) = 1.d0
      !adhs = 0.5
      ! initializing simulation box

      if(tcell.eq.2) then
        r(1) = lxyz_inv(0,-7,0) !(-6,10)
        r(2) = lxyz_inv(0,7,0) !(8,-10)
      elseif(tcell.eq.1) then

        ncell(1) = 1
        !rinit(1) = -Lsize(1)+int(cell_radius)+3

        r(1) = lxyz_inv(rinit(1),rinit(2),rinit(3))
        r1(1:3) =  rinit(1:3)
      elseif(tcell.eq.0) then

        do j = -Lsize(2)+int(cell_radius), Lsize(2)-1
          do k = -Lsize(3)+int(cell_radius), Lsize(3)-1-int(cell_radius)
            superposition = .false.
            do i=1,tcell
              distance(1:3) = (/ 0, j-lxyz(r(i),2), k-lxyz(r(i),3)  /)
              distance(1:3) = min(abs(distance(1:3)), 2*Lsize(1:3)-1-abs(distance(1:3)))
              if(sqrt(sum(1.0*distance(1:3)*distance(1:3))) .le. 2.d0*cell_radius) superposition = .true.
            end do
            if( .not.superposition) then

              tcell = tcell +  1
              r(tcell) = lxyz_inv(-Lsize(1)+int(cell_radius)+1,j,k)

            end if

          end do
        end do

      end if
      !r(tcell+1) = r(5)
      !r(5) = r(tcell)
      !r(tcell) = r(tcell+1)

      !ncell(1) = tcell-1
      !ncell(2) = 1

      ALLOCATE(cell(0:np_part,tcell))
      ALLOCATE(adhesion(0:np_part,tcell))
      ALLOCATE(hfield(0:np_part,tcell))
      ALLOCATE(hfield_lapl(0:np_part,tcell))
      ALLOCATE(r_cm(1:tcell,1:3))
      ALLOCATE(r_cm_part(1:tcell,1:3))
      ALLOCATE(r_cmg(1:tcell,1:3))
      ALLOCATE(volume(1:tcell))

      cell(0,:)%phi = 0.d0
      cell(0,:)%mu = 0.d0
      cell(0,:)%lapl_mu = 0.d0
      cell(0,:)%lapl_phi = 0.d0
      cell(0,:)%lapl_h = 0.d0
      write(*,*) np_sphere
      call single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first=.true.)
      call cahn_hilliard(cell(0:np_part,1), 100, np_part, 1.0, lxyz_part, lxyz_inv_part, dr)
      if(tcell.gt.1) then
         call single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first=.false.)
      end if


      call hfield_calc(cell, aux, r, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np_part )

      !call dderivatives_lapl(cell(0:np,1)%phi, cell(1:np_part,1)%lapl_phi ,&
      !     np_part, lxyz_part, lxyz_inv_part, dr)

      ! printing header
      call print_header(Lsize, tcell, ntype, ncell, dir_name, periodic)
      if(tcell.eq.2) then
        !call cm_calc(r_cm, cell, tcell, np, np_part, r, lxyz, lxyz_inv, lxyz_inv_part)
        write(*,'(A,F10.2,F10.2,F10.2)') "Cell 1 - Initial Position", lxyz(INT(r(1)),1:3)*dr(3)
        write(*,'(A,F10.2,F10.2,F10.2)') "Cell 2 - Initial Position", lxyz(INT(r(2)),1:3)*dr(3)
        write(*,'(A,F10.2)') "Distance between the cells", REAL(lxyz(r(1),2))-REAL(lxyz(r(2),2))
      end if

      call gen_cell_points(1.0,porous,np_porous)
      call substrate_init(porosity, s, np, porous, np_porous, Lsize, lxyz, lxyz_inv, iseed)
  !     cs = 0.d0
  !     call cross_section(s, cs, lxyz_inv, Lsize)
  !     write(*,*) cs
  !  s(:) = 0
  !  seq(:) = 0
  !      do ip=1, np
  !        if(abs(lxyz(ip,2)) .eq. abs(lxyz(r(1),2)+2*cell_radius-1) .and. lxyz(ip,3) .eq. lxyz(r(1),3)) then
  !        !if(lxyz(ip,2) .eq. lxyz(r(1),2)+2*cell_radius .and. lxyz(ip,3) .eq. lxyz(r(1),3)) then
  !          s(ip) = 1.d0
  !          !seq(ip) = 1.d0
  !          do k=1,np_porous
  !            s(lxyz_inv( lxyz(ip,1)+porous(k,1),lxyz(ip,2)+porous(k,2),lxyz(ip,3)+porous(k,3) ) ) = 1.d0
  !            if(abs(lxyz(ip,1)).eq.0 .or. abs(lxyz(ip,1)).eq.1 .or. abs(lxyz(ip,1)).eq.2 &
  !            .or. abs(lxyz(ip,1)).eq.3.or. abs(lxyz(ip,1)).eq.4.or. abs(lxyz(ip,1)).eq.5) then
  !              seq(lxyz_inv( lxyz(ip,1)+porous(k,1),lxyz(ip,2)+2+porous(k,2),lxyz(ip,3)+porous(k,3) ) ) = 1.d0
  !            else
  !              seq(lxyz_inv( lxyz(ip,1)+porous(k,1),lxyz(ip,2)+porous(k,2),lxyz(ip,3)+porous(k,3) ) ) = 1.d0
  !            end if
  !          end do
  !        end if
   !
  !      end do

      call sch(s(0:np), 100, np, 0.25, lxyz, lxyz_inv, dr)

      !do icell=1,tcell
      !  do i=1, np_sphere
      !    s(lxyz_inv(lxyz(r(icell),1)+sphere(i,1),lxyz(r(icell),2)+sphere(i,2),lxyz(r(icell),3)+sphere(i,3) ) ) = 0.d0
      !  end do
      !end do
      call density_calc(s, Lsize, lxyz_inv, dens, deviation)
      write(*,'(A,F10.4,A,F10.4)') "  Fiber density:", dens, "  Deviation:", deviation
!!!! buraco
!do ip=1, np
!  call vec_global2local(ip_part, r(1), ip, lxyz, lxyz_inv, lxyz_inv_part)
!  if(cell(ip_part,1)%phi.gt.0) then
!    s(ip) = 0.d0
!  end if
!end do
      !do ip=1, np
      !  if(lxyz(ip,1)<=-Lsize(1)+2*cell_radius+1) then
      !    s(ip) = 0.d0
      !  end if
      !end do
      ! calculating the h(s) function
      do ip=1, np
        shfield(ip) = hfunc(s(ip))
      end do
      ! calculating the laplacian of the h(s) field substrate
      call dderivatives_lapl(shfield(0:np), shfield_lapl(1:np), np, lxyz, lxyz_inv, dr)


      call output(s, 300, 7, 'subs0000',dir_name, 1, ntype, np, lxyz,dr)

      ! saving the initial condition
      !call output_aux(cell, 100, 7, 'phi0000',dir_name, 1, 1, np_part, lxyz_part)
      call output_aux(aux, 200, 7, 'aux0000',dir_name, 1, ntype, np, lxyz,dr)

      call output(chem, 300, 7, 'chem0000',dir_name, 1, ntype, np, lxyz,dr)
      cell(0,:)%phi = 0.d0
      cell(0,:)%mu = 0.d0

      nstep = 0
      counter = 0
      cm_calc_counter = 0
      vcounter = 0
      path(:) = 0.d0

      OPEN (UNIT=700,FILE=dir_name//'/vinst.dat')
      write(700,'(A)') "# time (min)    velocity (mu m/min)"
      write(*,'(A)') "Initiating the core program... "


      call dderivatives_grad(chem, gchem, np, lxyz, lxyz_inv, dr)
      !gchem(:,1) = s(:)*gchem(:,1)
      !gchem(:,2) = s(:)*gchem(:,2)
      !gchem(:,3) = s(:)*gchem(:,3)
      !gchem(:,1) = gchem(:,1)*s(:)
      !gchem(:,2) = gchem(:,2)*s(:)
      !gchem(:,3) = gchem(:,3)*s(:)

      do while(nstep<=tstep)
         nstep = nstep + 1

         call CPU_TIME(time_init)

         ! calculating gradient of vegf
         !call dderivatives_grad(cell, gg, np, lxyz, lxyz_inv, dr)

         call hfield_calc(cell, aux, r, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np_part )

         ! calculating the gradient of chemical concentration
         !call dderivatives_grad(chem, gchem, np, lxyz, lxyz_inv, dr) ! como os campos nao se alteram, vou por essa derivada la fora
        ! call dderivatives_lapl(chem(0:np), chem_lapl(1:np), np, lxyz, lxyz_inv, dr)
         ! chemical diffusion
         !do ip=1,np
        !    if(sources(ip)<1.d0) then
        !      chem(ip) = chem(ip) + dt*chem_lapl(ip)
        !    end if
        !  end do

        ! calculating the laplacian of MMPs
        call dderivatives_lapl(mmps(0:np), mmps_lapl(1:np), np, lxyz, lxyz_inv, dr)

        do ip=1,np
          call vec_global2local(ip_part, r(1), ip, lxyz, lxyz_inv, lxyz_inv_part)
          mmps(ip) = mmps(ip) + dt*(mmps_lapl(ip)+0.01*heaviside(2.d0*cell(ip_part,1)%phi-1.d0))!&
          !  -0.005*mmps(ip)*heaviside(-(2.d0*cell(ip_part,1)%phi-1.d0)))
          if(mmps(ip).gt.0.09) s(ip) = s(ip) - dt*(10.d0*s(ip)*mmps(ip))
        end do
         ! calculating laplacian of phi

         do icell = 1, tcell
            call dderivatives_lapl(cell(0:np_part,icell)%phi, cell(1:np_part,icell)%lapl_phi ,&
                 np_part, lxyz_part, lxyz_inv_part, dr)
         end do

         call volume_calc(volume, cell, tcell, np, np_part, lxyz, lxyz_inv, lxyz_inv_part)

         ! Chemical potential:  phi**3 - phi - epsilon*laplacian(phi )
         adhesion(:,:) = 0.d0
         do ip = 1, np_part
            do icell=1, tcell

              call vec_local2global(ip_global, r(icell), ip, lxyz, lxyz_inv, lxyz_part)

              !ip_global = lxyz_inv( lxyz contribution for the chemical energy
              ! the adhesion ter_part(ip,1) + lxyz(r(icell),1), lxyz_part(ip,2) + lxyz(r(icell),2))
              fnu = 0.d0
              do itype=1, ntype
                 ! functiona f(u_m,s,phi) - > f(cell(:,icell)%phi,aux(:,itype)%phi)
                 fnu = fnu + aux(ip_global,itype)%phi - hfunc( cell(ip,icell)%phi )*deltak(cell(ip,icell)%itype,itype)
                 ! function g_int(phi_m,aux) - calculate the adhesion term
                 adhesion(ip,icell) = adhesion(ip,icell) + adh1(cell(ip,icell)%itype,itype)* &
                 (aux(ip_global,itype)%phi - hfunc( cell(ip,icell)%phi )*deltak(cell(ip,icell)%itype,itype))
                 !eta(cell(ip,icell)%itype,itype)
              end do
              ! summing the first contribution for the chemical energy
              ! the adhesion term will be summed after
          !    hfield(ip,icell) = hfunc(cell(ip,icell)%phi)
              if(cell(ip,icell)%itype.eq.1) then
                chemresponse(1) = &
                  cell(lxyz_inv_part(lxyz_part(ip,1)+1,lxyz_part(ip,2),lxyz_part(ip,3)),icell)%phi*&
                  gchem(lxyz_inv(lxyz(ip_global,1)+1,lxyz(ip_global,2),lxyz(ip_global,3)),1) - &
                  cell(lxyz_inv_part(lxyz_part(ip,1)-1,lxyz_part(ip,2),lxyz_part(ip,3)),icell)%phi*&
                  gchem(lxyz_inv(lxyz(ip_global,1)-1,lxyz(ip_global,2),lxyz(ip_global,3)),1) + &
                  cell(lxyz_inv_part(lxyz_part(ip,1),lxyz_part(ip,2)+1,lxyz_part(ip,3)),icell)%phi*&
                  gchem(lxyz_inv(lxyz(ip_global,1),lxyz(ip_global,2)+1,lxyz(ip_global,3)),2) - &
                  cell(lxyz_inv_part(lxyz_part(ip,1),lxyz_part(ip,2)-1,lxyz_part(ip,3)),icell)%phi*&
                  gchem(lxyz_inv(lxyz(ip_global,1),lxyz(ip_global,2)-1,lxyz(ip_global,3)),2) + &
                  cell(lxyz_inv_part(lxyz_part(ip,1),lxyz_part(ip,2),lxyz_part(ip,3)+1),icell)%phi*&
                  gchem(lxyz_inv(lxyz(ip_global,1),lxyz(ip_global,2),lxyz(ip_global,3)+1),3) - &
                  cell(lxyz_inv_part(lxyz_part(ip,1),lxyz_part(ip,2),lxyz_part(ip,3)-1),icell)%phi*&
                  gchem(lxyz_inv(lxyz(ip_global,1),lxyz(ip_global,2),lxyz(ip_global,3)-1),3)
              end if
              cell(ip,icell)%mu = interface_width*cell(ip,icell)%lapl_phi +&
                   cell(ip,icell)%phi*(1.d0-cell(ip,icell)%phi)*(cell(ip,icell)%phi - 0.50 + &
                   vol_lagrangian*(volume_target-volume(icell)) - depletion_weight*fnu - scoef*hfunc(s(ip_global))) -&
                   chi(cell(ip,icell)%itype)*chemresponse(cell(ip,icell)%itype)/2.d0 + metcoef*(8.0-16.0*ran2(iseed) )*&
                   cell(ip,icell)%phi*(1.d0-cell(ip,icell)%phi)


            end do
         end do


         ! calculating laplacian(Gamma_l - hfunc(phi_m)*delta_k)

         do icell=1, tcell
           call dderivatives_lapl(hfield(0:np_part,icell), hfield_lapl(1:np_part,icell), &
           np_part, lxyz_part, lxyz_inv_part, dr)
           call dderivatives_lapl(adhesion(0:np_part,icell), cell(1:np_part,icell)%lapl_h, &
           np_part, lxyz_part, lxyz_inv_part, dr)
         end do

         ! including the cellular adhesion in the chemical potential

         do ip=1, np_part
           do icell=1, tcell
             call vec_local2global(ip_global, r(icell), ip, lxyz, lxyz_inv, lxyz_part)
             cell(ip,icell)%mu =  cell(ip,icell)%mu + cell(ip,icell)%phi*(1.0-cell(ip,icell)%phi)*&
                                  (cell(ip,icell)%lapl_h + adh2(cell(ip,icell)%itype)*hfield_lapl(ip,icell)) + &
                                  adhs*cell(ip,icell)%phi*(1.d0 -cell(ip,icell)%phi)*shfield_lapl(ip_global)
           end do
         end do

         ! Calculating laplacian of mu
         !do icell = 1, tcell
        !    call dderivatives_lapl(cell(0:np_part,icell)%mu, cell(1:np_part,icell)%lapl_mu, np_part, lxyz_part, lxyz_inv_part, dr)
         !end do

         cell(1:np_part, 1:tcell)%phi = cell(1:np_part,1:tcell)%phi + dt*(cell(1:np_part,1:tcell)%mu)

         cm_calc_counter = cm_calc_counter + 1
         if(cm_calc_counter.eq.10) then

           rt(1:3) = r_cmg(1,1:3)
           cm_calc_counter = 0
           call cm_calc_local(r_cm_part, cell, tcell, np_part, r, lxyz_part, lxyz_inv_part)
           do icell=1, tcell
             ip = lxyz_inv_part(int(r_cm_part(icell,1)),int(r_cm_part(icell,2)),int(r_cm_part(icell,3)))
             call vec_local2global(ip_global, r(icell), ip, lxyz, lxyz_inv, lxyz_part)
             r_cm(icell,1:3) = real(lxyz(ip_global,1:3)) + FRACTION(r_cm_part(icell,1:3))
           end do
           path(ip_global) = 1.0
           r2(:) = r_cm(1,1:3)
           vcounter = vcounter + 1
           if(vcounter.eq.20) then
             vinst(np_vinst) = sqrt(sum( (r2(1:3) - r1(1:3))**2))*(1.25/(20.0*dt*0.26)) ! 0.26
             write(700,'(F10.4, F10.4)') nstep*dt*0.26, vinst(np_vinst)
             r1(:) = r2(:)
             vcounter = 0
           end if

           call move(cell, r_cm, r, np, np_part, tcell, lxyz, lxyz_part, lxyz_inv, lxyz_inv_part, Lsize)


        !   write(*,*) r_cm_part(1,1:3), nstep
        !print*, nstep
           !write(*,*) "cell 2", r_cm(2,1:3)
        end if



         call CPU_TIME(time_end)

         ctime = ctime + (time_end - time_init)

         if(nstep.eq.100) then
            ctime = (ctime*(tstep-nstep) )/6000.d0

            if( ctime>60.d0) then
               write(*,'(A,F10.2)') "Estimated time (hour): ",ctime/60.d0
            else
               write(*,'(A,F10.2)') "Estimated time (min): ",ctime
            end if
         end if



         ! output

         counter = counter + 1


         if(counter.eq.output_period) then
           call sch(s(0:np), 100, np, 0.25, lxyz, lxyz_inv, dr)
            write(*,*) nstep
            counter = 0

          !  write(*,*) volume(1:tcell)
            write(file_name,'(I6)') nstep
          !  call output_aux(cell, nstep+1, 6, trim(file_name),dir_name, 1, 1, np_part, lxyz_part)

            OPEN (UNIT=100,FILE=dir_name//'/phi'//trim(file_name)//'.xyz')
            OPEN (UNIT=333,FILE=dir_name//'/mmps'//trim(file_name)//'.xyz')
            OPEN (UNIT=444, FILE=dir_name//'/s'//trim(file_name)//'.xyz')
            !OPEN(UNIT=1000, FILE=dir_name//'/path'//trim(file_name)//'.xyz')
            do ip=1, np
            !  if(path(ip).gt.0)  write(1000,'(F10.2,F10.2,F10.2,F10.4)') lxyz(ip,1:3)*dr(1:3), path(ip)
              if(mmps(ip).gt.0.001) write(333, '(I10,I10,I10,F10.2)') lxyz(ip,1:3), mmps(ip)
              if(s(ip).gt.0.001) write(444, '(I10,I10,I10,F10.2)') lxyz(ip,1:3), s(ip)
              do itype = 1, ntype
                 if(aux(ip,itype)%phi>0.0) then
                    write(100,'(F10.2,F10.2,F10.2, F10.2,I10)') lxyz(ip,1:3)*dr(1:3),aux(ip,itype)%phi, itype
                !    if(lxyz(ip,1).eq.lxyz(r(1),1)) then

                !      call vec_global2local(ip_part, r(1), ip, lxyz, lxyz_inv, lxyz_inv_part)
                !      call vec_global2local(ip_part2, r(2), ip, lxyz, lxyz_inv, lxyz_inv_part)

            !          write(nstep+2,'(I10,F10.2,F10.2)') lxyz(ip,2), cell(ip_part,1)%phi, cell(ip_part2,2)%phi

                      !write(nstep+2,'(I10,F10.2,F10.2,F10.2)') lxyz(ip,2), aux(ip,itype)%phi
                !    end if

                 end if
              end do
            end do
            close(333)
            close(444)
            close(100)
            !close(1000)
            !close(nstep+2)

         end if
         ! end of the output


      end do

      if(tcell.eq.2) then
        call cm_calc(r_cm, cell, tcell, np, np_part, r, lxyz, lxyz_inv, lxyz_inv_part)
        write(*,'(A,F10.2,F10.2)') "Cell 1 - End Position", r_cm(1,2)
        write(*,'(A,F10.2,F10.2)') "Cell 2 - End Position", r_cm(2,2)
        write(*,'(A,F10.2)') "Distance between the cells", r_cm(1,2)-r_cm(2,2)
      end if


      velocity(1) = (sum(path(1:np))*1.25)/(nstep*dt*0.26)
      write(*,'(A,F10.4)') "mean velocity: ", velocity(1)
      !DEALLOCATE(ncell)

      !DEALLOCATE(lxyz)
      !DEALLOCATE(lxyz_inv)
      !DEALLOCATE(lxyz_part)
      !DEALLOCATE(lxyz_inv_part)

      !DEALLOCATE(cell)
      !DEALLOCATE(aux)
      !DEALLOCATE(adhesion)
      !DEALLOCATE(hfield)
      !DEALLOCATE(hfield_lapl)
      !DEALLOCATE(gg)
      !DEALLOCATE(chem)
      !DEALLOCATE(gchem)
      !DEALLOCATE(r)
      !DEALLOCATE(r_cm_part)
      !DEALLOCATE(r_cm)
      !DEALLOCATE(volume)

    end subroutine run_cells


    subroutine cm_calc_local(r_cm_part, f, tcell, np_part, r, lxyz_part, lxyz_inv_part)

      implicit none

      integer, intent(in) :: tcell, np_part
      integer, allocatable, intent(in) :: lxyz_part(:,:), lxyz_inv_part(:,:,:)
      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      real, allocatable, intent(inout) :: r_cm_part(:,:)
      integer, allocatable, intent(inout) :: r(:)

      integer :: ip, icell, ncell, ip_part, ip_global_m
      real :: volume(1:tcell)

      volume(:) = 0.d0
      r_cm_part(:,:) = 0.d0
      do ip_part=1, np_part

         do icell=1, tcell
            !ip_global_m = r(icell)

            !call vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

            if(f(ip_part,icell)%phi.gt.0.d0) then
               volume(icell) = volume(icell) + f(ip_part,icell)%phi
               r_cm_part(icell,1:3) = r_cm_part(icell,1:3) + f(ip_part,icell)%phi*lxyz_part(ip_part,1:3)
            end if
         end do
      end do
      ! Volume = sum phi_i
      ! (sum phi_i r_i )/Volume
      r_cm_part(1:tcell,1) = r_cm_part(1:tcell,1)/volume(1:tcell)
      r_cm_part(1:tcell,2) = r_cm_part(1:tcell,2)/volume(1:tcell)
      r_cm_part(1:tcell,3) = r_cm_part(1:tcell,3)/volume(1:tcell)


    end subroutine cm_calc_local


    subroutine volume_calc(volume, f, tcell, np, np_part, lxyz, lxyz_inv, lxyz_inv_part)

      implicit none

      real, allocatable, intent(inout) :: volume(:)
      integer, intent(in) :: np, tcell, np_part
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), lxyz_inv_part(:,:,:)
      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      integer :: ip, icell, ncell, ip_part, ip_global_m

      volume(:) = 0.d0
      do ip=1, np

         do icell=1, tcell
            ip_global_m = r(icell)

            call vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

            if(f(ip_part,icell)%phi.gt.0.d0) then
               volume(icell) = volume(icell) + f(ip_part,icell)%phi
            end if
         end do
      end do
      ! Volume = sum phi_i


    end subroutine volume_calc

    subroutine cm_calc(r_cm, f, tcell, np, np_part, r, lxyz, lxyz_inv, lxyz_inv_part)

      implicit none

      integer, intent(in) :: np, tcell, np_part
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), lxyz_inv_part(:,:,:)
      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      real, allocatable, intent(inout) :: r_cm(:,:)
      integer, allocatable, intent(inout) :: r(:)

      integer :: ip, icell, ncell, ip_part, ip_global_m
      real :: volume(1:tcell)

      volume = 0.d0
      r_cm(tcell,1:3) = 0.d0
      do ip=1, np

         do icell=1, tcell
            ip_global_m = r(icell)

            call vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

            if(f(ip_part,icell)%phi.gt.0.d0) then
               volume(icell) = volume(icell) + f(ip_part,icell)%phi
               r_cm(icell,1:3) = r_cm(icell,1:3) + f(ip_part,icell)%phi*lxyz(ip,1:3)
            end if
         end do
      end do
      ! Volume = sum phi_i
      ! (sum phi_i r_i )/Volume
      r_cm(1:tcell,1) = r_cm(1:tcell,1)/volume(1:tcell)
      r_cm(1:tcell,2) = r_cm(1:tcell,2)/volume(1:tcell)
      r_cm(1:tcell,3) = r_cm(1:tcell,3)/volume(1:tcell)

    end subroutine cm_calc

    subroutine move(f, r_cm, r, np, np_part, tcell, lxyz, lxyz_part, lxyz_inv, lxyz_inv_part,Lsize)

      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      real, allocatable, intent(in) :: r_cm(:,:)
      integer, allocatable, intent(inout) :: r(:)
      integer, intent(in) :: np, np_part, tcell, Lsize(3)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_part(:,:), lxyz_inv(:,:,:), lxyz_inv_part(:,:,:)
      integer :: ip, icell, ncell, ip_part, ip_global_m
      real ::  delta_r(1:tcell,1:3), ftemp(1:np_part), dimg(1:tcell,1:3)

      do icell=1, tcell
        ! the geringonca is working, so ignore my previous comments Jan 12 2016
        ! I think the problem is the round error in the delta_r
        delta_r(icell,1:3) = r_cm(icell,1:3)-lxyz(r(icell),1:3)

        !dimg(icell,1) = min(abs(delta_r(icell,1)), 2*Lsize(1)-1-abs(delta_r(icell,1)))
        !dimg(icell,2) = min(abs(delta_r(icell,2)), 2*Lsize(2)-1-abs(delta_r(icell,2)))

        !if( sqrt(delta_r(icell,1)*delta_r(icell,1)+delta_r(icell,2)*delta_r(icell,2)) .ge. 1.d0 ) then

          do ip_part=1, np_part
              ! min image method is needed here! maybe.. i don't know
              ftemp(ip_part) = f(lxyz_inv_part( lxyz_part(ip_part,1) + int(anint(delta_r(icell,1))), &
                                                lxyz_part(ip_part,2) + int(anint(delta_r(icell,2))), &
                                                lxyz_part(ip_part,3) + int(anint(delta_r(icell,3)))), icell)%phi
          end do

          f(1:np_part,icell)%phi = ftemp(1:np_part)
          r(icell) = lxyz_inv(int(anint(r_cm(icell,1))),&
                              int(anint(r_cm(icell,2))),&
                              int(anint(r_cm(icell,3))) )
        !end if
      end do

    end subroutine move



end module run_cells_m
