module misc_m

  use global_m

  implicit none

  private

  public ::  cross_section, str2int, density_calc, hfield_calc, hfunc, deltak, heaviside, ran2, output_aux,&
   spherical_surface, gen_cell_points, vec_local2global, vec_global2local, output

  contains


  subroutine cross_section(s, cs, lxyz_inv, Lsize)
    implicit none
    real, allocatable, intent(in) :: s(:)
    real, intent(out) :: cs
    integer, allocatable, intent(in) :: lxyz_inv(:,:,:)
    integer, intent(in) :: Lsize(3)
    integer :: cluster(1:size(s)), dev
    integer, allocatable :: connections(:,:),cs_corrected(:), ct(:)
    integer :: i, j, k, l, m, n, x,y,z, n_corrected, att_new, att_old, ip2, cs_count(1000)
    logical :: out, out2
    n = 1
    cluster(:) = 0
    cs_count(:) = 0
    i = 0
    out = .false.
    out2 = .false.
    att_new = 1

    do i=-Lsize(1), Lsize(1)-1
      att_new = 1
    do while(.true.)
      att_old = 1
      out = .false.
      out2 = .false.
      do j=-Lsize(2),Lsize(2)-1
        do k=-Lsize(3), Lsize(3)-1
          ip = lxyz_inv(i,j,k)

          if(cs_count(n).eq.0) then
            if(s(ip).lt.1.and.cluster(ip).eq.0) then
              cs_count(n) = cs_count(n) + 1
              cluster(ip) = n
              att_new = 0
            end if
          else
            do y=-1,1
              do z = -1,1
                ip2 = lxyz_inv(i,lxyz(ip,2)+y,lxyz(ip,3)+z)
                if(cluster(ip2).eq.n .and. abs(y)+abs(z).gt.0 .and. s(ip).lt.1.and.cluster(ip).eq.0 ) then
                    cluster(ip) = n
                    att_old = 0
                    out2 = .true.
                end if
                if(out2) EXIT
              end do
              if(out2) EXIT
            end do

          end if
        end do
      end do
      if (att_new.ne.0.and.cs_count(n).eq.0) EXIT
      if (att_old.ne.0) then
        n = n + 1
        att_new = 1
      end if
    end do
    end do

    !connections(:,:) = 0.d0

    OPEN(UNIT=9988,file='cs')
    !cs_count(:) = 0.d0
    !ct(:) = 0.d0
    do ip=1, np
      if(cluster(ip).gt.0) then
         !cs_count(cluster(ip)) = cs_count(cluster(ip))  + 1
         !ct(cluster(ip)) = 1.d0
         write(9988,*) lxyz(ip,1:3), cluster(ip)
      end if
    end do
    CLOSE(9988)
    !n = sum(ct(:))
    !print*,"number of domains:", n

    cs =  1!sum(cs_count(1:n))/n
    !dev = sqrt(sum( (cs_count(1:n)-cs)**2/n) )
    !print*, dev
    !print*, n
  end subroutine cross_section

  subroutine str2int(str,int)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    read(str,*)  int
  end subroutine str2int

    subroutine density_calc(subs, Lsize, lxyz_inv, dens, deviation)

      implicit none

      real, allocatable, intent(in) :: subs(:)
      real, intent(out) :: dens, deviation
      integer, allocatable, intent(in) :: lxyz_inv(:,:,:)
      integer, intent(in) :: Lsize(3)
      integer :: i, j, k, ii, jj, kk, n, a
      real :: density(1000), mean, dev, mean2
      a = 20
      n = 0
      mean = 0.d0
      mean2 = 0.d0
      density(:) = 0.d0
      do i = -Lsize(1), Lsize(1)-a, a
        do j = -Lsize(2), Lsize(2)-a, a
          do k = -Lsize(3), Lsize(3)-a, a
            n = n + 1

            do ii = i, i+a
              do jj = j, j+a
                do kk = k, k+a

                  if(subs(lxyz_inv(ii,jj,kk))>0) density(n) = density(n) + subs(lxyz_inv(ii,jj,kk))/(real(a)**3)

                end do
              end do
            end do
          end do
        end do
      end do

      mean = sum(density(1:n))/real(n)
      dev = sqrt(sum((density(1:n)-mean)**2)/real(n))
      dens = mean
      deviation = dev

    end subroutine density_calc


    subroutine hfield_calc(cell, aux, r, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np_part )

      type(mesh_t), allocatable, intent(in) :: cell(:,:)
      integer, allocatable, intent(in) :: r(:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), lxyz_part(:,:), lxyz_inv_part(:,:,:), ncell(:)
      type(mesh_t), allocatable, intent(inout) :: aux(:,:)
      integer, intent(in) ::  ntype, np_part, tcell
      integer :: icell, itype, ip_part, ip

      aux(:,:)%phi = 0.d0
      do itype = 1, ntype
        do icell=1+ (itype-1)*(ncell(itype-1)), ncell(itype) + (itype-1)*(ncell(itype-1))
          do ip_part = 1, np_part
            call vec_local2global(ip, r(icell), ip_part, lxyz, lxyz_inv, lxyz_part)
            aux(ip,itype)%phi = aux(ip,itype)%phi + hfunc(cell(ip_part,icell)%phi)*deltak(cell(ip_part,icell)%itype,itype)
          end do
        end do
      end do

    end subroutine hfield_calc

    function hfunc(x)
      real :: hfunc, x,y
	    !y = 0.5*(x+1.0)
      hfunc = x**2*(3.0-2.0*x)
    end function hfunc

    function deltak(l,m)

      implicit none

      real :: deltak
      integer :: l, m

      if(l.eq.m) then
         deltak = 1.d0
      else
         deltak = 0.d0
      end if

    end function deltak

    function heaviside(x)

      implicit none

      real :: x, heaviside

      if ( x<0.d0) then
         heaviside = 0.d0
      else
         heaviside = 1.d0
      end if
    end function heaviside

    function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2,AM,EPS,RNMX
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
           NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer :: idum2,j,k,iv(NTAB),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1

            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
11          continue
            iy=iv(1)
         endif
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
         k=idum2/IQ2
         idum2=IA2*(idum2-k*IQ2)-k*IR2
         if (idum2.lt.0) idum2=idum2+IM2
         j=1+iy/NDIV
         iy=iv(j)-idum2
         iv(j)=idum
         if(iy.lt.1)iy=iy+IMM1
         ran2=min(AM*iy,RNMX)
         return
       end function ran2


    subroutine output_aux(aux, id, char_length, file_name, dir_name, itype, ntypes, np, lxyz, dr)

      implicit none
      type(mesh_t), intent(in) :: aux(:,:)
      integer, intent(in) :: id, char_length, itype, ntypes, np, lxyz(:,:)
      real, intent(in) :: dr(3)
      character(3), intent(in) :: dir_name
      character(len=char_length), intent(in) :: file_name
      integer :: itypes, ip


      OPEN (UNIT=id,FILE=trim(dir_name//'/sc'//file_name//'.xyz'))
      do ip=1, np
         do itypes = itype, ntypes
            !if(aux(ip,itypes)%phi>=0.d0) then
               write(id,'(F10.2,F10.2,F10.2,F10.2,I10)') lxyz(ip,1:3)*dr(1:3), aux(ip,itypes)%phi, itypes
            !end if
         end do
      end do
      close(id)
    end subroutine output_aux

    subroutine output(field, id, char_length, file_name, dir_name, itype, ntypes, np, lxyz,dr)

      implicit none
      real, allocatable, intent(in) :: field(:)
      real, intent(in) :: dr(3)
      integer, intent(in) :: id, char_length, itype, ntypes, np, lxyz(:,:)
      character(3), intent(in) :: dir_name
      character(len=char_length), intent(in) :: file_name
      integer :: itypes, ip


      OPEN (UNIT=id,FILE=trim(dir_name//'/'//file_name//'.xyz'))
      do ip=1, np
              !if (lxyz(ip,2).eq.0) then
                write(id,'(F10.2,F10.2,F10.2,F10.2)') lxyz(ip,1:3)*dr(1:3), field(ip)
            !  end if
               !write(id,'(I10,I10,F10.2,I10)') lxyz(ip,1:2), field(ip)
      end do
      close(id)
    end subroutine output

    subroutine spherical_surface(Rc, R, ndim, delta)

      implicit none

      real, intent(in) :: Rc, delta
      integer, intent(inout) :: ndim
      integer, allocatable, intent(inout) :: R(:,:)
      real :: theta, phi, M_Pi
      character(len=10) :: dr, ds, sdx, sdy, sdz
      integer :: i, dx, dy, dz
      M_Pi = 3.14159265359
      ! Note:
      ! No intrinsic exists to convert between a numeric value and a formatted character
      ! string representation  for instance, given the CHARACTER value '154',
      ! obtaining an INTEGER or REAL value with the value 154, or vice versa.
      ! Instead, this functionality is provided by internal-file I/O,
      ! as in the following example:

      ! Convert a string to a numeric value
      ! read (string,'(I10)') value
      ! print *, value

      ! Convert a value to a formatted string
      ! write (string2,'(I10)') value
      ! print *, string2


      i = 0
      theta = 0.0
      phi = 0.0

      do while (theta<=M_Pi)
         phi = 0.0
         do while (phi<=2.d0*M_Pi)
            dx =  int(anint(Rc * sin(theta) * cos(phi)) )
            dy =  int(anint(Rc * sin(theta) * sin(phi) ))
            dz =  int(anint(Rc * cos(theta)) )
            write(sdx,'(I2)') abs(dx)
            write(sdy,'(I2)') abs(dy)
            write(sdz,'(I2)') abs(dz)

            dr = trim(sdx)//trim(sdy)//trim(sdz)

            if (dr .ne. ds) then
               i = i + 1

               R(i,1) = dx
               R(i,2) = dy
               R(i,3) = dz

               ds = dr
            end if
            phi = phi + delta
         end do
         theta = theta + delta
      end do

      i = i +1
      R(i,1) = 0
      R(i,2) = 0
      R(i,3) = -Rc
      ndim = i

    end subroutine spherical_surface


    subroutine gen_cell_points(Rc,R,ndim)

      implicit none

      real, intent(in) :: Rc
      integer, intent(inout) :: ndim
      integer, intent(inout) :: R(:,:)
      real :: dx, i, j, k, M_Pi
      integer :: counter
      M_Pi = 3.14159265359
      counter = 0
      i = -Rc
      do while(i<=Rc)
         j = -Rc
         do while (j<=Rc)
           k = -Rc
           do while(k<=Rc)
               dx = sqrt(i**2 + j**2 + k**2)
               if(dx <= Rc) then
                  counter = counter + 1
                  R(counter,1) = int(anint(i))
                  R(counter,2) = int(anint(j))
                  R(counter,3) = int(anint(k))


               end if
               k = k + 1.d0
             end do
            j = j +1.d0
         end do
         i = i + 1.d0
      end do

      ndim = counter

    end subroutine gen_cell_points


    subroutine vec_local2global(ip, ip_global_m, ip_local, lxyz, lxyz_inv, lxyz_part)
      !  give the ip_local that you want to map in the global space and ip_global_m
      ! which represents the position of the small box
      ! the routine will return the ip (global) associated with the ip_local...
      ! hopefully... :D
      implicit none

      integer, allocatable, intent(in) :: lxyz(:,:),  lxyz_inv(:,:,:), lxyz_part(:,:)
      integer, intent(in) :: ip_global_m, ip_local
      integer, intent(out) :: ip

      ip =  lxyz_inv(lxyz_part(ip_local,1) + lxyz(ip_global_m,1), &
                     lxyz_part(ip_local,2) + lxyz(ip_global_m,2), &
                     lxyz_part(ip_local,3) + lxyz(ip_global_m,3) )

    end subroutine vec_local2global


    subroutine vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

      implicit none

      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), lxyz_inv_part(:,:,:)
      integer, intent(in) :: ip_global_m, ip
      integer, intent(out) :: ip_part

      ip_part =  lxyz_inv_part( lxyz(ip,1) - lxyz(ip_global_m,1), &
                                lxyz(ip,2) - lxyz(ip_global_m,2), &
                                lxyz(ip,3) - lxyz(ip_global_m,3))


    end subroutine  vec_global2local


end module misc_m
