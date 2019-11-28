program main
  !Numerical simulation of the Olami-Feder-Christensen model
  !developed by SO OZAWA (Master Student, Earthquake Research Institute, UTokyo)
  implicit none
  real(8),allocatable::s(:,:),c(:,:)
  real(8)::time,r,rate,dt,d
  integer::i,imax,j,jmax,k,kmax,a(2),seedsize
  integer,allocatable::mag(:),seed(:)
  !read(*,*) imax
  imax=10
  jmax=imax
  kmax=10000
  rate=1.d0 !stressing rate
  d=0.9d0 !dissipation rate

  allocate(s(0:imax+1,0:jmax+1),c(0:imax+1,0:jmax+1),mag(kmax))

  !get seed of random_number
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do i = 1, seedsize
      call system_clock(count=seed(i))
  end do
  call random_seed(put=seed(:))

  open(12,file='event.dat')
  open(13,file='stress.dat')

  !random initial condition
  s=0.d0
  do i=1,imax
    do j=1,jmax
      call random_number(r)
      s(i,j)=r
    end do
  end do

  !numerical simulation
  time=0.d0


  mag=1
  do k=1,kmax
    do i=1,imax
      do j=1,jmax
        write(13,'(2i5,f7.3)') i,j,s(i,j)
      end do
      write(13,*)
    end do

    dt=(1.d0-maxval(s(1:imax,1:jmax)))/rate
    s(1:imax,1:jmax)=s(1:imax,1:jmax)+dt*rate
    a=maxloc(s(1:imax,1:jmax))
    i=a(1)
    j=a(2)
    write(*,*) i,j,s(i,j)
    s(i,j)=0.d0
    s(i+1,j)=s(i+1,j)+0.25d0*d
    s(i-1,j)=s(i-1,j)+0.25d0*d
    s(i,j+1)=s(i,j+1)+0.25d0*d
    s(i,j-1)=s(i,j-1)+0.25d0*d

    do while(.true.)
      if(maxval(s(1:imax,1:jmax)).lt.1.d0) exit
      c=0.d0

      do i=1,imax
        do j=1,jmax
          if(s(i,j).ge.1.d0) then
            write(*,*) i,j,s(i,j)
            mag(k)=mag(k)+1
            !call topple(i,j,s,c)
            c(i,j)=-s(i,j)
            c(i+1,j)=c(i+1,j)+0.25d0*d*s(i,j)
            c(i-1,j)=c(i-1,j)+0.25d0*d*s(i,j)
            c(i,j+1)=c(i,j+1)+0.25d0*d*s(i,j)
            c(i,j-1)=c(i,j-1)+0.25d0*d*s(i,j)
          end if
        end do
      end do
      s=s+c


    end do
    time=time+dt
    !write(*,*) 'k=',k,mag(k)
    write(12,*)k,time,mag(k),sum(s(1:imax,1:jmax))/dble(imax*jmax)

    write(*,*)
    write(*,*)
    write(13,*)

  end do
   close(12)
   stop
end program main
