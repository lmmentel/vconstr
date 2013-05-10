module ioModule
 use varModule 

 private

 interface matprint 
    module procedure print_eigenproblem
    module procedure print_eigenproblem_with_header
    module procedure print_vector
    module procedure print_vector_with_header
    module procedure print_matrix
    module procedure print_matrix_with_header
    module procedure print_triangular_matrix_as_square
 end interface 

 interface write_gamess_vecs
    module procedure print_gamess_vec
 end interface

 interface header
     module procedure print_header
 end interface
 interface fmt_to_len
     module procedure adjust_string_to_length
 end interface

 public :: matprint, write_gamess_vecs, header, fmt_to_len

 interface EnergySummary
     module procedure print_final_energies
 end interface

 public :: EnergySummary

 contains 
 
 function adjust_string_to_length(string, length) result(frmtstring)
   character(len=*), intent(in)       :: string
   integer(ik),      intent(in)       :: length
   character(len=length) :: frmtstring

   if (len(string) < length) then 
     frmtstring = string//repeat(' ',(length-len(string)))
   else
     frmtstring = string(1:length)
   endif
   return     
 end function adjust_string_to_length

 subroutine print_header(header)
  character(len=*), intent(in) :: header

  write(*,'(/5x,a)') repeat('=',len(header))
  write(*,'(5x,a)') header
  write(*,'(5x,a/)') repeat('=',len(header))
 end subroutine print_header

 subroutine print_eigenproblem(M, V)
  integer(ik), parameter :: vec_per_page = 5 
  integer(ik)          :: i, j, iel, k, l
  real(dp), intent(in) :: M(:,:)
  real(dp), intent(in) :: V(:)

 if (mod(size(M,2),vec_per_page) == 0 ) then 
   iel = size(M,2)/vec_per_page
 else
   iel = size(M,2)/vec_per_page+1
 endif

 do i = 1, iel 
   if(i < iel) then
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5)
     write(*,'(/7x,5(f10.5)/)') (V(j+(i-1)*5),j=1,5)
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5)
       enddo
       write(*,*)
   else   
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5-(iel*5-size(M,2)))
     write(*,'(/7x,5(f10.5)/)') (V(j+(i-1)*5),j=1,5-(iel*5-size(M,2)))
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5-(iel*5-size(M,2)))
       enddo
   endif
 enddo    
 write(*,*)
 end subroutine print_eigenproblem

 subroutine print_eigenproblem_with_header(M, V, header)
  real(dp),     intent(in) :: M(:,:)
  real(dp),     intent(in) :: V(:)
  character(*), intent(in) :: header
  integer(ik),  parameter  :: vec_per_page = 5 
  integer(ik)              :: i, j, iel, k, l, lspace
  character(len=30)        :: headerfmt

  lspace = 7 + int((50-len(header))/2)
  write(headerfmt,'("(",i2,"x,a)")') lspace  
  write(*,*)
  write(*,headerfmt) repeat('=',len(trim(header))) 
  write(*,headerfmt) trim(header) 
  write(*,headerfmt) repeat('=',len(trim(header))) 

 if (mod(size(M,2),vec_per_page) == 0 ) then 
   iel = size(M,2)/vec_per_page
 else
   iel = size(M,2)/vec_per_page+1
 endif

 do i = 1, iel 
   if(i < iel) then
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5)
     write(*,'(/7x,5(f10.5)/)') (V(j+(i-1)*5),j=1,5)
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5)
       enddo
       write(*,*)
   else   
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5-(iel*5-size(M,2)))
     write(*,'(/7x,5(f10.5)/)') (V(j+(i-1)*5),j=1,5-(iel*5-size(M,2)))
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5-(iel*5-size(M,2)))
       enddo
   endif
 enddo    
 write(*,*)
 end subroutine print_eigenproblem_with_header


 subroutine print_vector(V)
  real(dp), intent(in) :: V(:)
  integer(ik)          :: i

  do i = 1, size(V)
    write(*,'(1x,i3,f10.5)') i, V(i)
  enddo
 end subroutine print_vector

 subroutine print_vector_with_header(V, header)
  real(dp),     intent(in) :: V(:)
  character(*), intent(in) :: header
  integer(ik)              :: i

  write(*,'(/4x,a)') repeat('=',len(trim(header))) 
  write(*,'(4x,a)') trim(header) 
  write(*,'(4x,a/)') repeat('=',len(trim(header))) 

  do i = 1, size(V)
    write(*,'(1x,i3,f10.5)') i, V(i)
  enddo
 end subroutine print_vector_with_header

 subroutine print_matrix(M)
  real(dp), intent(in) :: M(:,:)
  integer(ik)          :: i,j,iel,k,l

 if (mod(size(M,2),5) == 0 ) then 
   iel = size(M,2)/5
 else
   iel = size(M,2)/5+1
 endif

 do i = 1, iel 
   if(i < iel) then
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5)
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5)
       enddo
       write(*,*)
   else   
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5-(iel*5-size(M,2)))
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5-(iel*5-size(M,2)))
       enddo
   endif
 enddo    
 write(*,*)
 end subroutine print_matrix

 subroutine print_matrix_with_header(M, header)
  real(dp),     intent(in) :: M(:,:)
  character(*), intent(in) :: header
  integer(ik)              :: i,j,iel,k,l, lspace
  character(len=30)        :: headerfmt

  lspace = 7 + int((50-len(header))/2)
  write(headerfmt,'("(",i2,"x,a)")') lspace  
  write(*,*)
  write(*,headerfmt) repeat('=',len(trim(header))) 
  write(*,headerfmt) trim(header) 
  write(*,headerfmt) repeat('=',len(trim(header))) 

 if (mod(size(M,2),5) == 0 ) then 
   iel = size(M,2)/5
 else
   iel = size(M,2)/5+1
 endif

 do i = 1, iel 
   if(i < iel) then
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5)
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5)
       enddo
       write(*,*)
   else   
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5-(iel*5-size(M,2)))
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5-(iel*5-size(M,2)))
       enddo
   endif
 enddo    
 write(*,*)
 end subroutine print_matrix_with_header
 
 subroutine print_gamess_vec(coeffs)
  use commons
  real(dp),           intent(in) :: coeffs(:,:)
  integer(ik) :: i, j, l, ilab, llab, nlines
  integer(ik), parameter :: gunit=16

  if (mod(size(coeffs,1),5) == 0) then
    nlines = int(size(coeffs,1)/5)
  else
    nlines = int(size(coeffs,1)/5,ik)+1  
  endif
 
!  print *, 'nlines = ', nlines

  open(gunit, file=trim(GmsNOsFile), status="replace")
 ! write header and then orbitals
  write(gunit,'(/a,i3,a/)') ' $guess guess=moread norb=',size(coeffs,2),' punmo=.true. prtmo=.true. $end' 
  write(gunit,'(a)') ' $vec' 
  do i = 1, size(coeffs,2)
    if (i >= 100) then 
      ilab = mod(i, 100)
    else 
      ilab = i
    endif
    do l = 1, nlines
      if (l > 1000) then 
        llab = mod(l, 1000)
      else
        llab = l
      endif
      if (l < nlines) then
        write(gunit,gmsvfmt) ilab,llab,(coeffs(j+(l-1)*5,i),j=1,5)
      else
        write(gunit,gmsvfmt) ilab,llab,(coeffs(j+(l-1)*5,i),j=1,5-(l*5-size(coeffs,1)))
      endif
    enddo
  enddo
  write(gunit,'(a)') ' $end' 
  close(gunit) 

 end subroutine print_gamess_vec
 
 subroutine print_triangular_matrix_gamess(M, nb)
  real(dp), intent(in) :: M(:)
  integer,  intent(in) :: nb

  interface
      subroutine prtri(D, N)
        use varModule
        real(DPgamess),    intent(in) :: D(:)
        integer(IPgamess), intent(in) :: N
      end subroutine
  end interface

  call prtri(M, int(nb, kind=IPgamess))
 end subroutine print_triangular_matrix_gamess

 subroutine print_triangular_matrix_as_square(M, nb, header)
  real(dp),         intent(in) :: M(:)
  integer(ik),      intent(in) :: nb
  character(len=*), intent(in) :: header 
  real(dp), allocatable :: Mtemp(:,:)
  integer(ik) :: i, j, ij

  allocate(Mtemp(nb, nb))
  ij = 0
  do i = 1, nb 
    do j = 1, i-1
      ij = ij + 1
      Mtemp(i,j) = M(ij)
      Mtemp(j,i) = M(ij)
    enddo
    ij = ij + 1
    Mtemp(i,i) = M(ij)
  enddo
  call print_matrix_with_header(Mtemp, header) 
  deallocate(Mtemp)
 end subroutine print_triangular_matrix_as_square

 subroutine print_final_energies(dmft_Etotel, dmft_Ekin, dmft_Epot, dmft_Eee)
  use commons
  real(DP), intent(in) :: dmft_Etotel, dmft_Ekin, dmft_Epot, dmft_Eee 

  call print_header("FINAL ENERGIES")

  write(*,'(35x,3a16)') 'DMFT   ', 'CI    ', 'HF-SCF  '
  write(*,'(5x,a30,f16.8)')   "Nuclear Repulsion Energy      ", &
 & nuclear_repulsion_energy
  write(*,'(5x,a30,2f16.8)')  "Kinetic Energy                ", &
 & dmft_Ekin, kinetic_energy
  write(*,'(5x,a30,2f16.8)')  "Electron-Nuclei Energy        ", &
 & dmft_Epot, el_nuc_energy
  write(*,'(5x,a30,2f16.8)')  "One Electron Energy           ", &
 & dmft_Ekin+dmft_Epot, kinetic_energy+el_nuc_energy
  write(*,'(5x,a30,2f16.8)')  "Total Electronic Energy       ", &
 & dmft_Etotel, electronic_energy
  write(*,'(5x,a30,2f16.8)')  "Electron Interaction Energy   ", &
 & dmft_Eee, el_el_energy
  write(*,'(/5x,a30,3f16.8)') "Total Energy                  ", &
 & dmft_Etotel+nuclear_repulsion_energy, total_energy, scf_energy 
  write(*,'(/a)') '.....end of FINAL ENERGIES'
 end subroutine print_final_energies

end module ioModule
