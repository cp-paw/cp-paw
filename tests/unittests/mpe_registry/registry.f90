program mpe_registry
  use mpe_module
  implicit none
  integer :: nt, me, n, rank_, value, expected, i
  integer, allocatable :: color(:)

  call MPE$INIT
  call MPE$QUERY('~', nt, me)
  allocate(color(nt))
  color = 1
  call MPE$NEW('~', 'MONOMER', nt, color)
  call MPE$NEW('MONOMER', 'K', nt, color)
  color = [(mod(i-1, 2)+1, i=1,nt)]
  call MPE$NEW('MONOMER', 'SKALA_POINTS', nt, color)
  call MPE$QUERY('K', n, rank_)
  if (n /= nt .or. rank_ /= me) error stop 'pre-existing communicator changed'
  value = me
  call MPE$COMBINE('K', '+', value)
  if (value /= nt*(nt+1)/2) error stop 'pre-existing reduction changed'
  call MPE$QUERY('SKALA_POINTS', n, rank_)
  expected = count(color == color(me))
  if (n /= expected .or. rank_ /= (me+1)/2) error stop 'new communicator layout changed'
  value = 1
  call MPE$COMBINE('SKALA_POINTS', '+', value)
  if (value /= expected) error stop 'new reduction changed'
  call MPE$QUERY('MONOMER', n, rank_)
  if (n /= nt .or. rank_ /= me) error stop 'parent communicator changed'
  if (me == 1) print *, 'PASS: communicator insertion preserves existing groups'
  call MPE$EXIT
end program mpe_registry
