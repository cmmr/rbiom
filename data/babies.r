delayedAssign(
  x     = 'babies', 
  value = local({
    f <- system.file(package = 'rbiom', 'extdata', 'babies.bz2')
    if (nzchar(f)) rbiom::as_rbiom(f) else NULL
  })
)
