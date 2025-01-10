delayedAssign(
  x     = 'gems', 
  value = local({
    f <- system.file(package = 'rbiom', 'extdata', 'gems.bz2')
    if (nzchar(f)) rbiom::as_rbiom(f) else NULL
  })
)
