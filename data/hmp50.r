delayedAssign(
  x     = 'hmp50', 
  value = local({
    f <- system.file(package = 'rbiom', 'extdata', 'hmp50.bz2')
    if (nzchar(f)) rbiom::as_rbiom(f) else NULL
  })
)
