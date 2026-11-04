run_analysis_example <- function(script,args) {
  file <- system.file("examples",script,package="geoCensor")
  if(!nzchar(file))file <- testthat::test_path("..","..","inst","examples",script)
  log <- tempfile(fileext=".log")
  previous <- Sys.getenv("R_LIBS",unset=NA_character_)
  on.exit(if(is.na(previous))Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS=previous),add=TRUE)
  Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
  executable <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
  status <- suppressWarnings(system2(executable,
    args=shQuote(c(file,args)),stdout=log,stderr=log))
  list(status=status,log=readLines(log,warn=FALSE))
}
