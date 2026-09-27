Sys.setenv(PATH = paste("C:/rtools45/usr/bin", "C:/rtools45/ucrt64/bin", "C:/rtools45/mingw64/bin", Sys.getenv("PATH"), sep = ";"))
remotes::install_local("C:/Users/19473/AppData/Local/Temp/presto-master.tar.gz", dependencies = FALSE, upgrade = "never")
cat("PRESTO DONE:", requireNamespace("presto", quietly = TRUE), "\n")
