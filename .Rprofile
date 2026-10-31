if (identical(Sys.info()[["sysname"]], "Darwin") &&
    identical(R.version$arch, "aarch64")) {
  # Avoid source builds that link data.table against an incompatible libomp.
  options(
    pkgType = "binary",
    # The shared renv sandbox lock can stall indefinitely on this setup.
    renv.config.sandbox.enabled = FALSE
  )
}

source("renv/activate.R")
