
ldpred_coord <- function(
  ldref.file, 
  summary.stats.file, 
  out.file,
  id.col, 
  a1.col, a2.col, freq.col,
  chr.col, pos.col, 
  eff.col, pval.col, 
  n.col, beta, 
  exec
) {
  
  log <- system2(
    command = exec,
    args = c(
      "coord",
      "--gf", ldref.file,
      "--ssf", summary.stats.file,
      "--ssf-format=CUSTOM",
      "--out", out.file,
      "--rs", id.col,
      "--A1", a1.col, "--A2", a2.col, "--reffreq", freq.col,
      "--chr", chr.col, "--pos", pos.col,
      "--eff", eff.col, "--pval", pval.col,
      "--ncol", n.col,
      ifelse(beta, "--beta", "")
    ),
    stdout = TRUE, stderr = TRUE
  )
  
  return(
    list(
      coord_file = out.file
    )
  )
  
}