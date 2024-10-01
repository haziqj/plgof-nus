file_names <- c(
  "res_srs_type1.rda",
  "res_srs_power.rda",
  "res_complex_type1.rda",
  "res_complex_power.rda"
)
file_urls  <- c(
  "https://osf.io/ndu73/download",
  "https://osf.io/2d97y/download",
  "https://osf.io/638nw/download",
  "https://osf.io/gqde2/download"
)

for (i in seq_along(file_names)) {
  if (!file.exists(paste0("R/", file_names[i]))) {
    # If the file doesn't exist, download it
    download.file(file_urls[i], destfile = paste0("R/", file_names[i]), mode = "wb")
    message("File downloaded successfully.")
  } else {
    message("File already exists.")
  }
}
