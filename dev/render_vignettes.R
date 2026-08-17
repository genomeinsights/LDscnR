#!/usr/bin/env Rscript
# Render all package vignettes to ./rendered/ for sharing (HTML + best-effort PDF).
# `rendered/` is git-ignored (build artifacts); `dev/` is .Rbuildignore'd so this
# script never ships in the package. PDF needs LaTeX
# (install.packages("tinytex"); tinytex::install_tinytex()).
#
# Usage, from the package root:  Rscript dev/render_vignettes.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
dir.create("rendered", showWarnings = FALSE)
out <- normalizePath("rendered")   # absolute path (render() may change the working dir)

for (rmd in list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)) {
  base <- tools::file_path_sans_ext(basename(rmd))
  rmarkdown::render(rmd, output_file = file.path(out, paste0(base, ".html")), quiet = TRUE)
  message("HTML  -> rendered/", base, ".html")
  ok <- tryCatch({
    rmarkdown::render(rmd, output_format = rmarkdown::pdf_document(latex_engine = "xelatex"),
                      output_file = file.path(out, paste0(base, ".pdf")), quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  message(if (ok) paste0("PDF   -> rendered/", base, ".pdf")
          else "PDF   skipped (needs LaTeX: install.packages('tinytex'); tinytex::install_tinytex())")
}
message("\nShareable files in: ", out)
