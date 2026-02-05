#!/usr/bin/env Rscript

lints <- lintr::lint_dir("scripts")
if (length(lints) > 0) {
  print(lints)
  quit(status = 1)
}

message("No lints found.")
