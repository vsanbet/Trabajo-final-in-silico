required_packages <- c("dplyr", "ggplot2", "tidyverse")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
