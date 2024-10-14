# this top-level script may be used to run (knit) the individual
# Rmd files

# set up output path for the reports
config <- yaml::read_yaml("config.yaml")
report_dir <- file.path(config$out_root, 'html')
dir.create(report_dir)

rmarkdown::render(input = file.path(config$project_root, 'Rmd', '01_load_and_filter.Rmd'),
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

rmarkdown::render(input = file.path(config$project_root, 'Rmd', '02_integration.Rmd'),
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

rmarkdown::render(input = file.path(config$project_root, 'Rmd', '03_azimuth_annotation.Rmd'),
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

rmarkdown::render(input = file.path(config$project_root, 'Rmd', '04_clonotypes.Rmd'),
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

rmarkdown::render(input = file.path(config$project_root, 'Rmd', '05_diff_expr.Rmd'),
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

rmarkdown::render(input = file.path(config$project_root, 'Rmd', '06_figures.Rmd'),
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())
