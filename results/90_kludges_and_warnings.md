Kludges and Warnings
================
Keith Baggerly
2018-11-12

At present, the table of contents (toc) option is not working for
`github_document`. This is documented more fully as rmarkdown
[issue 1311](https://github.com/rstudio/rmarkdown/issues/1311), and
seems to involve the transition to PanDoc 2.0. Unfortunately, the
issue’s been open since March.

Since I like having tables of contents, I currently run things twice:
once, as an `md_document` to generate the table of contents, which I
then cut and paste to the Rmd file, and then again as a
`github_document` so the links will work when posted to GitHub.
