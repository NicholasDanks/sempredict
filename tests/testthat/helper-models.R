# Shared model specifications for the tests.
pd <- lavaan::PoliticalDemocracy
pd_syntax <- "
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
"
pd_y <- paste0("y", 5:8)
pd_x_da <- c(paste0("x", 1:3), paste0("y", 1:4))
pd_x_ea <- paste0("x", 1:3)
pd_sm <- seminr::relationships(seminr::paths(from = "ind60", to = "dem60"),
                               seminr::paths(from = c("ind60", "dem60"), to = "dem65"))
pd_mm_pls <- seminr::constructs(
  seminr::composite("ind60", paste0("x", 1:3), weights = seminr::mode_A),
  seminr::composite("dem60", paste0("y", 1:4), weights = seminr::mode_A),
  seminr::composite("dem65", paste0("y", 5:8), weights = seminr::mode_A))
pd_mm_plsc <- seminr::constructs(
  seminr::reflective("ind60", paste0("x", 1:3)),
  seminr::reflective("dem60", paste0("y", 1:4)),
  seminr::reflective("dem65", paste0("y", 5:8)))

mobi <- seminr::mobi
mobi_y <- paste0("CUSA", 1:3)
mobi_x <- paste0("CUEX", 1:3)
mobi_sm <- seminr::relationships(seminr::paths(from = "EXP", to = "SAT"))
mobi_mm_pls <- seminr::constructs(
  seminr::composite("EXP", mobi_x, weights = seminr::mode_A),
  seminr::composite("SAT", mobi_y, weights = seminr::mode_A))
mobi_mm_plsc <- seminr::constructs(
  seminr::reflective("EXP", mobi_x),
  seminr::reflective("SAT", mobi_y))

quiet_pls <- function(...) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(seminr::estimate_pls(...)))
  out
}
