# Contributing to sempredict

Thank you for considering a contribution.

## Reporting a problem

Open an issue at <https://github.com/NicholasDanks/sempredict/issues>. For a
wrong number or a crash, please include a minimal reproducible example: the
model syntax or `seminr` specification, the data (or a way to simulate it), the
call, the seed, and the output you got against the output you expected, with
`sessionInfo()`.

## Proposing a change

1. Fork the repository and create a branch from `main`.
2. Every estimator, construction or diagnostic must come with a test against an
   external oracle where one exists (for example `lavaan::lavPredictY()`, a
   hand-computed t-test, or a published number from the companion tutorial),
   not only against the package's own previous output.
3. Run `devtools::document()`, `devtools::test()` and
   `rcmdcheck::rcmdcheck(args = "--as-cran")`; the check should be clean.
4. Open a pull request describing what the change makes possible or what it
   now gets right.

## Scope

`sempredict` does not estimate models. Estimation belongs in `lavaan`,
`seminr` or another fitting package; this package supplies the prediction
constructions, the identical-fold cross-validation harness and the assessment
layer. Support for a further fitting package is welcome as a new
`sem_params()` method.

## Seeking support

Questions about usage can be raised as issues. Questions about the
statistical methodology are best directed to the references in the package
documentation and the companion tutorial.

## Code of conduct

Be respectful. Contributors are expected to follow the
[Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
