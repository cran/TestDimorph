# TestDimorph 0.5.8

-   Updated citation info.

# TestDimorph 0.5.7

-   Minor corrections.

# TestDimorph 0.5.5

-   Removed function `accu_model`.
-   Various enhancements and bug fixes.

# TestDimorph 0.4.1

-   Fixed a bug in multivariate raw data generation.

# TestDimorph 0.4.0

-   A major update.
-   New functions:`van_vark`, `Hedges_g`,`MI_index` and `D_index`.
-   `univariate` and `multivariate` functions report different types of ANOVA and MANOVA.
-   effect sizes are reported with confidence intervals.
-   added 3 new datasets.
-   added a list `models` of the supported models for `accu_model` function.

# TestDimorph 0.3.6

-   Minor fixes-no user visible changes.

# TestDimorph 0.3.5

-   Maintainer email changed.

# TestDimorph 0.3.3

-   The package has been re-written and optimized with fewer dependencies.
-   New function names with snake case letters are introduced
-   Function names written in camel case letters are deprecated
-   `accu_model` function now can do cross-validation using different models supported by `caret` package.

# TestDimorph 0.3.1

-   `mda` is removed from dependencies and mixture and flexible discriminant analyses are no longer options in the `AccuModel` function.

# TestDimorph 0.3.0

-   Fixed some issues with corrplots in `Tg` function.

# TestDimorph 0.2.9

-   Raw data generation now can be either by uni or multivariate lognormal or truncated distribution.
-   Added the package to github.
-   Output of most functions is in the form of tidy tibbles.
-   Effect size for uni and multivariate analyses.
-   `pMatrix` function is deprecated and included in `Tg` function.
-   `AccuModel` function can generate roc curves.
-   added random forest to methods of `AccuModel`
-   `aovSS` function can do different types of post hoc tests.
-   Pairwise comparisons can be expressed by means of different alphabetical letters.
-   Informative error messages in case of wrong type of input.

# TestDimorph 0.2.1

-   Minor fixes-no user visible changes.

# TestDimorph 0.2.0

-   Added a `NEWS.md` file to track changes to the package.
-   Removed `biotools` from the dependencies to allow the package to work on Mac OS X.
-   Added 2 new functions `RawGen` and `AccuModel`.
-   Now `pMatrix` function can generate corrplots.
-   Added a `README.md`file.
-   updated the the package's description and help files.
