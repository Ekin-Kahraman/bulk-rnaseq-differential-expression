# Contributing

Contributions are welcome when they improve reproducibility, statistical correctness, documentation, or output validation.

## Development setup

```bash
Rscript 000_install_dependencies.R
Rscript scripts/12_output_manifest.R
Rscript -e 'testthat::test_dir("tests/testthat")'
Rscript dev/lint.R
```

The full rebuild path is checked in CI. Run the narrow tests locally before opening a pull request.

## Pull request checklist

- Tests pass with `testthat`.
- `Rscript dev/lint.R` reports no lints.
- Changed analysis outputs are regenerated intentionally.
- `results/tables/output_manifest.csv` is updated when tracked outputs change.
- Any changed statistical threshold, design formula, contrast, or filtering rule is documented in the README.

## Scientific correctness

Do not change reported results silently. If a change affects differential-expression counts, enrichment results, or figures, include the before and after numbers in the pull request.

## Licence

By contributing, you agree that your contributions are licensed under the MIT licence.
