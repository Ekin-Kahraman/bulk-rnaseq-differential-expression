# Security policy

## Supported versions

Security and correctness fixes are made on the default branch and latest tagged release.

## Reporting a vulnerability

This repository processes public GEO data and local tabular files. If you find a vulnerability such as unsafe file handling, credential leakage, or a supply-chain issue in dependency restoration, do not open a public issue.

Email: evk23umu@uea.ac.uk
Subject: `bulk-rnaseq security: <short description>`

Please include:
- Reproduction steps
- R version and operating system
- Whether `renv` restore was used
- Expected impact

I aim to respond within 7 days.

## Data handling

Do not attach private patient metadata, raw clinical files, or access credentials to public issues. Use public accession IDs or redacted synthetic examples.
