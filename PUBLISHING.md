# Publication notes

The prepared package was integrated into the existing repository after checking
its clean working tree, configured remote, and fetched branch state. Existing
tracked files and both original research sources were preserved. Generated
build output is excluded from version control.

Saptarshi explicitly authorized committing and pushing this package. The
publication commit uses the message:

    Package Bose-polaron GPE solver with CLI, restart support, and CI

Local GCC/Make integration checks passed without installing dependencies.
CMake/CTest and Clang are checked by the hosted GitHub Actions workflow; see
[VALIDATION.md](VALIDATION.md) for local coverage and scientific limitations.

A tagged release remains pending review of the interacting particle-number
drift and confirmation that the scientific behavior is ready for distribution.
