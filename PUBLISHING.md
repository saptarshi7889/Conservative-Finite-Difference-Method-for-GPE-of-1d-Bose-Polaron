# Publishing v0.2.0

This release distributes source, not platform-specific binaries. Read
RELEASE_NOTES.md and VALIDATION.md for requirements and numerical limitations.
This release does not claim exact interacting norm conservation.

## Local preparation

With an existing Eigen installation:

```sh
make EIGEN_INCLUDE=/path/to/directory/containing/Eigen test
python3 scripts/package_release.py
cd dist
sha256sum -c SHA256SUMS
```

The packaging script creates `bose-polaron-gpe-0.2.0.tar.gz`,
`bose-polaron-gpe-0.2.0.zip`, and `SHA256SUMS` under dist/. Archives include
build inputs, tests, documentation, citation/license, both research programs
under `legacy/`, the comparison runner, and the README animation. Local
simulation data, build output, run logs, and other machine-specific example
scripts are excluded. GitHub's automatic source
archives contain the full tagged repository; the named package assets are smaller.
No dependencies are downloaded by the local packaging script.

## Commit, push, and release

The maintainer retains control of commits. Review the prepared changes, then
run these commands from the repository root yourself:

```sh
git add -A -- README.md CMakeLists.txt CITATION.cff PUBLISHING.md RELEASE_NOTES.md \
  VALIDATION.md legacy examples/L500_v05 scripts/package_release.py \
  finite_difference_modified_impurity_Thomas_generalized_without_Matrix_hermitian_vel_contd_new_cluster_upwind_scheme.cpp \
  finite_difference_modified_impurity_Thomas_generalized_without_Matrix_hermitian_vel_new_cluster_upwind_scheme_corrected.cpp
git commit -m "release: prepare v0.2.0 source distribution"
git push origin master
git tag -a v0.2.0 -m "Bose-polaron GPE solver v0.2.0"
git push origin v0.2.0
```

Pushing the tag triggers `.github/workflows/release.yml`. GitHub-hosted jobs
install Eigen, validate the tag/version match, build and test with GCC and
Clang via CMake/CTest, package the source, and test the extracted tar archive.
Only after both validation jobs pass does the publishing job create the GitHub
Release and upload both package archives and SHA256SUMS. The hosted runner's
built-in GitHub CLI uses GITHUB_TOKEN; no personal token is stored in the repo.
GitHub Actions must be enabled and permitted to create releases.

Do not publish a release manually for the same tag while this workflow is running.
Monitor the Source release workflow under Actions, then open the Releases page:
https://github.com/saptarshi7889/Conservative-Finite-Difference-Method-for-GPE-of-1d-Bose-Polaron/releases

No commit, tag, push, or GitHub Release was created by local preparation alone.
