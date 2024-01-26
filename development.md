# Information for developers

## Release management
To create a new release for smesh, perform the following steps:
1) Make sure that all PRs and changes that you want to go into the release are merged to
   `main` and that the latest commit on `main` has passed all CI tests.
2) Determine the currently released version of Smesh.jl, e.g., on the
   [release page](https://github.com/trixi-framework/smesh/releases). For this manual,
   we will assume that the latest release was `v0.2.3`.
3) Decide on the next version number. We follow [semantic versioning](https://semver.org/),
   thus each version is of the form `vX.Y.Z` where `X` is the major version, `Y` the minor
   version, and `Z` the patch version. In this manual, we assume that the major version is
   always `0`, thus the decision process on the new version is as follows:
   * If the new release contains *breaking changes* (i.e., user code might not work as
     before without modifications), increase the *minor* version by one and set the
     *patch* version to zero. In our example, the new version should thus be `v0.3.0`.
   * If the new release only contains minor modifications and/or bug fixes, the *minor*
     version is kept as-is and the *patch* version is increased by one. In our example, the
     new version should thus be `v0.2.4`.
4) Edit the string in the
   [`VERSION`](https://github.com/trixi-framework/Smesh.jl/blob/main/VERSION)
   and set it to the new version. Push/merge this change to `main`.
5) Go to GitHub and create a new release, either by going to the
   [release page](https://github.com/trixi-framework/smesh/releases) and clicking on "Draft
   a new release" or by directly following
   [this link](https://github.com/trixi-framework/smesh/releases/new).
   * Click on "Choose a tag", enter the new version (e.g., `v0.2.4`, *including* the `v`
     prefix), and select "Create a new tag".
   * Click on "Generate release notes".
   * Click on "Publish release".
7) The new release is immediately usable and visible on the release page \o/
8) To make sure people do not mistake the latest state of `main` as the latest release, we
   set the version in the `Project.toml` to a *development* version. The development version
   should be the latest released version, with the patch version incremented by one, and the
   `-dev` suffix added. For example, if you just released `v0.3.0`, the new development
   version should be `v0.3.1-dev`. If you just released `v0.2.4`, the new development
   version should be `v0.2.5-dev`.

Among other uses, smesh is the backend for the Julia package
[Smesh.jl](https://github.com/trixi-framework/Smesh.jl). For this purpose, it is
automatically built and distributed as a Julia JLL package. Therefore, after each release of
mesh, you should also update the corresponding build recipe in the
[Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil) repository:
1) If not yet done, create a fork of the
   [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil/) repository, which contains all
   build recipes for registered JLL packages.
2) Create a new branch in your Yggdrasil fork and modify the
   [`build_tarballs.jl`](https://github.com/JuliaPackaging/Yggdrasil/blob/master/S/smesh/build_tarballs.jl)
   file for smesh:
   * Change the existing version number in `build_tarballs.jl` to the version you just released (e.g.,
     [here](https://github.com/JuliaPackaging/Yggdrasil/blob/414237372f5bac40fc3cd8045727def18388a1d7/S/smesh/build_tarballs.jl#L6)
     the existing smesh version is `v0.1.0`)
   * Get the commit hash of the current smesh release, e.g., by going to the
     [latest release](https://github.com/trixi-framework/smesh/releases/latest) and then
     clicking on the short commit hash, then copying it from the browser URL.
   * Change the existing commit hash in `build_tarballs.jl` to the commit hash you just
     obtained (e.g.,
     [here](https://github.com/JuliaPackaging/Yggdrasil/blob/414237372f5bac40fc3cd8045727def18388a1d7/S/smesh/build_tarballs.jl#L11)
     the hash was `db69bf884d10c52577aee6795202136cfcf77178`)
3) Commit and push the changes to your fork.
4) Create a pull request from your fork's branch to Yggdrasil's `master` branch and name it
   appropriately (e.g., `[smesh] update to version v0.2.4`)
5) Wait for all tests in Yggdrasil to pass. If they don't, find the error and fix it.
6) Wait for the Yggdrasil PR to be merged to `master`.. After an appropriate waiting period
   (1-2 business days), you can ask nicely in the `#binarybuilder` on the Julia Slack if
   someone could merge your PR.
7) Once the Yggdrasil CI has built the new artifacts from its updated `master` branch, it
   will upload them to the
   [smesh\_jll.jl package repository](https://github.com/JuliaBinaryWrappers/smesh_jll.jl)
   and create a PR to Julia's general registry with the new version. Once all checks pass
   there and a grace period of 15 minutes has passed, the new release will be available.
