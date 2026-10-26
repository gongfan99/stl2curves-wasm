# Third-party notices

The reconstruction implementation and native CAD fixtures are derived from and
informed by **stl2curves**, Copyright (c) 2026 Duvdivan, MIT license.
The full upstream MIT notice is preserved in `LICENSE`.
Reference revision: `9f71aaa16eb689e759435d2e1fbf2ba7468cadfb`.
Source: https://github.com/Duvdivan/stl2curves

Dependencies are resolved through the vcpkg registry at
`817892b0ecbc6af994dd6df999c73faf795adbf2`:

| Dependency | License | Source |
| --- | --- | --- |
| Open CASCADE Technology 8.0.1 | LGPL-2.1 with OCCT exception | https://github.com/Open-Cascade-SAS/OCCT/tree/V8_0_1 |
| Eigen | MPL-2.0 | https://gitlab.com/libeigen/eigen |
| nanoflann | BSD-2-Clause | https://github.com/jlblancoc/nanoflann |
| libzip | BSD-3-Clause | https://libzip.org/ |
| pugixml | MIT | https://github.com/zeux/pugixml |
| zlib (libzip dependency) | Zlib | https://zlib.net/ |
| vcpkg port definitions and patches | MIT | https://github.com/microsoft/vcpkg |

The OpenCascade overlay retains vcpkg's pinned source hash and upstream patches;
local changes select headless toolkits and retain exception handling. Build
scripts and manifests in this repository identify the dependency sources and
configuration needed to rebuild the distributed modules. Packaging copies
vcpkg's installed dependency copyright files into `dist/licenses/`.
