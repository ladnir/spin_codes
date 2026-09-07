# libOTe integration

This implementation is an overlay rooted exactly like a libOTe checkout:

- `libOTe/Tools/RiffleCode/` contains the installed public header, the compiled
  wrapper, the fixed-width kernel, and its generated schedule;
- `libOTe_Tests/` contains the native `CLP` test adapter and independent test
  support;
- `integration/libote-root.patch` adds the compiled wrapper and native test to
  their existing targets while keeping AVX2/VPCLMUL flags source-local.

Copy `libOTe/` and `libOTe_Tests/` into the root of a libOTe checkout, then
apply the patch. The whitespace flags tolerate upstream's mixed line endings:

```bash
git apply --ignore-space-change --ignore-whitespace \
  integration/libote-root.patch
```

The patch registers the test in `libOTe_Tests/UnitTests.cpp`:

```cpp
#include "libOTe_Tests/StructuredSpin_Tests.h"
// Inside tests_libOTe::Tests initialization:
tc.add("StructuredSpin_correctness_test ", StructuredSpin_correctness_test);
```

The existing libOTe install rule recursively installs `*.h` below `libOTe/`,
so no install-rule change is required. A downstream consumer uses the normal
package target:

```cmake
find_package(libOTe REQUIRED)
target_link_libraries(my_target PRIVATE oc::libOTe)
```

and includes:

```cpp
#include <libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.h>
```

The construction requires an AVX2 and VPCLMULQDQ-capable processor. The patch
applies those flags only to the implementation and native test source; they do
not become usage requirements of `oc::libOTe`.
