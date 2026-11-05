"""Refuses to load an x86-64-v3 build on a CPU without the instructions it uses."""

import platform
import sys

# Set by CMake when it copies the package: "ON" when the extension is built for x86-64-v3.
_X86_64_V3 = "@CHARGEFW_X86_64_V3@"

# The features of x86-64-v2 and x86-64-v3, as Linux names them in /proc/cpuinfo (pni is SSE3, abm
# is LZCNT).
_X86_64_V3_FLAGS = frozenset(
    {
        "cx16", "lahf_lm", "pni", "popcnt", "sse4_1", "sse4_2", "ssse3",
        "avx", "avx2", "bmi1", "bmi2", "f16c", "fma", "abm", "movbe", "xsave",
    }
)  # fmt: skip


def _supports_x86_64_v3() -> bool:
    if sys.platform == "win32":
        # Windows reports AVX2 only from Windows 10 version 2004 (build 19041); before, the answer
        # is always no.
        if sys.getwindowsversion().build < 19041:
            return True
        import ctypes

        pf_avx2_instructions_available = 40
        return bool(
            ctypes.windll.kernel32.IsProcessorFeaturePresent(pf_avx2_instructions_available)
        )
    try:
        with open("/proc/cpuinfo", encoding="ascii", errors="replace") as cpuinfo:
            for line in cpuinfo:
                if line.startswith("flags"):
                    return _X86_64_V3_FLAGS.issubset(line.partition(":")[2].split())
    except OSError:
        pass
    return True  # Unknown: let the extension load.


def _check_cpu() -> None:
    if (
        _X86_64_V3 == "ON"
        and platform.machine().lower() in {"x86_64", "amd64"}
        and not _supports_x86_64_v3()
    ):
        raise ImportError(
            "ChargeFW needs an x86-64 CPU with AVX2 and FMA (x86-64-v3: Intel Haswell, AMD "
            "Excavator, or newer). On this CPU, build the library and command line from source "
            "with CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON."
        )


_check_cpu()
