#pragma once
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>

namespace stl2curves::detail {
inline const char *failure_message(const Standard_Failure &error,
                                   const char *fallback = "OpenCascade failure") {
#if OCC_VERSION_HEX >= 0x080000
    const char *message = error.what();
#else
    // Older native development kernels do not expose std::exception::what().
    const char *message = error.GetMessageString();
#endif
    return message ? message : fallback;
}
} // namespace stl2curves::detail
