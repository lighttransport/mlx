@echo off

rem This script generates a C++ function that provides the CPU
rem code for use with kernel generation.
rem Copyright @ 2024 Light Transport Entertainment Inc.
rem Copyright @ 2023-24 Apple Inc.


set OUTPUT_FILE=%1
set MSVCBIN=%2
set SRCDIR=%3

> %OUTPUT_FILE% (
@echo.const char* get_kernel_preamble(^) {
@echo.return R"preamble(
)

call %MSVCBIN% /std:c++17 /E /I%SRCDIR% /Tp %SRCDIR%/mlx/backend/common/compiled_preamble.h >> %OUTPUT_FILE%

>> %OUTPUT_FILE% (
@echo.using namespace mlx::core;
@echo.using namespace mlx::core::detail;
@echo.^)preamble";
@echo.}
)
