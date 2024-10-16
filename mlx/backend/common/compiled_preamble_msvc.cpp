// Since MSVC cannot embed large string literal to source file,
// Use resource file.
//
// Reference:
// https://blog.kowalczyk.info/article/zy/embedding-binary-resources-on-windows.html
//

#include <string>

#ifdef _MSC_VER
#define COMPILED_PREAMBLE_SOURCE_ID 21264
#include <windows.h>

static std::string get_compiled_preamble_source_str() {

    HMODULE module_handle = GetModuleHandle(nullptr);
    HGLOBAL     res_handle = nullptr;
    HRSRC       res;
    char *      res_data;
    DWORD       res_size;

    res = FindResource(module_handle, MAKEINTRESOURCE(COMPILED_PREAMBLE_SOURCE_ID), RT_RCDATA);
    if (!res)
        return std::string();
    res_handle = LoadResource(nullptr, res);
    if (!res_handle)
        return std::string();
    res_data = (char*)LockResource(res_handle);
    res_size = SizeofResource(nullptr, res);
    /* you can now use the resource data */

    return std::string(res_data, res_size);
}

const char* get_kernel_preamble() {
  static std::string code;
  if (code.empty()) {
    code = get_compiled_preamble_source_str();
  }

  return code.c_str();
}

#endif
