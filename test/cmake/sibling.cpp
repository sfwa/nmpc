#ifdef _MSVC_LANG
static_assert(_MSVC_LANG == 201402L, "NMPC changed the parent's C++ standard");
#else
static_assert(__cplusplus == 201402L, "NMPC changed the parent's C++ standard");
#endif

// This should compile even when NMPC's own targets enable -Wall -Wextra
// -Werror.
int main(int argc, char**) { return 0; }
