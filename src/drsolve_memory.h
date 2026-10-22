/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_MEMORY_H
#define DRSOLVE_MEMORY_H

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

/* Process lifetime peak resident memory, including all worker threads. */
static inline double drsolve_peak_memory_mib(void)
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS counters;
    /* Resolve dynamically to avoid adding a PSAPI link dependency. */
    typedef BOOL (WINAPI *memory_info_fn)(HANDLE, PPROCESS_MEMORY_COUNTERS, DWORD);
    HMODULE psapi = LoadLibraryA("psapi.dll");
    if (!psapi) return -1.0;
    memory_info_fn memory_info = (memory_info_fn) GetProcAddress(psapi, "GetProcessMemoryInfo");
    BOOL ok = memory_info && memory_info(GetCurrentProcess(), &counters, sizeof(counters));
    FreeLibrary(psapi);
    if (!ok) return -1.0;
    return (double) counters.PeakWorkingSetSize / (1024.0 * 1024.0);
#else
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1.0;
#ifdef __APPLE__
    /* macOS reports bytes; Linux and other Unix targets report KiB. */
    return (double) usage.ru_maxrss / (1024.0 * 1024.0);
#else
    return (double) usage.ru_maxrss / 1024.0;
#endif
#endif
}

#endif
