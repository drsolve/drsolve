/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_TEST_ENV_H
#define DRSOLVE_TEST_ENV_H

#include <stdlib.h>

#ifdef _WIN32
#include <errno.h>

static inline int dr_test_setenv(const char *name, const char *value, int overwrite)
{
    if (!overwrite && getenv(name) != NULL)
        return 0;
    int error = _putenv_s(name, value);
    if (error) {
        errno = error;
        return -1;
    }
    return 0;
}

static inline int dr_test_unsetenv(const char *name)
{
    /* The Windows CRT removes a variable when its value is empty. */
    return dr_test_setenv(name, "", 1);
}
#else
#define dr_test_setenv setenv
#define dr_test_unsetenv unsetenv
#endif

#endif
