/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <windows.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "dixon_wrapper.h"
#include "../test_env.h"
int main(int argc, char **argv)
{
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) printf("<%s>\n", argv[i]);
        return strcmp(argv[1], "--exit-7") == 0 ? 7 : 0;
    }
    assert(dr_test_setenv("DRSOLVE_PML_MUL", "fft", 1) == 0);
    assert(strcmp(getenv("DRSOLVE_PML_MUL"), "fft") == 0);
    assert(dr_test_setenv("DRSOLVE_PML_MUL", "other", 0) == 0);
    assert(strcmp(getenv("DRSOLVE_PML_MUL"), "fft") == 0);
    char path[MAX_PATH], error[512];
    assert(GetModuleFileNameA(NULL, path, sizeof(path)) > 0);
    for (int mode = 0; mode < 4; ++mode) {
        for (int submode = 0; submode < (mode == 3 ? 3 : 1); ++submode) {
            for (int method = -1; method <= 5; ++method) {
                dixon_gui_request_t request = {0};
                dixon_gui_response_t response;
                dixon_gui_response_init(&response);
                request.mode = mode;
                request.dixon_path = path;
                request.polynomials = "x+y\r\nx*y";
                request.eliminate_vars = "x";
                request.field = "257";
                request.omega = "2.81";
                request.random_degrees = "[3]*3";
                request.random_mode = submode;
                request.step1_method = request.step4_method = method;
                if (!dixon_gui_run_request(&request, &response, error, sizeof(error))) {
                    fprintf(stderr, "%s\n", error); return 1;
                }
                assert(response.exit_code == 0);
                const char *out = response.stdout_text;
                assert(strstr(out, "<257>"));
                assert(!strstr(out, "--ideal"));
                assert(!!strstr(out, "<--step1>") == (method >= 0));
                assert(!!strstr(out, "<--step4>") == (method >= 0));
                assert(!!strstr(out, "<-s>") == (mode == 0 || (mode == 3 && submode == 1)));
                assert(!!strstr(out, "<--comp>") == (mode == 2 || (mode == 3 && submode == 2)));
                assert(!!strstr(out, "<--random>") == (mode == 3));
                assert(strstr(out, mode == 3 ? "<[3]*3>" : "<x+y, x*y>"));
                dixon_gui_response_clear(&response);
            }
        }
    }
    /* Raw CLI input needs no field, preserves quoting, and ignores GUI methods. */
    const char *commands[] = {
        "--help",
        "-s \"x^2-1, y-x\" 257",
        "-r \"[2]*3\" 65537 --seed 12345",
        "--step1 0 --step4 5 -o \"result with spaces.dr\"",
        "\"literal & argument\" \"\"",
        "--exit-7"
    };
    const char *expected[] = {
        "<--help>\n",
        "<-s>\n<x^2-1, y-x>\n<257>\n",
        "<-r>\n<[2]*3>\n<65537>\n<--seed>\n<12345>\n",
        "<--step1>\n<0>\n<--step4>\n<5>\n<-o>\n<result with spaces.dr>\n",
        "<literal & argument>\n<>\n",
        "<--exit-7>\n"
    };
    for (unsigned i = 0; i < sizeof(commands) / sizeof(*commands); ++i) {
        dixon_gui_request_t request = {0};
        dixon_gui_response_t response;
        dixon_gui_response_init(&response);
        request.mode = DIXON_GUI_MODE_COMMAND;
        request.dixon_path = path;
        request.command_args = commands[i];
        request.step1_method = request.step4_method = 4;
        assert(dixon_gui_run_request(&request, &response, error, sizeof(error)));
        assert(response.exit_code == (i == 5 ? 7 : 0));
        /* Windows text stdout uses CRLF; compare the actual argv after normalization. */
        char *src = response.stdout_text, *dst = src;
        while (*src) {
            if (*src != '\r') *dst++ = *src;
            ++src;
        }
        *dst = '\0';
        assert(strcmp(response.stdout_text, expected[i]) == 0);
        dixon_gui_response_clear(&response);
    }
    dixon_gui_request_t empty = {0};
    dixon_gui_response_t response;
    dixon_gui_response_init(&response);
    empty.mode = DIXON_GUI_MODE_COMMAND;
    empty.dixon_path = path;
    empty.command_args = " \t ";
    assert(!dixon_gui_run_request(&empty, &response, error, sizeof(error)));
    assert(strstr(error, "Enter CLI arguments"));
    dixon_gui_response_clear(&response);
    puts("PASS: Windows environment, 42 form requests, raw CLI quoting and exit codes");
    return 0;
}
