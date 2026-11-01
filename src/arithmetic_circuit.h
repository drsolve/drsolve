/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_ARITHMETIC_CIRCUIT_H
#define DRSOLVE_ARITHMETIC_CIRCUIT_H
#include <stdio.h>
#include "bezout_block.h"
typedef struct dr_circuit dr_circuit;
/* DRCIRCUIT 1, prime P, eliminate x|y (default y); nodes 0=x, 1=y.
 * Each following const/add/sub/mul/pow/lin line appends one node.
 * References must point backwards. End with outputs A B.
 * lin C c0 id0 c1 id1 ... denotes C + sum(ci * node[id_i]).
 * pow A E requires an unsigned integer exponent. Comments start with #.
 * This experimental solver certifies equal positive total degrees and
 * constant nonzero leading coefficients in the elimination variable.
 * Other circuits are rejected; no expansion or verification bypass occurs. */
dr_circuit *dr_circuit_read(FILE *input,char *error,size_t error_size);
void dr_circuit_clear(dr_circuit *circuit);
void dr_circuit_source(dr_bezout_source *source,const dr_circuit *circuit);
ulong dr_circuit_prime(const dr_circuit *circuit);
slong dr_circuit_degree(const dr_circuit *circuit);
int dr_circuit_cli_main(int argc,char **argv);
#endif
