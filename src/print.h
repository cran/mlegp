
#ifndef __print__
#define __print__

#include <stdarg.h>

#ifdef __useR__
#include <R_ext/Print.h>
#include <R.h>
#else 
#include <stdio.h>
#endif



void printout(char *str, ...) {
	va_list argp;
	va_start(argp, str);
	#ifdef __useR__
	Rvprintf(str, argp);
	R_FlushConsole();
	#else
	vfprintf(stdout, str, argp);
	fflush(stdout);
	#endif
 	va_end(argp);
}

void printerr(char *str, ...) {
	va_list argp;
	va_start(argp, str);
	#ifdef __useR__
	REvprintf(str, argp);
	R_FlushConsole();
	#else
	vfprintf(stderr, str, argp);
	fflush(stderr);
	#endif
 	va_end(argp);
}

#endif


