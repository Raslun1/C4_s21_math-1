#include <float.h>
#include <stdio.h>
#define ld long double
#define s21_EPS 1e-14L
#define s21_EXP 2.71828182845904L
#define s21_M_PI 3.141592653589793238L
#define s21_2PI (2 * s21_M_PI)
#define s21_M_E 2.7182818284590L
#define s21_INFINITY 1.0 / 0.0
#define s21_NAN 0.0 / 0.0

int s21_abs(int x);
ld s21_acos(double x);
ld s21_asin(double x);
ld s21_atan(double x);
ld s21_sin(double x);
ld s21_cos(double x);
ld s21_tan(double x);
ld s21_exp(double x);
ld s21_fabs(double x);
ld s21_log(double x);
ld s21_sqrt(double x);
ld s21_floor(double x);
ld s21_ceil(double x);
ld s21_fmax(double first, double second);
ld s21_pow(double base, double exp);
ld s21_fmod(double x, double y);
int s21_isnan(double x);
int s21_isinf(double x);