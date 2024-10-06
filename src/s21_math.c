#include "s21_math.h"

int s21_isnan(double x) { return x != x; }

// int s21_isinf(double x) { return !s21_isnan(x) && s21_isnan(x - x); }
int s21_isinf(double x) { return !s21_isnan(x) && x == x + 1 && x != 0; }

ld s21_atan(double x) {
  ld result = 0.0;
  if (x > 1)
    result = (s21_M_PI / 2) - s21_atan(1 / x);
  else if (x < -1)
    result = -(s21_M_PI / 2) - s21_atan(1 / x);
  else if (x == 1)
    result = s21_M_PI / 4;
  else if (x == -1)
    result = -s21_M_PI / 4;
  else if (s21_isnan(x))
    result = s21_NAN;
  else if (s21_isinf(x))
    result = (x > 0) ? s21_M_PI / 2 : -s21_M_PI / 2;
  else {
    ld sign = -1.0, numerator = x, shag = 0.0;

    for (ld i = 0.0; i < 1000.0; i++) {
      sign *= -1.0;
      shag = sign / (2.0 * i + 1);
      if (x > -1.0 && x < 1.0)
        shag *= numerator;
      else
        shag /= numerator;
      result += shag;
      numerator *= x * x;
    }
  }
  return result;
}

ld s21_asin(double x) {
  ld result = 0.0;
  if (x == 1.0)
    result = s21_M_PI / 2;
  else if (x == -1.0)
    result = -s21_M_PI / 2;
  else {
    if (x > -1 && x < 1) {
      result = s21_atan(x / (s21_sqrt(1 - x * x)));
    } else
      result = s21_NAN;
  }
  return result;
}

ld s21_acos(double x) {
  ld result = 0.0;
  if (x == 0.0)
    result = s21_M_PI / 2;
  else if (x == -1)
    result = s21_M_PI;
  else if (x == 1)
    result = 0.0;
  else {
    if (x > -1 && x < 1) {
      if (x > 0.0)
        result = s21_atan(s21_sqrt(1 - x * x) / x);
      else
        result = s21_M_PI + s21_atan(s21_sqrt(1 - x * x) / x);
    } else
      result = s21_NAN;
  }
  return result;
}

ld s21_cos(double x) {
  ld sumCos = 0.0;
  if (s21_isnan(x)) {
    sumCos = s21_NAN;
  } else if (s21_isinf(x)) {
    sumCos = s21_NAN;
  } else {
    for (; x < -s21_2PI || x > s21_2PI;) {
      if (x < -s21_2PI)
        x += s21_2PI;
      else
        x -= s21_2PI;
    }

    ld sign = -1.0, numerator = 1.0, denominator = 1.0;
    for (ld i = 1.0; i < 1000.0; i++) {
      sign *= -1.0;
      sumCos += sign * numerator / denominator;
      numerator *= x * x;
      denominator *=
          (2.0 * i) * (2.0 * i - 1.0);  // берем нынешний и прошлый. (факториал)
    }
  }
  return sumCos;
}

ld s21_sin(double x) {
  ld sumSin = 0.0;
  if (s21_isnan(x)) {
    sumSin = s21_NAN;
  } else if (s21_isinf(x)) {
    sumSin = s21_NAN;
  } else {
    for (; x < -s21_2PI || x > s21_2PI;) {
      if (x < -s21_2PI)
        x += s21_2PI;
      else
        x -= s21_2PI;
    }
    ld sign = -1.0, numerator = x, denominator = 1.0;
    for (ld i = 1.0; i < 1000.0; i++) {
      sign *= -1.0;
      sumSin += sign * numerator / denominator;
      numerator *= x * x;
      denominator *= (2.0 * i + 1) * (2 * i);
    }
  }
  return sumSin;
}

ld s21_tan(double x) {
  ld result = 0.0;
  if (s21_isnan(x) || s21_isinf(x))
    result = s21_NAN;
  else {
    result = s21_sin(x) / s21_cos(x);
  }
  return result;
}

ld s21_fabs(double x) {
  ld res;
  if (s21_isnan(x)) {
    res = s21_NAN;
  } else if (s21_isinf(x)) {
    res = s21_INFINITY;
  } else if (x < 0) {
    res = -x;
  } else
    res = x;
  return res;
}

int s21_abs(int x) {
  int res;
  if (s21_isnan(x)) {
    res = s21_NAN;
  } else if (s21_isinf(x)) {
    res = s21_INFINITY;
  } else if (x < 0) {
    res = -x;
  } else
    res = x;
  return res;
}

ld s21_exp(double x) {
  ld sum_exp = 0;
  if (!s21_isnan(x) && !s21_isinf(x)) {
    sum_exp = 1;
    ld add_value = 1;
    ld i = 1;
    while (s21_fabs(add_value) > s21_EPS) {
      add_value *= x / i;
      i++;
      sum_exp += add_value;
      if (sum_exp > DBL_MAX) {
        sum_exp = s21_INFINITY;
        break;
      }
    }
  } else if (s21_isnan(x))
    sum_exp = s21_NAN;
  else if (s21_isinf(x) && x > 0)
    sum_exp = s21_INFINITY;
  return sum_exp;
}

ld s21_ceil(double x) {
  ld res = 0;
  if (!s21_isnan(x) && !s21_isinf(x)) {
    res = (long long int)x;
    if (x == res)
      res = x;
    else if (x > 0)
      res++;
  } else
    res = x;
  return res;
}

ld s21_floor(double x) {
  ld res = 0;
  if (!s21_isnan(x) && !s21_isinf(x)) {
    res = (long long int)x;
    if (x == res)
      res = x;
    else if (x < 0)
      res--;
  } else
    res = x;
  return res;
}

ld s21_log(double x) {
  ld result = 0.0;
  ld compare = 0.0;
  if (s21_isnan(x) || x < 0.0) {
    result = s21_NAN;
  } else if (s21_isinf(x) && x > 0.0) {
    result = s21_INFINITY;
  } else if (x == 0.0) {
    result = -s21_INFINITY;
  } else {
    for (int i = 0; i < 100; i++) {
      compare = result;
      result = compare + 2 * (x - s21_exp(compare)) / (x + s21_exp(compare));
    }
  }
  return result;
}

ld s21_fmax(double first, double second) {
  ld max = 0.0;
  if (s21_isnan(second))
    max = first;
  else if (s21_isnan(first))
    max = second;
  else if (s21_isinf(first) || s21_isinf(second))
    max = s21_INFINITY;
  else if (first >= second)
    max = first;
  else if (second > first)
    max = second;
  return max;
}

ld s21_sqrt(double x) {
  ld start = 0.0, end = s21_fmax(1.0, x), middle = 0.0;
  if (x >= 0.0) {
    if (s21_isnan(x))
      middle = s21_NAN;
    else if (s21_isinf(x))
      middle = s21_INFINITY;
    else {
      middle = (start + end) / 2;
      while ((middle - start) > s21_EPS) {
        if (middle * middle > x)
          end = middle;
        else
          start = middle;
        middle = (start + end) / 2;
      }
    }
  } else
    middle = s21_NAN;
  return middle;
}

ld s21_pow(double base, double exp) {
  ld result = 0;
  if ((base == 0 && exp == 0) || exp == 0) {
    result = 1;
  } else if (base == 0 && exp < 0) {
    result = s21_INFINITY;
  } else if (base == 0) {
    if (!s21_isnan(exp) && !s21_isinf(exp)) {
      result = 0;
    } else if (s21_isnan(exp)) {
      result = s21_NAN;
    } else if (s21_isinf(exp)) {
      result = 0;
    }
  } else if (base == -s21_INFINITY && exp == -s21_INFINITY)
    result = 0.0;
  else if (!s21_isnan(base) && exp == -s21_INFINITY) {
    if (base == 1 || base == -1)
      result = 1.0;
    else if (base > 1 || base < -1)
      result = 0.0;
    else
      result = s21_INFINITY;
  } else if (!s21_isnan(base) && exp == s21_INFINITY) {
    if (base == 1 || base == -1)
      result = 1.0;
    else if (base > -1 && base < 1)
      result = 0.0;
    else
      result = s21_INFINITY;
  } else if (base == 1 && s21_isnan(exp))
    result = 1.0;
  else if ((s21_isnan(base) && s21_isinf(exp)) ||
           (s21_isnan(exp) && s21_isinf(base)))
    result = s21_NAN;
  else {
    result = s21_exp(s21_fabs(exp) * s21_log(s21_fabs(base)));  // log
    if (base < 0 && s21_fmod(exp, 2.)) result = -result;        // fmod
    if (exp < 0) result = 1.0 / result;
  }
  return result;
}

ld s21_fmod(double x, double y) {
  ld result;
  if (y == 0) {
    result = s21_NAN;
  } else if (x == 0) {
    if (s21_isnan(y)) {
      result = s21_NAN;
    } else
      result = 0;
  } else if (s21_isnan(x) || s21_isnan(y)) {
    result = s21_NAN;
  } else if (s21_isinf(x)) {
    result = s21_NAN;
  } else if (y == s21_INFINITY) {
    result = 1;
  } else if (y == -s21_INFINITY) {
    result = -1;
  } else {
    long long int mod = (long long int)(x / y);
    result = x - y * mod;
  }
  return result;
}