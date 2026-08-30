/**
# params.h

Typed accessors for values loaded by `parse_params.h`. Missing keys use the
call-site default; malformed values emit a warning and use that same default.
*/

#ifndef HYPHAL_FLOW_PARAMS_H
#define HYPHAL_FLOW_PARAMS_H

#include "parse_params.h"

#include <errno.h>
#include <limits.h>
#include <math.h>

static inline bool param_text_equal_ignore_case (const char * left,
                                                 const char * right)
{
  while (*left && *right) {
    if (tolower ((unsigned char) *left) != tolower ((unsigned char) *right))
      return false;
    left++, right++;
  }
  return *left == '\0' && *right == '\0';
}

static inline const char * param_string (const char * key,
                                         const char * default_value)
{
  const char * value = params_raw (key);
  return value ? value : default_value;
}

static inline double param_double (const char * key, double default_value)
{
  const char * text = params_raw (key);
  if (!text)
    return default_value;

  errno = 0;
  char * end = NULL;
  double value = strtod (text, &end);
  while (end && *end && isspace ((unsigned char) *end))
    end++;
  if (errno || end == text || (end && *end) || !isfinite (value)) {
    fprintf (stderr, "WARNING: invalid double '%s=%s'; using %g\n",
             key, text, default_value);
    return default_value;
  }
  return value;
}

static inline int param_int (const char * key, int default_value)
{
  const char * text = params_raw (key);
  if (!text)
    return default_value;

  errno = 0;
  char * end = NULL;
  long value = strtol (text, &end, 10);
  while (end && *end && isspace ((unsigned char) *end))
    end++;
  if (errno || end == text || (end && *end) ||
      value < INT_MIN || value > INT_MAX) {
    fprintf (stderr, "WARNING: invalid integer '%s=%s'; using %d\n",
             key, text, default_value);
    return default_value;
  }
  return (int) value;
}

static inline bool param_bool (const char * key, bool default_value)
{
  const char * text = params_raw (key);
  if (!text)
    return default_value;
  if (param_text_equal_ignore_case (text, "1") ||
      param_text_equal_ignore_case (text, "true") ||
      param_text_equal_ignore_case (text, "yes") ||
      param_text_equal_ignore_case (text, "on"))
    return true;
  if (param_text_equal_ignore_case (text, "0") ||
      param_text_equal_ignore_case (text, "false") ||
      param_text_equal_ignore_case (text, "no") ||
      param_text_equal_ignore_case (text, "off"))
    return false;

  fprintf (stderr, "WARNING: invalid boolean '%s=%s'; using %s\n",
           key, text, default_value ? "true" : "false");
  return default_value;
}

#endif // HYPHAL_FLOW_PARAMS_H
