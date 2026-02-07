/**
# params.h

Lightweight runtime parameter reader for hyphal-flow simulation cases.

Parameter files are plain `key=value` text files with optional `#` comments.
Typical usage in a case file:

```c
#include "params.h"
...
params_init_from_argv(argc, argv);
Ec_h = param_double("Ec_h", 0.0);
```
*/

#ifndef HYPHAL_FLOW_PARAMS_H
#define HYPHAL_FLOW_PARAMS_H

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PARAMS_MAX_ENTRIES
#define PARAMS_MAX_ENTRIES 256
#endif

#ifndef PARAMS_KEY_LEN
#define PARAMS_KEY_LEN 128
#endif

#ifndef PARAMS_VALUE_LEN
#define PARAMS_VALUE_LEN 256
#endif

typedef struct {
  char key[PARAMS_KEY_LEN];
  char value[PARAMS_VALUE_LEN];
} ParamEntry;

static ParamEntry _params_entries[PARAMS_MAX_ENTRIES];
static int _params_count = 0;
static bool _params_loaded = false;
static bool _params_warned_missing = false;
static char _params_file[PARAMS_VALUE_LEN] = "case.params";

static inline char *_params_trim (char *s)
{
  if (!s)
    return s;
  while (*s && isspace((unsigned char)*s))
    s++;
  size_t n = strlen(s);
  while (n > 0 && isspace((unsigned char)s[n - 1]))
    s[--n] = '\0';
  return s;
}

static inline int _params_find_key (const char *key)
{
  for (int i = 0; i < _params_count; i++)
    if (!strcmp(_params_entries[i].key, key))
      return i;
  return -1;
}

static inline void _params_set_value (const char *key, const char *value)
{
  int idx = _params_find_key(key);
  if (idx < 0) {
    if (_params_count >= PARAMS_MAX_ENTRIES) {
      fprintf(stderr,
              "WARNING: params.h entry limit reached (%d), skipping '%s'\n",
              PARAMS_MAX_ENTRIES, key);
      return;
    }
    idx = _params_count++;
    strncpy(_params_entries[idx].key, key, PARAMS_KEY_LEN - 1);
    _params_entries[idx].key[PARAMS_KEY_LEN - 1] = '\0';
  }
  strncpy(_params_entries[idx].value, value, PARAMS_VALUE_LEN - 1);
  _params_entries[idx].value[PARAMS_VALUE_LEN - 1] = '\0';
}

static inline int params_load (const char *filename)
{
  _params_count = 0;
  _params_loaded = true;

  FILE *fp = fopen(filename, "r");
  if (!fp) {
    if (!_params_warned_missing) {
      fprintf(stderr,
              "WARNING: Parameter file '%s' not found. Using defaults.\n",
              filename);
      _params_warned_missing = true;
    }
    return -1;
  }

  char line[PARAMS_KEY_LEN + PARAMS_VALUE_LEN + 64];
  while (fgets(line, sizeof(line), fp)) {
    char *comment = strchr(line, '#');
    if (comment)
      *comment = '\0';

    char *eq = strchr(line, '=');
    if (!eq)
      continue;

    *eq = '\0';
    char *key = _params_trim(line);
    char *value = _params_trim(eq + 1);

    if (!key || !value || !*key || !*value)
      continue;

    _params_set_value(key, value);
  }

  fclose(fp);
  return 0;
}

static inline void params_init_from_argv (int argc, const char *argv[])
{
  if (argc > 1 && argv[1] && argv[1][0]) {
    strncpy(_params_file, argv[1], PARAMS_VALUE_LEN - 1);
    _params_file[PARAMS_VALUE_LEN - 1] = '\0';
  } else {
    strncpy(_params_file, "case.params", PARAMS_VALUE_LEN - 1);
    _params_file[PARAMS_VALUE_LEN - 1] = '\0';
  }
  (void) params_load(_params_file);
}

static inline void _params_ensure_loaded (void)
{
  if (!_params_loaded)
    (void) params_load(_params_file);
}

static inline const char *param_string (const char *key,
                                        const char *default_value)
{
  _params_ensure_loaded();
  int idx = _params_find_key(key);
  return idx >= 0 ? _params_entries[idx].value : default_value;
}

static inline double param_double (const char *key, double default_value)
{
  const char *s = param_string(key, NULL);
  if (!s)
    return default_value;

  errno = 0;
  char *end = NULL;
  double v = strtod(s, &end);
  while (end && *end && isspace((unsigned char)*end))
    end++;

  if (errno != 0 || end == s || (end && *end != '\0')) {
    fprintf(stderr,
            "WARNING: Invalid double for '%s' ('%s'), using default %g\n",
            key, s, default_value);
    return default_value;
  }

  return v;
}

static inline int param_int (const char *key, int default_value)
{
  const char *s = param_string(key, NULL);
  if (!s)
    return default_value;

  errno = 0;
  char *end = NULL;
  long v = strtol(s, &end, 10);
  while (end && *end && isspace((unsigned char)*end))
    end++;

  if (errno != 0 || end == s || (end && *end != '\0') ||
      v < INT_MIN || v > INT_MAX) {
    fprintf(stderr,
            "WARNING: Invalid int for '%s' ('%s'), using default %d\n",
            key, s, default_value);
    return default_value;
  }

  return (int) v;
}

static inline bool param_bool (const char *key, bool default_value)
{
  const char *s = param_string(key, NULL);
  if (!s)
    return default_value;

  if (!strcasecmp(s, "1") || !strcasecmp(s, "true") ||
      !strcasecmp(s, "yes") || !strcasecmp(s, "on"))
    return true;
  if (!strcasecmp(s, "0") || !strcasecmp(s, "false") ||
      !strcasecmp(s, "no") || !strcasecmp(s, "off"))
    return false;

  fprintf(stderr,
          "WARNING: Invalid bool for '%s' ('%s'), using default %d\n",
          key, s, default_value ? 1 : 0);
  return default_value;
}

#endif // HYPHAL_FLOW_PARAMS_H
