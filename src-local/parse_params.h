/**
# parse_params.h

Low-level parser for the project `key=value` parameter files. It strips
whitespace, ignores `#` comments, lets the last duplicate key win and reads
the file passed as `argv[1]` (falling back to `case.params`). Typed accessors
live in `params.h`.
*/

#ifndef HYPHAL_FLOW_PARSE_PARAMS_H
#define HYPHAL_FLOW_PARSE_PARAMS_H

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PARAMS_MAX_ENTRIES
# define PARAMS_MAX_ENTRIES 256
#endif
#ifndef PARAMS_KEY_LEN
# define PARAMS_KEY_LEN 128
#endif
#ifndef PARAMS_VALUE_LEN
# define PARAMS_VALUE_LEN 256
#endif

typedef struct {
  char key[PARAMS_KEY_LEN];
  char value[PARAMS_VALUE_LEN];
} ParamEntry;

static ParamEntry params_entries[PARAMS_MAX_ENTRIES];
static int params_count = 0;
static bool params_loaded = false;
static bool params_warned_missing = false;
static char params_file[PARAMS_VALUE_LEN] = "case.params";

static inline char * params_trim (char * text)
{
  if (!text)
    return text;
  while (*text && isspace ((unsigned char) *text))
    text++;
  size_t length = strlen (text);
  while (length > 0 && isspace ((unsigned char) text[length - 1]))
    text[--length] = '\0';
  return text;
}

static inline int params_find_key (const char * key)
{
  for (int entry = 0; entry < params_count; entry++)
    if (!strcmp (params_entries[entry].key, key))
      return entry;
  return -1;
}

static inline void params_set_value (const char * key, const char * value)
{
  int entry = params_find_key (key);
  if (entry < 0) {
    if (params_count >= PARAMS_MAX_ENTRIES) {
      fprintf (stderr, "WARNING: parameter limit (%d) reached; ignoring '%s'\n",
               PARAMS_MAX_ENTRIES, key);
      return;
    }
    entry = params_count++;
    strncpy (params_entries[entry].key, key, PARAMS_KEY_LEN - 1);
    params_entries[entry].key[PARAMS_KEY_LEN - 1] = '\0';
  }
  strncpy (params_entries[entry].value, value, PARAMS_VALUE_LEN - 1);
  params_entries[entry].value[PARAMS_VALUE_LEN - 1] = '\0';
}

static inline int params_load (const char * filename)
{
  params_count = 0;
  params_loaded = true;

  FILE * stream = fopen (filename, "r");
  if (!stream) {
    if (!params_warned_missing) {
      fprintf (stderr,
               "WARNING: parameter file '%s' not found; using defaults\n",
               filename);
      params_warned_missing = true;
    }
    return -1;
  }

  char line[PARAMS_KEY_LEN + PARAMS_VALUE_LEN + 64];
  while (fgets (line, sizeof line, stream)) {
    size_t length = strlen (line);
    if (length > 0 && line[length - 1] != '\n') {
      int next = fgetc (stream);
      if (next != '\n' && next != EOF) {
        while ((next = fgetc (stream)) != '\n' && next != EOF)
          ;
        fprintf (stderr, "WARNING: overlong parameter line ignored\n");
        continue;
      }
    }
    char * comment = strchr (line, '#');
    if (comment)
      *comment = '\0';

    char * separator = strchr (line, '=');
    if (!separator)
      continue;

    *separator = '\0';
    char * key = params_trim (line);
    char * value = params_trim (separator + 1);
    if (key && value && *key && *value)
      params_set_value (key, value);
  }

  fclose (stream);
  return 0;
}

static inline void params_init_from_argv (int argc, const char * argv[])
{
  const char * selected = argc > 1 && argv[1] && argv[1][0]
    ? argv[1] : "case.params";
  strncpy (params_file, selected, PARAMS_VALUE_LEN - 1);
  params_file[PARAMS_VALUE_LEN - 1] = '\0';
  (void) params_load (params_file);
}

static inline void params_ensure_loaded (void)
{
  if (!params_loaded)
    (void) params_load (params_file);
}

static inline const char * params_raw (const char * key)
{
  params_ensure_loaded ();
  int entry = params_find_key (key);
  return entry >= 0 ? params_entries[entry].value : NULL;
}

#endif // HYPHAL_FLOW_PARSE_PARAMS_H
