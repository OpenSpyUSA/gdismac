#pragma once

#include <glib.h>

G_BEGIN_DECLS

typedef enum
{
  GDIS_ELEMENT_FAMILY_UNKNOWN = 0,
  GDIS_ELEMENT_FAMILY_NONMETAL,
  GDIS_ELEMENT_FAMILY_NOBLE_GAS,
  GDIS_ELEMENT_FAMILY_ALKALI_METAL,
  GDIS_ELEMENT_FAMILY_ALKALINE_EARTH_METAL,
  GDIS_ELEMENT_FAMILY_TRANSITION_METAL,
  GDIS_ELEMENT_FAMILY_POST_TRANSITION_METAL,
  GDIS_ELEMENT_FAMILY_METALLOID,
  GDIS_ELEMENT_FAMILY_HALOGEN,
  GDIS_ELEMENT_FAMILY_LANTHANIDE,
  GDIS_ELEMENT_FAMILY_ACTINIDE
} GdisElementFamily;

typedef struct
{
  guint atomic_number;
  const gchar *symbol;
  const gchar *name;
  guint row;
  guint column;
  GdisElementFamily family;
  gdouble covalent_radius;
  gdouble vdw_radius;
  gdouble draw_radius;
  gdouble color_rgb[3];
} GdisElementInfo;

void gdis_element_normalize_symbol(const char *input, char output[4]);
const GdisElementInfo *gdis_element_lookup(const char *symbol);
const GdisElementInfo *gdis_element_lookup_atomic_number(guint atomic_number);
const GdisElementInfo *gdis_element_fallback(void);
const char *gdis_element_family_label(GdisElementFamily family);
guint gdis_element_count(void);

G_END_DECLS
