#include "gdis_elements.h"

#include <string.h>

typedef struct
{
  const gchar *symbol;
  const gchar *name;
} GdisElementSeed;

static GdisElementInfo gdis_elements[119];
static gboolean gdis_elements_initialized = FALSE;

static const GdisElementSeed gdis_element_seeds[119] = {
  { "", "Unknown" },
  { "H", "Hydrogen" },
  { "He", "Helium" },
  { "Li", "Lithium" },
  { "Be", "Beryllium" },
  { "B", "Boron" },
  { "C", "Carbon" },
  { "N", "Nitrogen" },
  { "O", "Oxygen" },
  { "F", "Fluorine" },
  { "Ne", "Neon" },
  { "Na", "Sodium" },
  { "Mg", "Magnesium" },
  { "Al", "Aluminum" },
  { "Si", "Silicon" },
  { "P", "Phosphorus" },
  { "S", "Sulfur" },
  { "Cl", "Chlorine" },
  { "Ar", "Argon" },
  { "K", "Potassium" },
  { "Ca", "Calcium" },
  { "Sc", "Scandium" },
  { "Ti", "Titanium" },
  { "V", "Vanadium" },
  { "Cr", "Chromium" },
  { "Mn", "Manganese" },
  { "Fe", "Iron" },
  { "Co", "Cobalt" },
  { "Ni", "Nickel" },
  { "Cu", "Copper" },
  { "Zn", "Zinc" },
  { "Ga", "Gallium" },
  { "Ge", "Germanium" },
  { "As", "Arsenic" },
  { "Se", "Selenium" },
  { "Br", "Bromine" },
  { "Kr", "Krypton" },
  { "Rb", "Rubidium" },
  { "Sr", "Strontium" },
  { "Y", "Yttrium" },
  { "Zr", "Zirconium" },
  { "Nb", "Niobium" },
  { "Mo", "Molybdenum" },
  { "Tc", "Technetium" },
  { "Ru", "Ruthenium" },
  { "Rh", "Rhodium" },
  { "Pd", "Palladium" },
  { "Ag", "Silver" },
  { "Cd", "Cadmium" },
  { "In", "Indium" },
  { "Sn", "Tin" },
  { "Sb", "Antimony" },
  { "Te", "Tellurium" },
  { "I", "Iodine" },
  { "Xe", "Xenon" },
  { "Cs", "Cesium" },
  { "Ba", "Barium" },
  { "La", "Lanthanum" },
  { "Ce", "Cerium" },
  { "Pr", "Praseodymium" },
  { "Nd", "Neodymium" },
  { "Pm", "Promethium" },
  { "Sm", "Samarium" },
  { "Eu", "Europium" },
  { "Gd", "Gadolinium" },
  { "Tb", "Terbium" },
  { "Dy", "Dysprosium" },
  { "Ho", "Holmium" },
  { "Er", "Erbium" },
  { "Tm", "Thulium" },
  { "Yb", "Ytterbium" },
  { "Lu", "Lutetium" },
  { "Hf", "Hafnium" },
  { "Ta", "Tantalum" },
  { "W", "Tungsten" },
  { "Re", "Rhenium" },
  { "Os", "Osmium" },
  { "Ir", "Iridium" },
  { "Pt", "Platinum" },
  { "Au", "Gold" },
  { "Hg", "Mercury" },
  { "Tl", "Thallium" },
  { "Pb", "Lead" },
  { "Bi", "Bismuth" },
  { "Po", "Polonium" },
  { "At", "Astatine" },
  { "Rn", "Radon" },
  { "Fr", "Francium" },
  { "Ra", "Radium" },
  { "Ac", "Actinium" },
  { "Th", "Thorium" },
  { "Pa", "Protactinium" },
  { "U", "Uranium" },
  { "Np", "Neptunium" },
  { "Pu", "Plutonium" },
  { "Am", "Americium" },
  { "Cm", "Curium" },
  { "Bk", "Berkelium" },
  { "Cf", "Californium" },
  { "Es", "Einsteinium" },
  { "Fm", "Fermium" },
  { "Md", "Mendelevium" },
  { "No", "Nobelium" },
  { "Lr", "Lawrencium" },
  { "Rf", "Rutherfordium" },
  { "Db", "Dubnium" },
  { "Sg", "Seaborgium" },
  { "Bh", "Bohrium" },
  { "Hs", "Hassium" },
  { "Mt", "Meitnerium" },
  { "Ds", "Darmstadtium" },
  { "Rg", "Roentgenium" },
  { "Cn", "Copernicium" },
  { "Nh", "Nihonium" },
  { "Fl", "Flerovium" },
  { "Mc", "Moscovium" },
  { "Lv", "Livermorium" },
  { "Ts", "Tennessine" },
  { "Og", "Oganesson" }
};

static const guint8 gdis_element_positions[119][2] = {
  { 0, 0 },
  { 1, 1 }, { 18, 1 }, { 1, 2 }, { 2, 2 }, { 13, 2 }, { 14, 2 }, { 15, 2 }, { 16, 2 }, { 17, 2 }, { 18, 2 },
  { 1, 3 }, { 2, 3 }, { 13, 3 }, { 14, 3 }, { 15, 3 }, { 16, 3 }, { 17, 3 }, { 18, 3 },
  { 1, 4 }, { 2, 4 }, { 3, 4 }, { 4, 4 }, { 5, 4 }, { 6, 4 }, { 7, 4 }, { 8, 4 }, { 9, 4 }, { 10, 4 }, { 11, 4 }, { 12, 4 }, { 13, 4 }, { 14, 4 }, { 15, 4 }, { 16, 4 }, { 17, 4 }, { 18, 4 },
  { 1, 5 }, { 2, 5 }, { 3, 5 }, { 4, 5 }, { 5, 5 }, { 6, 5 }, { 7, 5 }, { 8, 5 }, { 9, 5 }, { 10, 5 }, { 11, 5 }, { 12, 5 }, { 13, 5 }, { 14, 5 }, { 15, 5 }, { 16, 5 }, { 17, 5 }, { 18, 5 },
  { 1, 6 }, { 2, 6 }, { 3, 6 }, { 4, 9 }, { 5, 9 }, { 6, 9 }, { 7, 9 }, { 8, 9 }, { 9, 9 }, { 10, 9 }, { 11, 9 }, { 12, 9 }, { 13, 9 }, { 14, 9 }, { 15, 9 }, { 16, 9 }, { 17, 9 }, { 4, 6 }, { 5, 6 }, { 6, 6 }, { 7, 6 }, { 8, 6 }, { 9, 6 }, { 10, 6 }, { 11, 6 }, { 12, 6 }, { 13, 6 }, { 14, 6 }, { 15, 6 }, { 16, 6 }, { 17, 6 }, { 18, 6 },
  { 1, 7 }, { 2, 7 }, { 3, 7 }, { 4, 10 }, { 5, 10 }, { 6, 10 }, { 7, 10 }, { 8, 10 }, { 9, 10 }, { 10, 10 }, { 11, 10 }, { 12, 10 }, { 13, 10 }, { 14, 10 }, { 15, 10 }, { 16, 10 }, { 17, 10 }, { 4, 7 }, { 5, 7 }, { 6, 7 }, { 7, 7 }, { 8, 7 }, { 9, 7 }, { 10, 7 }, { 11, 7 }, { 12, 7 }, { 13, 7 }, { 14, 7 }, { 15, 7 }, { 16, 7 }, { 17, 7 }, { 18, 7 }
};

static gboolean gdis_element_is_atomic_number(guint atomic_number, const guint *values, gsize count);
static GdisElementFamily gdis_element_guess_family(guint atomic_number, guint row, guint column);
static guint gdis_element_guess_period(guint atomic_number, guint row);
static void gdis_element_default_color(GdisElementFamily family, gdouble rgb[3]);
static void gdis_element_apply_common_override(GdisElementInfo *element);
static void gdis_elements_init_once(void);

void
gdis_element_normalize_symbol(const char *input, char output[4])
{
  g_return_if_fail(output != NULL);

  output[0] = '\0';
  output[1] = '\0';
  output[2] = '\0';
  output[3] = '\0';

  if (!input || !*input)
    return;

  output[0] = g_ascii_toupper(input[0]);
  if (input[1] && g_ascii_isalpha(input[1]))
    output[1] = g_ascii_tolower(input[1]);
  if (input[2] && g_ascii_isalpha(input[2]))
    output[2] = g_ascii_tolower(input[2]);
}

const GdisElementInfo *
gdis_element_lookup(const char *symbol)
{
  char normalized[4];

  gdis_elements_init_once();
  gdis_element_normalize_symbol(symbol, normalized);
  if (!normalized[0])
    return gdis_element_fallback();

  for (guint i = 1; i < G_N_ELEMENTS(gdis_elements); i++)
    {
      if (g_strcmp0(normalized, gdis_elements[i].symbol) == 0)
        return &gdis_elements[i];
    }

  return gdis_element_fallback();
}

const GdisElementInfo *
gdis_element_lookup_atomic_number(guint atomic_number)
{
  gdis_elements_init_once();
  if (atomic_number == 0 || atomic_number >= G_N_ELEMENTS(gdis_elements))
    return gdis_element_fallback();
  return &gdis_elements[atomic_number];
}

const GdisElementInfo *
gdis_element_fallback(void)
{
  gdis_elements_init_once();
  return &gdis_elements[6];
}

const char *
gdis_element_family_label(GdisElementFamily family)
{
  switch (family)
    {
    case GDIS_ELEMENT_FAMILY_NONMETAL:
      return "Nonmetal";
    case GDIS_ELEMENT_FAMILY_NOBLE_GAS:
      return "Noble gas";
    case GDIS_ELEMENT_FAMILY_ALKALI_METAL:
      return "Alkali metal";
    case GDIS_ELEMENT_FAMILY_ALKALINE_EARTH_METAL:
      return "Alkaline earth metal";
    case GDIS_ELEMENT_FAMILY_TRANSITION_METAL:
      return "Transition metal";
    case GDIS_ELEMENT_FAMILY_POST_TRANSITION_METAL:
      return "Post-transition metal";
    case GDIS_ELEMENT_FAMILY_METALLOID:
      return "Metalloid";
    case GDIS_ELEMENT_FAMILY_HALOGEN:
      return "Halogen";
    case GDIS_ELEMENT_FAMILY_LANTHANIDE:
      return "Lanthanide";
    case GDIS_ELEMENT_FAMILY_ACTINIDE:
      return "Actinide";
    case GDIS_ELEMENT_FAMILY_UNKNOWN:
    default:
      return "Unknown";
    }
}

guint
gdis_element_count(void)
{
  return 118;
}

static gboolean
gdis_element_is_atomic_number(guint atomic_number, const guint *values, gsize count)
{
  for (gsize i = 0; i < count; i++)
    {
      if (values[i] == atomic_number)
        return TRUE;
    }
  return FALSE;
}

static GdisElementFamily
gdis_element_guess_family(guint atomic_number, guint row, guint column)
{
  static const guint metalloids[] = {5u, 14u, 32u, 33u, 51u, 52u};
  static const guint nonmetals[] = {1u, 6u, 7u, 8u, 15u, 16u, 34u};
  static const guint post_transition[] = {13u, 31u, 49u, 50u, 81u, 82u, 83u, 84u, 113u, 114u, 115u, 116u};

  if (atomic_number >= 57u && atomic_number <= 71u)
    return GDIS_ELEMENT_FAMILY_LANTHANIDE;
  if (atomic_number >= 89u && atomic_number <= 103u)
    return GDIS_ELEMENT_FAMILY_ACTINIDE;
  if (gdis_element_is_atomic_number(atomic_number, nonmetals, G_N_ELEMENTS(nonmetals)))
    return GDIS_ELEMENT_FAMILY_NONMETAL;
  if (column == 18u)
    return GDIS_ELEMENT_FAMILY_NOBLE_GAS;
  if (column == 17u)
    return GDIS_ELEMENT_FAMILY_HALOGEN;
  if (column == 1u)
    return GDIS_ELEMENT_FAMILY_ALKALI_METAL;
  if (column == 2u)
    return GDIS_ELEMENT_FAMILY_ALKALINE_EARTH_METAL;
  if (gdis_element_is_atomic_number(atomic_number, metalloids, G_N_ELEMENTS(metalloids)))
    return GDIS_ELEMENT_FAMILY_METALLOID;
  if ((row >= 4u && row <= 7u && column >= 3u && column <= 12u) || atomic_number == 57u || atomic_number == 89u)
    return GDIS_ELEMENT_FAMILY_TRANSITION_METAL;
  if (gdis_element_is_atomic_number(atomic_number, post_transition, G_N_ELEMENTS(post_transition)))
    return GDIS_ELEMENT_FAMILY_POST_TRANSITION_METAL;
  if (column >= 13u && column <= 16u)
    return GDIS_ELEMENT_FAMILY_POST_TRANSITION_METAL;
  return GDIS_ELEMENT_FAMILY_UNKNOWN;
}

static guint
gdis_element_guess_period(guint atomic_number, guint row)
{
  if (atomic_number >= 57u && atomic_number <= 71u)
    return 6u;
  if (atomic_number >= 89u && atomic_number <= 103u)
    return 7u;
  return CLAMP(row, 1u, 7u);
}

static void
gdis_element_default_color(GdisElementFamily family, gdouble rgb[3])
{
  switch (family)
    {
    case GDIS_ELEMENT_FAMILY_NONMETAL:
      rgb[0] = 0.62;
      rgb[1] = 0.72;
      rgb[2] = 0.78;
      break;
    case GDIS_ELEMENT_FAMILY_NOBLE_GAS:
      rgb[0] = 0.44;
      rgb[1] = 0.84;
      rgb[2] = 0.92;
      break;
    case GDIS_ELEMENT_FAMILY_ALKALI_METAL:
      rgb[0] = 0.60;
      rgb[1] = 0.52;
      rgb[2] = 0.90;
      break;
    case GDIS_ELEMENT_FAMILY_ALKALINE_EARTH_METAL:
      rgb[0] = 0.42;
      rgb[1] = 0.72;
      rgb[2] = 0.76;
      break;
    case GDIS_ELEMENT_FAMILY_TRANSITION_METAL:
      rgb[0] = 0.54;
      rgb[1] = 0.68;
      rgb[2] = 0.82;
      break;
    case GDIS_ELEMENT_FAMILY_POST_TRANSITION_METAL:
      rgb[0] = 0.68;
      rgb[1] = 0.70;
      rgb[2] = 0.74;
      break;
    case GDIS_ELEMENT_FAMILY_METALLOID:
      rgb[0] = 0.84;
      rgb[1] = 0.70;
      rgb[2] = 0.36;
      break;
    case GDIS_ELEMENT_FAMILY_HALOGEN:
      rgb[0] = 0.32;
      rgb[1] = 0.80;
      rgb[2] = 0.44;
      break;
    case GDIS_ELEMENT_FAMILY_LANTHANIDE:
      rgb[0] = 0.74;
      rgb[1] = 0.54;
      rgb[2] = 0.88;
      break;
    case GDIS_ELEMENT_FAMILY_ACTINIDE:
      rgb[0] = 0.90;
      rgb[1] = 0.54;
      rgb[2] = 0.72;
      break;
    case GDIS_ELEMENT_FAMILY_UNKNOWN:
    default:
      rgb[0] = 0.64;
      rgb[1] = 0.74;
      rgb[2] = 0.82;
      break;
    }
}

static void
gdis_element_apply_common_override(GdisElementInfo *element)
{
  g_return_if_fail(element != NULL);

  if (g_strcmp0(element->symbol, "H") == 0)
    {
      element->covalent_radius = 0.31;
      element->vdw_radius = 1.20;
      element->draw_radius = 0.32;
      element->color_rgb[0] = 0.96;
      element->color_rgb[1] = 0.96;
      element->color_rgb[2] = 0.96;
      return;
    }
  if (g_strcmp0(element->symbol, "C") == 0)
    {
      element->covalent_radius = 0.76;
      element->vdw_radius = 1.70;
      element->draw_radius = 0.72;
      element->color_rgb[0] = 0.35;
      element->color_rgb[1] = 0.37;
      element->color_rgb[2] = 0.40;
      return;
    }
  if (g_strcmp0(element->symbol, "N") == 0)
    {
      element->covalent_radius = 0.71;
      element->vdw_radius = 1.55;
      element->draw_radius = 0.68;
      element->color_rgb[0] = 0.20;
      element->color_rgb[1] = 0.45;
      element->color_rgb[2] = 0.92;
      return;
    }
  if (g_strcmp0(element->symbol, "O") == 0)
    {
      element->covalent_radius = 0.66;
      element->vdw_radius = 1.52;
      element->draw_radius = 0.66;
      element->color_rgb[0] = 0.88;
      element->color_rgb[1] = 0.18;
      element->color_rgb[2] = 0.22;
      return;
    }
  if (g_strcmp0(element->symbol, "F") == 0)
    {
      element->covalent_radius = 0.57;
      element->vdw_radius = 1.47;
      element->draw_radius = 0.60;
      return;
    }
  if (g_strcmp0(element->symbol, "Na") == 0)
    {
      element->covalent_radius = 1.66;
      element->vdw_radius = 2.27;
      element->draw_radius = 1.02;
      return;
    }
  if (g_strcmp0(element->symbol, "Mg") == 0)
    {
      element->covalent_radius = 1.41;
      element->vdw_radius = 1.73;
      element->draw_radius = 0.96;
      return;
    }
  if (g_strcmp0(element->symbol, "Al") == 0)
    {
      element->covalent_radius = 1.21;
      element->vdw_radius = 1.84;
      element->draw_radius = 0.95;
      return;
    }
  if (g_strcmp0(element->symbol, "Si") == 0)
    {
      element->covalent_radius = 1.11;
      element->vdw_radius = 2.10;
      element->draw_radius = 0.92;
      return;
    }
  if (g_strcmp0(element->symbol, "P") == 0)
    {
      element->covalent_radius = 1.07;
      element->vdw_radius = 1.80;
      element->draw_radius = 1.06;
      element->color_rgb[0] = 0.95;
      element->color_rgb[1] = 0.56;
      element->color_rgb[2] = 0.20;
      return;
    }
  if (g_strcmp0(element->symbol, "S") == 0)
    {
      element->covalent_radius = 1.05;
      element->vdw_radius = 1.80;
      element->draw_radius = 1.02;
      element->color_rgb[0] = 0.93;
      element->color_rgb[1] = 0.79;
      element->color_rgb[2] = 0.24;
      return;
    }
  if (g_strcmp0(element->symbol, "Cl") == 0)
    {
      element->covalent_radius = 1.02;
      element->vdw_radius = 1.75;
      element->draw_radius = 0.99;
      element->color_rgb[0] = 0.19;
      element->color_rgb[1] = 0.74;
      element->color_rgb[2] = 0.30;
      return;
    }
  if (g_strcmp0(element->symbol, "Ca") == 0)
    {
      element->covalent_radius = 1.76;
      element->vdw_radius = 2.31;
      element->draw_radius = 1.08;
      return;
    }
  if (g_strcmp0(element->symbol, "Br") == 0)
    {
      element->covalent_radius = 1.20;
      element->vdw_radius = 1.85;
      element->draw_radius = 1.05;
      return;
    }
  if (g_strcmp0(element->symbol, "I") == 0)
    {
      element->covalent_radius = 1.39;
      element->vdw_radius = 1.98;
      element->draw_radius = 1.10;
      return;
    }
  if (g_strcmp0(element->symbol, "Au") == 0)
    {
      element->color_rgb[0] = 0.86;
      element->color_rgb[1] = 0.74;
      element->color_rgb[2] = 0.26;
      return;
    }
}

static void
gdis_elements_init_once(void)
{
  if (gdis_elements_initialized)
    return;

  for (guint atomic_number = 1; atomic_number < G_N_ELEMENTS(gdis_elements); atomic_number++)
    {
      GdisElementInfo *element;
      guint row;
      guint column;
      guint period;
      gdouble covalent_radius;

      element = &gdis_elements[atomic_number];
      row = gdis_element_positions[atomic_number][1];
      column = gdis_element_positions[atomic_number][0];

      element->atomic_number = atomic_number;
      element->symbol = gdis_element_seeds[atomic_number].symbol;
      element->name = gdis_element_seeds[atomic_number].name;
      element->row = row;
      element->column = column;
      element->family = gdis_element_guess_family(atomic_number, row, column);

      period = gdis_element_guess_period(atomic_number, row);
      covalent_radius = 0.35 + (gdouble) period * 0.18;
      switch (element->family)
        {
        case GDIS_ELEMENT_FAMILY_ALKALI_METAL:
          covalent_radius += 0.38;
          break;
        case GDIS_ELEMENT_FAMILY_ALKALINE_EARTH_METAL:
          covalent_radius += 0.22;
          break;
        case GDIS_ELEMENT_FAMILY_TRANSITION_METAL:
          covalent_radius += 0.02;
          break;
        case GDIS_ELEMENT_FAMILY_POST_TRANSITION_METAL:
          covalent_radius += 0.10;
          break;
        case GDIS_ELEMENT_FAMILY_METALLOID:
          covalent_radius += 0.05;
          break;
        case GDIS_ELEMENT_FAMILY_HALOGEN:
          covalent_radius -= 0.10;
          break;
        case GDIS_ELEMENT_FAMILY_NOBLE_GAS:
          covalent_radius += 0.04;
          break;
        case GDIS_ELEMENT_FAMILY_NONMETAL:
          covalent_radius -= 0.06;
          break;
        case GDIS_ELEMENT_FAMILY_LANTHANIDE:
        case GDIS_ELEMENT_FAMILY_ACTINIDE:
          covalent_radius += 0.12;
          break;
        case GDIS_ELEMENT_FAMILY_UNKNOWN:
        default:
          break;
        }

      element->covalent_radius = MAX(covalent_radius, 0.28);
      element->vdw_radius = MAX(element->covalent_radius + 0.90, 1.20);
      element->draw_radius = MAX(element->covalent_radius * 0.82, 0.30);
      gdis_element_default_color(element->family, element->color_rgb);
      gdis_element_apply_common_override(element);
    }

  gdis_elements_initialized = TRUE;
}
