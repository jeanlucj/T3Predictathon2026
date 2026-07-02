import pandas as pd
import re
from datetime import datetime

# Input and output
in_path = "data/processed/historical_env_metadata.csv"
out_path = "data/processed/historical_env_metadata_completed.csv"

df = pd.read_csv(in_path)


# 1. Curated lat/lon lookup table
latlon = {
    "Stillwater, OK": (36.1156, -97.0584),
    "Brookings, SD": (44.3114, -96.7983),
    "Manhattan, KS": (39.1836, -96.5717),
    "Garden City, KS": (37.9717, -100.8727),
    "Hays, KS": (38.8790, -99.3268),
    "Fort Collins, CO": (40.5853, -105.0844),
    "Akron, CO": (40.1603, -103.2144),
    "Julesburg, CO": (40.9889, -102.2630),
    "Burlington, CO": (39.3061, -102.2696),
    "Bushland, TX": (35.1920, -102.0649),
    "Prosper, TX": (33.2360, -96.8011),
    "Chillicothe, TX": (34.2562, -99.5140),
    "Clovis, NM": (34.4048, -103.2052),
    "Farmington, NM": (36.7281, -108.2187),
    "Sidney, NE": (41.1420, -102.9770),
    "Alliance, NE": (42.1013, -102.8727),
    "North Platte, NE": (41.1403, -100.7601),
    "Lincoln, NE": (40.8136, -96.7026),
    "Colby, KS": (39.3950, -101.0524),
    "Lahoma, OK": (36.3870, -98.1287),
    "Winfield, KS": (37.2398, -96.9953),
    "Fort Cobb, OK": (35.0987, -98.4353),
    "Archer, WY": (41.1400, -104.6400),
    "Bozeman, MT": (45.6770, -111.0429),
    "Moccasin, MT": (47.0591, -109.9507),
    "Minot, ND": (48.2325, -101.2963),
    "Carrington, ND": (47.4492, -99.1268),
    "Hettinger, ND": (46.0017, -102.6388),
    "Dickinson, ND": (46.8792, -102.7896),
    "Langdon, ND": (48.7600, -98.3687),
    "Prosper, ND": (46.9747, -97.1231),
    "Roseau, MN": (48.8469, -95.7603),
    "Crookston, MN": (47.7744, -96.6081),
    "Lamberton, MN": (44.2311, -95.2639),
    "Kimball, MN": (45.3130, -94.2972),
    "St. Paul, MN": (44.9537, -93.0900),
    "Ithaca, NY": (42.4430, -76.5019),
    "Urbana, IL": (40.1106, -88.2073),
}


# 2. Classify winter vs spring wheat
def classify_crop(studyName, locationName):
    name = studyName.upper()

    if name.startswith(("SRPN", "ARS-SRPN", "NRPN", "HWWPANEL", "CORNELLMASTER", "SDK-WHEAT")):
        return "winter"
    if name.startswith(("NDK-WHEAT", "UMN-WHEAT")):
        return "spring"

    if ", SD" in locationName or ", ND" in locationName or ", MN" in locationName:
        return "spring"
    if ", OK" in locationName or ", KS" in locationName or ", TX" in locationName or ", CO" in locationName or ", NE" in locationName or ", WY" in locationName:
        return "winter"

    return "winter"


# 3. Extract study year (with 2-digit → 20XX rule)
def extract_year(studyName):
    # First try 4-digit year
    m4 = re.search(r"(\d{4})", studyName)
    if m4:
        return int(m4.group(1))

    # Then try 2-digit year
    m2 = re.search(r"(\d{2})(?!\d)", studyName)
    if m2:
        return 2000 + int(m2.group(1))

    return None

df["study_year"] = df["studyName"].apply(extract_year)


# 4. Infer planting/harvest dates
def infer_dates(row):
    year = row["study_year"]
    if year is None:
        return None, None

    crop = row["crop_type"]

    if crop == "winter":
        return datetime(year - 1, 10, 1).date(), datetime(year, 7, 1).date()
    else:
        return datetime(year, 5, 1).date(), datetime(year, 8, 15).date()

df["crop_type"] = df.apply(lambda r: classify_crop(r["studyName"], r["locationName"]), axis=1)
df[["plantingDate", "harvestDate"]] = df.apply(lambda r: pd.Series(infer_dates(r)), axis=1)


# 5. Assign lat/lon
df["latitude"] = df["locationName"].apply(lambda x: latlon.get(x, (None, None))[0])
df["longitude"] = df["locationName"].apply(lambda x: latlon.get(x, (None, None))[1])


# 6. Save
df.to_csv(out_path, index=False)
print(f"✓ Wrote {out_path} with {len(df)} rows.")