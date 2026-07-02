import pandas as pd
from pathlib import Path

INPUT = Path("data/raw/metadata.csv")
OUTPUT = Path("data/processed/metadata_predictathon_clean.csv")

def main():
    df = pd.read_csv(INPUT)

    # Keep only the fields needed for weather extraction
    df = df[["studyName", "locationName", "latitude", "longitude", "plantingDate", "harvestDate"]]

    # Convert planting/harvest dates to YYYYMMDD
    df["plantingDate"] = pd.to_datetime(df["plantingDate"]).dt.strftime("%Y%m%d")
    df["harvestDate"] = pd.to_datetime(df["harvestDate"]).dt.strftime("%Y%m%d")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT, index=False)

    print(f"Wrote cleaned metadata to {OUTPUT}")

if __name__ == "__main__":
    main()