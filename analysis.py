import pandas as pd
import yaml
def load_data(path="egfr_inhibitors.csv"):
return pd.read_csv(path)
def main():
df = load_data()
print("Skeleton loaded successfully.")
print(df.head())
if_name_ == "_main_":
main()
