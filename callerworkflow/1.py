import os
import pandas as pd

# Directory where the CSV files are stored
directory_path = '/data/home/std_12/ShiHe/sarcomabam_csv/csv_files'

# List to hold dataframes
df_list = []

# Iterate through all CSV files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".csv"):
        file_path = os.path.join(directory_path, filename)
        df = pd.read_csv(file_path)
        df_list.append(df)

# Combine all dataframes into one
combined_df = pd.concat(df_list, ignore_index=True)

# Save the combined dataframe to a new CSV file
output_file_path = '/data/home/std_12/ShiHe/sarcomabam_csv/combined_files.csv'
combined_df.to_csv(output_file_path, index=False)

# Output the path of the saved file
print(f"Combined CSV file saved at: {output_file_path}")
