import requests
import pandas as pd

def read_ko_ids_from_file(ko_ids_file):
    """
    Read KO IDs from a file with two columns (gene ID and KO ID) and return the KO IDs as a list.

    Args:
        ko_ids_file (str): The path to the file containing gene IDs and KO IDs.

    Returns:
        list: A list of KO IDs.
    """
    ko_ids = []
    with open(ko_ids_file, "r") as file:
        for line in file:
            # Split the line into two columns: gene_id and ko_id
            columns = line.split()
            if len(columns) > 1:  # Check if there is a KO ID in the second column
                ko_id = columns[1].strip()  # Get the KO ID (second column)
                if ko_id:  # Only add if the KO ID is not empty
                    ko_ids.append(ko_id)
    return ko_ids

def fetch_kegg_info(ko_id):
    """
    Fetch detailed information about a KEGG Orthology (KO) ID using the KEGG REST API.

    Args:
        ko_id (str): The KEGG Orthology ID (e.g., "K01077").

    Returns:
        str: The detailed information as plain text, or None if an error occurred.
    """
    url = f"http://rest.kegg.jp/get/ko:{ko_id}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            print(f"Error: Unable to fetch data for KO ID {ko_id}. Status code: {response.status_code}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Error: An exception occurred while fetching data for KO ID {ko_id}: {e}")
        return None

def parse_kegg_info(raw_data):
    """
    Parse KEGG data into a structured dictionary.

    Args:
        raw_data (str): The raw KEGG response text.

    Returns:
        dict: A dictionary containing parsed KEGG data.
    """
    data = {}
    lines = raw_data.splitlines()
    current_key = None
    for line in lines:
        
        if line.startswith("GENES"):  # Stop processing when 'GENES' is encountered
            break
        
        if line.startswith(" "):  # Continuation of the previous field
            if current_key:
                data[current_key] += f" {line.strip()}"
        else:
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                current_key = parts[0]
                data[current_key] = parts[1]
            elif len(parts) == 1:  # A field with no value
                current_key = parts[0]
                data[current_key] = ""
    return data

def kegg_to_dataframe(ko_ids):
    """
    Fetch and parse KEGG data for multiple KO IDs and create a pandas DataFrame.

    Args:
        ko_ids (list): A list of KEGG Orthology IDs (e.g., ["K01077", "K00001"]).

    Returns:
        pd.DataFrame: A DataFrame containing parsed KEGG data.
    """
    all_data = []
    for ko_id in ko_ids:
        print(f"Processing KO ID: {ko_id}")
        raw_data = fetch_kegg_info(ko_id)
        if raw_data:
            parsed_data = parse_kegg_info(raw_data)
            parsed_data["KO_ID"] = ko_id  # Add the KO ID to the data
            all_data.append(parsed_data)
    return pd.DataFrame(all_data)

if __name__ == "__main__":
    # Path to the file containing the gene IDs and KO IDs
    ko_ids_file = "data/user_ko.txt"

    # Read KO IDs from the file
    ko_ids = read_ko_ids_from_file(ko_ids_file)

    # Fetch data and create the DataFrame
    df = kegg_to_dataframe(ko_ids)

    # Save the DataFrame to a CSV file
    output_file = "kegg_data.csv"
    df.to_csv(output_file, index=False)

    print(f"KEGG data saved to {output_file}")