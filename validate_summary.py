import csv
import pandas as pd
from datetime import datetime

# Define the column names based on the TSV structure
columns = [
    'MeasureSetID', 'MeasureSetAcc', 'MeasureSetVersion', 'AlleleID', 'MeasureType', 'CanonicalSPDI',
    'GRCh37_Chromosome', 'GRCh37_PositionVCF', 'GRCh37_ReferenceAlleleVCF', 'GRCh37_AlternateAlleleVCF', 'GRCh37_Start', 'GRCh37_Stop',
    'GRCh38_Chromosome', 'GRCh38_PositionVCF', 'GRCh38_ReferenceAlleleVCF', 'GRCh38_AlternateAlleleVCF', 'GRCh38_Start', 'GRCh38_Stop',
    'MolecularConsequence', 'FunctionalConsequence',
    'ModeOfInheritance', 'PreferredValues', 'Citations', 'Comments', 'FamilyData',
    'RecordStatus', 'ClinicalSignificance', 'ReviewStatus', 'Description', 'DateLastEvaluated'
]

def check_date_format(date_str):
    """Check if a date is in the correct format (YYYY-MM-DD)."""
    try:
        datetime.strptime(date_str, '%Y-%m-%d')
        return True, ""  # Valid date
    except ValueError:
        return False, "Invalid date format (expected YYYY-MM-DD)"  # Invalid date format

def validate_missing_and_invalid(file_path):
    """Validate the first 50,000 rows for missing values and invalid data."""
    # Read the first 50,000 rows into a pandas DataFrame
    with open(file_path, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = [row for _, row in zip(range(50000), reader)]  # Read up to 50,000 rows
    df = pd.DataFrame(rows)

    # Check for missing values
    missing_data = df.isnull().sum()
    print("Missing Data Summary:")
    print(missing_data[missing_data > 0])  # Print columns with missing data
    print("\n")

    # Check for invalid data
    invalid_data = {}

    def validate_numeric(column_name, value):
        """Check if the value is numeric."""
        try:
            float(value)
            return True, ""
        except ValueError:
            return False, f"Invalid numeric value in {column_name}"

    def validate_position(column_name, value):
        """Check if the position fields are valid numbers."""
        if value.isdigit():
            return True, ""
        else:
            return False, f"Invalid position value in {column_name}, expected numeric."

    # Prepare the invalid data results for specific columns
    for column in ['GRCh37_PositionVCF', 'GRCh38_PositionVCF', 'GRCh37_Start', 'GRCh37_Stop', 'GRCh38_Start', 'GRCh38_Stop']:
        invalid_data[column] = {'invalid': [], 'valid': []}
        for value in df[column]:
            is_valid, reason = validate_position(column, value)
            if not is_valid:
                invalid_data[column]['invalid'].append((value, reason))
            else:
                invalid_data[column]['valid'].append(value)

    # Check invalid date format for 'DateLastEvaluated'
    invalid_data['DateLastEvaluated'] = {'invalid': [], 'valid': []}
    for value in df['DateLastEvaluated']:
        is_valid, reason = check_date_format(value)
        if not is_valid:
            invalid_data['DateLastEvaluated']['invalid'].append((value, reason))
        else:
            invalid_data['DateLastEvaluated']['valid'].append(value)

    # Print invalid data summary with reasons and valid examples
    print("Invalid and Valid Data Summary:")

    for column, data in invalid_data.items():
        total_checked = len(df[column])  # Total number of entries checked
        total_invalid = len(data['invalid'])  # Number of invalid entries
        total_valid = len(data['valid'])  # Number of valid entries

        if total_invalid > 0 or total_valid > 0:
            print(f"\nColumn: {column}")

            # Show up to 10 invalid examples
            print(f"Invalid examples (up to 10):")
            for i, (value, reason) in enumerate(data['invalid'][:10]):
                print(f"Invalid value '{value}': {reason}")

            # Show up to 10 valid examples
            print(f"\nValid examples (up to 10):")
            for i, value in enumerate(data['valid'][:10]):
                print(f"Valid value '{value}'")

            # Calculate and print the percentage of invalid entries
            invalid_percentage = (total_invalid / total_checked) * 100 if total_checked > 0 else 0
            print(f"Invalid percentage: {invalid_percentage:.2f}%")
    
    print("\nValidation complete.")

def main():
    file_path = 'output.tsv'
    try:
        validate_missing_and_invalid(file_path)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found. Please check the file path.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    main()
