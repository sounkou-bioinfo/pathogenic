import lxml.etree as ET
import csv
import os
from tqdm import tqdm

# Define TSV column headers
columns = [
    'MeasureSetID', 'MeasureSetAcc', 'MeasureSetVersion', 'AlleleID', 'MeasureType', 'CanonicalSPDI',
    'GRCh37_Chromosome', 'GRCh37_PositionVCF', 'GRCh37_ReferenceAlleleVCF', 'GRCh37_AlternateAlleleVCF', 'GRCh37_Start', 'GRCh37_Stop',
    'GRCh38_Chromosome', 'GRCh38_PositionVCF', 'GRCh38_ReferenceAlleleVCF', 'GRCh38_AlternateAlleleVCF', 'GRCh38_Start', 'GRCh38_Stop',
    'MolecularConsequence', 'FunctionalConsequence',
    'ModeOfInheritance', 'PreferredValues', 'Citations', 'Comments', 'FamilyData',
    'RecordStatus', 'ClinicalSignificance', 'ReviewStatus', 'Description', 'DateLastEvaluated'
]

def extract_variant_data(clinvar_set):
    """Extract all required fields from a ClinVarSet element into a dictionary."""
    ref_assertion = clinvar_set.find('ReferenceClinVarAssertion')
    if ref_assertion is None:
        print("Warning: No ReferenceClinVarAssertion found in ClinVarSet")
        return None

    # Basic fields from ReferenceClinVarAssertion
    record_status = ref_assertion.findtext('RecordStatus', default='')

    # Classification details
    germline_classification = ref_assertion.find('Classifications/GermlineClassification')
    if germline_classification is not None:
        review_status = germline_classification.findtext('ReviewStatus', default='')
        description_elem = germline_classification.find('Description')
        description = description_elem.text if description_elem is not None else ''
        date_last_evaluated = description_elem.get('DateLastEvaluated', '') if description_elem is not None else ''
        clinical_significance = description
    else:
        review_status, description, date_last_evaluated, clinical_significance = '', '', '', ''

    # MeasureSet details
    measure_set = ref_assertion.find('MeasureSet[@Type="Variant"]')
    if measure_set is None:
        print("Warning: No MeasureSet with Type='Variant' found")
        return None

    measure_set_id = measure_set.get('ID', '')
    measure_set_acc = measure_set.get('Acc', '')
    measure_set_version = measure_set.get('Version', '')

    measure = measure_set.find('Measure')
    if measure is None:
        print("Warning: No Measure found in MeasureSet")
        return None

    allele_id = measure.get('ID', '')
    measure_type = measure.get('Type', '')
    canonical_spdi = measure.findtext('CanonicalSPDI', default='')

    # Sequence Locations for GRCh37 and GRCh38
    sequence_locations = {}
    for seq_loc in measure.findall('SequenceLocation'):
        assembly = seq_loc.get('Assembly')
        if assembly in ['GRCh37', 'GRCh38']:
            sequence_locations[assembly] = {
                'Chromosome': seq_loc.get('Chr', ''),
                'PositionVCF': seq_loc.get('positionVCF', ''),
                'ReferenceAlleleVCF': seq_loc.get('referenceAlleleVCF', ''),
                'AlternateAlleleVCF': seq_loc.get('alternateAlleleVCF', ''),
                'Start': seq_loc.get('start', ''),
                'Stop': seq_loc.get('stop', '')
            }

    # Attributes (MolecularConsequence, FunctionalConsequence)
    attributes = {}
    for attr in measure.findall('AttributeSet/Attribute'):
        attr_type = attr.get('Type')
        if attr_type in ['MolecularConsequence', 'FunctionalConsequence']:
            attributes[attr_type] = attr.text or ''

    # Lists joined with semicolons
    mode_of_inheritance = [attr.text for attr in ref_assertion.findall('AttributeSet/Attribute[@Type="ModeOfInheritance"]') if attr.text]
    preferred_values = [elem.text for elem in measure.findall('Name/ElementValue[@Type="Preferred"]') if elem.text]
    citations = [cit.findtext('ID') for cit in clinvar_set.findall('.//Citation') if cit.find('ID') is not None]
    comments = [com.text for com in clinvar_set.findall('.//Comment') if com.text]
    family_data = [fd.get('NumFamilies', '') for fd in clinvar_set.findall('.//FamilyData') if fd.get('NumFamilies')]

    # Construct the row dictionary
    row = {
        'MeasureSetID': measure_set_id,
        'MeasureSetAcc': measure_set_acc,
        'MeasureSetVersion': measure_set_version,
        'AlleleID': allele_id,
        'MeasureType': measure_type,
        'CanonicalSPDI': canonical_spdi,
        'GRCh37_Chromosome': sequence_locations.get('GRCh37', {}).get('Chromosome', ''),
        'GRCh37_PositionVCF': sequence_locations.get('GRCh37', {}).get('PositionVCF', ''),
        'GRCh37_ReferenceAlleleVCF': sequence_locations.get('GRCh37', {}).get('ReferenceAlleleVCF', ''),
        'GRCh37_AlternateAlleleVCF': sequence_locations.get('GRCh37', {}).get('AlternateAlleleVCF', ''),
        'GRCh37_Start': sequence_locations.get('GRCh37', {}).get('Start', ''),
        'GRCh37_Stop': sequence_locations.get('GRCh37', {}).get('Stop', ''),
        'GRCh38_Chromosome': sequence_locations.get('GRCh38', {}).get('Chromosome', ''),
        'GRCh38_PositionVCF': sequence_locations.get('GRCh38', {}).get('PositionVCF', ''),
        'GRCh38_ReferenceAlleleVCF': sequence_locations.get('GRCh38', {}).get('ReferenceAlleleVCF', ''),
        'GRCh38_AlternateAlleleVCF': sequence_locations.get('GRCh38', {}).get('AlternateAlleleVCF', ''),
        'GRCh38_Start': sequence_locations.get('GRCh38', {}).get('Start', ''),
        'GRCh38_Stop': sequence_locations.get('GRCh38', {}).get('Stop', ''),
        'MolecularConsequence': attributes.get('MolecularConsequence', ''),
        'FunctionalConsequence': attributes.get('FunctionalConsequence', ''),
        'ModeOfInheritance': ';'.join(mode_of_inheritance),
        'PreferredValues': ';'.join(preferred_values),
        'Citations': ';'.join(citations),
        'Comments': ';'.join(comments),
        'FamilyData': ';'.join(family_data),
        'RecordStatus': record_status,
        'ClinicalSignificance': clinical_significance,
        'ReviewStatus': review_status,
        'Description': description,
        'DateLastEvaluated': date_last_evaluated
    }
    return row

def parse_variant_xml(xml_file, tsv_file):
    """Parse the XML file and write variants to a TSV with progress tracking."""
    print(f"Starting to parse XML file: {xml_file}")
    print(f"File size: {os.path.getsize(xml_file) / (1024**3):.2f} GB")

    with open(xml_file, 'rb') as f, open(tsv_file, 'w', newline='') as tsvfile:
        # Set up TSV writer
        writer = csv.DictWriter(tsvfile, fieldnames=columns, delimiter='\t')
        writer.writeheader()

        # Get total file size for progress tracking
        total_size = os.path.getsize(xml_file)
        context = ET.iterparse(f, events=('end',), tag='ClinVarSet')
        count = 0

        # Wrap the parsing loop with tqdm for a progress bar
        with tqdm(total=total_size, unit='B', unit_scale=True, desc="Parsing ClinVar XML") as pbar:
            for event, elem in context:
                if elem.tag == 'ClinVarSet':
                    count += 1
                    if count % 10000 == 0:
                        print(f"Processed {count} ClinVarSet elements")
                    row = extract_variant_data(elem)
                    if row:
                        writer.writerow(row)
                    # Clear memory
                    elem.clear()
                    while elem.getprevious() is not None:
                        del elem.getparent()[0]
                # Update progress bar based on current file position
                current_pos = f.tell()
                pbar.n = current_pos  # Update to current position
                pbar.refresh()

        print(f"\nFinished parsing. Total ClinVarSet elements processed: {count}")

def main():
    xml_file = 'ClinVarRCVRelease_00-latest.xml'
    tsv_file = 'output.tsv'
    try:
        parse_variant_xml(xml_file, tsv_file)
        print(f"Parsing completed successfully. Output written to {tsv_file}")
    except FileNotFoundError:
        print(f"Error: The file '{xml_file}' was not found. Please check the file path.")
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}. The file may be corrupted or not valid XML.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    main()
