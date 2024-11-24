from pyteomics import mgf
import re
import os

def mgf_extractScans(mgf_file_path, target_scan):
    """
    Extract a specific scan from an MGF file.

    Parameters:
    - mgf_file_path: str, path to the MGF file.
    - target_scan: int, the scan number to extract.

    Returns:
    - spectrum: dict, containing information about the extracted scan.
    """

    # mgf.read iterates through scans in the MGF file
    with mgf.read(mgf_file_path) as mgf_reader:
        for spectrum in mgf_reader:
            try:
                scan_number = int(spectrum['params']['scans'])
                if scan_number == int(target_scan):
                    # Found the target scan, store its information and break the loop
                    return spectrum
            #sometimes, 'SCANS=xxx' line is missing. check for the scan in the 'TITLE' line.
            except KeyError:
                pass

            title = spectrum['params'].get('title', '')
            try:
                #Different tools export mgf title format in a different way. Try to capture the scan in any of the cases.
                #scan_number = int(re.split('scan=|_?scan:|Scan |Index: ', title)[1].rstrip('"'))
                pattern = r"(scan=|scan:|Scan|Index:)\s*(\d+)"
                match = re.search(pattern, title)
                scan_number = int(match.group(2))

                if scan_number == int(target_scan):
                    # Found the target scan, store its information and break the loop
                    spectrum['params']['scans'] = scan_number
                    return spectrum

            except (IndexError, ValueError):
                pass

    # If the loop completes without finding the target scan, print a message
    print(f"Scan {target_scan} is not found in {os.path.basename(mgf_file_path)}.")
    return None  # Return None if the scan is not found