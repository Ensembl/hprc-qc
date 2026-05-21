#!/usr/bin/env python3
"""
Fetch NCBI assembly report for a given assembly accession.
"""
import re
import subprocess
import sys
from html.parser import HTMLParser


class NCBIDirectoryParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.directories = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            for attr, value in attrs:
                if attr == 'href' and value.endswith('/') and not value.startswith('?'):
                    self.directories.append(value.rstrip('/'))


def get_ncbi_ftp_path(accession):
    """Construct NCBI FTP path for assembly accession."""
    parts = accession.split('_')
    prefix = parts[0]  # GCA or GCF
    num = parts[1].split('.')[0]

    # Split into groups of 3 digits: e.g., 041900255 -> 041/900/255
    group1 = num[0:3]
    group2 = num[3:6]
    group3 = num[6:9]

    base_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{group1}/{group2}/{group3}/"

    # List directory to find the full assembly name
    result = subprocess.run(
        ['curl', '-s', base_url],
        capture_output=True,
        text=True,
        check=True
    )

    parser = NCBIDirectoryParser()
    parser.feed(result.stdout)

    # Find directory matching our accession
    matching_dirs = [d for d in parser.directories if d.startswith(accession)]

    if not matching_dirs:
        raise ValueError(f"No assembly directory found for {accession}")

    if len(matching_dirs) > 1:
        print(f"Warning: Multiple directories found for {accession}, using first one", file=sys.stderr)

    assembly_dir = matching_dirs[0]
    return f"{base_url}{assembly_dir}/{assembly_dir}_assembly_report.txt"


def main():
    if len(sys.argv) != 3:
        print("Usage: fetch_ncbi_assembly_report.py <accession> <output_file>", file=sys.stderr)
        sys.exit(1)

    accession = sys.argv[1]
    output_file = sys.argv[2]

    try:
        report_url = get_ncbi_ftp_path(accession)
        print(f"Downloading assembly report from: {report_url}", file=sys.stderr)

        # Download the assembly report
        result = subprocess.run(
            ['curl', '-f', '-s', report_url],
            capture_output=True,
            text=True,
            check=True
        )

        # Write to output file
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        print(f"Successfully downloaded assembly report for {accession}", file=sys.stderr)

    except subprocess.CalledProcessError as e:
        print(f"Error downloading assembly report: {e}", file=sys.stderr)
        print(f"Creating empty placeholder file", file=sys.stderr)
        # Create empty file so pipeline doesn't fail
        with open(output_file, 'w') as f:
            f.write(f"# Assembly report not found for {accession}\n")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        # Create empty file so pipeline doesn't fail
        with open(output_file, 'w') as f:
            f.write(f"# Assembly report not found for {accession}\n")


if __name__ == '__main__':
    main()
