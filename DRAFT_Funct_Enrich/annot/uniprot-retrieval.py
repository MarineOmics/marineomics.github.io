import zlib
import gzip
import requests
from requests.adapters import HTTPAdapter, Retry
import sys
import shutil

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


def process_accessions(accessions):
    accession_batches = [accessions[i:i+500] for i in range(0, len(accessions), 500)]
    all_lines = []

    for accession_batch in accession_batches:
        accession_query = '%29%20OR%20%28accession%3A'.join(accession_batch)
        url = f"https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id&format=tsv&query=%28%28accession%3A{accession_query}%29%29&size=500"

        progress = 0
        lines = []
        for batch, total in get_batch(url):
            # Decompress each batch as we want to extract the header
            decompressed = zlib.decompress(batch.content, 16 + zlib.MAX_WBITS)
            batch_lines = [line for line in decompressed.decode("utf-8").split("\n") if line]
            if not progress:
                # First line so print TSV header
                lines = [batch_lines[0]]
            lines += batch_lines[1:]
            progress = len(lines) - 1
            print(f"{progress} / {total}")

        all_lines.extend(lines)

    return all_lines


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python uniprot-retrieval.py <input_file>")
        sys.exit(1)

    accession_file = sys.argv[1]
    with open(accession_file, 'r') as f:
        accessions = f.read().splitlines()

    retrieved_data = process_accessions(accessions)

    # Write to a temporary gzip file
    temp_filename = "uniprot-retrieval-temp.tsv.gz"
    with gzip.open(temp_filename, "wt", encoding="utf-8") as f:
        f.write('\n'.join(retrieved_data))

    # Merge the temporary file with the existing output file (if it exists)
    try:
        with gzip.open("uniprot-retrieval.tsv.gz", "rb") as f_existing, open(temp_filename, "rb") as f_temp:
            with gzip.open("uniprot-retrieval-merged.tsv.gz", "wb") as f_merged:
                shutil.copyfileobj(f_existing, f_merged)
                shutil.copyfileobj(f_temp, f_merged)
        
        # Replace the original output file with the merged file
        shutil.move("uniprot-retrieval-merged.tsv.gz", "uniprot-retrieval.tsv.gz")
    except FileNotFoundError:
        # If the existing output file doesn't exist, rename the temporary file
        shutil.move(temp_filename, "uniprot-retrieval.tsv.gz")
