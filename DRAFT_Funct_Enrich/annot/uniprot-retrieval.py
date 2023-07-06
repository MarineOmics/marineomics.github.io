import re
import zlib
import gzip
import requests
from requests.adapters import HTTPAdapter, Retry
import sys


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


def main(accession_file):
    with open(accession_file, 'r') as f:
        accessions = f.read().splitlines()

    accession_query = '%29%20OR%20%28accession%3A'.join(accessions)
    url = f"https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_p%2Cgo_c%2Cgo_f%2Cgo%2Cgo_id&format=tsv&query=%28%28accession%3A{accession_query}%29%29&size=500"

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

    # Save lines to a gzip file
    with gzip.open("uniprot-retrieval.tsv.gz", "wt", encoding="utf-8") as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    main(sys.argv[1])
