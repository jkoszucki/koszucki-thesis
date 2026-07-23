#!/usr/bin/env python3

"""
Script to scrape K-PAM K-antigen pages and download all PDB files linked on them.

Usage:
  python download_k_antigen_pdbs.py
"""

import os
import requests
from bs4 import BeautifulSoup
import re

BASE_URL = 'https://www.iith.ac.in/K-PAM/'
ANTIGEN_INDEX = 'k_antigen.html'
OUTPUT_DIR = 'pdb_files'

def get_kk_links(base_url: str, antigen_index: str) -> list:
    """
    Fetch links for all KK pages from the main antigen index.
    """
    resp = requests.get(base_url + antigen_index)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, 'html.parser')

    # Collect all anchors whose text starts with "KK" and href ends with ".html"
    kk_anchors = soup.find_all('a', href=True, text=re.compile(r'^KK', re.IGNORECASE))
    return [a['href'] for a in kk_anchors if a['href'].lower().endswith('.html')]


def extract_pdb_links(soup: BeautifulSoup, base_url: str) -> list:
    """
    From a KK page soup, find all <a href="*.pdb"> links.
    Returns full URLs.
    """
    pdb_links = []
    for a in soup.find_all('a', href=True):
        href = a['href'].strip()
        if href.lower().endswith('.pdb'):
            # make absolute
            if href.startswith('http'):
                full = href
            else:
                full = base_url + href.lstrip('/')
            pdb_links.append(full)
    return pdb_links


def download_file(url: str, out_dir: str):
    """
    Download a file from `url` into `out_dir`, preserving filename.
    """
    local_filename = os.path.join(out_dir, os.path.basename(url))
    # skip if already downloaded
    if os.path.exists(local_filename):
        print(f"  → Skipping (already exists): {local_filename}")
        return

    print(f"  ↓ Downloading {url} ...")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f"    saved to {local_filename}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    kk_pages = get_kk_links(BASE_URL, ANTIGEN_INDEX)
    kk_pages = sorted(kk_pages, key=lambda x: int(re.search(r'\d+', x).group()))

    print(f"Found {len(kk_pages)} KK pages. Starting download of PDB files...\n")

    for href in kk_pages:
        kk_id = href.replace('.html', '').upper()
        page_url = BASE_URL + href
        print(f"Processing {kk_id} ({page_url})")

        try:
            resp = requests.get(page_url)
            resp.raise_for_status()
            soup = BeautifulSoup(resp.text, 'html.parser')

            pdb_urls = extract_pdb_links(soup, BASE_URL)
            if not pdb_urls:
                print("  ! No PDB links found on this page.")
                continue

            for pdb_url in pdb_urls:
                download_file(pdb_url, OUTPUT_DIR)

        except requests.HTTPError as e:
            print(f"  ! Failed to fetch {kk_id}: {e}")
        except Exception as e:
            print(f"  ! Unexpected error on {kk_id}: {e}")

    print("\nDone. All available PDBs are in the ‘pdb_files/’ folder.")


if __name__ == '__main__':
    main()
