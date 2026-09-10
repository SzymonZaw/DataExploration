"""Download prospectively selected Stage 2.11C H3/H2 control inputs.

Acquisition only. No H3 similarity statistics are computed here.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from urllib.request import Request, urlopen

CONTROLS = {
    "GSE263713": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263713/suppl/GSE263713_raw_counts.tsv.gz",
        "filename": "GSE263713_raw_counts.tsv.gz",
        "role": "H3 orthogonal temporal candidate",
    },
    "GSE129486_gene_tpm": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129486/suppl/GSE129486_rnaseq-data-1_gene-tpm.tsv.gz",
        "filename": "GSE129486_rnaseq-data-1_gene-tpm.tsv.gz",
        "role": "H2 inflammatory nuisance control",
    },
    "GSE129486_metadata": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129486/suppl/GSE129486_rnaseq-data-1_metadata.tsv.gz",
        "filename": "GSE129486_rnaseq-data-1_metadata.tsv.gz",
        "role": "H2 inflammatory nuisance control metadata",
    },
}


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = Request(url, headers={"User-Agent": "DataExploration-H3/1.0"})
    with urlopen(request, timeout=60) as response, temporary.open("wb") as handle:
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)
    temporary.replace(destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    checksum_lines = []
    for name, item in CONTROLS.items():
        destination = args.data_dir / item["filename"]
        if destination.exists() and not args.force:
            print(f"EXISTS {name}: {destination}")
        else:
            print(f"DOWNLOAD {name}: {item['url']}")
            download(item["url"], destination)
        digest = sha256(destination)
        checksum_lines.append(f"{digest}  {destination.as_posix()}")
        print(f"SHA256 {name}: {digest}  ROLE={item['role']}")

    checksum_path = args.data_dir / "H3_CONTROL_INPUTS.sha256"
    checksum_path.write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    print(f"WROTE {checksum_path}")
    print("Acquisition complete. No H3 similarity analysis was performed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
