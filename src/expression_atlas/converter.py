"""Client for the Expression Atlas RData Converter AWS service."""

from __future__ import annotations

import io
import json
import logging
import os
import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class ConverterError(Exception):
    """Error from the converter service."""

    def __init__(self, message: str, status_code: int | None = None):
        super().__init__(message)
        self.status_code = status_code


@dataclass
class ConvertedBundle:
    """Container for converted Expression Atlas data."""

    # Expression matrix (genes x samples)
    matrix: np.ndarray | None = None

    # Gene annotations
    genes: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Sample annotations
    samples: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Metadata from conversion
    meta: dict[str, Any] = field(default_factory=dict)

    # Row names (gene IDs)
    rownames: list[str] = field(default_factory=list)

    # Column names (sample IDs)
    colnames: list[str] = field(default_factory=list)

    @property
    def shape(self) -> tuple[int, int]:
        """Return (n_genes, n_samples) shape."""
        if self.matrix is not None:
            return self.matrix.shape
        return (len(self.rownames), len(self.colnames))


class ConverterClient:
    """
    Client for the Expression Atlas RData Converter service (AWS).

    This client calls the AWS App Runner/ECS service to convert .RData files
    to portable formats that Python can read without R.

    Parameters
    ----------
    service_url : str, optional
        URL of the converter service. Defaults to CONVERTER_URL env var.
    use_iam_auth : bool
        If True, use AWS IAM authentication (SigV4).
        If False, use API key from CONVERTER_API_KEY env var.
    cache_dir : Path, optional
        Directory to cache downloaded bundles. Defaults to temp dir.
    timeout : int
        Request timeout in seconds.

    Examples
    --------
    >>> client = ConverterClient(use_iam_auth=False)
    >>> bundle = client.convert_and_load(
    ...     "ftp://ftp.ebi.ac.uk/.../E-MTAB-7841-atlasExperimentSummary.Rdata",
    ...     "E-MTAB-7841"
    ... )
    >>> print(bundle.matrix.shape)
    (58735, 48)
    """

    def __init__(
        self,
        service_url: str | None = None,
        use_iam_auth: bool = False,
        cache_dir: Path | None = None,
        timeout: int = 600,
    ):
        self.service_url = service_url or os.environ.get("CONVERTER_URL", "")
        self.use_iam_auth = use_iam_auth
        self.cache_dir = cache_dir or Path(tempfile.gettempdir()) / "atlas_converter_cache"
        self.timeout = timeout

        if not self.service_url:
            logger.warning(
                "CONVERTER_URL not set. Converter client will not work. "
                "Set CONVERTER_URL environment variable or pass service_url parameter."
            )

    def _get_auth_headers(self) -> dict[str, str]:
        """Get authentication headers for the request."""
        headers = {"Content-Type": "application/json"}

        if self.use_iam_auth:
            # Use AWS SigV4 signing for IAM auth
            try:
                from botocore.auth import SigV4Auth
                from botocore.awsrequest import AWSRequest
                from botocore.session import Session

                session = Session()
                credentials = session.get_credentials()
                if credentials:
                    # Create a request to sign
                    aws_request = AWSRequest(
                        method="POST",
                        url=f"{self.service_url.rstrip('/')}/convert",
                        headers=headers,
                    )
                    SigV4Auth(credentials, "execute-api", os.environ.get("AWS_REGION", "us-east-1")).add_auth(aws_request)
                    headers.update(dict(aws_request.headers))
            except ImportError:
                logger.warning(
                    "botocore not installed. Cannot use IAM auth. "
                    "Install with: pip install botocore"
                )
            except Exception as e:
                logger.warning(f"Failed to sign request: {e}")
        else:
            # Use API key
            api_key = os.environ.get("CONVERTER_API_KEY", "")
            if api_key:
                headers["X-API-Key"] = api_key

        return headers

    def convert(
        self,
        rdata_url: str,
        accession: str,
        output_format: str = "mtx_bundle",
        assay_name: str | None = None,
        force: bool = False,
    ) -> dict[str, Any]:
        """
        Request conversion of an .RData file.

        Parameters
        ----------
        rdata_url : str
            URL to the .RData file.
        accession : str
            Experiment accession (e.g., E-MTAB-7841).
        output_format : str
            Output format (mtx_bundle or tsv_bundle).
        assay_name : str, optional
            Specific assay to extract.
        force : bool
            Force re-conversion even if cached.

        Returns
        -------
        dict
            Response from the converter service including signed_url.
        """
        if not self.service_url:
            raise ConverterError("CONVERTER_URL not configured")

        endpoint = f"{self.service_url.rstrip('/')}/convert"

        payload = {
            "rdata_url": rdata_url,
            "accession": accession,
            "output_format": output_format,
            "force": force,
        }
        if assay_name:
            payload["assay_name"] = assay_name

        headers = self._get_auth_headers()

        logger.info(f"Requesting conversion for {accession}")

        try:
            req = Request(
                endpoint,
                data=json.dumps(payload).encode("utf-8"),
                headers=headers,
                method="POST",
            )
            with urlopen(req, timeout=self.timeout) as response:
                result = json.loads(response.read().decode("utf-8"))

            if result.get("status") == "error":
                raise ConverterError(
                    result.get("detail", result.get("error", "Unknown error"))
                )

            logger.info(
                f"Conversion {'cache hit' if result.get('cache_hit') else 'complete'} "
                f"for {accession}"
            )
            return result

        except Exception as e:
            if isinstance(e, ConverterError):
                raise
            raise ConverterError(f"Request failed: {e}") from e

    def download_bundle(self, signed_url: str, accession: str) -> Path:
        """
        Download converted bundle from signed URL.

        Parameters
        ----------
        signed_url : str
            Signed URL from convert() response.
        accession : str
            Experiment accession (for cache path).

        Returns
        -------
        Path
            Path to extracted bundle directory.
        """
        # Create cache directory
        bundle_dir = self.cache_dir / accession
        bundle_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Downloading bundle for {accession}")

        try:
            with urlopen(signed_url, timeout=self.timeout) as response:
                zip_data = response.read()

            # Extract zip
            with zipfile.ZipFile(io.BytesIO(zip_data)) as zf:
                zf.extractall(bundle_dir)

            logger.info(f"Bundle extracted to {bundle_dir}")
            return bundle_dir

        except Exception as e:
            raise ConverterError(f"Failed to download bundle: {e}") from e

    def load_bundle(self, bundle_dir: Path) -> dict[str, ConvertedBundle]:
        """
        Load converted data from bundle directory.

        Parameters
        ----------
        bundle_dir : Path
            Path to extracted bundle directory.

        Returns
        -------
        dict[str, ConvertedBundle]
            Dict mapping dataset name to ConvertedBundle.
        """
        results = {}

        # Find all dataset directories
        for item in bundle_dir.iterdir():
            if item.is_dir() and item.name.startswith("dataset_"):
                dataset_name = item.name.replace("dataset_", "")
                results[dataset_name] = self._load_dataset(item)

        # Load metadata
        meta_path = bundle_dir / "meta.json"
        if meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
            # Attach to each bundle
            for bundle in results.values():
                bundle.meta = meta

        return results

    def _load_dataset(self, dataset_dir: Path) -> ConvertedBundle:
        """Load a single dataset from directory."""
        bundle = ConvertedBundle()

        # Load matrix
        mtx_path = dataset_dir / "matrix.mtx"
        tsv_path = dataset_dir / "counts.tsv.gz"

        if mtx_path.exists():
            bundle.matrix = self._load_mtx(mtx_path)
        elif tsv_path.exists():
            df = pd.read_csv(tsv_path, sep="\t", index_col=0, compression="gzip")
            bundle.matrix = df.values
            bundle.rownames = df.index.tolist()
            bundle.colnames = df.columns.tolist()

        # Load row/column names from separate files if MTX format
        barcodes_path = dataset_dir / "barcodes.tsv"
        features_path = dataset_dir / "features.tsv"

        if barcodes_path.exists():
            bundle.colnames = pd.read_csv(barcodes_path, header=None)[0].tolist()
        if features_path.exists():
            bundle.rownames = pd.read_csv(features_path, header=None)[0].tolist()

        # Load genes (rowData)
        genes_path = dataset_dir / "genes.csv"
        if genes_path.exists():
            bundle.genes = pd.read_csv(genes_path, index_col=0)
            if not bundle.rownames:
                bundle.rownames = bundle.genes.index.tolist()

        # Load samples (colData)
        samples_path = dataset_dir / "samples.csv"
        if samples_path.exists():
            bundle.samples = pd.read_csv(samples_path, index_col=0)
            if not bundle.colnames:
                bundle.colnames = bundle.samples.index.tolist()

        return bundle

    def _load_mtx(self, mtx_path: Path) -> np.ndarray:
        """Load Matrix Market file."""
        try:
            from scipy.io import mmread

            sparse_matrix = mmread(str(mtx_path))
            return sparse_matrix.toarray()
        except ImportError:
            logger.warning("scipy not installed, cannot load MTX files efficiently")
            # Fallback: manual parsing (slow)
            return self._parse_mtx_manual(mtx_path)

    def _parse_mtx_manual(self, mtx_path: Path) -> np.ndarray:
        """Parse MTX file manually (fallback if scipy not available)."""
        with open(mtx_path) as f:
            # Skip comments
            line = f.readline()
            while line.startswith("%"):
                line = f.readline()

            # Read dimensions
            parts = line.strip().split()
            nrows, ncols, _ = int(parts[0]), int(parts[1]), int(parts[2])

            # Create dense matrix
            matrix = np.zeros((nrows, ncols))

            # Read entries
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    i, j, val = int(parts[0]) - 1, int(parts[1]) - 1, float(parts[2])
                    matrix[i, j] = val

        return matrix

    def convert_and_load(
        self,
        rdata_url: str,
        accession: str,
        output_format: str = "mtx_bundle",
        assay_name: str | None = None,
        force: bool = False,
    ) -> dict[str, ConvertedBundle]:
        """
        Convert .RData and load the result in one call.

        Parameters
        ----------
        rdata_url : str
            URL to the .RData file.
        accession : str
            Experiment accession.
        output_format : str
            Output format.
        assay_name : str, optional
            Specific assay to extract.
        force : bool
            Force re-conversion.

        Returns
        -------
        dict[str, ConvertedBundle]
            Dict mapping dataset name to ConvertedBundle.
        """
        # Check local cache first
        bundle_dir = self.cache_dir / accession
        meta_path = bundle_dir / "meta.json"

        if not force and meta_path.exists():
            logger.info(f"Loading from local cache: {bundle_dir}")
            return self.load_bundle(bundle_dir)

        # Request conversion
        result = self.convert(rdata_url, accession, output_format, assay_name, force)

        # Download and extract
        bundle_dir = self.download_bundle(result["signed_url"], accession)

        # Load and return
        return self.load_bundle(bundle_dir)

    def is_configured(self) -> bool:
        """Check if the converter client is properly configured."""
        return bool(self.service_url)
