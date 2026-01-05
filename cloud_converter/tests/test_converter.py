"""Tests for the Cloud Converter service."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

# Import from parent directory
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from server.schema import ConvertRequest, ConvertResponse, DatasetInfo
from server.security import URLValidator


class TestURLValidator:
    """Tests for URL security validation."""

    def test_valid_ftp_url(self) -> None:
        """Valid Expression Atlas FTP URLs should pass."""
        url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7841/E-MTAB-7841-atlasExperimentSummary.Rdata"
        assert URLValidator.validate_url(url) is True

    def test_valid_https_url(self) -> None:
        """Valid HTTPS URLs should pass."""
        url = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7841/E-MTAB-7841-atlasExperimentSummary.Rdata"
        assert URLValidator.validate_url(url) is True

    def test_blocked_localhost(self) -> None:
        """Localhost URLs should be blocked (SSRF prevention)."""
        urls = [
            "http://localhost/test.Rdata",
            "http://127.0.0.1/test.Rdata",
            "ftp://127.0.0.1/test.Rdata",
            "http://[::1]/test.Rdata",
        ]
        for url in urls:
            assert URLValidator.validate_url(url) is False, f"Should block: {url}"

    def test_blocked_private_ip(self) -> None:
        """Private IP ranges should be blocked (SSRF prevention)."""
        urls = [
            "http://192.168.1.1/test.Rdata",
            "http://10.0.0.1/test.Rdata",
            "http://172.16.0.1/test.Rdata",
            "ftp://192.168.0.1/test.Rdata",
        ]
        for url in urls:
            assert URLValidator.validate_url(url) is False, f"Should block: {url}"

    def test_blocked_metadata_endpoint(self) -> None:
        """Cloud metadata endpoints should be blocked (SSRF prevention)."""
        urls = [
            "http://169.254.169.254/latest/meta-data/",
            "http://metadata.google.internal/",
        ]
        for url in urls:
            assert URLValidator.validate_url(url) is False, f"Should block: {url}"

    def test_blocked_non_rdata_extension(self) -> None:
        """Non-.Rdata files should be blocked."""
        urls = [
            "ftp://ftp.ebi.ac.uk/file.txt",
            "https://ftp.ebi.ac.uk/file.zip",
            "ftp://ftp.ebi.ac.uk/file.Rda",  # Wrong extension
        ]
        for url in urls:
            assert URLValidator.validate_url(url) is False, f"Should block: {url}"


class TestConvertRequest:
    """Tests for request validation."""

    def test_valid_request(self) -> None:
        """Valid request should parse correctly."""
        req = ConvertRequest(
            rdata_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7841/E-MTAB-7841-atlasExperimentSummary.Rdata",
            accession="E-MTAB-7841",
        )
        assert req.accession == "E-MTAB-7841"
        assert req.output_format == "mtx_bundle"
        assert req.force is False

    def test_invalid_accession_format(self) -> None:
        """Invalid accession format should be rejected."""
        with pytest.raises(ValueError):
            ConvertRequest(
                rdata_url="ftp://ftp.ebi.ac.uk/test.Rdata",
                accession="invalid",
            )

    def test_invalid_output_format(self) -> None:
        """Invalid output format should be rejected."""
        with pytest.raises(ValueError):
            ConvertRequest(
                rdata_url="ftp://ftp.ebi.ac.uk/test.Rdata",
                accession="E-MTAB-7841",
                output_format="invalid_format",
            )


class TestConvertResponse:
    """Tests for response serialization."""

    def test_success_response(self) -> None:
        """Success response should include all fields."""
        resp = ConvertResponse(
            status="success",
            accession="E-MTAB-7841",
            signed_url="https://expression-atlas-converter.s3.amazonaws.com/converted/E-MTAB-7841/abc123/bundle.zip?X-Amz-Signature=...",
            datasets=[
                DatasetInfo(
                    name="rnaseq",
                    type="SummarizedExperiment",
                    n_genes=58735,
                    n_samples=48,
                    assay_names=["counts"],
                )
            ],
            cache_hit=False,
        )
        data = resp.model_dump()
        assert data["status"] == "success"
        assert len(data["datasets"]) == 1
        assert data["datasets"][0]["name"] == "rnaseq"

    def test_error_response(self) -> None:
        """Error response should include error message."""
        resp = ConvertResponse(
            status="error",
            accession="E-MTAB-7841",
            error="Failed to parse .RData file",
        )
        assert resp.error == "Failed to parse .RData file"


class TestCacheKey:
    """Tests for cache key generation."""

    def test_cache_key_deterministic(self) -> None:
        """Cache key should be deterministic for same inputs."""
        from server.storage import compute_cache_key

        key1 = compute_cache_key(
            rdata_url="ftp://ftp.ebi.ac.uk/test.Rdata",
            output_format="mtx_bundle",
            assay_name=None,
        )
        key2 = compute_cache_key(
            rdata_url="ftp://ftp.ebi.ac.uk/test.Rdata",
            output_format="mtx_bundle",
            assay_name=None,
        )
        assert key1 == key2

    def test_cache_key_unique_for_different_formats(self) -> None:
        """Different formats should produce different cache keys."""
        from server.storage import compute_cache_key

        key1 = compute_cache_key(
            rdata_url="ftp://ftp.ebi.ac.uk/test.Rdata",
            output_format="mtx_bundle",
            assay_name=None,
        )
        key2 = compute_cache_key(
            rdata_url="ftp://ftp.ebi.ac.uk/test.Rdata",
            output_format="tsv_bundle",
            assay_name=None,
        )
        assert key1 != key2


# Integration tests (require running service)
@pytest.mark.integration
class TestConverterIntegration:
    """Integration tests for the full conversion pipeline."""

    def test_convert_rnaseq_experiment(self) -> None:
        """Test converting an RNA-seq experiment."""
        from client.converter_client import ConverterClient

        client = ConverterClient()
        if not client.is_configured():
            pytest.skip("CONVERTER_URL not configured")

        bundles = client.convert_and_load(
            "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7841/E-MTAB-7841-atlasExperimentSummary.Rdata",
            "E-MTAB-7841",
        )

        assert "rnaseq" in bundles
        bundle = bundles["rnaseq"]
        assert bundle.matrix is not None
        assert bundle.matrix.shape[0] > 0  # Has genes
        assert bundle.matrix.shape[1] > 0  # Has samples


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
