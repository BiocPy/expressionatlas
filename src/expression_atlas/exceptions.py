"""Custom exceptions for Expression Atlas client."""

from __future__ import annotations


class ExpressionAtlasError(Exception):
    """Base exception for Expression Atlas errors."""

    pass


class InvalidAccessionError(ExpressionAtlasError):
    """Raised when an experiment accession is invalid."""

    def __init__(self, accession: str) -> None:
        self.accession = accession
        super().__init__(
            f'"{accession}" does not look like an ArrayExpress/BioStudies '
            "experiment accession. Expected format: E-XXXX-#### (e.g., E-MTAB-1624)"
        )


class DownloadError(ExpressionAtlasError):
    """Raised when experiment download fails."""

    def __init__(self, accession: str, reason: str) -> None:
        self.accession = accession
        self.reason = reason
        super().__init__(
            f"Error downloading experiment {accession}: {reason}. "
            f"Check https://www.ebi.ac.uk/gxa/experiments/{accession} or "
            "contact https://www.ebi.ac.uk/about/contact/support/gxa"
        )


class APIError(ExpressionAtlasError):
    """Raised when BioStudies API request fails."""

    def __init__(self, status_code: int, message: str | None = None) -> None:
        self.status_code = status_code
        self.message = message
        msg = f"BioStudies API error (HTTP {status_code})"
        if message:
            msg += f": {message}"
        msg += ". Please try again later or contact https://www.ebi.ac.uk/about/contact/support/gxa"
        super().__init__(msg)
