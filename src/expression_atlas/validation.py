"""Validation utilities for Expression Atlas."""

import re
from collections.abc import Sequence

from expression_atlas.exceptions import InvalidAccessionError

# Pattern: E-XXXX-#### where XXXX is 4 word characters and #### is one or more digits
ACCESSION_PATTERN = re.compile(r"^E-\w{4}-\d+$")


def is_valid_accession(accession: str) -> bool:
    """
    Check if experiment accession matches expected ArrayExpress/BioStudies format.

    Valid format: E-XXXX-#### (e.g., E-MTAB-1624, E-GEOD-11175)

    Parameters
    ----------
    accession : str
        The experiment accession to validate.

    Returns
    -------
    bool
        True if valid, False otherwise.

    Examples
    --------
    >>> is_valid_accession("E-MTAB-1624")
    True
    >>> is_valid_accession("E-GEOD-11175")
    True
    >>> is_valid_accession("DRP000391")
    False
    >>> is_valid_accession("invalid")
    False
    """
    if not accession:
        return False
    return bool(ACCESSION_PATTERN.match(accession))


def validate_accession(accession: str) -> str:
    """
    Validate accession and raise error if invalid.

    Parameters
    ----------
    accession : str
        The experiment accession to validate.

    Returns
    -------
    str
        The validated accession (unchanged if valid).

    Raises
    ------
    InvalidAccessionError
        If the accession format is invalid.
    """
    if not is_valid_accession(accession):
        raise InvalidAccessionError(accession)
    return accession


def filter_valid_accessions(
    accessions: Sequence[str],
    raise_on_invalid: bool = False,
) -> list[str]:
    """
    Filter a list of accessions to only include valid ones.

    Parameters
    ----------
    accessions : Sequence[str]
        List of experiment accessions to filter.
    raise_on_invalid : bool, optional
        If True, raise error on first invalid accession.
        If False (default), silently skip invalid accessions.

    Returns
    -------
    list[str]
        List containing only valid accessions.

    Raises
    ------
    InvalidAccessionError
        If raise_on_invalid is True and an invalid accession is found.
    """
    valid = []
    for acc in accessions:
        if is_valid_accession(acc):
            valid.append(acc)
        elif raise_on_invalid:
            raise InvalidAccessionError(acc)
    return valid
