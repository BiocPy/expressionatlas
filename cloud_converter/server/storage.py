"""S3 storage operations for upload and presigned URL generation."""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path

import boto3
from botocore.exceptions import ClientError

from config import get_settings

logger = logging.getLogger(__name__)


def get_s3_client():
    """Get S3 client (uses default credentials from IAM role or env vars)."""
    settings = get_settings()
    return boto3.client("s3", region_name=settings.aws_region)


def compute_cache_key(
    rdata_url: str,
    output_format: str,
    assay_name: str | None,
) -> str:
    """
    Compute deterministic cache key for a conversion request.

    Parameters
    ----------
    rdata_url : str
        URL to the .RData file.
    output_format : str
        Output format (mtx_bundle or tsv_bundle).
    assay_name : str or None
        Specific assay name or None for first.

    Returns
    -------
    str
        SHA256 hash as cache key.
    """
    key_parts = f"{rdata_url}|{output_format}|{assay_name or 'default'}"
    return hashlib.sha256(key_parts.encode()).hexdigest()[:32]


def get_s3_key(accession: str, cache_key: str) -> str:
    """Get S3 key for a conversion bundle."""
    return f"converted/{accession}/{cache_key}/bundle.zip"


def check_cache(accession: str, cache_key: str) -> str | None:
    """
    Check if a converted bundle exists in cache.

    Parameters
    ----------
    accession : str
        Experiment accession.
    cache_key : str
        Cache key for the conversion.

    Returns
    -------
    str or None
        S3 key if exists, None otherwise.
    """
    settings = get_settings()
    s3 = get_s3_client()

    s3_key = get_s3_key(accession, cache_key)

    try:
        s3.head_object(Bucket=settings.s3_bucket_name, Key=s3_key)
        logger.info(f"Cache hit for {accession}/{cache_key}")
        return s3_key
    except ClientError as e:
        if e.response["Error"]["Code"] == "404":
            logger.info(f"Cache miss for {accession}/{cache_key}")
            return None
        raise


def upload_bundle(
    local_zip_path: Path,
    meta_dict: dict,
    accession: str,
    cache_key: str,
) -> str:
    """
    Upload conversion bundle to S3.

    Parameters
    ----------
    local_zip_path : Path
        Path to the local zip file.
    meta_dict : dict
        Metadata dictionary to store alongside.
    accession : str
        Experiment accession.
    cache_key : str
        Cache key for the conversion.

    Returns
    -------
    str
        S3 key path.
    """
    settings = get_settings()
    s3 = get_s3_client()

    # Upload zip
    zip_key = get_s3_key(accession, cache_key)
    s3.upload_file(
        str(local_zip_path),
        settings.s3_bucket_name,
        zip_key,
        ExtraArgs={"ContentType": "application/zip"},
    )
    logger.info(f"Uploaded bundle to s3://{settings.s3_bucket_name}/{zip_key}")

    # Upload meta.json separately for quick inspection
    meta_key = f"converted/{accession}/{cache_key}/meta.json"
    s3.put_object(
        Bucket=settings.s3_bucket_name,
        Key=meta_key,
        Body=json.dumps(meta_dict, indent=2, default=str),
        ContentType="application/json",
    )
    logger.info(f"Uploaded metadata to s3://{settings.s3_bucket_name}/{meta_key}")

    return zip_key


def generate_signed_url(s3_key: str) -> tuple[str, datetime]:
    """
    Generate a presigned URL for downloading an S3 object.

    Parameters
    ----------
    s3_key : str
        S3 object key.

    Returns
    -------
    tuple[str, datetime]
        Presigned URL and expiration time.
    """
    settings = get_settings()
    s3 = get_s3_client()

    expiry_seconds = settings.signed_url_expiry_minutes * 60
    expires_at = datetime.now(timezone.utc) + timedelta(seconds=expiry_seconds)

    url = s3.generate_presigned_url(
        "get_object",
        Params={
            "Bucket": settings.s3_bucket_name,
            "Key": s3_key,
        },
        ExpiresIn=expiry_seconds,
    )

    logger.info(f"Generated presigned URL for {s3_key}, expires at {expires_at}")
    return url, expires_at


def get_cached_meta(accession: str, cache_key: str) -> dict | None:
    """
    Get cached metadata without downloading the full bundle.

    Parameters
    ----------
    accession : str
        Experiment accession.
    cache_key : str
        Cache key.

    Returns
    -------
    dict or None
        Metadata dict if exists.
    """
    settings = get_settings()
    s3 = get_s3_client()

    meta_key = f"converted/{accession}/{cache_key}/meta.json"

    try:
        response = s3.get_object(Bucket=settings.s3_bucket_name, Key=meta_key)
        content = response["Body"].read().decode("utf-8")
        return json.loads(content)
    except ClientError as e:
        if e.response["Error"]["Code"] == "NoSuchKey":
            return None
        raise
