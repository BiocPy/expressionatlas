"""FastAPI server for Expression Atlas RData conversion."""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import tempfile
import time
import urllib.request
import zipfile
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path

from fastapi import Depends, FastAPI, Header, HTTPException, Request
from fastapi.responses import JSONResponse

from config import Settings, get_settings
from schema import (
    ConversionMeta,
    ConvertRequest,
    ConvertResponse,
    DatasetInfo,
    ErrorResponse,
)
from security import SecurityError, validate_url, verify_api_key
from storage import (
    check_cache,
    compute_cache_key,
    generate_signed_url,
    get_cached_meta,
    get_s3_key,
    upload_bundle,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan handler."""
    logger.info("Starting Expression Atlas Converter Service")
    yield
    logger.info("Shutting down Expression Atlas Converter Service")


app = FastAPI(
    title="Expression Atlas RData Converter",
    description="Converts Expression Atlas .RData files to portable formats",
    version="1.0.0",
    lifespan=lifespan,
)


async def verify_auth(
    x_api_key: str | None = Header(default=None),
    settings: Settings = Depends(get_settings),
) -> bool:
    """Verify request authentication."""
    try:
        return verify_api_key(x_api_key)
    except SecurityError as e:
        raise HTTPException(status_code=401, detail=str(e))


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "timestamp": datetime.now(timezone.utc).isoformat()}


@app.post("/convert", response_model=ConvertResponse)
async def convert_rdata(
    request: ConvertRequest,
    settings: Settings = Depends(get_settings),
    _auth: bool = Depends(verify_auth),
):
    """
    Convert an Expression Atlas .RData file to portable formats.

    Returns a signed URL to download the converted bundle.
    """
    start_time = time.time()

    # Validate URL for security
    try:
        validate_url(request.rdata_url)
    except SecurityError as e:
        logger.warning(f"URL validation failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))

    # Compute cache key
    cache_key = compute_cache_key(
        request.rdata_url,
        request.output_format.value,
        request.assay_name,
    )
    logger.info(
        f"Processing {request.accession}, cache_key={cache_key}, force={request.force}"
    )

    # Check cache
    if not request.force:
        cached_path = check_cache(request.accession, cache_key)
        if cached_path:
            # Get cached metadata
            meta_dict = get_cached_meta(request.accession, cache_key)
            if meta_dict:
                signed_url, expires_at = generate_signed_url(cached_path)
                logger.info(f"Cache hit for {request.accession} in {time.time() - start_time:.2f}s")

                return ConvertResponse(
                    signed_url=signed_url,
                    cache_hit=True,
                    meta=ConversionMeta(
                        accession=request.accession,
                        source_url=request.rdata_url,
                        cache_key=cache_key,
                        datasets=[
                            DatasetInfo(**ds) for ds in meta_dict.get("datasets", {}).values()
                        ],
                        r_version=meta_dict.get("r_version"),
                        bioconductor_version=meta_dict.get("bioconductor_version"),
                        converted_at=meta_dict.get("converted_at", datetime.now(timezone.utc)),
                    ),
                    expires_at=expires_at,
                )

    # Create temp directory for this conversion
    with tempfile.TemporaryDirectory(prefix=f"convert_{request.accession}_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        rdata_path = tmpdir_path / "input.RData"
        outdir_path = tmpdir_path / "output"
        zip_path = tmpdir_path / "bundle.zip"

        # Download .RData file
        download_start = time.time()
        try:
            logger.info(f"Downloading {request.rdata_url}")
            download_rdata(request.rdata_url, rdata_path, settings)
            logger.info(f"Download completed in {time.time() - download_start:.2f}s")
        except Exception as e:
            logger.error(f"Download failed: {e}")
            raise HTTPException(
                status_code=502,
                detail=f"Failed to download .RData file: {e}",
            )

        # Run R extraction
        r_start = time.time()
        try:
            logger.info("Running R extraction")
            run_r_extraction(
                rdata_path,
                outdir_path,
                request.output_format.value,
                request.assay_name,
                settings,
            )
            logger.info(f"R extraction completed in {time.time() - r_start:.2f}s")
        except Exception as e:
            logger.error(f"R extraction failed: {e}")
            raise HTTPException(
                status_code=500,
                detail=f"R extraction failed: {e}",
            )

        # Read meta.json
        meta_path = outdir_path / "meta.json"
        if not meta_path.exists():
            raise HTTPException(
                status_code=500,
                detail="R extraction did not produce meta.json",
            )

        with open(meta_path) as f:
            meta_dict = json.load(f)

        # Add conversion metadata
        meta_dict["accession"] = request.accession
        meta_dict["source_url"] = request.rdata_url
        meta_dict["cache_key"] = cache_key
        meta_dict["converted_at"] = datetime.now(timezone.utc).isoformat()

        # Create zip
        logger.info("Creating zip bundle")
        create_zip_bundle(outdir_path, zip_path)

        # Upload to S3
        upload_start = time.time()
        logger.info("Uploading to S3")
        s3_key = upload_bundle(zip_path, meta_dict, request.accession, cache_key)
        logger.info(f"Upload completed in {time.time() - upload_start:.2f}s")

        # Generate presigned URL
        signed_url, expires_at = generate_signed_url(s3_key)

        total_time = time.time() - start_time
        logger.info(
            f"Conversion complete for {request.accession} in {total_time:.2f}s "
            f"(download={time.time() - download_start:.2f}s, "
            f"R={time.time() - r_start:.2f}s)"
        )

        return ConvertResponse(
            signed_url=signed_url,
            cache_hit=False,
            meta=ConversionMeta(
                accession=request.accession,
                source_url=request.rdata_url,
                cache_key=cache_key,
                datasets=[
                    DatasetInfo(**ds) for ds in meta_dict.get("datasets", {}).values()
                ],
                r_version=meta_dict.get("r_version"),
                bioconductor_version=meta_dict.get("bioconductor_version"),
                converted_at=datetime.now(timezone.utc),
            ),
            expires_at=expires_at,
        )


def download_rdata(url: str, dest_path: Path, settings: Settings) -> None:
    """Download .RData file with size limits and timeout."""
    # Use urllib for both FTP and HTTPS
    try:
        with urllib.request.urlopen(url, timeout=settings.download_timeout_seconds) as response:
            # Check content length if available
            content_length = response.headers.get("Content-Length")
            if content_length:
                size_mb = int(content_length) / (1024 * 1024)
                if size_mb > settings.max_rdata_size_mb:
                    raise ValueError(
                        f"File too large: {size_mb:.1f} MB > {settings.max_rdata_size_mb} MB limit"
                    )

            # Download in chunks
            total_size = 0
            max_bytes = settings.max_rdata_size_mb * 1024 * 1024

            with open(dest_path, "wb") as f:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    total_size += len(chunk)
                    if total_size > max_bytes:
                        raise ValueError(
                            f"File exceeded size limit of {settings.max_rdata_size_mb} MB"
                        )
                    f.write(chunk)

    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to download: {e}")


def run_r_extraction(
    rdata_path: Path,
    outdir_path: Path,
    output_format: str,
    assay_name: str | None,
    settings: Settings,
) -> None:
    """Run R extraction script."""
    # Find extract.R script
    script_dir = Path(__file__).parent / "r"
    extract_script = script_dir / "extract.R"

    if not extract_script.exists():
        raise FileNotFoundError(f"R extraction script not found: {extract_script}")

    # Build command
    cmd = [
        "Rscript",
        str(extract_script),
        "--input",
        str(rdata_path),
        "--outdir",
        str(outdir_path),
        "--format",
        output_format,
    ]

    if assay_name:
        cmd.extend(["--assay_name", assay_name])

    # Run with timeout
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=settings.r_timeout_seconds,
            check=False,
        )

        # Log R output
        if result.stdout:
            for line in result.stdout.strip().split("\n"):
                logger.info(f"R: {line}")
        if result.stderr:
            for line in result.stderr.strip().split("\n"):
                logger.warning(f"R stderr: {line}")

        if result.returncode != 0:
            raise RuntimeError(
                f"R script failed with exit code {result.returncode}: {result.stderr}"
            )

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"R extraction timed out after {settings.r_timeout_seconds}s")


def create_zip_bundle(source_dir: Path, zip_path: Path) -> None:
    """Create zip file from output directory."""
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for file_path in source_dir.rglob("*"):
            if file_path.is_file():
                arcname = file_path.relative_to(source_dir)
                zf.write(file_path, arcname)


@app.exception_handler(Exception)
async def generic_exception_handler(request: Request, exc: Exception):
    """Handle unexpected exceptions."""
    logger.exception(f"Unexpected error: {exc}")
    return JSONResponse(
        status_code=500,
        content=ErrorResponse(
            error="Internal server error",
            detail=str(exc),
        ).model_dump(),
    )


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
