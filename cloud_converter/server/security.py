"""Security utilities for SSRF prevention and URL validation."""

from __future__ import annotations

import ipaddress
import socket
from urllib.parse import urlparse

from config import get_settings


class SecurityError(Exception):
    """Raised when a security check fails."""

    pass


# Internal/private IP ranges that should never be accessed
BLOCKED_IP_RANGES = [
    ipaddress.ip_network("10.0.0.0/8"),  # Private
    ipaddress.ip_network("172.16.0.0/12"),  # Private
    ipaddress.ip_network("192.168.0.0/16"),  # Private
    ipaddress.ip_network("127.0.0.0/8"),  # Loopback
    ipaddress.ip_network("169.254.0.0/16"),  # Link-local
    ipaddress.ip_network("100.64.0.0/10"),  # Carrier-grade NAT
    ipaddress.ip_network("::1/128"),  # IPv6 loopback
    ipaddress.ip_network("fc00::/7"),  # IPv6 private
    ipaddress.ip_network("fe80::/10"),  # IPv6 link-local
]

# Cloud metadata endpoints (must block for SSRF prevention)
BLOCKED_HOSTS = [
    # AWS
    "169.254.169.254",  # AWS EC2 metadata
    "169.254.170.2",    # AWS ECS task metadata
    # GCP
    "metadata.google.internal",
    "metadata",
    # Azure
    "169.254.169.254",  # Azure instance metadata
]


def validate_url(url: str) -> str:
    """
    Validate URL for security.

    Checks:
    - URL scheme is allowed (https or ftp)
    - Domain is in allowlist
    - Resolved IP is not internal/private
    - Not a cloud metadata endpoint

    Parameters
    ----------
    url : str
        URL to validate.

    Returns
    -------
    str
        The validated URL.

    Raises
    ------
    SecurityError
        If URL fails any security check.
    """
    settings = get_settings()

    # Parse URL
    try:
        parsed = urlparse(url)
    except Exception as e:
        raise SecurityError(f"Invalid URL format: {e}")

    # Check scheme
    if parsed.scheme not in ("https", "ftp"):
        raise SecurityError(f"URL scheme '{parsed.scheme}' not allowed. Use https:// or ftp://")

    # Check hostname exists
    hostname = parsed.hostname
    if not hostname:
        raise SecurityError("URL must have a hostname")

    # Check against blocked hosts
    hostname_lower = hostname.lower()
    if hostname_lower in BLOCKED_HOSTS:
        raise SecurityError(f"Access to '{hostname}' is blocked")

    # Check domain allowlist
    domain_allowed = False
    for allowed in settings.allowed_domains:
        if hostname_lower == allowed or hostname_lower.endswith(f".{allowed}"):
            domain_allowed = True
            break

    if not domain_allowed:
        raise SecurityError(
            f"Domain '{hostname}' not in allowlist. "
            f"Allowed domains: {settings.allowed_domains}"
        )

    # Resolve hostname and check IP
    try:
        # Get all IPs for the hostname
        addr_info = socket.getaddrinfo(hostname, None, socket.AF_UNSPEC, socket.SOCK_STREAM)
        ips = [info[4][0] for info in addr_info]
    except socket.gaierror as e:
        raise SecurityError(f"Failed to resolve hostname '{hostname}': {e}")

    for ip_str in ips:
        try:
            ip = ipaddress.ip_address(ip_str)
            for blocked_range in BLOCKED_IP_RANGES:
                if ip in blocked_range:
                    raise SecurityError(
                        f"Resolved IP '{ip}' is in blocked range {blocked_range}"
                    )
        except ValueError:
            # Not a valid IP, skip
            continue

    # All checks passed
    return url


def verify_api_key(provided_key: str | None) -> bool:
    """
    Verify API key if API key auth is enabled.

    Parameters
    ----------
    provided_key : str or None
        The API key from the request header.

    Returns
    -------
    bool
        True if valid or API key auth is disabled.

    Raises
    ------
    SecurityError
        If API key is required but missing or invalid.
    """
    settings = get_settings()

    # If no API key configured, skip check (assume IAM auth)
    if not settings.converter_api_key:
        return True

    if not provided_key:
        raise SecurityError("API key required but not provided")

    # Constant-time comparison to prevent timing attacks
    import secrets

    if not secrets.compare_digest(provided_key, settings.converter_api_key):
        raise SecurityError("Invalid API key")

    return True
