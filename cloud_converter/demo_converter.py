#!/usr/bin/env python
"""
Demo script showing how to use the AWS Converter service.

This demonstrates the "escape hatch" for loading .RData files without
having R installed locally, by using an AWS App Runner service that has R.

Prerequisites:
1. Deploy the App Runner service (see README.md)
2. Set environment variables:
   - CONVERTER_URL: App Runner service URL
   - CONVERTER_API_KEY: API key for authentication (if not using IAM)

Usage:
    python demo_converter.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from client.converter_client import ConvertedBundle, ConverterClient


def main() -> None:
    """Demo the AWS Converter client."""
    
    print("=" * 60)
    print("Expression Atlas AWS Converter Demo")
    print("=" * 60)
    
    # Check configuration
    converter_url = os.environ.get("CONVERTER_URL", "")
    if not converter_url:
        print("\n⚠️  CONVERTER_URL not set!")
        print("\nTo use the AWS converter:")
        print("1. Deploy the App Runner service (see cloud_converter/README.md)")
        print("2. Set CONVERTER_URL environment variable to the service URL")
        print("\nExample:")
        print("  export CONVERTER_URL=https://your-service.us-east-1.awsapprunner.com")
        print("\nFor local testing, run the service locally:")
        print("  cd cloud_converter/server")
        print("  pip install -r requirements.txt")
        print("  uvicorn main:app --port 8080")
        print("  export CONVERTER_URL=http://localhost:8080")
        return
    
    print(f"\n✓ Converter URL: {converter_url}")
    
    # Create client
    # use_iam_auth=True requires botocore and AWS credentials
    # use_iam_auth=False requires CONVERTER_API_KEY to be set
    client = ConverterClient(
        use_iam_auth=False,  # Use API key auth
    )
    
    # Example: Convert an RNA-seq experiment
    accession = "E-MTAB-7841"  # Small-ish experiment for demo
    rdata_url = (
        f"ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/"
        f"{accession}/{accession}-atlasExperimentSummary.Rdata"
    )
    
    print(f"\n📥 Converting experiment: {accession}")
    print(f"   Source: {rdata_url}")
    
    try:
        # Convert and load in one call
        bundles = client.convert_and_load(rdata_url, accession)
        
        print(f"\n✓ Conversion successful!")
        print(f"   Datasets found: {list(bundles.keys())}")
        
        # Explore the data
        for name, bundle in bundles.items():
            print(f"\n📊 Dataset: {name}")
            print(f"   Matrix shape: {bundle.shape} (genes × samples)")
            
            if not bundle.genes.empty:
                print(f"   Gene annotations: {list(bundle.genes.columns)}")
            
            if not bundle.samples.empty:
                print(f"   Sample annotations: {list(bundle.samples.columns)}")
            
            # Show first few gene names
            if bundle.rownames:
                print(f"   First 5 genes: {bundle.rownames[:5]}")
            
            # Show first few sample names
            if bundle.colnames:
                print(f"   First 5 samples: {bundle.colnames[:5]}")
            
            # Show metadata
            if bundle.meta:
                print(f"   Metadata: {bundle.meta.get('accession', 'N/A')}")
        
        # Example: Convert the bundle to a SummarizedExperiment for compatibility
        print("\n" + "=" * 60)
        print("Converting to expression-atlas format...")
        
        if "rnaseq" in bundles:
            bundle = bundles["rnaseq"]
            
            # Import the rcompat classes
            from expression_atlas.rcompat import SummarizedExperiment
            
            sumexp = SummarizedExperiment()
            sumexp.assays["counts"] = bundle.matrix
            sumexp.rownames = bundle.rownames
            sumexp.colnames = bundle.colnames
            sumexp.rowData = bundle.genes
            sumexp.colData = bundle.samples
            sumexp.metadata = bundle.meta
            
            print(f"\n✓ Created SummarizedExperiment!")
            print(f"   {sumexp}")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 60)
    print("Demo complete!")


if __name__ == "__main__":
    main()
