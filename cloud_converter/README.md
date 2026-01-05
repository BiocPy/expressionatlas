# Expression Atlas RData Converter - AWS Service

This is an **optional** AWS service that converts Expression Atlas `.RData` files 
to portable formats (Matrix Market, CSV) that Python can read without R installed.

## Architecture

```
┌─────────────────┐     POST /convert      ┌──────────────────────┐
│  Python Client  │ ───────────────────────▶│  App Runner / ECS    │
│  (your package) │                         │  (FastAPI + R)       │
└─────────────────┘                         └──────────┬───────────┘
        │                                              │
        │  Download bundle.zip                         │ Upload
        │  via presigned URL                           ▼
        │                                   ┌──────────────────────┐
        └───────────────────────────────────│      S3 Bucket       │
                                            │  converted/<acc>/    │
                                            └──────────────────────┘
```

## Output Bundle Structure

```
bundle.zip
├── dataset_rnaseq/           # or dataset_<array_design> for microarray
│   ├── matrix.mtx            # Matrix Market format (sparse) or counts.tsv.gz
│   ├── genes.csv             # rowData (gene annotations)
│   ├── samples.csv           # colData (sample annotations)
│   └── barcodes.tsv          # column names (for MTX compatibility)
├── meta.json                 # Provenance, dimensions, assay names
└── features.tsv              # row names (for MTX compatibility)
```

## Setup

### 1. AWS Prerequisites

```bash
# Set your region
export AWS_REGION=us-east-1

# Create ECR repository
aws ecr create-repository --repository-name atlas-converter

# Create S3 bucket for outputs
aws s3 mb s3://expression-atlas-converter --region $AWS_REGION

# Create IAM role for App Runner / ECS (example policy)
aws iam create-role --role-name atlas-converter-role \
    --assume-role-policy-document file://trust-policy.json

# Attach S3 permissions
aws iam put-role-policy --role-name atlas-converter-role \
    --policy-name s3-access \
    --policy-document '{
        "Version": "2012-10-17",
        "Statement": [{
            "Effect": "Allow",
            "Action": ["s3:GetObject", "s3:PutObject", "s3:ListBucket"],
            "Resource": [
                "arn:aws:s3:::expression-atlas-converter",
                "arn:aws:s3:::expression-atlas-converter/*"
            ]
        }]
    }'
```

### 2. Create `.env` file

```bash
cp .env.example .env
# Edit .env with your AWS settings
```

### 3. Build and Push to ECR (using CodeBuild)

Your repository should have:
- `Dockerfile`
- `buildspec.yml`
- `server/` folder with FastAPI code

**Option A: AWS CodeBuild (CI/CD)**

1. Go to AWS Console → CodeBuild → Create build project
2. Source: Connect to GitHub and select your repo
3. Environment: 
   - Managed image → Amazon Linux 2 → Standard → aws/codebuild/amazonlinux2-x86_64-standard:5.0
   - Privileged: ✓ (required for Docker builds)
4. Service role: Create new or use existing with ECR permissions
5. Buildspec: Use buildspec.yml from repository
6. Environment variables:
   - `AWS_ACCOUNT_ID` = your 12-digit account ID
   - `AWS_DEFAULT_REGION` = us-east-1
   - `IMAGE_REPO_NAME` = atlas-converter

**Option B: Manual Build**

```bash
# Login to ECR
aws ecr get-login-password --region $AWS_REGION | docker login --username AWS --password-stdin $AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com

# Build and push
cd cloud_converter
docker build -t atlas-converter .
docker tag atlas-converter:latest $AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com/atlas-converter:latest
docker push $AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com/atlas-converter:latest
```

### 4. Deploy to AWS App Runner

```bash
# Create App Runner service (via CLI)
aws apprunner create-service \
    --service-name atlas-converter \
    --source-configuration '{
        "ImageRepository": {
            "ImageIdentifier": "'$AWS_ACCOUNT_ID'.dkr.ecr.'$AWS_REGION'.amazonaws.com/atlas-converter:latest",
            "ImageRepositoryType": "ECR",
            "ImageConfiguration": {
                "Port": "8080",
                "RuntimeEnvironmentVariables": {
                    "AWS_REGION": "'$AWS_REGION'",
                    "S3_BUCKET_NAME": "expression-atlas-converter"
                }
            }
        },
        "AuthenticationConfiguration": {
            "AccessRoleArn": "arn:aws:iam::'$AWS_ACCOUNT_ID':role/atlas-converter-role"
        }
    }' \
    --instance-configuration '{
        "Cpu": "2 vCPU",
        "Memory": "4 GB"
    }'
```

Or use the AWS Console for a visual setup.

## Local Development

### Install R + Bioconductor (one-time)

```bash
# On Ubuntu/Debian
sudo apt-get install r-base

# In R console
install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "S4Vectors", "Matrix"))
```

### Run locally

```bash
cd cloud_converter/server
pip install -r requirements.txt
uvicorn main:app --reload --port 8080
```

### Test the endpoint

```bash
curl -X POST http://localhost:8080/convert \
    -H "Content-Type: application/json" \
    -d '{"rdata_url": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7841/E-MTAB-7841-atlasExperimentSummary.Rdata", "accession": "E-MTAB-7841"}'
```

## API Reference

### POST /convert

**Request:**
```json
{
    "rdata_url": "ftp://ftp.ebi.ac.uk/.../E-MTAB-7841-atlasExperimentSummary.Rdata",
    "accession": "E-MTAB-7841",
    "output_format": "mtx_bundle",  // optional, default: mtx_bundle
    "assay_name": null,             // optional, uses first assay if null
    "force": false                  // optional, bypass cache if true
}
```

**Response (success):**
```json
{
    "status": "success",
    "signed_url": "https://expression-atlas-converter.s3.amazonaws.com/...",
    "cache_hit": false,
    "meta": {
        "accession": "E-MTAB-7841",
        "datasets": ["rnaseq"],
        "dimensions": {"rnaseq": [58735, 48]},
        "assay_names": {"rnaseq": ["counts"]}
    },
    "expires_at": "2026-01-05T15:30:00Z"
}
```

**Response (error):**
```json
{
    "status": "error",
    "error": "Conversion failed",
    "detail": "R extraction error: unsupported object type"
}
```

## Python Client Usage

```python
from cloud_converter.client import ConverterClient

# Initialize client
client = ConverterClient(
    service_url="https://your-service.us-east-1.awsapprunner.com",
    use_iam_auth=False  # Use API key instead
)

# Convert and download
bundle_path = client.convert_and_download(
    rdata_url="ftp://ftp.ebi.ac.uk/.../E-MTAB-7841-atlasExperimentSummary.Rdata",
    accession="E-MTAB-7841"
)

# Load the converted data
data = client.load_bundle(bundle_path)
# data.matrix - scipy sparse matrix
# data.genes - pandas DataFrame
# data.samples - pandas DataFrame
# data.meta - dict
```

## Security Features

- **URL Allowlist**: Only EBI/Expression Atlas FTP/HTTPS URLs allowed
- **SSRF Prevention**: Blocks internal IPs, localhost, metadata endpoints
- **Size Limits**: Rejects .RData files over configured limit
- **Authentication**: Requires API key (or IAM for internal services)
- **Presigned URLs**: Time-limited access to converted bundles in S3

## Cost Optimization

- **Caching**: Converted bundles cached in S3 by deterministic hash
- **App Runner**: Auto-scales to zero when not in use
- **S3 Lifecycle**: Consider adding lifecycle rules to expire old conversions
