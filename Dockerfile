# Track 2 test image. Offline by default: runs the unit tests only (no data download, no patient data).
FROM python:3.10-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONPATH=src

WORKDIR /app

COPY requirements-track2.lock ./
RUN python -m pip install --upgrade pip && pip install -r requirements-track2.lock

COPY pytest.ini README.md THIRD_PARTY_LICENSES.md ./
COPY config ./config
COPY protocol ./protocol
COPY schemas ./schemas
COPY src ./src
COPY tests ./tests
COPY out ./out

# Default: run the offline test suite. To reproduce the analysis, run
#   docker run --rm -v "$PWD/data:/app/data" -v "$PWD/out:/app/out" IMAGE python -m track2.run_pipeline ...
CMD ["python", "-m", "pytest", "-q"]
