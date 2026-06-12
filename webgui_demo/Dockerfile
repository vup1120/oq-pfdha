# syntax=docker/dockerfile:1
# Dockerfile for the oq-pfdha web GUI (real engine runs).
#
# IMPORTANT: build from the REPOSITORY ROOT, because the image installs the
# full oq-pfdha package:
#
#     docker build -f webgui_demo/Dockerfile -t pfdha-gui .
#     docker run -p 7860:7860 pfdha-gui
#
# Layout matches Hugging Face Spaces (Docker SDK): port 7860, user 1000.
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential git \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 user
USER user
ENV PATH="/home/user/.local/bin:$PATH"
WORKDIR /app

# Engine first (resolves its own numpy/pandas stack), then GUI libraries.
COPY --chown=user . .
RUN pip install --no-cache-dir --user -e . \
    && pip install --no-cache-dir --user -r webgui_demo/requirements.txt

EXPOSE 7860
CMD ["streamlit", "run", "webgui_demo/app.py", \
     "--server.port", "7860", \
     "--server.address", "0.0.0.0", \
     "--server.headless", "true", \
     "--browser.gatherUsageStats", "false"]
