# Installation

This section provides instructions for installing the `pfdha` software and its dependencies.

## Prerequisites

Before installing `pfdha`, ensure you have the following prerequisites met:

-   **Python**: `pfdha` requires Python 3.11 or newer. You can check your Python version by running:
    ```bash
    python3 --version
    ```
-   **Git**: Required to clone the project repository.

## Installation Steps

!!! tip "Use a Virtual Environment"
    We strongly recommend installing `pfdha` in a dedicated Python virtual environment to avoid conflicts with other packages or system-level Python installations.

1.  **Clone the Repository**:
    First, clone the `pfdha` repository from GitHub to your local machine:
    ```bash
    git clone https://github.com/vup1120/oq-pfdha.git
    cd pfdha
    ```

2.  **Create and Activate a Virtual Environment**:
    From the root of the `pfdha` directory, create and activate a new virtual environment:
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```
    You will see `(venv)` prefixed to your shell prompt, indicating the environment is active.

3.  **Install the Package**:
    Install `pfdha` and its Python dependencies using `pip`. The `-e` flag installs the package in "editable" mode, which is convenient for development as it links the `fdha` command to your source code directory.
    ```bash
    pip install -e .
    ```

## Verifying the Installation

After the installation is complete, you can verify it by checking the help message of the `fdha` command:

```bash
fdha --help
```

This should print the main help message for the tool, confirming that the installation was successful and the command-line interface is accessible on your path.

With the installation complete, you are now ready to run your first calculation. Proceed to the [Quick Start](02-QuickStart.md) guide.
