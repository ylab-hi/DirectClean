"""Allow running DirectClean as ``python -m directclean``."""

from directclean.cli.app import app

if __name__ == "__main__":
    app()