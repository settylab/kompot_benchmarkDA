"""
Logging utility for benchmarkDA.
Provides structured, colorized logging with clear severity levels.
Built on Python's standard logging module for better interoperability.
"""

import sys
import logging
from pathlib import Path
from datetime import datetime


# ANSI color codes
class Colors:
    RESET = "\033[0m"
    BOLD = "\033[1m"

    # Standard colors
    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"

    # Bright colors
    BRIGHT_BLACK = "\033[90m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    BRIGHT_WHITE = "\033[97m"


class ColoredFormatter(logging.Formatter):
    """Formatter that adds colors to terminal output."""

    COLORS = {
        "DEBUG": Colors.BRIGHT_BLACK,
        "INFO": Colors.BLUE,
        "STEP": f"{Colors.BOLD}{Colors.CYAN}",
        "SUCCESS": Colors.GREEN,
        "WARNING": Colors.YELLOW,
        "ERROR": Colors.RED,
        "CRITICAL": f"{Colors.BOLD}{Colors.RED}",
    }

    def __init__(self):
        super().__init__(
            fmt="[%(asctime)s] %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )

    def format(self, record):
        # Add color if terminal supports it
        if sys.stdout.isatty():
            # Save original levelname
            original_levelname = record.levelname

            # Color the levelname
            if original_levelname in self.COLORS:
                record.levelname = f"{self.COLORS[original_levelname]}{original_levelname}{Colors.RESET}"

            # Format the record
            formatted = super().format(record)

            # Color the timestamp
            formatted = formatted.replace("[", f"{Colors.BRIGHT_BLACK}[", 1)
            formatted = formatted.replace("]", f"]{Colors.RESET}", 1)

            # Restore original levelname for other handlers
            record.levelname = original_levelname

            return formatted
        return super().format(record)


class PlainFormatter(logging.Formatter):
    """Plain formatter without colors for file output."""

    def __init__(self):
        super().__init__(
            fmt="[%(asctime)s] %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )


class BenchmarkLogger:
    """
    Simple, clear logger for benchmarking workflows.
    Built on Python's standard logging module.

    Levels:
    - DEBUG: Detailed diagnostic information
    - INFO: General informational messages
    - STEP: Major workflow steps (e.g., "Generating labels", "Running method X")
    - SUCCESS: Successful completion messages
    - WARNING: Warning messages
    - ERROR: Error messages
    - CRITICAL: Critical failures
    """

    # Custom log levels
    STEP = 25  # Between INFO (20) and WARNING (30)
    SUCCESS = 22  # Between INFO (20) and STEP (25)

    def __init__(self, name: str, log_file: Path = None, verbose: bool = True):
        """
        Initialize logger.

        Args:
            name: Logger name (usually script or module name)
            log_file: Optional file path to write logs to
            verbose: If True, show DEBUG messages
        """
        self.name = name
        self.verbose = verbose
        self.log_file = log_file

        # Create underlying Python logger
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.DEBUG if verbose else logging.INFO)
        self.logger.propagate = False  # Don't propagate to root logger

        # Add custom log levels
        logging.addLevelName(self.STEP, "STEP")
        logging.addLevelName(self.SUCCESS, "SUCCESS")

        # Clear any existing handlers
        self.logger.handlers.clear()

        # Add console handler with color formatter
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
        console_handler.setFormatter(ColoredFormatter())
        self.logger.addHandler(console_handler)

        # Add file handler if requested
        if log_file:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(PlainFormatter())
            self.logger.addHandler(file_handler)
            self.file_handler = file_handler
        else:
            self.file_handler = None

    def debug(self, message: str):
        """Log debug message (only if verbose=True)."""
        self.logger.debug(message)

    def info(self, message: str):
        """Log informational message."""
        self.logger.info(message)

    def step(self, message: str):
        """Log major workflow step."""
        self.logger.log(self.STEP, message)

    def success(self, message: str):
        """Log success message."""
        self.logger.log(self.SUCCESS, message)

    def warning(self, message: str):
        """Log warning message."""
        self.logger.warning(message)

    def error(self, message: str):
        """Log error message."""
        self.logger.error(message)

    def critical(self, message: str):
        """Log critical failure."""
        self.logger.critical(message)

    def progress(self, current: int, total: int, item: str, status: str = None):
        """
        Log progress for iterative operations.

        Args:
            current: Current iteration number
            total: Total iterations
            item: Description of current item
            status: Optional status (OK, FAIL, SKIP, etc.)
        """
        # Map status to log level
        status_levels = {
            "OK": self.SUCCESS,
            "SUCCESS": self.SUCCESS,
            "FAIL": logging.ERROR,
            "ERROR": logging.ERROR,
            "SKIP": logging.WARNING,
            "TIMEOUT": logging.WARNING,
        }

        if status:
            level = status_levels.get(status, logging.INFO)
            message = f"[{status} {current}/{total}] {item}"
            self.logger.log(level, message)
        else:
            message = f"[{current}/{total}] {item}"
            self.logger.info(message)

    def section(self, title: str, width: int = 60):
        """Print a section header."""
        separator = "=" * width

        # Use plain print for sections to maintain formatting
        if sys.stdout.isatty():
            terminal_msg = f"{Colors.BOLD}{Colors.BLUE}{separator}{Colors.RESET}\n{Colors.BOLD}{Colors.BLUE}{title}{Colors.RESET}\n{Colors.BOLD}{Colors.BLUE}{separator}{Colors.RESET}"
            print(terminal_msg, flush=True)
        else:
            print(f"{separator}\n{title}\n{separator}", flush=True)

        # Write to file handler if exists
        if self.file_handler:
            file_msg = f"{separator}\n{title}\n{separator}"
            self.file_handler.stream.write(file_msg + "\n")
            self.file_handler.stream.flush()

    def close(self):
        """Close log handlers."""
        handlers = self.logger.handlers[:]
        for handler in handlers:
            handler.close()
            self.logger.removeHandler(handler)

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - close handlers."""
        self.close()


def get_logger(
    name: str, log_file: Path = None, verbose: bool = True
) -> BenchmarkLogger:
    """
    Get a logger instance.

    Args:
        name: Logger name
        log_file: Optional log file path
        verbose: Show debug messages

    Returns:
        BenchmarkLogger instance
    """
    return BenchmarkLogger(name, log_file, verbose)
