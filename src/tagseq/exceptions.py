"""Custom exceptions."""


class TagseqError(Exception):
    """Base exception."""


class ConfigError(TagseqError):
    """Invalid or missing configuration."""


class ExternalToolError(TagseqError):
    """External tool returned non-zero."""
